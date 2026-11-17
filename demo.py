"""A real persistent round trip for mass and flux, in a disposable directory.

Demonstrates partial adjacency, double-origin faces, quadratic integration and
a persistent reopen.  Everything is exact rational arithmetic.
"""
import json
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory

from remap import Archive, Cell, Target
from remap.__main__ import json_value

UNIT = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def payload(value):
    return json.loads(json.dumps(json_value(value)))


def main():
    root = Path(__file__).resolve().parent
    source = json.loads((root / "samples/source.json").read_text(encoding="utf-8"))
    plan = json.loads((root / "samples/plan.json").read_text(encoding="utf-8"))
    with TemporaryDirectory(prefix="remap-demo-") as directory:
        path = str(Path(directory) / "observations.sqlite")
        with Archive(path) as archive:
            revision = archive.put_source(
                source["name"], [Cell(c["id"], c["vertices"], c["density"])
                                 for c in source["cells"]])
            plan_id = archive.save_plan(
                plan["name"], revision, [Target(t["id"], t["lower"], t["upper"])
                                        for t in plan["targets"]])
            old = archive.execute_plan(plan_id, "demo-west")
            assert old["transfers"][0]["mass"] == F(11, 384)
            assert old["source_ledger"][0]["uncovered_mass"] == F(5, 384)

            # --- flux: partial adjacency (left touches r0 and r1 only) -----
            partial = archive.put_source("demo-partial", [
                Cell("s", UNIT, (1, 0, 0, 0))])
            partial_plan = archive.save_plan("demo-partial-adjacency", partial, [
                Target("left", (0, 0, 0), (F(1, 2), 1, 1)),
                Target("r0", (F(1, 2), 0, 0), (1, F(1, 4), 1)),
                Target("r1", (F(1, 2), F(1, 4), 0), (1, 1, 1))])
            partial_field = archive.put_velocity("demo-partial-velocity", partial, {
                "s": ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))})
            partial_result = archive.execute_flux(
                partial_plan, partial_field, "demo-partial")
            exchange = {(e["a"], e["b"]): e["flux"]
                        for e in partial_result["exchanges"]}
            assert exchange[("left", "r0")] == F(3, 32)
            assert exchange[("left", "r1")] == F(1, 32)
            assert exchange[("r0", "r1")] == 0

            # --- flux: double-origin faces and a quadratic integrand ------
            split = archive.put_source("demo-split", [
                Cell("s", UNIT, (1, 0, 0, 0))])
            split_plan = archive.save_plan("demo-split-adjacency", split, [
                Target("left", (0, 0, 0), (F(1, 2), 1, 1)),
                Target("right", (F(1, 2), 0, 0), (1, 1, 1))])
            split_field = archive.put_velocity("demo-split-velocity", split, {
                "s": ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))})
            split_result = archive.execute_flux(
                split_plan, split_field, "demo-split")
            x_zero = [face for face in split_result["faces"]
                      if all(p[0] == 0 for p in face["vertices"])]
            assert sum((f["flux"] for f in x_zero), F(0)) == -F(1, 2)
            assert all(f["origins"]["source_faces"] and
                       f["origins"]["target_faces"] for f in x_zero)
            assert sum((e["flux"] for e in split_result["exchanges"]), F(0)) \
                == F(1, 8)

            quadratic = archive.put_source("demo-quadratic", [
                Cell("s", UNIT, (0, 1, 0, 0))])
            quadratic_plan = archive.save_plan("demo-quadratic-plan", quadratic,
                                               [Target("all", (0, 0, 0), (1, 1, 1))])
            quadratic_field = archive.put_velocity(
                "demo-quadratic-velocity", quadratic,
                {"s": ((0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))})
            quadratic_result = archive.execute_flux(
                quadratic_plan, quadratic_field, "demo-quadratic")
            assert quadratic_result["source_ledger"][0]["external_net"] == F(1, 12)
            assert quadratic_result["targets"][0]["net_outward"] == F(1, 12)

            # --- replacements must not disturb the bound snapshots -------
            archive.put_source("demo-split", [
                Cell("s", UNIT, (2, 0, 0, 0))])
            archive.put_velocity("demo-split-velocity", split, {
                "s": ((0, 0, 0, 0),) * 3})
            replacement = archive.put_source(source["name"], [
                Cell("cell-A", source["cells"][0]["vertices"], (1, 0, 0, 0))])
            new_plan = archive.save_plan("new-observation", replacement,
                                        [Target("west", (0, 0, 0), ("1/2", 1, 1))])
            new = archive.execute_plan(new_plan, "demo-replacement")
            flux_snapshot = payload(split_result)
            partial_snapshot = payload(partial_result)
            split_result["faces"].clear()

        # --- reopen: history, bindings and replays are unchanged ---------
        with Archive(path) as reopened:
            assert reopened.read_result(old["id"]) == old
            assert reopened.execute_plan(plan_id, "demo-west") == old
            assert new["transfers"][0]["mass"] == F(7, 48)
            assert payload(reopened.read_flux(split_result["id"])) == flux_snapshot
            assert reopened.execute_flux(
                split_plan, split_field, "demo-split") == reopened.read_flux(
                    split_result["id"])
            assert payload(reopened.read_flux(partial_result["id"])) \
                == partial_snapshot
            print(json.dumps(json_value({
                "historical": old,
                "replacement": new,
                "flux": {
                    "split": {"faces": len(flux_snapshot["faces"]),
                              "exchange": str(sum(
                                  (F(e["flux"]) for e in flux_snapshot["exchanges"]),
                                  F(0))),
                              "x_zero_flux": -F(1, 2)},
                    "partial_adjacency": {
                        "left_r0": str(exchange[("left", "r0")]),
                        "left_r1": str(exchange[("left", "r1")]),
                        "r0_r1": str(exchange[("r0", "r1")])},
                    "quadratic": str(quadratic_result["targets"][0]["net_outward"])},
                "reopened": True}), ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
