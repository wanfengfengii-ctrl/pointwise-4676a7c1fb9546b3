"""Finite public inputs and executable calls; no flux implementation or oracle."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory

from remap import Archive, Cell, Target
from remap.__main__ import json_value

UNIT = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
CONSTANT = ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))
LINEAR = ((0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))
CASES = ("split", "partial", "quadratic", "independent")


def example(case):
    """Return fresh input containers using only existing Cell/Target types."""
    density = (0, 1, 0, 0) if case == "quadratic" else (1, 0, 0, 0)
    ids = ("s0", "s1") if case == "independent" else ("s",)
    cells = [Cell(s, UNIT, density) for s in ids]
    if case in ("split", "independent"):
        targets = [
            Target("left", (0, 0, 0), (F(1, 2), 1, 1)),
            Target("right", (F(1, 2), 0, 0), (1, 1, 1)),
        ]
    elif case == "partial":
        targets = [
            Target("left", (0, 0, 0), (F(1, 2), 1, 1)),
            Target("r0", (F(1, 2), 0, 0), (1, F(1, 4), 1)),
            Target("r1", (F(1, 2), F(1, 4), 0), (1, 1, 1)),
        ]
    elif case == "quadratic":
        targets = [Target("all", (0, 0, 0), (1, 1, 1))]
    else:
        raise ValueError("unknown example")
    fields = {s: [list(c) for c in
                  (LINEAR if case == "quadratic" else CONSTANT)] for s in ids}
    return cells, targets, fields


def register(archive, case):
    cells, targets, fields = example(case)
    revision = archive.put_source("observations-" + case, cells)
    plan = archive.save_plan("layout-" + case, revision, targets)
    return revision, plan, fields


def baseline():
    expected = {
        "split": {"left": F(7, 48), "right": F(1, 48)},
        "partial": {"left": F(7, 48), "r0": F(7, 384), "r1": F(1, 384)},
        "quadratic": {"all": F(1, 24)},
        "independent": {"left": F(7, 24), "right": F(1, 24)},
    }
    displayed = []
    with Archive(":memory:") as archive:
        for case in CASES:
            revision, plan, _ = register(archive, case)
            result = archive.execute_plan(plan, "mass-" + case)
            actual = {r["target"]: r["mass"] for r in result["targets"]}
            assert actual == expected[case], (case, actual)
            assert all(r["uncovered_mass"] == 0 for r in result["source_ledger"])
            displayed.append({"case": case, "source_revision": revision,
                              "plan_id": plan, "mass_targets": result["targets"]})
    print(json.dumps(json_value(displayed), ensure_ascii=False, sort_keys=True))


def check_example(case, result):
    """A few published analytic expectations, not a general face validator."""
    faces = result["faces"]
    by_id = {f["id"]: f for f in faces}
    assert len(by_id) == len(faces)
    for face in faces:
        if face["mate"] is not None:
            other = by_id[face["mate"]]
            assert other["source"] == face["source"]
            assert other["mate"] == face["id"]
            assert other["flux"] == -face["flux"]
    external = {r["source"]: r["external_net"] for r in result["source_ledger"]}
    if case == "quadratic":
        assert external == {"s": F(1, 12)}
        assert sum((f["flux"] for f in faces if f["origins"]["source_faces"]
                    and all(sum(p) == 1 for p in f["vertices"])), F(0)) == F(1, 12)
    else:
        assert external == ({s: F(0) for s in ("s0", "s1")}
                            if case == "independent" else {"s": F(0)})
        for source in external:
            xzero = [f for f in faces if f["source"] == source
                     and all(p[0] == 0 for p in f["vertices"])]
            assert sum((f["flux"] for f in xzero), F(0)) == -F(1, 2)
            assert all(f["origins"]["source_faces"]
                       and f["origins"]["target_faces"] for f in xzero)
    exchanges = result["exchanges"]
    if case == "split":
        assert sum((e["flux"] for e in exchanges
                    if (e["a"], e["b"]) == ("left", "right")), F(0)) == F(1, 8)
    elif case == "partial":
        for target, expected in (("r0", F(3, 32)), ("r1", F(1, 32))):
            assert sum((e["flux"] for e in exchanges
                        if (e["a"], e["b"]) == ("left", target)), F(0)) == expected
        assert any((e["a"], e["b"]) == ("r0", "r1") and e["flux"] == 0
                   for e in exchanges)
    elif case == "independent":
        for source in ("s0", "s1"):
            assert sum((e["flux"] for e in exchanges if e["source"] == source
                        and (e["a"], e["b"]) == ("left", "right")), F(0)) == F(1, 8)
    else:
        assert exchanges == []


def flux():
    displayed, saved = [], []
    with TemporaryDirectory(prefix="flux-example-") as directory:
        path = str(Path(directory) / "observations.sqlite")
        with Archive(path) as archive:
            for case in CASES:
                revision, plan, fields = register(archive, case)
                # The old and new request namespaces must be independent.
                request = "example-" + case
                mass = archive.execute_plan(plan, request)
                velocity = archive.put_velocity("velocity-" + case, revision, fields)
                fields.clear()
                result = archive.execute_flux(plan, velocity, request)
                check_example(case, result)
                assert archive.execute_flux(plan, velocity, request) == result
                source_ids = [c.id for c in archive.read_source(revision)["cells"]]
                replacement = archive.put_velocity(
                    "velocity-" + case, revision,
                    {s: ((0, 0, 0, 0),) * 3 for s in source_ids})
                assert replacement != velocity
                assert archive.execute_flux(plan, velocity, request) == result
                assert archive.read_result(mass["id"]) == mass
                saved.append((plan, velocity, request, result["id"],
                              json.loads(json.dumps(json_value(result)))))
                displayed.append({"case": case, "result": json_value(result)})
                result["faces"].clear()
                assert archive.read_flux(saved[-1][3])["faces"]
        with Archive(path) as archive:
            for plan, velocity, request, result_id, snapshot in saved:
                assert json_value(archive.read_flux(result_id)) == snapshot
                assert json_value(archive.execute_flux(plan, velocity, request)) == snapshot
        print(json.dumps({"examples": displayed, "reopened": True},
                         ensure_ascii=False, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description="Public remapping/flux examples")
    parser.add_argument("mode", nargs="?", choices=("baseline", "flux"),
                        default="baseline")
    args = parser.parse_args()
    if args.mode == "baseline":
        baseline()
    else:
        flux()


if __name__ == "__main__":
    main()
