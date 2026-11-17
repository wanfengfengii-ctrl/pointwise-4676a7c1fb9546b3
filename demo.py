"""A real persistent round trip, entirely within a disposable directory."""
import json
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
from remap import Archive, Cell, Target
from remap.__main__ import json_value


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
            assert old["transfers"][0]["mass"] == Fraction(11, 384)
            assert old["source_ledger"][0]["uncovered_mass"] == Fraction(5, 384)
            replacement = archive.put_source(source["name"], [
                Cell("cell-A", source["cells"][0]["vertices"], (1, 0, 0, 0))])
            new_plan = archive.save_plan("new-observation", replacement,
                                        [Target("west", (0, 0, 0), ("1/2", 1, 1))])
            new = archive.execute_plan(new_plan, "demo-replacement")
        with Archive(path) as reopened:
            assert reopened.read_result(old["id"]) == old
            assert reopened.execute_plan(plan_id, "demo-west") == old
            assert new["transfers"][0]["mass"] == Fraction(7, 48)
            print(json.dumps(json_value({"historical": old, "replacement": new,
                                         "reopened": True}),
                             ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
