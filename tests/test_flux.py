import json
import sqlite3
import unittest
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory

from remap import (Archive, Cell, NotFoundError, RequestConflict, Target,
                   ValidationError)
from tests.common import UNIT, cell, half, whole


def field(v=(1, 0, 0, 0)):
    return {"s": [v, (0, 0, 0, 0), (0, 0, 0, 0)]}


def flat(id="s", vertices=UNIT):
    """Unit tetrahedron carrying rho = 1 (tests.common.cell defaults to rho=x)."""
    return Cell(id, vertices, (1, 0, 0, 0))


class FluxCase(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory(prefix="remap-flux-")
        self.addCleanup(self.directory.cleanup)
        self.path = str(Path(self.directory.name) / "flux.sqlite")
        self.archive = Archive(self.path)
        self.addCleanup(self.archive.close)


class TestVelocityRevision(FluxCase):
    def test_immutable_versions_and_input_isolation(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        fields = {"s": [list(x) for x in field()["s"]]}
        first = a.put_velocity("field", revision, fields)
        fields.clear()
        second = a.put_velocity("field", revision, field((0, 0, 1, 0)))
        self.assertNotEqual(first, second)
        plan = a.save_plan("layout", revision, [whole()])
        old = a.execute_flux(plan, first, "old")
        self.assertEqual(old["velocity_revision"], first)
        self.assertEqual(old["receipt"]["velocity_revision"], first)
        # A later same-name version never replaces the bound one.
        fields = {"s": [list(x) for x in field((5, 0, 0, 0))["s"]]}
        third = a.put_velocity("field", revision, fields)
        self.assertNotIn(third, (first, second))
        self.assertEqual(a.execute_flux(plan, first, "old"), old)

    def test_reject_invalid_fields(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        invalid = [
            {},
            {"s": [(1, 0), (0, 0, 0, 0), (0, 0, 0, 0)]},
            {"s": [(1, 0, 0, 0), (0, 0, 0, 0)]},
            {"s": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0),
                   (0, 0, 0, 0)]},
            {"s": [(True, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]},
            {"s": [(0.5, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]},
            {"s": [("1/0", 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]},
            {"s": [(1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)], "x":
             [(0, 0, 0, 0)] * 3},
            {"": [(0, 0, 0, 0)] * 3, "s": [(0, 0, 0, 0)] * 3},
            [1, 2, 3],
        ]
        for fields in invalid:
            with self.assertRaises(ValidationError):
                a.put_velocity("field", revision, fields)
        with self.assertRaises(ValidationError):
            a.put_velocity("", revision, field())
        with self.assertRaises(ValidationError):
            a.put_velocity("field", True, field())
        with self.assertRaises(ValidationError):
            a.put_velocity("field", 0, field())
        with self.assertRaises(NotFoundError):
            a.put_velocity("field", 999, field())

    def test_bound_to_source_cells(self):
        a = self.archive
        revision = a.put_source("obs", [cell("s"), cell("t")])
        with self.assertRaises(ValidationError):
            a.put_velocity("field", revision, field())
        with self.assertRaises(ValidationError):
            a.put_velocity("field", revision, {
                "s": field()["s"], "t": field()["s"], "x": field()["s"]})


class TestFluxPlanBinding(FluxCase):
    def test_source_binding_must_match(self):
        a = self.archive
        old = a.put_source("obs", [cell()])
        old_plan = a.save_plan("layout", old, [whole()])
        new = a.put_source("obs", [cell(density=(1, 0, 0, 0))])
        new_plan = a.save_plan("layout", new, [whole()])
        velocity = a.put_velocity("field", new, field())
        with self.assertRaises(ValidationError):
            a.execute_flux(old_plan, velocity, "cross")
        self.assertEqual(a.list_results(), [])
        # A same-named cell in another revision is not accepted.
        self.assertEqual(a.execute_flux(new_plan, velocity, "same")[
            "source_revision"], new)

    def test_request_namespaces_and_conflicts(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        plan = a.save_plan("layout", revision, [half()])
        other = a.save_plan("other", revision, [whole()])
        velocity = a.put_velocity("field", revision, field())
        second = a.put_velocity("field", revision, field((0, 1, 0, 0)))
        mass = a.execute_plan(plan, "shared")
        flux = a.execute_flux(plan, velocity, "shared")
        self.assertEqual(mass["receipt"]["request_id"], "shared")
        self.assertEqual(flux["receipt"]["request_id"], "shared")
        self.assertEqual(a.execute_flux(plan, velocity, "shared"), flux)
        with self.assertRaises(RequestConflict):
            a.execute_flux(plan, second, "shared")
        with self.assertRaises(RequestConflict):
            a.execute_flux(other, velocity, "shared")
        with self.assertRaises(NotFoundError):
            a.execute_flux(plan, 999, "unknown")
        with self.assertRaises(NotFoundError):
            a.execute_flux(999, velocity, "unknown")
        self.assertEqual(len(a.list_results()), 1)


class TestFluxEvidence(FluxCase):
    def execute(self, cells, targets, fields, request="r"):
        a = self.archive
        revision = a.put_source("obs", cells)
        plan = a.save_plan("layout", revision, targets)
        velocity = a.put_velocity("field", revision, fields)
        return a.execute_flux(plan, velocity, request)

    def test_split_signed_flux_and_double_origin(self):
        result = self.execute([flat()], [half("left"), Target(
            "right", (F(1, 2), 0, 0), (1, 1, 1))], field())
        self.assertEqual(result["source_ledger"],
                         [{"source": "s", "external_net": F(0)}])
        self.assertEqual([(r["target"], r["net_outward"])
                          for r in result["targets"]],
                         [("left", F(0)), ("right", F(0))])
        exchanges = [e for e in result["exchanges"]
                     if (e["a"], e["b"]) == ("left", "right")]
        self.assertEqual(sum((e["flux"] for e in exchanges), F(0)), F(1, 8))
        x_zero = [f for f in result["faces"]
                  if all(p[0] == 0 for p in f["vertices"])]
        self.assertEqual(sum((f["flux"] for f in x_zero), F(0)), -F(1, 2))
        for face in x_zero:
            self.assertTrue(face["origins"]["source_faces"])
            self.assertEqual(face["origins"]["target_faces"],
                             [["left", 0, "lower"]])
        # The documented moments for the unrefined x=0 triangle.
        triangle = [f for f in x_zero if len(f["vertices"]) == 3][0]
        self.assertEqual(triangle["moments"]["1"], [-F(1, 2), F(0), F(0)])
        self.assertEqual(triangle["moments"]["y"], [-F(1, 6), F(0), F(0)])
        self.assertEqual(triangle["moments"]["yy"], [-F(1, 12), F(0), F(0)])
        self.assertEqual(triangle["moments"]["yz"], [-F(1, 24), F(0), F(0)])

    def test_mate_reciprocity_and_zero_exchange(self):
        result = self.execute(
            [flat()],
            [half("left"), Target("r0", (F(1, 2), 0, 0), (1, F(1, 4), 1)),
             Target("r1", (F(1, 2), F(1, 4), 0), (1, 1, 1))], field())
        by_id = {f["id"]: f for f in result["faces"]}
        for face in result["faces"]:
            if face["mate"] is not None:
                other = by_id[face["mate"]]
                self.assertEqual(other["mate"], face["id"])
                self.assertEqual(other["flux"], -face["flux"])
                for name, value in face["moments"].items():
                    self.assertEqual(other["moments"][name],
                                     [-x for x in value])
        pair = [e for e in result["exchanges"]
                if (e["a"], e["b"]) == ("left", "r0")]
        self.assertEqual(sum((e["flux"] for e in pair), F(0)), F(3, 32))
        other = [e for e in result["exchanges"]
                 if (e["a"], e["b"]) == ("left", "r1")]
        self.assertEqual(sum((e["flux"] for e in other), F(0)), F(1, 32))
        zero = [e for e in result["exchanges"] if (e["a"], e["b"]) == ("r0", "r1")]
        self.assertTrue(zero)
        self.assertEqual(sum((e["flux"] for e in zero), F(0)), F(0))
        for exchange in result["exchanges"]:
            self.assertLess(exchange["a"], exchange["b"])

    def test_quadratic_flux_and_divergence(self):
        result = self.execute(
            [Cell("s", UNIT, (0, 1, 0, 0))], [whole("all")],
            {"s": [(0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0)]})
        self.assertEqual(result["source_ledger"],
                         [{"source": "s", "external_net": F(1, 12)}])
        self.assertEqual(result["targets"],
                         [{"target": "all", "net_outward": F(1, 12)}])
        self.assertEqual(result["exchanges"], [])
        slant = [f for f in result["faces"]
                 if all(sum(p) == 1 for p in f["vertices"])]
        self.assertEqual(sum((f["flux"] for f in slant), F(0)), F(1, 12))
        self.assertTrue(all(f["origins"]["source_faces"] for f in slant))

    def test_independent_sources_do_not_share(self):
        result = self.execute(
            [flat("s0"), flat("s1")],
            [half("left"), Target("right", (F(1, 2), 0, 0), (1, 1, 1))],
            {"s0": field()["s"], "s1": field()["s"]})
        self.assertEqual(result["source_ledger"],
                         [{"source": "s0", "external_net": F(0)},
                          {"source": "s1", "external_net": F(0)}])
        for source in ("s0", "s1"):
            selected = [e for e in result["exchanges"]
                        if e["source"] == source
                        and (e["a"], e["b"]) == ("left", "right")]
            self.assertEqual(sum((e["flux"] for e in selected), F(0)), F(1, 8))
            x_zero = [f for f in result["faces"] if f["source"] == source
                      and all(p[0] == 0 for p in f["vertices"])]
            self.assertEqual(sum((f["flux"] for f in x_zero), F(0)), -F(1, 2))
            self.assertTrue(all(f["source"] == source for f in x_zero))

    def test_zero_and_low_dimension_targets(self):
        result = self.execute(
            [flat()], [half("west"), Target("far", (2, 2, 2), (3, 3, 3)),
                       Target("touch", (-1, 0, 0), (0, 1, 1))], field())
        self.assertEqual([(r["target"], r["net_outward"])
                          for r in result["targets"]],
                         [("far", F(0)), ("touch", F(0)), ("west", F(0))])
        touches = [f for f in result["faces"] if f["target"] == "touch"]
        self.assertEqual(touches, [])

    def test_isolated_returns_and_reopen(self):
        result = self.execute([flat()], [half()], field())
        first_face = result["faces"][0]
        first_face["vertices"][0][0] = F(999)
        first_face["flux"] = F(999)
        result["exchanges"].clear()
        result["receipt"]["request_id"] = "tampered"
        with Archive(self.path) as reopened:
            again = reopened.read_flux(result["id"])
            self.assertNotEqual(again["faces"][0]["vertices"][0][0], F(999))
            self.assertEqual(again["receipt"]["request_id"], "r")
            self.assertEqual(reopened.read_flux(again["id"]), again)


class TestFluxCompatibility(FluxCase):
    def test_old_database_without_flux_tables(self):
        import sqlite3
        a = self.archive
        revision = a.put_source("obs", [flat()])
        plan = a.save_plan("layout", revision, [whole()])
        mass = a.execute_plan(plan, "old")
        a.close()
        with sqlite3.connect(self.path) as sql:
            for table in ("flux_receipts", "flux_results",
                          "velocity_revisions"):
                sql.execute("DROP TABLE " + table)
        with Archive(self.path) as reopened:
            self.assertEqual(reopened.read_result(mass["id"]), mass)
            velocity = reopened.put_velocity(
                "field", revision, {"s": field()["s"]})
            result = reopened.execute_flux(plan, velocity, "new")
            # Constant rho and constant v over a closed body: no net outflow.
            self.assertEqual(result["source_ledger"][0]["external_net"], F(0))
            self.assertEqual(reopened.read_flux(result["id"]), result)
        self.archive = Archive(self.path)
        self.addCleanup(self.archive.close)


class TestFluxFailures(FluxCase):
    def counts(self):
        with sqlite3.connect(self.path) as sql:
            return tuple(sql.execute("SELECT count(*) FROM " + table).fetchone()[0]
                         for table in ("flux_results", "flux_receipts"))

    def test_compute_failure_leaves_nothing(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        plan = a.save_plan("layout", revision, [half()])
        velocity = a.put_velocity("field", revision, field())
        before = self.counts()
        with sqlite3.connect(self.path) as sql:
            sql.execute("CREATE TRIGGER reject_flux_receipt BEFORE INSERT ON "
                        "flux_receipts BEGIN SELECT RAISE(ABORT,"
                        " 'receipt rejected'); END")
        with self.assertRaises(sqlite3.IntegrityError):
            a.execute_flux(plan, velocity, "publish-failed")
        self.assertEqual(self.counts(), before)
        with sqlite3.connect(self.path) as sql:
            sql.execute("DROP TRIGGER reject_flux_receipt")
        success = a.execute_flux(plan, velocity, "publish-failed")
        self.assertEqual(self.counts(), (before[0] + 1, before[1] + 1))
        self.assertEqual(a.read_flux(success["id"]), success)

    def test_old_results_survive_new_failures(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        plan = a.save_plan("layout", revision, [half()])
        velocity = a.put_velocity("field", revision, field())
        mass = a.execute_plan(plan, "old-mass")
        flux = a.execute_flux(plan, velocity, "old-flux")
        with sqlite3.connect(self.path) as sql:
            sql.execute("CREATE TRIGGER reject_flux BEFORE INSERT ON "
                        "flux_results BEGIN SELECT RAISE(ABORT,'no'); END")
        with self.assertRaises(sqlite3.IntegrityError):
            a.execute_flux(plan, velocity, "new-flux")
        with sqlite3.connect(self.path) as sql:
            sql.execute("DROP TRIGGER reject_flux")
        self.assertEqual(a.read_result(mass["id"]), mass)
        self.assertEqual(a.read_flux(flux["id"]), flux)
        self.assertEqual(a.execute_flux(plan, velocity, "old-flux"), flux)

    def test_new_namespace_does_not_disturb_requests(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        plan = a.save_plan("layout", revision, [half()])
        velocity = a.put_velocity("field", revision, field())
        mass = a.execute_plan(plan, "request")
        with self.assertRaises(RequestConflict):
            other = a.save_plan("other", revision, [whole()])
            a.execute_plan(other, "request")
        self.assertEqual(a.execute_plan(plan, "request"), mass)
        a.execute_flux(plan, velocity, "request")
        self.assertEqual(a.execute_plan(plan, "request"), mass)

    def test_read_flux_validation(self):
        with self.assertRaises(NotFoundError):
            self.archive.read_flux(999)
        with self.assertRaises(ValidationError):
            self.archive.read_flux(True)


if __name__ == "__main__":
    unittest.main()
