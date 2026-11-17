import unittest
from fractions import Fraction as F
from remap import Archive, Target
from tests.common import cell, configured, half, whole


class TestRemap(unittest.TestCase):
    def setUp(self):
        self.archive = Archive(":memory:")
        self.addCleanup(self.archive.close)

    def test_full_cover_constant(self):
        a = self.archive
        revision, plan = configured(a, [whole()], (6, 0, 0, 0))
        result = a.execute_plan(plan, "full")
        self.assertEqual(result["source_revision"], revision)
        self.assertEqual(result["transfers"][0]["mass"], F(1))
        self.assertEqual(result["source_ledger"],
                         [{"source": "s", "total_mass": F(1),
                           "covered_mass": F(1), "uncovered_mass": F(0)}])
        self.assertEqual(result["targets"], [{"target": "all", "mass": F(1)}])
        self.assertEqual(result["receipt"]["result_id"], result["id"])
        self.assertEqual(a.read_result(result["id"]), result)

    def test_affine_partial_cover(self):
        a = self.archive
        _, plan = configured(a)
        result = a.execute_plan(plan, "affine")
        row = result["transfers"][0]
        self.assertEqual(row["volume"], F(7, 48))
        self.assertEqual(row["first"], [F(11, 384), F(5, 128), F(5, 128)])
        self.assertEqual(row["mass"], F(11, 384))
        self.assertEqual(result["source_ledger"], [
            {"source": "s", "total_mass": F(1, 24), "covered_mass": F(11, 384),
             "uncovered_mass": F(5, 384)}])
        self.assertEqual(result["targets"][0]["mass"], F(11, 384))

    def test_affine_all_coefficients_and_translation(self):
        a = self.archive
        vertices = ((-3, 2, -5), (-1, 2, -5), (-3, 5, -5), (-3, 2, -1))
        # rho = 20 + x + 2*y + 3*z; values 6,8,12,18 at the vertices.
        revision = a.put_source("translated", [
            cell(density=(20, 1, 2, 3), vertices=vertices)])
        plan = a.save_plan("cut", revision, [Target("t", (-3, 2, -5), (-2, 5, -1))])
        result = a.execute_plan(plan, "translated")
        expected = 24*(6*F(7, 48)+2*F(11, 384)+6*F(5, 128)+12*F(5, 128))
        self.assertEqual(result["transfers"][0]["mass"], expected)
        self.assertEqual(result["source_ledger"][0]["total_mass"], F(44))
        self.assertEqual(result["source_ledger"][0]["uncovered_mass"], 44-expected)

    def test_subdivision_multisource_and_zero_entries(self):
        a = self.archive
        revision = a.put_source("many", [cell("z"), cell("a", (0, 0, 0, 0))])
        p = a.save_plan("one", revision, [half()])
        q = a.save_plan("split", revision, [
            Target("b", ("1/4", 0, 0), ("1/2", 1, 1)),
            Target("a", (0, 0, 0), ("1/4", 1, 1)),
            Target("outside", (2, 2, 2), (3, 3, 3))])
        one, split = a.execute_plan(p, "one"), a.execute_plan(q, "split")
        self.assertEqual(one["source_ledger"], split["source_ledger"])
        self.assertEqual(sum(t["mass"] for t in split["targets"]), F(11, 384))
        self.assertEqual(len(split["transfers"]), 6)
        self.assertEqual([(r["source"], r["target"]) for r in split["transfers"]],
                         [(s, t) for s in ("a", "z") for t in ("a", "b", "outside")])
        for row in split["source_ledger"]:
            self.assertEqual(sum(t["mass"] for t in split["transfers"]
                                 if t["source"] == row["source"]) + row["uncovered_mass"],
                             row["total_mass"])
        self.assertEqual(split["targets"][-1]["mass"], F(0))
