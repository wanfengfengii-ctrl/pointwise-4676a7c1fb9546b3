import unittest
from remap import Archive, PlanError, Target, NotFoundError
from tests.common import cell, whole


class TestPlans(unittest.TestCase):
    def setUp(self):
        self.archive = Archive(":memory:")
        self.addCleanup(self.archive.close)
        self.revision = self.archive.put_source("obs", [cell()])

    def test_overlap_and_touch(self):
        a, revision = self.archive, self.revision
        left = Target("left", (0, 0, 0), (1, 1, 1))
        for targets in (
                [left, Target("overlap", ("1/2", 0, 0), (2, 1, 1))],
                [left, left], [], [Target("flat", (0, 0, 0), (0, 1, 1))]):
            with self.assertRaises(PlanError):
                a.save_plan("bad", revision, targets)
            self.assertEqual(a.list_plans(), [])
        right = Target("right", (1, 0, 0), (2, 1, 1))
        p = a.save_plan("touching", revision, [right, left])
        self.assertEqual([t.id for t in a.read_plan(p)["targets"]], ["left", "right"])
        self.assertEqual(a.execute_plan(p, "touch")["source_ledger"][0]["uncovered_mass"], 0)

    def test_plan_freezes_input_and_source(self):
        a = self.archive
        lower, upper = [0, 0, 0], [1, 1, 1]
        p = a.save_plan("same-name", self.revision, [Target("t", lower, upper)])
        lower[0], upper[0] = 10, 11
        new = a.put_source("obs", [cell(density=(9, 0, 0, 0))])
        q = a.save_plan("same-name", new, [whole()])
        self.assertNotEqual(a.read_plan(p)["target_revision"],
                            a.read_plan(q)["target_revision"])
        self.assertEqual(a.read_plan(p)["source_revision"], self.revision)
        self.assertEqual(a.read_plan(p)["targets"][0].lower, (0, 0, 0))
        before = a.list_plans()
        with self.assertRaises(PlanError):
            a.save_plan("unknown", 999, [whole()])
        self.assertEqual(a.list_plans(), before)
        with self.assertRaises(NotFoundError):
            a.read_plan(999)
