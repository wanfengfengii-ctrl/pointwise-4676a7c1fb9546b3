import sqlite3
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from remap import Archive, Target
from tests.common import configured


class TestFailures(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory(prefix="remap-failures-")
        self.addCleanup(self.directory.cleanup)
        self.path = str(Path(self.directory.name) / "archive.sqlite")
        self.archive = Archive(self.path)
        self.addCleanup(self.archive.close)
        _, self.plan = configured(self.archive, [
            Target("left", (0, 0, 0), ("1/2", 1, 1)),
            Target("right", ("1/2", 0, 0), (1, 1, 1))])
        self.old = self.archive.execute_plan(self.plan, "old")

    def counts(self):
        with sqlite3.connect(self.path) as sql:
            return tuple(sql.execute("SELECT count(*) FROM " + table).fetchone()[0]
                         for table in ("results", "receipts"))

    def test_compute_abort(self):
        observed = []
        problem = RuntimeError("consumer stopped after second overlap")

        def observe(event):
            observed.append(event)
            if len(observed) == 2:
                raise problem

        before = self.counts()
        with self.assertRaises(RuntimeError) as caught:
            self.archive.execute_plan(self.plan, "compute-failed", observer=observe)
        self.assertIs(caught.exception, problem)
        self.assertEqual([e["target"] for e in observed], ["left", "right"])
        self.assertTrue(all(e["moments"].volume > 0 for e in observed))
        self.assertEqual(self.counts(), before)
        self.assertEqual(self.archive.read_result(self.old["id"]), self.old)
        success = self.archive.execute_plan(self.plan, "compute-failed")
        self.assertEqual(self.counts(), (before[0]+1, before[1]+1))
        self.assertEqual(self.archive.execute_plan(self.plan, "compute-failed"), success)

    def test_publish_rollback_retry(self):
        before = self.counts()
        # Genuine SQLite failure after results INSERT, at receipts INSERT.
        with sqlite3.connect(self.path) as sql:
            sql.execute("CREATE TRIGGER reject_receipt BEFORE INSERT ON receipts "
                        "BEGIN SELECT RAISE(ABORT, 'receipt rejected'); END")
        with self.assertRaises(sqlite3.IntegrityError):
            self.archive.execute_plan(self.plan, "publish-failed")
        self.assertEqual(self.counts(), before)
        with Archive(self.path) as reopened:
            self.assertEqual(reopened.read_result(self.old["id"]), self.old)
        with sqlite3.connect(self.path) as sql:
            sql.execute("DROP TRIGGER reject_receipt")
        success = self.archive.execute_plan(self.plan, "publish-failed")
        self.assertEqual(self.counts(), (before[0]+1, before[1]+1))
        events = []
        self.assertEqual(self.archive.execute_plan(
            self.plan, "publish-failed", observer=events.append), success)
        self.assertEqual(events, [])
        with Archive(self.path) as reopened:
            self.assertEqual(reopened.execute_plan(self.plan, "publish-failed"), success)
