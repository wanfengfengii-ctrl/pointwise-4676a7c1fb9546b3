import sqlite3
import unittest
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
from remap import Archive, NotFoundError, RequestConflict
from tests.common import cell, configured, half


class TestArchive(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory(prefix="remap-archive-")
        self.addCleanup(self.directory.cleanup)
        self.path = str(Path(self.directory.name) / "archive.sqlite")

    def test_frozen_revision_reopen(self):
        with Archive(self.path) as a:
            old_revision, old_plan = configured(a)
            old = a.execute_plan(old_plan, "old")
            new_revision = a.put_source("observations", [cell(density=(1, 0, 0, 0))])
            self.assertNotEqual(old_revision, new_revision)
            self.assertEqual(a.source_head("observations"), new_revision)
            new_plan = a.save_plan("layout", new_revision, [half()])
            new = a.execute_plan(new_plan, "new")
            self.assertEqual(new["transfers"][0]["mass"], F(7, 48))
            self.assertEqual(a.execute_plan(old_plan, "old-again")["transfers"],
                             old["transfers"])
        with Archive(self.path) as a:
            self.assertEqual(a.read_result(old["id"]), old)
            self.assertEqual(a.execute_plan(old_plan, "old"), old)
            self.assertEqual(a.read_plan(old_plan)["source_revision"], old_revision)
            self.assertEqual(a.read_result(new["id"]), new)
            old["transfers"][0]["first"][0] = F(999)
            old["receipt"]["request_id"] = "corrupt-return"
            self.assertEqual(a.read_result(old["id"])["transfers"][0]["mass"],
                             F(11, 384))
            self.assertEqual(a.read_result(old["id"])["receipt"]["request_id"], "old")

    def test_replay_conflict_and_not_found(self):
        with Archive(self.path) as a:
            revision, p = configured(a)
            q = a.save_plan("another", revision, [half()])
            first = a.execute_plan(p, "request")
            self.assertEqual(a.execute_plan(p, "request"), first)
            with self.assertRaises(RequestConflict):
                a.execute_plan(q, "request")
            with self.assertRaises(NotFoundError):
                a.execute_plan(999, "unknown")
            with self.assertRaises(NotFoundError):
                a.read_result(999)
            self.assertEqual(len(a.list_results()), 1)

    def test_catalog_write_rollback(self):
        with Archive(self.path) as a:
            revision, _ = configured(a)
            before_sources, before_plans = a.list_sources(), a.list_plans()
            with sqlite3.connect(self.path) as sql:
                sql.execute("CREATE TRIGGER reject_head BEFORE UPDATE ON source_heads "
                            "BEGIN SELECT RAISE(ABORT, 'head rejected'); END")
            with self.assertRaises(sqlite3.IntegrityError):
                a.put_source("observations", [cell(density=(2, 0, 0, 0))])
            self.assertEqual(a.list_sources(), before_sources)
            self.assertEqual(a.source_head("observations"), revision)
            with sqlite3.connect(self.path) as sql:
                sql.execute("DROP TRIGGER reject_head")
                sql.execute("CREATE TRIGGER reject_plan BEFORE INSERT ON plans "
                            "BEGIN SELECT RAISE(ABORT, 'plan rejected'); END")
                count = sql.execute("SELECT count(*) FROM target_revisions").fetchone()[0]
            with self.assertRaises(sqlite3.IntegrityError):
                a.save_plan("cannot-save", revision, [half()])
            self.assertEqual(a.list_plans(), before_plans)
            with sqlite3.connect(self.path) as sql:
                self.assertEqual(sql.execute(
                    "SELECT count(*) FROM target_revisions").fetchone()[0], count)
                sql.execute("DROP TRIGGER reject_plan")
            a.save_plan("retry", revision, [half()])
