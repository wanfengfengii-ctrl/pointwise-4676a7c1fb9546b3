import unittest
from fractions import Fraction as F
from remap import Archive, Cell, NotFoundError, ValidationError
from tests.common import UNIT, cell


class TestSources(unittest.TestCase):
    def setUp(self):
        self.archive = Archive(":memory:")
        self.addCleanup(self.archive.close)

    def test_invalid_source_atomic(self):
        a = self.archive
        revision = a.put_source("obs", [cell()])
        before = a.list_sources()
        invalid = [
            [cell(), Cell("bad", ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)),
                          (1, 0, 0, 0))],
            [cell(), cell("bad", (-1, 0, 0, 0))],
            [cell(), cell()],
            [],
            [Cell("bad", UNIT, (True, 0, 0, 0))],
            [Cell("bad", UNIT, (0.5, 0, 0, 0))],
            [Cell("bad", UNIT, ("1/0", 0, 0, 0))],
            [Cell("", UNIT, (1, 0, 0, 0))],
        ]
        for cells in invalid:
            with self.assertRaises(ValidationError):
                a.put_source("obs", cells)
            self.assertEqual(a.list_sources(), before)
            self.assertEqual(a.source_head("obs"), revision)

    def test_exact_revisions_and_detached_reads(self):
        a = self.archive
        points = [list(p) for p in UNIT]
        density = ["1/3", "2/3", 0, 0]
        old = a.put_source("obs", [Cell("z", points, density), cell("a")])
        points[0][0] = 900
        density[0] = "99"
        saved = a.read_source(old)
        self.assertEqual([c.id for c in saved["cells"]], ["a", "z"])
        self.assertEqual(saved["cells"][1].density[0], F(1, 3))
        self.assertEqual(saved["cells"][1].vertices[0], (F(0),)*3)
        new = a.put_source("obs", [cell(density=(2, 0, 0, 0))])
        self.assertEqual(a.source_head("obs"), new)
        self.assertEqual(a.read_source(old), saved)
        saved["name"] = "changed"
        self.assertEqual(a.read_source(old)["name"], "obs")
        self.assertEqual([r["is_head"] for r in a.list_sources()], [False, True])
        with self.assertRaises(NotFoundError):
            a.read_source(999)
        with self.assertRaises(NotFoundError):
            a.source_head("absent")
