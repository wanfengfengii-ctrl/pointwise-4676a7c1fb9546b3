import unittest
from fractions import Fraction as F
from itertools import permutations
from remap import ValidationError
from remap.geometry import intersection_moments, tetrahedron_moments
from remap.models import Moments, ZERO
from tests.common import UNIT


def analytic_x(t):
    # Integrate the triangular section 0<=y,z, y+z<=1-x.
    return Moments((1-(1-t)**3)/6,
                   (t*t/4-t**3/3+t**4/8,
                    (1-(1-t)**4)/24, (1-(1-t)**4)/24))


def add(parts):
    return Moments(sum((p.volume for p in parts), F(0)),
                   tuple(sum((p.first[i] for p in parts), F(0)) for i in range(3)))


class TestGeometry(unittest.TestCase):
    def test_analytic_moments(self):
        self.assertEqual(tetrahedron_moments(UNIT),
                         Moments(F(1, 6), (F(1, 24),)*3))
        for axis in range(3):
            for t in (F(1, 4), F(1, 2), F(2, 3), F(1)):
                upper = [F(1)]*3
                upper[axis] = t
                expected = analytic_x(t)
                first = [expected.first[1]]*3
                first[axis] = expected.first[0]
                self.assertEqual(intersection_moments(UNIT, (0, 0, 0), upper),
                                 Moments(expected.volume, tuple(first)))
        self.assertEqual(intersection_moments(
            UNIT, (-2, -2, -2), (-1, -1, -1)), ZERO)

    def test_degeneracy_and_subdivision(self):
        # Closed face, edge, point contacts contribute no three-dimensional measure.
        contacts = [((-1, -1, -1), (0, 2, 2)),
                    ((-1, -1, -1), (0, 0, 2)),
                    ((-1, -1, -1), (0, 0, 0)),
                    ((1, 0, 0), (2, 1, 1)),
                    ((F(1, 2), F(1, 2), 0), (1, 1, 1))]
        expected = analytic_x(F(1, 2))
        for vertices in permutations(UNIT):
            for lower, upper in contacts:
                self.assertEqual(intersection_moments(vertices, lower, upper), ZERO)
            self.assertEqual(intersection_moments(
                vertices, (0, 0, 0), (1, 1, 1)), tetrahedron_moments(UNIT))
            self.assertEqual(intersection_moments(
                vertices, (0, 0, 0), (F(1, 2), 1, 1)), expected)
        pieces = []
        for x in range(2):
            for y in range(2):
                for z in range(2):
                    pieces.append(intersection_moments(
                        UNIT, (F(x, 2), F(y, 2), F(z, 2)),
                        (F(x+1, 2), F(y+1, 2), F(z+1, 2))))
        self.assertEqual(add(pieces), tetrahedron_moments(UNIT))
        # A box plane contains a non-axis-aligned source edge.
        skew = ((0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1))
        full = tetrahedron_moments(skew)
        self.assertEqual(intersection_moments(skew, (0, 0, 0), (1, 1, 1)), full)
        self.assertEqual(add([
            intersection_moments(skew, (0, 0, 0), (F(1, 2), 1, 1)),
            intersection_moments(skew, (F(1, 2), 0, 0), (1, 1, 1))]), full)

    def test_translation_and_global_coordinates(self):
        offset, scale = (-3, 2, -5), (2, 3, 4)
        vertices = tuple(tuple(F(offset[i])+scale[i]*p[i] for i in range(3))
                         for p in UNIT)
        base = analytic_x(F(1, 2))
        expected = Moments(24*base.volume,
                           tuple(24*(offset[i]*base.volume+scale[i]*base.first[i])
                                 for i in range(3)))
        self.assertEqual(intersection_moments(vertices, offset, (-2, 5, -1)),
                         expected)
        # The small cube is wholly inside a translated, sloping tetrahedron.
        lo, hi = (F(-3), F(2), F(-5)), (F(-11, 4), F(9, 4), F(-19, 4))
        volume = F(1, 64)
        self.assertEqual(intersection_moments(vertices, lo, hi),
                         Moments(volume, tuple(volume*(a+b)/2 for a, b in zip(lo, hi))))

    def test_oblique_homothetic_cap(self):
        vertices = ((0, 0, 0), (2, 1, 0), (2, 0, 1), (2, 1, 1))
        # x<=1 scales this entire tetrahedron by 1/2 about the origin.
        expected = Moments(F(1, 24), (F(1, 32), F(1, 96), F(1, 96)))
        for ordering in permutations(vertices):
            self.assertEqual(intersection_moments(
                ordering, (-1, -1, -1), (1, 2, 2)), expected)

    def test_reject_invalid_geometry(self):
        with self.assertRaises(ValidationError):
            tetrahedron_moments(((0, 0, 0),)*4)
        with self.assertRaises(ValidationError):
            intersection_moments(UNIT, (0, 0, 0), (0, 1, 1))
        with self.assertRaises(ValidationError):
            intersection_moments(UNIT, (0.0, 0, 0), (1, 1, 1))
