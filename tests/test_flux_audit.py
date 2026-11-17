import unittest
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory

import verify_flux


class TestIndependentAudit(unittest.TestCase):
    """Runs the independent analytic auditor over every published scenario."""

    def test_all_scenarios_reconcile(self):
        with TemporaryDirectory(prefix="remap-audit-") as directory:
            report = verify_flux.run_scenarios(
                str(Path(directory) / "audit.sqlite"))
        self.assertEqual(
            set(report),
            {"split", "partial", "quadratic", "independent", "skew", "grid",
             "uncovered", "permutation", "reopened", "random"})
        self.assertEqual(report["split"]["external_net"], {"s": "0"})
        self.assertEqual(report["partial"]["exchanges"], 3)
        self.assertEqual(report["quadratic"]["external_net"], {"s": "1/12"})
        self.assertEqual(report["independent"]["external_net"],
                         {"s0": "0", "s1": "0"})
        self.assertEqual(report["uncovered"]["external_net"],
                         {"s0": "0", "s1": "0"})

    def test_analytic_solid_moments(self):
        unit = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        solid = verify_flux.intersection_solid(
            unit, (0, 0, 0), (F(1, 2), 1, 1))
        self.assertEqual(verify_flux.solid_volume(solid), F(7, 48))
        self.assertEqual(verify_flux.solid_first_moments(solid),
                         (F(11, 384), F(5, 128), F(5, 128)))

    def test_simplex_integral_matches_flux(self):
        tri = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        self.assertEqual(verify_flux.triangle_flux(
            tri, (0, 1, 0, 0), ((0, 1, 0, 0), (0, 0, 0, 0), (0, 0, 0, 0))),
            F(1, 12))
        # Constant density, constant velocity on a unit square sheet.
        square = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
        self.assertEqual(verify_flux.triangle_flux(
            square, (2, 0, 0, 0), ((0, 0, 0, 0), (0, 0, 0, 0),
                                   (1, 0, 0, 0))), F(2))
        # The same sheet with the opposite winding flips the sign.
        self.assertEqual(verify_flux.triangle_flux(
            list(reversed(square)), (2, 0, 0, 0),
            ((0, 0, 0, 0), (0, 0, 0, 0), (1, 0, 0, 0))), F(-2))


if __name__ == "__main__":
    unittest.main()
