"""Regression tests for stable, public MultiFischer behavior."""

import math
import unittest

from multifischer.pss import AbsorbanceData, multifischer, qyratio_array


class ScientificCoreTests(unittest.TestCase):
    def test_qy_ratio_grid_includes_endpoints(self):
        values = qyratio_array((0.5, 2.0))

        self.assertEqual(len(values), 200)
        self.assertEqual(values[0], 0.5)
        self.assertEqual(values[-1], 2.0)

    def test_standard_fischer_expression_at_x_one(self):
        absorbance = AbsorbanceData(
            dark_irr1=0.4,
            dark_irr2=0.2,
            dark_max=0.6,
            irr1_max=0.3,
            irr2_max=0.45,
            irr1_at_irr1=0.2,
            irr2_at_irr2=0.18,
        )

        d1 = absorbance.irr1_at_irr1 - absorbance.dark_irr1
        d2 = absorbance.irr2_at_irr2 - absorbance.dark_irr2
        n_value = (
            (absorbance.irr2_max - absorbance.dark_max)
            / (absorbance.irr1_max - absorbance.dark_max)
        )
        expected = (
            d2 / absorbance.dark_irr2 - d1 / absorbance.dark_irr1
        ) / (
            1.0
            + d2 / absorbance.dark_irr2
            - n_value * (1.0 + d1 / absorbance.dark_irr1)
        )

        self.assertTrue(
            math.isclose(
                multifischer(absorbance, X=1.0, target="irr1"),
                expected,
                rel_tol=1e-12,
                abs_tol=1e-12,
            )
        )


if __name__ == "__main__":
    unittest.main()
