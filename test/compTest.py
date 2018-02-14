import unittest

import numpy as np

from component.component import Component, ComponentList
from property_set.comp_factor_helper import calculate_comp_factor


class PipeDefinitionTest(unittest.TestCase):
    """
        data for this tests is retrived from: ISO12213-2-2005 natural gas calculation of compressibility factor part 2.
    """
    def setUp(self):
        self.component_list = ComponentList(
            Component.from_db('CO2'),
            Component.from_db('N2'),
            Component.from_db('H2'),
            Component.from_db('CO'),
            Component.from_db('C1'),
            Component.from_db('C2'),
            Component.from_db('C3'),
            Component.from_db('i-C4'),
            Component.from_db('n-C4'),
            Component.from_db('i-C5'),
            Component.from_db('n-C5'),
            Component.from_db('C6'),
            Component.from_db('C7'),
            Component.from_db('C8'),
        )
        self.condition = np.array([
            # first column is Pressure in bar
            # second column is Temperature in degree centigerad
            # other columns are compersibility factor for diffrent gases 
            [ 60, -3.15, 0.84053, 0.83348, 0.79380, 0.88550, 0.82609, 0.85380],  # noqa: E201, E241
            [ 60,  6.85, 0.86199, 0.85596, 0.82206, 0.90144, 0.84969, 0.87370],  # noqa: E201, E241
            [ 60, 16.85, 0.88006, 0.87484, 0.84544, 0.91501, 0.86944, 0.89052],  # noqa: E201, E241
            [ 60, 36.85, 0.90867, 0.90446, 0.88183, 0.93674, 0.90052, 0.91723],  # noqa: E201, E241
            [ 60, 56.85, 0.93011, 0.92696, 0.90868, 0.95318, 0.92368, 0.93730],  # noqa: E201, E241
            [120, -3.15, 0.72133, 0.71044, 0.64145, 0.81024, 0.69540, 0.75074],  # noqa: E201, E241
            [120,  6.85, 0.76025, 0.75066, 0.68971, 0.83782, 0.73780, 0.78586],  # noqa: E201, E241
            [120, 16.85, 0.79317, 0.78475, 0.73123, 0.86137, 0.77369, 0.81569],  # noqa: E201, E241
            [120, 36.85, 0.84515, 0.83863, 0.79698, 0.89913, 0.83022, 0.86311],  # noqa: E201, E241
            [120, 56.85, 0.88383, 0.87870, 0.84553, 0.92766, 0.87211, 0.89862],  # noqa: E201, E241
        ])
        self.compositions = np.array([
            # gas molar composition:
            # rows are:
            # CO2, N2, H2, CO, C1, C2, c3, i-c, n-c, i-c, n-c, c6, c7, c8
            [0.006,  0.005,  0.015,  0.016,  0.076,  0.011],   # noqa: E201, E241
            [0.003,  0.031,  0.010,  0.100,  0.057,  0.117],   # noqa: E201, E241
            [0.00,   0.00,   0.00,   0.095,  0.00,   0.00],    # noqa: E201, E241
            [0.00,   0.00,   0.00,   0.010,  0.00,   0.00],    # noqa: E201, E241
            [0.965,  0.907,  0.859,  0.735,  0.812,  0.826],   # noqa: E201, E241
            [0.018,  0.0450, 0.085,  0.033,  0.043,  0.035],   # noqa: E201, E241
            [0.0045, 0.0084, 0.023,  0.0074, 0.009,  0.0075],  # noqa: E201, E241
            [0.001,  0.0010, 0.0035, 0.0012, 0.0015, 0.0012],  # noqa: E201, E241
            [0.001,  0.0015, 0.0035, 0.0012, 0.0015, 0.0012],  # noqa: E201, E241
            [0.0005, 0.0003, 0.0005, 0.0004, 0.00,   0.0004],  # noqa: E201, E241
            [0.0003, 0.0004, 0.0005, 0.0004, 0.00,   0.0004],  # noqa: E201, E241
            [0.0007, 0.0004, 0.00,   0.0002, 0.00,   0.0002],  # noqa: E201, E241
            [0.00,   0.00,   0.00,   0.0001, 0.00,   0.0001],  # noqa: E201, E241
            [0.00,   0.00,   0.00,   0.0001, 0.00,   0.00],    # noqa: E201, E241
        ])

    def testCompositions(self):
        self.assertEqual(len(self.component_list.MW), 14)

    def testGas(self):
        _, compositions_num = self.compositions.shape
        for i in range(compositions_num):
            composition = self.compositions[:, i]
            z_calculated = calculate_comp_factor(self.condition[:, 0], self.condition[:, 1], composition, self.component_list)
            z_iso = self.condition[:, 2 + i]
            result = abs(z_calculated - z_iso) / z_iso
            for r in result:
                self.assertLessEqual(r, 0.001)
