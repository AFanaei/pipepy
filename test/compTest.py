import json
import os
import unittest
import random

from property_set.comp_factor import CompFactorInterpolator
from property_set.comp_factor_helper import calculate_comp_factor


class PipeDefinitionTest(unittest.TestCase):
    def setUp(self):
        self.z = CompFactorInterpolator()
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "z.json")

        with open(address) as fp:
            points = json.load(fp)

        # testing unit conversion in pipe pressure is absulute.
        self.point = random.choice(points)

    def testDefinition(self):
        result = calculate_comp_factor(self.point, self.point, [0.5, 0.5])

        self.assertLessEqual(abs(result[0] - self.point[2]), 0.001 * self.point[2])
