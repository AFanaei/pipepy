import json
import os
import random
import unittest

import numpy as np

from property_set.property_set import PropertySet


class ZTest(unittest.TestCase):
    def setUp(self):
        self.ps = PropertySet()
        self.ps.MW = 16.04 / 1000
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "z.json")

        with open(address) as fp:
            points = json.load(fp)

        self.points = [random.choice(points) for _ in range(20)]
        self.interp = []
        for i in range(20):
            index = random.randint(0, len(points) - 1)
            if index + 1 == len(points):
                continue
            self.interp.append([points[index], points[index + 1]])

    def testZ(self):
        for i in self.points:
            self.ps.P = i[0] * 10**5
            self.ps.T = i[1] + 273.15
            self.ps.invalidate_cache()
            self.assertLessEqual(abs(self.ps.Z - i[2]), i[2] * 0.001)

    def testZInterp(self):
        for i in self.interp:
            p = [np.average([i[0][j], i[1][j]]) for j in range(3)]
            self.ps.P = p[0] * 10**5
            self.ps.T = p[1] + 273.15
            self.ps.invalidate_cache()
            # we cant check exactly just check the interpolated value is between the two bounds.
            # check to make shure calculated z is between the two points.
            one = i[0][2] - self.ps.Z
            two = self.ps.Z - i[1][2]
            self.assertGreater(one * two, 0,
                               "z:{},P:{},T:{} not between {k[0][0]}-{k[1][0]},{k[0][1]}-{k[1][1]}"
                               .format(self.ps.Z, p[0], p[1], k=i))
