import json
import os
import random
import unittest

import numpy as np

from property_set.comp_factor import CompFactorInterpolator


class ZTest(unittest.TestCase):
    def setUp(self):
        self.z = CompFactorInterpolator()
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
            res = self.z.get_z(i[0], i[1])
            self.assertLessEqual(abs(res - i[2]), i[2] * 0.001, "diffrence in result")

    def testZInterp(self):
        for i in self.interp:
            p = [np.average([i[0][j], i[1][j]]) for j in range(3)]
            res = self.z.get_z(p[0], p[1])
            # we cant check exactly just check the interpolated value is between the two bounds.
            # check to make shure calculated z is between the two points.
            one = i[0][2] - res
            two = res - i[1][2]
            self.assertGreater(one * two, 0,
                               "z:{},P:{},T:{} not between {k[0][0]}-{k[1][0]},{k[0][1]}-{k[1][1]}"
                               .format(res, p[0], p[1], k=i))
