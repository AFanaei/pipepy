import json
import os
import unittest
import numpy as np
import random

from sympy.physics.units import *

from pipeline.pipe import Pipe, Zcalculator


class ZTest(unittest.TestCase):
    def setUp(self):
        self.z = Zcalculator.instance()
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "z.json")

        with open(address) as fp:
            points = json.load(fp)

        self.points = [random.choice(points) for i in range(20)]
        self.interp = []
        for i in range(20):
            index = random.randint(0, len(points) - 1)
            if index + 1 == len(points):
                continue
            self.interp.append([points[index], points[index + 1]])

        # testing unit conversion in pipe pressure is absulute.
        self.pipeInlet = random.choice(points)
        self.p = Pipe(num_nodes=8, length=1 * km, teta=0, diameter=10 * inch, molar_mass=16.04 * g / mol,
                      inlet={'P': self.pipeInlet[0] * bar,
                             'T': (self.pipeInlet[1] + 273.15) * K,
                             'm': 22.28 * kg / s})

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
                    "z:{},P:{},T:{} not between {k[0][0]}-{k[1][0]},{k[0][1]}-{k[1][1]}".format(res, p[0], p[1], k=i))

    def testDefinition(self):
        self.assertIsNotNone(self.p)
        self.assertLessEqual(abs(self.p.nodes[0].Z - self.pipeInlet[2]), 0.001 * self.pipeInlet[2])


class PipeTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.p = Pipe(num_nodes=8, length=1 * km, teta=0, diameter=0.254 * m, molar_mass=16.0428 * g / mol,
                      inlet={'P': 1761580 * pa, 'T': 322.7366 * K, 'm': 22.28 * kg / s})
        # this test is only for ch4 in case of other components update the test.
        self.properties = {'Z': 0.9712604, 'ro': 10.84355, 'v': 40.5526}

    def testDefinition(self):
        self.assertIsNotNone(self.p)
        self.assertAlmostEqual(self.p.nodes[0].Z, self.properties['Z'], delta=0.001*self.properties['Z'])
        self.assertAlmostEqual(self.p.nodes[0].ro, self.properties['ro'], delta=0.001*self.properties['ro'])
        self.assertAlmostEqual(self.p.nodes[0].v, self.properties['v'], delta=0.001 * self.properties['v'])


if __name__ == '__main__':
    unittest.main()
