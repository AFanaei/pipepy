import json
import os
import unittest
import random
from sympy.physics import units as U

from equipment.pipe import Pipe, Zcalculator


class PipeDefinitionTest(unittest.TestCase):
    def setUp(self):
        self.z = Zcalculator.instance()
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "z.json")

        with open(address) as fp:
            points = json.load(fp)

        # testing unit conversion in pipe pressure is absulute.
        self.pipeInlet = random.choice(points)
        self.p = Pipe(num_nodes=8, length=1 * U.km, teta=0, diameter=10 * U.inch, molar_mass=16.04 * U.g / U.mol,
                      inlet={'P': self.pipeInlet[0] * U.bar,
                             'T': (self.pipeInlet[1] + 273.15) * U.K,
                             'm': 22.28 * U.kg / U.s})

    def testDefinition(self):
        self.assertIsNotNone(self.p)
        self.assertLessEqual(abs(self.p.nodes[0].Z - self.pipeInlet[2]), 0.001 * self.pipeInlet[2])


class PipeIsothermTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.p = Pipe(num_nodes=5, length=1 * U.km, teta=0, diameter=0.254 * U.m, molar_mass=16.0428 * U.g / U.mol,
                      inlet={'P': 1761580 * U.pa, 'T': 322.737 * U.K, 'm': 22.2816 * U.kg / U.s}, epsilon=4.572e-05 * U.m)
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.9712604, 'ro': 10.84355, 'v': 40.5526, 'v_w': 459.8984}
        self.propertiesOutlet = {'Z': 0.980507, 'ro': 7.18855, 'v': 61.1714, 'v_w': 461.51}

    def testDefinition(self):
        self.assertIsNotNone(self.p)
        self.assertAlmostEqual(self.p.nodes[0].Z, self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.p.nodes[0].ro, self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.p.nodes[0].v, self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])
        self.assertAlmostEqual(self.p.nodes[0].v_w, self.propertiesInlet['v_w'],
                               delta=0.01 * self.propertiesInlet['v_w'])

    def testSolve(self):
        self.assertIsNotNone(self.p)

        self.p.solve_steady_state()

        self.assertAlmostEqual(self.p.nodes[-1].Z, self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.p.nodes[-1].ro, self.propertiesOutlet['ro'], delta=0.01 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.p.nodes[-1].v, self.propertiesOutlet['v'], delta=0.01 * self.propertiesOutlet['v'])
        self.assertAlmostEqual(self.p.nodes[-1].v_w, self.propertiesOutlet['v_w'],
                               delta=0.01 * self.propertiesOutlet['v_w'])


class PipeTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.p = Pipe(num_nodes=7, length=1 * U.km, teta=0, diameter=0.254 * U.m, molar_mass=16.0428 * U.g / U.mol,
                      inlet={'P': 1761580 * U.pa, 'T': 322.737 * U.K, 'm': 22.2816 * U.kg / U.s}, epsilon=4.572e-05 * U.m,
                      isotherm=False, ambient_t=(10 + 273.15) * U.K, heat_transfer_coef=25 * U.W / ((U.m ** 2) * U.K))
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.97126, 'ro': 10.8436, 'v': 40.5526}
        self.propertiesOutlet = {'Z': 0.97669, 'ro': 7.68319, 'v': 57.2332}

    def testDefinition(self):
        self.assertIsNotNone(self.p)
        self.assertAlmostEqual(self.p.nodes[0].Z, self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.p.nodes[0].ro, self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.p.nodes[0].v, self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])

    def testSolve(self):
        self.assertIsNotNone(self.p)

        self.p.solve_steady_state()

        self.assertAlmostEqual(self.p.nodes[-1].Z, self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.p.nodes[-1].ro, self.propertiesOutlet['ro'], delta=0.01 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.p.nodes[-1].v, self.propertiesOutlet['v'], delta=0.01 * self.propertiesOutlet['v'])
