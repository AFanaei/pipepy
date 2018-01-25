import json
import os
import unittest
import random
from sympy.physics import units as U

from equipment.pipe import Pipe
from property_set.comp_factor import CompFactorInterpolator
from stream.material_stream import MaterialStream


class PipeDefinitionTest(unittest.TestCase):
    def setUp(self):
        self.z = CompFactorInterpolator()
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "z.json")

        with open(address) as fp:
            points = json.load(fp)

        # testing unit conversion in pipe pressure is absulute.
        self.pipeInlet = random.choice(points)
        self.pipe = Pipe(num_nodes=8, length=1 * U.km, teta=0, diameter=10 * U.inch, inlet=MaterialStream(
            P=self.pipeInlet[0] * U.bar,
            T=(self.pipeInlet[1] + 273.15) * U.K,
            m=22.28 * U.kg / U.s,
            MW=16.04 * U.g / U.mol
        ))

    def testDefinition(self):
        self.assertIsNotNone(self.pipe)
        self.assertLessEqual(abs(self.pipe.ps.Z[0] - self.pipeInlet[2]), 0.001 * self.pipeInlet[2])


class PipeIsothermTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.pipe = Pipe(num_nodes=5, length=1 * U.km, teta=0, diameter=0.254 * U.m, epsilon=4.572e-05 * U.m, inlet=MaterialStream(
            P=1761580 * U.pa,
            T=322.737 * U.K,
            m=22.2816 * U.kg / U.s,
            MW=16.0428 * U.g / U.mol
        ))
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.9712604, 'ro': 10.84355, 'v': 40.5526, 'v_w': 459.8984}
        self.propertiesOutlet = {'Z': 0.980507, 'ro': 7.18855, 'v': 61.1714, 'v_w': 461.51}

    def testDefinition(self):
        self.assertIsNotNone(self.pipe)
        self.assertAlmostEqual(self.pipe.ps.Z[0], self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[0], self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[0], self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])
        self.assertAlmostEqual(self.pipe.ps.v_w[0], self.propertiesInlet['v_w'],
                               delta=0.01 * self.propertiesInlet['v_w'])

    def testSolve(self):
        self.pipe.solve_steady_state()

        self.assertAlmostEqual(self.pipe.ps.Z[-1], self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[-1], self.propertiesOutlet['ro'], delta=0.01 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[-1], self.propertiesOutlet['v'], delta=0.01 * self.propertiesOutlet['v'])
        self.assertAlmostEqual(self.pipe.ps.v_w[-1], self.propertiesOutlet['v_w'],
                               delta=0.01 * self.propertiesOutlet['v_w'])


class PipeTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.pipe = Pipe(num_nodes=7, length=1 * U.km, teta=0, diameter=0.254 * U.m, epsilon=4.572e-05 * U.m, inlet=MaterialStream(
            P=1761580 * U.pa,
            T=322.737 * U.K,
            m=22.2816 * U.kg / U.s,
            MW=16.0428 * U.g / U.mol
        ), isotherm=False, ambient_t=(10 + 273.15) * U.K, heat_transfer_coef=25 * U.W / ((U.m ** 2) * U.K))
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.97126, 'ro': 10.8436, 'v': 40.5526}
        self.propertiesOutlet = {'Z': 0.97669, 'ro': 7.68319, 'v': 57.2332}

    def testDefinition(self):
        self.assertIsNotNone(self.pipe)
        self.assertAlmostEqual(self.pipe.ps.Z[0], self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[0], self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[0], self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])

    def testSolve(self):
        self.pipe.solve_steady_state()

        self.assertAlmostEqual(self.pipe.ps.Z[-1], self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[-1], self.propertiesOutlet['ro'], delta=0.01 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[-1], self.propertiesOutlet['v'], delta=0.01 * self.propertiesOutlet['v'])
