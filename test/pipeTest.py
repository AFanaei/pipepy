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
            self.points = json.load(fp)

    def testDefinition(self):
        pipeInlet = random.choice(self.points)
        pipe = Pipe(num_nodes=8, length=1 * U.km, teta=0, diameter=10 * U.inch, inlet=MaterialStream(
            P=pipeInlet[0] * U.bar,
            T=(pipeInlet[1] + 273.15) * U.K,
            m=22.28 * U.kg / U.s,
            MW=16.04 * U.g / U.mol
        ))
        self.assertIsNotNone(pipe)
        self.assertLessEqual(abs(pipe.ps.Z[0] - pipeInlet[2]), 0.001 * pipeInlet[2])

    def testInterpolation(self):
        interpolator = CompFactorInterpolator()
        for p in self.points:
            self.assertLessEqual(abs(interpolator.get_z(p[0], p[1]) - p[2]), 0.001 * p[2])


class PipeIsothermTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.pipe = Pipe(num_nodes=5, length=10 * U.km, teta=0, diameter=0.4 * U.m, epsilon=4.572e-05 * U.m, inlet=MaterialStream(
            P=2000000 * U.pa,
            T=293.15 * U.K,
            m=22.2816 * U.kg / U.s,
            MW=16.0428 * U.g / U.mol
        ))
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.954547, 'ro': 13.790961, 'v': 12.8571}
        self.propertiesOutlet = {'Z': 0.963394, 'ro': 10.930067, 'v': 16.2223}
        self.outlet_p = 1599796.4

    def testDefinition(self):
        self.assertIsNotNone(self.pipe)
        self.assertAlmostEqual(self.pipe.ps.Z[0], self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[0], self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[0], self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])

    def testSolve(self):
        self.pipe.solve_steady_state()

        self.assertAlmostEqual(self.pipe.ps.Z[-1], self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.pipe.domain.P[-1], self.outlet_p, delta=0.001 * self.outlet_p)
        self.assertAlmostEqual(self.pipe.ps.ro[-1], self.propertiesOutlet['ro'], delta=0.001 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[-1], self.propertiesOutlet['v'], delta=0.001 * self.propertiesOutlet['v'])


class PipeTest(unittest.TestCase):
    def setUp(self):
        # use absulute pressure.
        self.pipe = Pipe(num_nodes=5, length=10 * U.km, teta=0, diameter=0.4 * U.m, epsilon=4.572e-05 * U.m, inlet=MaterialStream(
            P=2000000 * U.pa,
            T=293.15 * U.K,
            m=22.2816 * U.kg / U.s,
            MW=16.0428 * U.g / U.mol
        ), isotherm=False, ambient_t=(10 + 273.15) * U.K, heat_transfer_coef=25 * U.W / ((U.m ** 2) * U.K))
        # this test is only for ch4 in case of other components update the test.
        self.propertiesInlet = {'Z': 0.954547, 'ro': 13.790961, 'v': 12.8571}
        self.propertiesOutlet = {'Z': 0.958328, 'ro': 11.502703, 'v': 15.4148}
        self.outlet_p = 1615544.6

    def testDefinition(self):
        self.assertIsNotNone(self.pipe)
        self.assertAlmostEqual(self.pipe.ps.Z[0], self.propertiesInlet['Z'], delta=0.001 * self.propertiesInlet['Z'])
        self.assertAlmostEqual(self.pipe.ps.ro[0], self.propertiesInlet['ro'], delta=0.001 * self.propertiesInlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[0], self.propertiesInlet['v'], delta=0.001 * self.propertiesInlet['v'])

    def testSolve(self):
        self.pipe.solve_steady_state()

        self.assertAlmostEqual(self.pipe.ps.Z[-1], self.propertiesOutlet['Z'], delta=0.001 * self.propertiesOutlet['Z'])
        self.assertAlmostEqual(self.pipe.domain.P[-1], self.outlet_p, delta=0.001 * self.outlet_p)
        self.assertAlmostEqual(self.pipe.ps.ro[-1], self.propertiesOutlet['ro'], delta=0.001 * self.propertiesOutlet['ro'])
        self.assertAlmostEqual(self.pipe.pc.v[-1], self.propertiesOutlet['v'], delta=0.001 * self.propertiesOutlet['v'])
