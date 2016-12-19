import json
import os
import numpy as np

from scipy import interpolate, constants


class InputError(Exception):
    pass


class Zcalculator:
    m_instance = None
    address = os.path.join(
        os.path.dirname(__file__), os.pardir, 'databases', "z.json")

    def __init__(self):
        with open(self.address) as fp:
            points = np.array(json.load(fp))
        self.tck = interpolate.bisplrep(points[:, 0], points[:, 1], points[:, 2])

    @staticmethod
    def instance():
        if Zcalculator.m_instance is None:
            Zcalculator.m_instance = Zcalculator()
        return Zcalculator.m_instance

    def get_z(self, P, T):
        return interpolate.bisplev(float(P), float(T), self.tck)


class Node:
    def __init__(self, pipe=None, state=None):
        self.P = 0
        self.T = 0
        self.m = 0
        self.pipe = pipe
        if state is not None:
            self.P = state['P'].as_two_terms()[0]
            self.T = state['T'].as_two_terms()[0]
            self.m = state['m'].as_two_terms()[0]

    @property
    def Z(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return Zcalculator.instance().get_z((self.P-101325)/10**5, self.T-273.15)

    @property
    def v(self):
        return self.m/(self.ro*self.pipe.A)

    @property
    def ro(self):
        return self.P*self.pipe.M/(constants.R*self.T)

    def equations(self, pre, next):
        dx = 2 * self.pipe.dx
        res = []
        # equation for dp/dt
        res[0] = -self.Z * constants.R * self.T / self.pipe.A * (next.m - pre.m) / dx
        # equation for dm/dt
        res[1] = -(next.m ** 2 * next.Z * constants.R * next.T / next.P / self.pipe.A) / dx \
                 +(pre.m ** 2 * pre.Z * constants.R * pre.T / pre.P / self.pipe.A) / dx \
                - self.pipe.A*(next.P-pre.P)/dx \
                - self.pipe.f_r*self.Z*constants.R*self.T/(2*self.pipe.D*self.pipe.A) * self.m*abs(self.m)/self.P \
                - self.pipe.A*self.P*constants.g*np.sin(self.pipe.teta)/(self.Z*constants.R*self.T)


class Pipe:
    def __init__(self, inlet=None, outlet=None, num_nodes=None, teta=None, length=None, diameter=None, molar_mass=None):
        self.inlet = inlet
        self.outlet = outlet
        self.num_nodes = num_nodes
        self.length = length.as_two_terms()[0]
        self.dx = length / (num_nodes + 1)
        self.D = diameter.as_two_terms()[0]
        self.A = np.pi*self.D**2/4
        self.teta = teta
        self.M = molar_mass.as_two_terms()[0]

        self.nodes = [Node(self) for i in range(0, num_nodes + 2)]
        self.nodes[0] = Node(self, self.inlet)
        self.nodes[-1] = Node(self, self.outlet)

        self.isFeasible()

    @property
    def f_r(self):
        return 0.0004

    def isFeasible(self):
        number_equations = self.num_nodes*3 + 3
        number_vars = self.num_nodes*3
        if self.inlet is not None:
            number_vars += 0 if self.inlet['P'] is None else 1
            number_vars += 0 if self.inlet['T'] is None else 1
            number_vars += 0 if self.inlet['m'] is None else 1
        if self.outlet is not None:
            number_vars += 0 if self.outlet['P'] is None else 1
            number_vars += 0 if self.outlet['T'] is None else 1
            number_vars += 0 if self.outlet['m'] is None else 1

        if number_vars!=number_equations:
            raise InputError("Bad input")
