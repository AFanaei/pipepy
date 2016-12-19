import json
import os
import numpy as np

from scipy import interpolate, constants, optimize


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

    def get_z(self, p, t):
        return interpolate.bisplev(float(p), float(t), self.tck)


class Node:
    def __init__(self, pipe=None, state=None):
        self.P = 0
        self.T = 0
        self.m = 0
        self.pipe = pipe
        if state is not None:
            self.P = state['P'].as_two_terms()[0] if 'P' in state else 0
            self.T = state['T'].as_two_terms()[0] if 'T' in state else 0
            self.m = state['m'].as_two_terms()[0] if 'm' in state else 0

            self.is_boundry_p = state['P'].as_two_terms()[0] if 'P' in state else False
            self.is_boundry_T = state['T'].as_two_terms()[0] if 'T' in state else False
            self.is_boundry_m = state['m'].as_two_terms()[0] if 'm' in state else False

    @property
    def Z(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return Zcalculator.instance().get_z(self.P / 10 ** 5, self.T - 273.15)

    @property
    def v(self):
        return self.m / (self.ro * self.pipe.A)

    @property
    def ro(self):
        return self.P * self.pipe.M / (self.Z * constants.R * self.T)

    def equations(self, pre_node, next_node):
        dx = 2 * self.pipe.dx
        res = [0, 0]
        if pre_node is None:
            dx = self.pipe.dx
            pre_node = self

        if next_node is None:
            dx = self.pipe.dx
            next_node = self

        # equation for dp/dt
        if hasattr(self, 'is_boundry_p'):
            res[0] = self.is_boundry_p - self.P
        else:
            res[0] = -self.Z * constants.R * self.T / self.pipe.A * (next_node.m - pre_node.m) / dx
        # equation for dm/dt
        if hasattr(self, 'is_boundry_m'):
            res[1] = self.is_boundry_m - self.m
        else:
            res[1] = -(next_node.m ** 2 * next_node.Z * constants.R * next_node.T / next_node.P / self.pipe.A) / dx \
                     + (pre_node.m ** 2 * pre_node.Z * constants.R * pre_node.T / pre_node.P / self.pipe.A) / dx \
                     - self.pipe.A * (next_node.P - pre_node.P) / dx \
                     - self.pipe.f_r * self.Z * constants.R * self.T / (2 * self.pipe.D * self.pipe.A) * self.m * abs(self.m) / self.P \
                     - self.pipe.A * self.P * constants.g * np.sin(self.pipe.teta) / (self.Z * constants.R * self.T)
        return res


class Pipe:
    def __init__(self, inlet=None, outlet=None, num_nodes=None, teta=None, length=None, diameter=None, molar_mass=None):
        self.inlet = inlet
        self.outlet = outlet
        self.num_nodes = num_nodes
        self.length = length.as_two_terms()[0]
        self.dx = self.length / (num_nodes + 1)
        self.D = diameter.as_two_terms()[0]
        self.A = np.pi * self.D ** 2 / 4
        self.teta = teta
        self.M = molar_mass.as_two_terms()[0]

        self.nodes = [Node(self) for i in range(0, num_nodes + 2)]
        self.nodes[0] = Node(self, self.inlet)
        self.nodes[-1] = Node(self, self.outlet)

        self.is_feasible()

    @property
    def f_r(self):
        return 0.855

    def is_feasible(self):
        number_equations = self.num_nodes * 3 + 3
        number_vars = self.num_nodes * 3
        if self.inlet is not None:
            number_vars += 0 if self.inlet['P'] is None else 1
            number_vars += 0 if self.inlet['T'] is None else 1
            number_vars += 0 if self.inlet['m'] is None else 1
        if self.outlet is not None:
            number_vars += 0 if self.outlet['P'] is None else 1
            number_vars += 0 if self.outlet['T'] is None else 1
            number_vars += 0 if self.outlet['m'] is None else 1

        if number_vars != number_equations:
            raise InputError("Bad input")

    def solve_steady_state(self):
        if self.inlet is not None and 'P' in self.inlet:
            p = self.inlet['P'].as_two_terms()[0]
        else:
            p = self.outlet['P'].as_two_terms()[0]

        if self.inlet is not None and 'm' in self.inlet:
            m = self.inlet['m'].as_two_terms()[0]
        else:
            m = self.outlet['m'].as_two_terms()[0]

        x = [float(p)] * (self.num_nodes + 2) + [float(m)] * (self.num_nodes + 2)
        x = optimize.fsolve(self._system_equations, x)

        for i, (p, m) in enumerate(zip(x[:len(x) // 2], x[len(x) // 2:])):
            self.nodes[i].P = p
            self.nodes[i].m = m
            self.nodes[i].T = 322.7366

    def _system_equations(self, x):
        y = []
        for i, (p, m) in enumerate(zip(x[:len(x) // 2], x[len(x) // 2:])):
            self.nodes[i].P = p
            self.nodes[i].m = m
            self.nodes[i].T = 322.7366

        for i in range(len(x) // 2):
            try:
                pre_node = self.nodes[i - 1]
            except IndexError:
                pre_node = None
            try:
                next_node = self.nodes[i + 1]
            except IndexError:
                next_node = None
            y.extend(self.nodes[i].equations(pre_node, next_node))
        return y
