import json
import os
import numpy as np

from scipy import interpolate, constants, optimize, math


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

    def get_z(self, p, t, dp=False, dt=False):
        if dp:
            return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck, dx=1)
        if dt:
            return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck, dy=1)
        return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck)


class Node:
    def __init__(self, pipe=None, state=None):
        self.P = 0
        self.T = 0
        self.m = 0
        self.pipe = pipe
        self.f_r_old = 0.001
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
    def dz_dp(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return Zcalculator.instance().get_z(self.P / 10 ** 5, self.T - 273.15, dp=True)/10**5

    @property
    def dz_dt(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return Zcalculator.instance().get_z(self.P / 10 ** 5, self.T - 273.15, dt=True)

    @property
    def v_w(self):
        return math.sqrt(self.Z*constants.R*self.T/self.pipe.M /
                    (1-self.P/self.Z*self.dz_dp-self.P/(self.ro*self.c_p*self.T)*(1+self.T/self.Z*self.dz_dt)**2))

    @property
    def c_p(self):
        # TODO: add specific heat from correlation.
        return 2314

    @property
    def omega(self):
        return self.pipe.U*(math.pi*self.pipe.D)*(self.pipe.ambient_T-self.T)

    @property
    def tav_w(self):
        return self.f_r*self.ro*self.v*abs(self.v)/8

    @property
    def v(self):
        return self.m / (self.ro * self.pipe.A)

    @property
    def ro(self):
        return self.P * self.pipe.M / (self.Z * constants.R * self.T)

    @property
    def f_r(self):
        def func(x):
            t = math.log10(self.pipe.epsilon/(3.7*self.pipe.D)+2.51/(self.Re*math.sqrt(x)))
            return 1/math.sqrt(x)+2*t
        self.f_r_old = optimize.fsolve(func, np.array(self.f_r_old))
        return float(self.f_r_old)

    @property
    def Re(self):
        # TODO: add viscosity from correlation.
        return self.ro*self.v*self.pipe.D/1.15*10**5

    def num_equations(self):
        if self.pipe.isotherm:
            return 2
        return 3

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
            res[0] = -self.Z * constants.R * self.T / self.pipe.A / self.pipe.M * (next_node.m - pre_node.m) / dx
        # equation for dm/dt
        if hasattr(self, 'is_boundry_m'):
            res[1] = self.is_boundry_m - self.m
        else:
            res[1] = - (next_node.m**2*next_node.Z*constants.R*next_node.T/next_node.P/self.pipe.A/self.pipe.M)/dx \
                     + (pre_node.m**2*pre_node.Z*constants.R*pre_node.T/pre_node.P/self.pipe.A/self.pipe.M)/dx \
                     - self.pipe.A*(next_node.P-pre_node.P)/dx \
                     - self.f_r*self.Z*constants.R*self.T/(2*self.pipe.M*self.pipe.D*self.pipe.A)*self.m*abs(self.m)/self.P \
                     - self.pipe.M*self.pipe.A*self.P*constants.g*np.sin(self.pipe.teta)/(self.Z*constants.R*self.T)
        # equation for dT/dt
        if self.pipe.isotherm:
            return res

        if hasattr(self, 'is_boundry_T'):
            res.append(self.is_boundry_T - self.T)
        else:
            res.append(+ self.v*(next_node.T-pre_node.T)/dx
                       + self.v_w**2/self.c_p*(1+self.T/self.Z*self.dz_dt)*(next_node.v-pre_node.v)/dx
                       - self.v_w**2/(self.c_p*self.P)*(1-self.P/self.Z*self.dz_dp)
                       * (self.omega+self.tav_w*math.pi*self.pipe.D*self.v)/self.pipe.A)
        return res


class Pipe:
    def __init__(self, inlet=None, outlet=None, num_nodes=None, teta=None, length=None, diameter=None, molar_mass=None,
                 epsilon=None, ambient_t=None, isotherm=True, heat_transfer_coef=None, ):
        self.inlet = inlet
        self.outlet = outlet
        self.num_nodes = num_nodes
        self.length = length.as_two_terms()[0]
        self.dx = self.length / (num_nodes + 1)
        self.D = diameter.as_two_terms()[0]
        self.A = np.pi * self.D ** 2 / 4
        self.teta = teta
        self.M = molar_mass.as_two_terms()[0]
        self.epsilon = epsilon and epsilon.as_two_terms()[0]
        self.isotherm = isotherm
        self.ambient_T = ambient_t and ambient_t.as_two_terms()[0]
        self.U = heat_transfer_coef and heat_transfer_coef.as_two_terms()[0]

        self.nodes = [Node(self) for i in range(0, num_nodes + 2)]
        self.nodes[0] = Node(self, self.inlet)
        self.nodes[-1] = Node(self, self.outlet)

        self.is_feasible()

    def is_feasible(self):
        if not self.isotherm and (self.ambient_T is None or self.U is None):
            raise InputError("pipe is not isotherm and ambient T or U does not specified.")

        number_equations = (self.num_nodes+1) * self.nodes[0].num_equations()
        number_vars = self.num_nodes * self.nodes[0].num_equations()
        if self.inlet is not None:
            number_vars += 0 if self.inlet['P'] is None else 1
            number_vars += 0 if self.inlet['m'] is None else 1
            number_vars += 1 if self.inlet['T'] is not None and not self.isotherm else 0
        if self.outlet is not None:
            number_vars += 0 if self.outlet['P'] is None else 1
            number_vars += 0 if self.outlet['m'] is None else 1
            number_vars += 1 if self.outlet['T'] is not None and not self.isotherm else 0

        if number_vars != number_equations:
            raise InputError("number of equations {} != variables {}.".format(number_equations,number_vars))

    def solve_steady_state(self):
        p, t, m = self._initialize_by_boundry()
        x = [float(p)] * (self.num_nodes + 2) + [float(m)] * (self.num_nodes + 2)
        if not self.isotherm:
            x += [float(t)] * (self.num_nodes + 2)
        x = optimize.fsolve(self._system_equations, x)

        # for i, (p, m) in enumerate(zip(x[:len(x) // 2], x[len(x) // 2:])):
        #     self.nodes[i].P = p
        #     self.nodes[i].m = m

    def _initialize_by_boundry(self):
        if self.inlet is not None and 'P' in self.inlet:
            p = self.inlet['P'].as_two_terms()[0]
        else:
            p = self.outlet['P'].as_two_terms()[0]

        if self.inlet is not None and 'm' in self.inlet:
            m = self.inlet['m'].as_two_terms()[0]
        else:
            m = self.outlet['m'].as_two_terms()[0]

        if self.inlet is not None and 'T' in self.inlet:
            t = self.inlet['T'].as_two_terms()[0]
        else:
            t = self.outlet['T'].as_two_terms()[0]

        for node in self.nodes:
            node.P = p
            node.m = m
            node.T = t

        return p, t, m

    def _system_equations(self, x):
        y = []
        len_x = len(x)
        d = len_x//len(self.nodes)
        l = [x[i*len_x//d:(i+1)*len_x//d] for i in range(d)]
        for i, tup in enumerate(zip(*l)):
            self.nodes[i].P = tup[0]
            self.nodes[i].m = tup[1]
            try:
                self.nodes[i].T = tup[2]
            except IndexError:
                pass

        for i,_ in enumerate(self.nodes):
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
