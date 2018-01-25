import numpy as np

from scipy import constants, optimize, math
from sympy.physics import units as U
from sympy.physics.units import convert_to

from pde.domain import Domain, DomainVar
from property_set.property_set import PropertySet
from utils.functional import cached_property


class InputError(Exception):
    pass


class PipePropertyCalculator(object):
    def __init__(self, ps, pipe):
        self.ps = ps
        self.pipe = pipe

    def invalidate_cache(self):
        cached_property.invalidate_cache(self)    

    @cached_property
    def omega(self):
        return self.pipe.U * (math.pi * self.pipe.D) * (self.pipe.ambient_T - self.pipe.domain.T)

    @cached_property
    def tav_w(self):
        return self.f_r * self.ps.ro * self.v * abs(self.v) / 8

    @cached_property
    def v(self):
        return self.pipe.domain.m / (self.ps.ro * self.pipe.A)

    @cached_property
    def f_r(self):
        return 0.01356143
        # def func(x):
        #     t = math.log10(self.pipe.epsilon / (3.7 * self.pipe.D) + 2.51 / (self.Re * math.sqrt(x)))
        #     return 1 / math.sqrt(x) + 2 * t
        # self.f_r_old = optimize.fsolve(func, np.array(self.f_r_old))
        # return float(self.f_r_old)

    @cached_property
    def Re(self):
        # TODO: add viscosity from correlation.
        return self.ps.ro * self.v * self.pipe.D / self.ps.vis


class Pipe(object):
    def __init__(self, inlet=None, outlet=None, num_nodes=None, teta=0, length=None, diameter=None,
                 epsilon=0, ambient_t=None, isotherm=True, heat_transfer_coef=0):
        self.inlet_stream = inlet
        self.outlet_stream = outlet
        self.D = np.array(convert_to(diameter, U.m).args[0], dtype=np.float32)
        self.A = np.pi * self.D ** 2 / 4
        self.teta = teta
        self.epsilon = epsilon and np.array(convert_to(epsilon, U.m).args[0], dtype=np.float32)
        self.isotherm = isotherm
        self.ambient_T = ambient_t and np.array(convert_to(ambient_t, U.K).args[0], dtype=np.float32)
        self.U = heat_transfer_coef and np.array(convert_to(heat_transfer_coef, U.W / ((U.m ** 2) * U.K)).args[0], dtype=np.float32)

        self.domain = Domain(np.array(convert_to(length, U.m).args[0], dtype=np.float32), num_nodes + 2)
        # define domain var and their boundry conditions
        DomainVar('P', self.domain, (0, self.inlet_stream.P))
        DomainVar('m', self.domain, (0, self.inlet_stream.m))
        DomainVar('T', self.domain, (0, self.inlet_stream.T))

        # intial conditions
        self.domain.P[:] = self.inlet_stream.P
        self.domain.m[:] = self.inlet_stream.m
        self.domain.T[:] = self.inlet_stream.T

        self.ps = PropertySet()
        self.ps.P = self.domain.P
        self.ps.T = self.domain.T
        self.ps.MW = inlet.MW
        self.pc = PipePropertyCalculator(self.ps, self)

    def solve_steady_state(self):
        return optimize.fsolve(self._system_equations, self.domain.initialize_vars())

    def invalidate_cache(self):
        self.domain.invalidate_cache()
        self.ps.invalidate_cache()
        self.pc.invalidate_cache()

    def _system_equations(self, x):
        self.domain.import_array_to_vars(x)
        self.ps.P = self.domain.P
        self.ps.T = self.domain.T
        self.invalidate_cache()

        expr = (self.domain.m**2 * self.ps.Z * self.domain.T / self.domain.P / self.ps.MW) * (constants.R / self.A)
        self.domain.dPdt = -(self.ps.Z * self.domain.T / self.ps.MW) * (constants.R / self.A) * self.domain.dmdx
        self.domain.dmdt = \
            - np.gradient(expr) / self.domain.dx \
            - self.A * self.domain.dPdx \
            - self.pc.f_r * self.ps.Z * constants.R * self.domain.T / (2 * self.ps.MW * self.D * self.A) * self.domain.m * abs(self.domain.m) / self.domain.P \
            - self.ps.MW * self.A * self.domain.P * constants.g * np.sin(self.teta) / (self.ps.Z * constants.R * self.domain.T)

        # equation for dT/dt
        if self.isotherm:
            return self.domain.export_vars_to_array()

        self.domain.dTdt = \
            + self.pc.v * self.domain.dTdx \
            + self.ps.v_w**2 / self.ps.c_p * (1 + self.domain.T / self.ps.Z * self.ps.dZdT) * np.gradient(self.pc.v) / self.domain.dx \
            - self.ps.v_w**2 / (self.ps.c_p * self.domain.P) * (1 - self.domain.P / self.ps.Z * self.ps.dZdP) \
            * (self.pc.omega + self.pc.tav_w * math.pi * self.D * self.pc.v) / self.A

        return self.domain.export_vars_to_array()