from numpy import sqrt
from scipy import constants

from utils.functional import cached_property

from .comp_factor import CompFactorInterpolator


class PropertySet(object):
    def __init__(self, z_calculator=None):
        self.z_calculator = z_calculator if z_calculator is not None else CompFactorInterpolator()
        self.P = None   # unit is pa
        self.T = None   # unit is kelvin
        self.MW = None  # unit is g/mol

    def invalidate_cache(self):
        cached_property.invalidate_cache(self)

    @cached_property
    def Z(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return self.z_calculator.get_z(self.P / 10 ** 5, self.T - 273.15)

    @cached_property
    def dz_dp(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return self.z_calculator.get_z(self.P / 10 ** 5, self.T - 273.15, dp=True) / 10**5

    @cached_property
    def dz_dt(self):
        # in the given correlation p is in bar gauge and T is in centigrade but in our model every thing is in SI
        return self.z_calculator.get_z(self.P / 10 ** 5, self.T - 273.15, dt=True)

    @cached_property
    def c_p(self):
        # TODO: add specific heat from correlation.
        return 2314

    @cached_property
    def ro(self):
        return self.P * self.MW / (self.Z * constants.R * self.T)

    @cached_property
    def v_w(self):
        value = self.Z * constants.R * self.T / self.MW / (1 - self.P / self.Z * self.dz_dp - self.P / (self.ro * self.c_p * self.T) * (1 + self.T / self.Z * self.dz_dt)**2)
        return sqrt(value)

    @cached_property
    def vis(self):
        # TODO: add viscosity from correlation.
        return 1.15 * 10**5
