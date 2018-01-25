from sympy.physics import units as U
from sympy.physics.units import convert_to


class MaterialStream:
    def __init__(self, P, T, m, MW):
        self.P = convert_to(P, U.pa).args[0]
        self.T = convert_to(T, U.K).args[0]
        self.m = convert_to(m, U.kg / U.s).args[0]
        self.MW = convert_to(MW, U.kg / U.mol).args[0]
