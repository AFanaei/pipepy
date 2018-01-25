import numpy as np
from sympy.physics import units as U
from sympy.physics.units import convert_to


class MaterialStream:
    def __init__(self, P, T, m, MW):
        self.P = np.array(convert_to(P, U.pa).args[0], dtype=np.float32)
        self.T = np.array(convert_to(T, U.K).args[0], dtype=np.float32)
        self.m = np.array(convert_to(m, U.kg / U.s).args[0], dtype=np.float32)
        self.MW = np.array(convert_to(MW, U.kg / U.mol).args[0], dtype=np.float32)
