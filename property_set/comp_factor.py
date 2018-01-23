import os
import json
import numpy as np
from scipy import interpolate


class CompFactorInterpolator:
    address = os.path.join(
        os.path.dirname(__file__), os.pardir, 'databases', "z.json")

    def __init__(self):
        with open(self.address) as fp:
            points = np.array(json.load(fp))
        self.tck = interpolate.bisplrep(points[:, 0], points[:, 1], points[:, 2])

    def get_z(self, p, t, dp=False, dt=False):
        if dp:
            return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck, dx=1)
        if dt:
            return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck, dy=1)
        return interpolate.bisplev(np.array(float(p)), np.array(float(t)), self.tck)
