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

    def get_z(self, pres, temp, dp=False, dt=False):
        if hasattr(pres, "__iter__") is False:
            pres = [pres]
            temp = [temp]
        result = np.zeros(len(pres), np.float32)
        for i, (p, t) in enumerate(zip(np.array(pres, dtype=np.float32), np.array(temp, dtype=np.float32))):
            if dp:
                result[i] = interpolate.bisplev(p, t, self.tck, dx=1)
            elif dt:
                result[i] = interpolate.bisplev(p, t, self.tck, dy=1)
            else:
                result[i] = interpolate.bisplev(p, t, self.tck)
        return result
