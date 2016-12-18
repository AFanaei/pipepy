import json
import os
import numpy as np
from scipy import interpolate


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
        return interpolate.bisplev(P, T, self.tck)


class Node:
    def __init__(self, state=None):
        self.P = 0
        self.T = 0
        self.m = 0
        if state is not None:
            self.P = state.P
            self.T = state.T
            self.m = state.m

    def equations(self, pre, next):
        pass


class Pipe:
    def __init__(self, inlet=None, outlet=None, num_nodes=None, teta=None,  ):
        self.inlet = inlet
        self.outlet = outlet
        self.num_nodes = num_nodes
        self.isFeasible()

        self.nodes = [Node() for i in range(0, num_nodes+2)]
        self.nodes[0] = Node(self.inlet)
        self.nodes[-1] = Node(self.outlet)

    def isFeasible(self):
        num_boundry = 0 if self.inlet.P is None else 1
        num_boundry += 0 if self.inlet.T is None else 1
        num_boundry += 0 if self.inlet.m is None else 1

        num_boundry += 0 if self.outlet.P is None else 1
        num_boundry += 0 if self.outlet.T is None else 1
        num_boundry += 0 if self.outlet.m is None else 1

        if num_boundry != 3:
            raise InputError("Bad input")

