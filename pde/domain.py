import re
from collections import OrderedDict

import numpy as np

from utils.functional import cached_property


class DomainVar(object):
    def __init__(self, name, domain, boundry):
        self.domain = domain
        self.value = np.zeros(domain.num_nodes, dtype=np.float32)
        self._ddt = np.zeros(domain.num_nodes, dtype=np.float32)
        self.boundry = boundry
        domain.add_var_to_domain(name, self)

    @property
    def ddt(self):
        k, v = self.boundry
        self._ddt[k] = self.value[k] - v
        return self._ddt

    @ddt.setter
    def ddt(self, value):
        self._ddt = value

    @cached_property
    def ddx(self):
        return np.gradient(self.value) / self.domain.dx


class Domain(object):
    num_nodes = 0
    length = 0
    dx = 0
    domain_vars = OrderedDict()

    def __init__(self, length, num_nodes):
        self.num_nodes = num_nodes
        self.length = length
        self.dx = self.length / (num_nodes - 1)

    def __getattr__(self, name):
        if name in self.domain_vars:
            return self.domain_vars[name].value
        if re.match(r"d(.*?)dx", name):
            var_name = re.findall(r"d(.*?)dx", name)[0]
            return self.domain_vars[var_name].ddx
        if re.match(r"d(.*?)dt", name):
            var_name = re.findall(r"d(.*?)dt", name)[0]
            return self.domain_vars[var_name].ddt
        raise AttributeError("Attribute {} not found".format(name))

    def __setattr__(self, name, value):
        if re.match(r"d(.*?)dt", name):
            var_name = re.findall(r"d(.*?)dt", name)[0]
            self.domain_vars[var_name].ddt = value
        object.__setattr__(self, name, value)

    def add_var_to_domain(self, name, var):
        self.domain_vars.update({name: var})

    def import_array_to_vars(self, input_array):
        # TODO: vectorize this
        in_arr = np.array(input_array, dtype=np.float32)
        in_reshaped = in_arr.reshape(len(self.domain_vars), self.num_nodes)
        for i, (k, v) in enumerate(self.domain_vars.iteritems()):
            v.value = in_reshaped[i, :]
        self.invalidate_cache()

    def export_vars_to_array(self):
        # TODO: vectorize this
        out_arr = np.zeros((len(self.domain_vars), self.num_nodes), dtype=np.float32)
        for i, (k, v) in enumerate(self.domain_vars.iteritems()):
            out_arr[i, :] = v.ddt
        return out_arr.reshape(len(self.domain_vars) * (self.num_nodes))

    def initialize_vars(self):
        out_arr = np.zeros((len(self.domain_vars), self.num_nodes), dtype=np.float32)
        for i, (k, v) in enumerate(self.domain_vars.iteritems()):
            out_arr[i, :] = v.value
        return out_arr.reshape(len(self.domain_vars) * (self.num_nodes))

    def invalidate_cache(self):
        for (k, v) in self.domain_vars.iteritems():
            cached_property.invalidate_cache(v)
