import re
from collections import OrderedDict

import numpy as np

from utils.functional import cached_property


class DomainVar(object):
    def __init__(self, name, domain):
        self.domain = domain
        self.value = np.zeros(domain.num_nodes)
        self.ddt = np.zeros(domain.num_nodes - 1)

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
            var_name = re.findall(r"d(.*?)dx", name)[0]
            return self.domain_vars[var_name].ddt
        raise AttributeError("Attribute {} not found".format(name))

    def __setattr__(self, name, value):
        if re.match(r"d(.*?)dt", name):
            var_name = re.findall(r"d(.*?)dx", name)[0]
            self.domain_vars[var_name].ddt = value
        object.__setattr__(self, name, value)

    def add_var_to_domain(self, name, var):
        self.domain_vars.append((name, var))

    def import_array_to_vars(self, input_array):
        in_arr = np.array(input_array)
        in_reshaped = in_arr.reshape(len(self.domain_vars), self.num_nodes - 1)
        for i, (k, v) in enumerate(self.domain_vars.iteritems()):
            v.value = in_reshaped[i, :]
        self.invalidate_cache()

    def export_vars_to_array(self):
        out_arr = np.zeros(len(self.domain_vars), self.num_nodes - 1)
        for i, (k, v) in enumerate(self.domain_vars.iteritems()):
            out_arr[i, :] = v.value
        return out_arr

    def invalidate_cache(self):
        for (k, v) in self.domain_vars.iteritems():
            cached_property.invalidate_cache(v)
