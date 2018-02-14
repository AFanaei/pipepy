import json
import os

component_list = None


class Component(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @staticmethod
    def from_db(name):
        address = os.path.join(
            os.path.dirname(__file__), os.pardir, 'databases', "c.json")

        with open(address) as fp:
            data = json.load(fp)
            component_data = filter(lambda x: x['name'] == name, data)[0]

        return Component(**component_data)


class ComponentList(object):
    components = None

    def __init__(self, *args):
        self.components = args

    def __getattr__(self, name):
        return [getattr(c, name) for c in self.components]
