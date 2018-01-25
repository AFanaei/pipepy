class cached_property(object):
    def __init__(self, func, name=None):
        self.func = func
        self.__doc__ = getattr(func, '__doc__')
        self.name = name or func.__name__

    def __get__(self, instance, cls=None):
        if instance is None:
            return self
        res = instance.__dict__[self.name] = self.func(instance)
        return res

    @staticmethod
    def invalidate_cache(obj):
        for property_name in dir(type(obj)):
            if isinstance(getattr(type(obj), property_name), cached_property):
                obj.__dict__.pop(property_name, None)