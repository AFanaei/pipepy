import unittest

from pipeline.pipe import Pipe, Zcalculator


class ZTest(unittest.TestCase):
    def setUp(self):
        self.z = Zcalculator.instance()
        self.points=[
            [1, -10, 0.9967102],
            [115, -10, 0.7071962],
            [40, 10, 0.9012377],
            [75, 20, 0.8507508],
            [65, 30, 0.8822946],
            [95, 40, 0.863464],
            [120, 50, 0.8638056],
        ]
        self.interp = [
            [2, 20, 0.9953242],
            [22, 23, 0.9518844],
            [33, 17, 0.924142],
            [42, 18, 0.9063737],
            [53, 18, 0.8848565],
        ]

    def testZ(self):
        for i in self.points:
            res = self.z.get_z(i[0], i[1])
            self.assertLessEqual(abs(res-i[2]), i[2]*0.001, "diffrence in result")

    def testZInterp(self):
        for i in self.interp:
            res = self.z.get_z(i[0], i[1])
            self.assertLessEqual(abs(res-i[2]), i[2]*0.001, "diffrence in result")


if __name__ == '__main__':
    unittest.main()
