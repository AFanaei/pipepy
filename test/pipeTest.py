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

    def testZ(self):
        for i in self.points:
            res = self.z.get_z(i[0], i[1])
            self.assertLessEqual(abs(res-i[2]), 0.01, "diffrence in result")


if __name__ == '__main__':
    unittest.main()
