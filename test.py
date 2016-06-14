import unittest
from math import isclose

from delaunay import DelaunayMap


class Test(unittest.TestCase):
    def setUp(self):
        self.points = [
            (0, 0, 0),
            (1, 0, 0.5),
            (0, 1, 0.5),
            (1, 1, 1),

        ]
        self.d = DelaunayMap(self.points)

    def test_usage(self):
        d = DelaunayMap([(0, 0, 0), (1, 0, 1), (0, 1, 1)])
        assert d[0.25, 0.25] == 0.5

    def test1(self):
        assert self.d[0.25, 0.25] == 0.25

    def test2(self):
        assert self.d[0.5, 0.5] == 0.5

    def test3(self):
        assert self.d[0.5, 0] == 0.25

    def test4(self):
        assert self.d[0, 0.5] == 0.25

    def test5(self):
        assert self.d[0.125, 0.125] == 0.125

    def test_fixpoints(self):
        for p in self.points:
            assert self.d[p[0], p[1]] == p[2]


if __name__ == '__main__':
    unittest.main()
