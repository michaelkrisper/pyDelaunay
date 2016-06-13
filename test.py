import unittest
from math import isclose

from delaunay import DelaunayMap


class Test(unittest.TestCase):
    def setUp(self):
        points = [(0, 0, 0),
                  (1, 0, 1),
                  (0, 1, 1),
                  (1, 1, 1)]
        self.d = DelaunayMap(points)

    def test_usage(self):
        d = DelaunayMap([(0, 0, 0), (1, 0, 1), (0, 1, 1)])
        result = d[0.25, 0.25]
        assert result == 0.5

    def test1(self):
        assert isclose(self.d[0.25, 0.25], 0.5)

    def test2(self):
        assert isclose(self.d[0.5, 0.5], 1)

    def test3(self):
        assert isclose(self.d[0.5, 0], 0.5)

    def test4(self):
        assert isclose(self.d[0, 0.5], 0.5)

    def test5(self):
        assert isclose(self.d[0.125, 0.125], 0.25)


if __name__ == '__main__':
    unittest.main()
