import unittest
from math import isclose
from delaunay import Point, DelaunayMap


class Test(unittest.TestCase):
    def setUp(self):
        points = [Point(0, 0, 0),
                  Point(1, 0, 1),
                  Point(0, 1, 1),
                  Point(1, 1, 1)]
        self.d = DelaunayMap(points)

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
