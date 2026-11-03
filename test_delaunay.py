import unittest
from delaunay import DelaunayMap


class Test(unittest.TestCase):
    def test_usage(self):
        d = DelaunayMap((0, 0, 0), (1, 0, 1), (0, 1, 1))
        assert d[0.25, 0.25] == 0.5

    def test0(self):
        d = DelaunayMap((0, 0, 0), (0, 1, 0), (1, 0, 0))
        assert d[0, 0] == 0
        assert d[1, 0] == 0
        assert d[0, 1] == 0
        assert d[0.5, 0.5] == 0
        assert d[0.25, 0.25] == 0

    def test1(self):
        d = DelaunayMap((0, 0, 1), (0, 1, 1), (1, 0, 1))
        assert d[0, 0] == 1
        assert d[1, 0] == 1
        assert d[0, 1] == 1
        assert d[0.5, 0.5] == 1
        assert d[0.25, 0.25] == 1

    def test2(self):
        d = DelaunayMap((0, 0, 0), (0, 1, 1), (1, 0, 1))
        assert d[0, 0] == 0
        assert d[1, 0] == 1
        assert d[0, 1] == 1
        assert d[0.5, 0.5] == 1
        assert d[0, 0.5] == 0.5
        assert d[0.5, 0] == 0.5
        assert d[0.25, 0.25] == 0.5
        assert d[0.125, 0.125] == 0.25
        assert d[0.25, 0] == 0.25
        assert d[0, 0.25] == 0.25

    def test_fixpoints(self):
        with open("map.csv") as f:
            points = [tuple(float(x) for x in line.strip().split(",")) for line in f.readlines()]
        d = DelaunayMap(*points)

        for x, y, z in points:
            assert z == d[x, y]

    def test_Simple_DelaunayMap(self):
        assert 0 == DelaunayMap([0, 0, 0], [1, 0, 0], [0, 1, 0])[0.25, 0.25]

    def test_DelaunayMapTriangle(self):
        map = DelaunayMap([0, 0, 0], [1, 0, 1], [0, 1, 2])

        # fixed points
        assert 0 == map[0, 0]
        assert 1 == map[1, 0]
        assert 2 == map[0, 1]

        # interpolations
        assert 0.5 == map[0.5, 0]
        assert 1 == map[0, 0.5]
        assert 1.5 == map[0.5, 0.5]
        assert 0.25 == map[0.25, 0]
        assert 0.5 == map[0, 0.25]
        assert 0.75 == map[0.25, 0.25]
        assert 0.75 == map[0.75, 0]
        assert 1.5 == map[0, 0.75]

    def test_DelaunayMapPlane(self):
        map = DelaunayMap([0, 0, 0, 3], [1, 0, 1], [0, 1, 2], [1, 1, 3])
        # fixed points
        assert 0 == map[0, 0]
        assert 1 == map[1, 0]
        assert 2 == map[0, 1]
        assert 3 == map[1, 1]

        # interpolations
        assert 0.5 == map[0.5, 0]
        assert 1 == map[0, 0.5]
        assert 2 == map[1, 0.5]
        assert 2.5 == map[0.5, 1]

        assert 1.5 == map[0.5, 0.5]

        assert 0.75 == map[0.25, 0.25]
        assert 2.25 == map[0.75, 0.75]

        assert 1.75 == map[0.25, 0.75]
        assert 1.25 == map[0.75, 0.25]

    def test_Delaunay_DuplicatePoints(self):
        try:
            DelaunayMap((0, 0, 0), (1, 0, 1), (1, 1, 3), (0, 1, 2), (1, 1, 5))
        except Exception:
            pass

    def test_call_with_list(self):
        DelaunayMap(*[(0, 0, 0, 4), (1, 0, 1), (1, 1, 3), (0, 1, 2)])


if __name__ == '__main__':
    unittest.main()
