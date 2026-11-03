import random
import unittest
from math import isclose

from delaunay import DelaunayMap, Point, Triangle


def random_points(n, seed=42, spread=100.0):
    rng = random.Random(seed)
    xy = set()
    while len(xy) < n:
        xy.add((round(rng.uniform(-spread, spread), 3), round(rng.uniform(-spread, spread), 3)))
    return [(x, y, 2 * x - 3 * y + 5) for x, y in sorted(xy)]


def convex_hull_size(points):
    pts = sorted({(p[0], p[1]) for p in points})

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def half(seq):
        hull = []
        for p in seq:
            while len(hull) >= 2 and cross(hull[-2], hull[-1], p) <= 0:
                hull.pop()
            hull.append(p)
        return hull

    return len(half(pts)) + len(half(reversed(pts))) - 2


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
        points = random_points(150)
        d = DelaunayMap(*points)

        for x, y, z in points:
            assert isclose(z, d[x, y], abs_tol=1e-9)

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
        with self.assertRaisesRegex(Exception, "must be unique"):
            DelaunayMap((0, 0, 0), (1, 0, 1), (1, 1, 3), (0, 1, 2), (1, 1, 5))

    def test_call_with_list(self):
        DelaunayMap(*[(0, 0, 0, 4), (1, 0, 1), (1, 1, 3), (0, 1, 2)])

    def test_too_few_points(self):
        with self.assertRaisesRegex(Exception, "at least 3 points"):
            DelaunayMap((0, 0, 0), (1, 0, 1))

    def test_negative_coordinates(self):
        d = DelaunayMap((-10, -10, 0), (-5, -10, 5), (-10, -5, 5))
        assert isclose(d[-8, -8], 4)

    def test_outside_convex_hull_returns_none(self):
        d = DelaunayMap((0, 0, 0), (1, 0, 1), (0, 1, 1))
        assert d[2, 2] is None
        assert d[-0.1, 0] is None

    def test_triangle_count_matches_euler_formula(self):
        points = random_points(100, seed=7)
        d = DelaunayMap(*points)
        triangles = list(d.triangle_indices)
        assert len(triangles) == 2 * len(points) - 2 - convex_hull_size(points)

    def test_legacy_alias(self):
        d = DelaunayMap((0, 0, 0), (1, 0, 1), (0, 1, 1))
        assert list(d.triangles_indizes) == list(d.triangle_indices)

    def test_triangles_cover_all_points(self):
        points = random_points(80, seed=3)
        d = DelaunayMap(*points)
        used = {i for tr in d.triangle_indices for i in tr}
        assert used == set(range(len(points)))

    def test_empty_circumcircle_property(self):
        points = random_points(60, seed=11)
        d = DelaunayMap(*points)
        vertices = [Point(p[0], p[1], p[2], i) for i, p in enumerate(points)]
        for a, b, c in d.triangle_indices:
            tr = Triangle(vertices[a], vertices[b], vertices[c])
            for v in vertices:
                if v.index not in (a, b, c):
                    assert not tr.contains_in_circumcircle(v)

    def test_linear_surface_is_reproduced(self):
        points = random_points(100, seed=5, spread=10)
        d = DelaunayMap(*points)
        rng = random.Random(1)
        for _ in range(200):
            x, y = rng.uniform(-2, 2), rng.uniform(-2, 2)
            z = d[x, y]
            if z is not None:
                assert isclose(z, 2 * x - 3 * y + 5, abs_tol=1e-6)

    def test_circumcircle_independent_of_orientation(self):
        p1, p2, p3 = Point(0, 0, 0, 0), Point(1, 0, 0, 1), Point(0, 1, 0, 2)
        inside, outside = Point(0.5, 0.5, 0, 3), Point(2, 2, 0, 4)
        for tr in (Triangle(p1, p2, p3), Triangle(p1, p3, p2)):
            assert tr.contains_in_circumcircle(inside)
            assert not tr.contains_in_circumcircle(outside)

    def test_coplanar_triangle_rejected(self):
        with self.assertRaisesRegex(Exception, "coplanar"):
            Triangle(Point(0, 0, 0, 0), Point(0, 1, 0, 1), Point(0, 2, 0, 2))


if __name__ == '__main__':
    unittest.main()
