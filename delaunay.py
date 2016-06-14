#! /usr/bin/env python3
# coding=utf-8
"""
A basic implementation of delaunay triangulation.

usage:
    d = DelaunayMap([(0, 0, 0), (1, 0, 1), (0, 1, 1)])
    d[0.25, 0.25]
 -> 0.5
"""

from collections import Counter
from typing import List, Tuple

import time

__author__ = "Michael Krisper"
__email__ = "michael.krisper@gmail.com"
__date__ = "2016-06-14"


class Point:
    def __init__(self, x: float, y: float, z: float):
        self.X = x
        self.Y = y
        self.Z = z

    def __repr__(self):
        return "({X},{Y},{Z})".format(**vars(self))

    def __sub__(self, other: 'Point') -> 'Point':
        return Point(self.X - other.X, self.Y - other.Y, self.Z - other.Z)

    def __eq__(self, other):
        return (self.X, self.Y, self.Z) == (other.X, other.Y, other.Z)

    def __hash__(self):
        return hash((self.X, self.Y, self.Z))


class Edge:
    def __init__(self, p1: Point, p2: Point):
        self.P1 = p1
        self.P2 = p2

    def __hash__(self):
        return hash(self.P1) ^ hash(self.P2)

    def __eq__(self, other):
        return (self.P1 == other.P1 and self.P2 == other.P2) or (self.P2 == other.P1 and self.P1 == other.P2)

    def __repr__(self):
        return "({P1}->{P2})".format(**vars(self))


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        if (p1.X == p2.X == p3.X) or (p1.Y == p2.Y == p3.Y):
            raise Exception("triangle is coplanar --> cannot extrapolated by a plane.")

        self.P1 = p1
        self.P2 = p2
        self.P3 = p3

    def __repr__(self):
        return "({P1},{P2},{P3})".format(**vars(self))

    def contains_in_circumcircle(self, p: Point) -> bool:
        p0 = self.P1 - p
        p1 = self.P2 - p
        p2 = self.P3 - p
        det01 = p0.X * p1.Y - p0.Y * p1.X
        det12 = p1.X * p2.Y - p1.Y * p2.X
        det20 = p2.X * p0.Y - p2.Y * p0.X
        return (p0.X ** 2 + p0.Y ** 2) * det12 + (p1.X ** 2 + p1.Y ** 2) * det20 + (p2.X ** 2 + p2.Y ** 2) * det01 > 0

    def shares_vertex(self, tr: 'Triangle') -> bool:
        return {self.P1, self.P2, self.P3} & {tr.P1, tr.P2, tr.P3}

    def is_inside(self, x: float, y: float) -> bool:
        v0 = self.P3 - self.P1
        v1 = self.P2 - self.P1
        v2 = Point(x, y, 0) - self.P1
        dot00 = v0.X ** 2 + v0.Y ** 2
        dot01 = v0.X * v1.X + v0.Y * v1.Y
        dot02 = v0.X * v2.X + v0.Y * v2.Y
        dot11 = v1.X ** 2 + v1.Y ** 2
        dot12 = v1.X * v2.X + v1.Y * v2.Y
        inv_denom = 1 / dot00 * dot11 - dot01 ** 2
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom
        return u >= 0 and v >= 0 and u + v <= 1


class Plane:
    def __init__(self, tr: Triangle):
        ab = tr.P2 - tr.P1
        ac = tr.P3 - tr.P1
        self.X = ab.Y * ac.Z - ab.Z * ac.Y
        self.Y = ab.Z * ac.X - ab.X * ac.Z
        self.Z = ab.X * ac.Y - ab.Y * ac.X
        self.W = tr.P1.X * self.X + tr.P1.Y * self.Y + tr.P1.Z * self.Z

    def __repr__(self):
        return "({X},{Y},{Z},{W})".format(**vars(self))


class DelaunayMap:
    def __init__(self, _points: List[Tuple[float, float, float]]):
        assert len(_points) >= 3
        points = [Point(*p) for p in _points]
        min_x, min_y = min(p.X for p in points), min(p.Y for p in points)
        max_x, max_y = max(p.X for p in points), max(p.Y for p in points)
        supertriangle = Triangle(Point(min_x - 1, min_y - 1, 0), Point(max_x * 3, 0, 0), Point(min_x - 1, max_y * 3, 0))
        triangles = [supertriangle]

        for p in points:
            container_triangles = [tr for tr
                                   in triangles
                                   if tr.contains_in_circumcircle(p)]
            edges = []  # type: List[Edge]
            for tr in container_triangles:
                triangles.remove(tr)
                edges += [Edge(tr.P1, tr.P2), Edge(tr.P2, tr.P3), Edge(tr.P3, tr.P1)]

            unique_edges = [edge
                            for edge, count in Counter(edges).items()
                            if count == 1]

            triangles += [Triangle(e.P1, e.P2, p) for e in unique_edges]

        self.triangles = [tr
                          for tr in triangles
                          if not tr.shares_vertex(supertriangle)]

    def __getitem__(self, xy: Tuple[float, float]) -> float:
        x, y = xy
        for tr in self.triangles:  # type: Triangle
            if tr.is_inside(x, y):
                pl = Plane(tr)
                return (pl.W - pl.X * x - pl.Y * y) / pl.Z


if __name__ == '__main__':
    start = time.clock()
    d = DelaunayMap([(0, 0, 0), (1, 0, 1), (0, 1, 1)])
    result = d[0.25, 0.25]
    print(time.clock() - start)
