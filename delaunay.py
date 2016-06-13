from collections import Counter
from itertools import groupby
from typing import List, Tuple
from math import isclose


class Point:
    def __init__(self, x: float, y: float, z: float):
        self.Z = z
        self.X = x
        self.Y = y


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.p3 = p3
        self.p2 = p2
        self.p1 = p1

    def contains_in_circumcircle(self, p: Point):
        p0X = self.p1.X - p.X
        p0Y = self.p1.Y - p.Y
        p1X = self.p2.X - p.X
        p1Y = self.p2.Y - p.Y
        p2X = self.p3.X - p.X
        p2Y = self.p3.Y - p.Y
        p0Square = p0X * p0X + p0Y * p0Y
        p1Square = p1X * p1X + p1Y * p1Y
        p2Square = p2X * p2X + p2Y * p2Y
        det01 = p0X * p1Y - p1X * p0Y
        det12 = p1X * p2Y - p2X * p1Y
        det20 = p2X * p0Y - p0X * p2Y
        result = p0Square * det12 + p1Square * det20 + p2Square * det01
        return result > 0

class Edge:
    def __init__(self, p1: Point, p2: Point):
        self.p2 = p2
        self.p1 = p1


class Plane:
    def __init__(self, tr: Triangle):
        # todo
        pass

class DelaunayMap:
    def __init__(self, points: List[Point]):
        assert len(points) > 3

        supertriangle = Triangle(Point(-1000, -1000, 0),
                                 Point(1000, 1000, 0),
                                 Point(1000, 0, 0))
        triangles = [supertriangle]

        for p in points:
            container_triangles = [tr for tr in triangles if tr.contains_in_circumcircle(p)]
            for tr in container_triangles:
                triangles.remove(tr)

            edges = [Edge(x, y)
                     for tr in container_triangles
                     for x, y in ((tr.P1, tr.P2), (tr.P2, tr.P3), (tr.P3, tr.P1))]
            c = Counter(edges)
            unique_edges = []
            for count, e in c.items():
                if count == 1:
                    unique_edges.append(e)

            newTriangles = []
            for e in unique_edges:
                newTriangles.append(Triangle(e.P1, e.P2, p))

            triangles += newTriangles
        self.triangles = [tr for tr in triangles if not tr.shares_vertex(supertriangle)]

    def __getitem__(self, xy: Tuple[float, float]):
        x, y = xy
        for tr in self.triangles:
            if tr.is_inside(x, y):
                pl = Plane(tr)
                return x, y, (pl.W - pl.X * x - pl.Y * y) / pl.Z

if __name__ == '__main__':
    points = [Point(0, 0, 0),
              Point(1, 0, 1),
              Point(0, 1, 1),
              Point(1, 1, 1)]
    d = DelaunayMap(points)

    print(d[0.5, 0.5])
    assert d[0.5, 0.5] == 0.5

    print("finished")
