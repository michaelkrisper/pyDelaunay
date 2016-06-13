from typing import List, Tuple
from math import isclose


class Point:
    def __init__(self, x: float, y: float, z: float):
        self.z = z
        self.x = x
        self.y = y


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.p3 = p3
        self.p2 = p2
        self.p1 = p1


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
        self.points = points
        supertriangle = Triangle(Point(-1000, -1000, 0),
                                 Point(1000, 1000, 0),
                                 Point(1000, 0, 0))

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
