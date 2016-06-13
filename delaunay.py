from typing import List, Tuple


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


class DelaunayMap:
    def __init__(self, points: List[Point]):
        self.points = points
        supertriangle = Triangle(Point(-1000, -1000, 0),
                                 Point(1000, 1000, 0),
                                 Point(1000, 0, 0))

        p = points

    def __getitem__(self, xy: Tuple[float, float]):
        x, y = xy
        return x, y, 0


if __name__ == '__main__':
    points = [Point(0, 0, 0),
              Point(1, 0, 1),
              Point(0, 1, 1),
              Point(1, 1, 1)]
    d = DelaunayMap(points)

    print(d[0.1, 0.1])
    print("finished")
