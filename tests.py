import numpy as np
import neighborhood_analysis as na
from neighborhood_analysis import CellCombs, get_bbox, get_point_neighbors, get_bbox_neighbors,comb_bootstrap

from time import time

types = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
points = np.random.randint(0, 1000, (10000, 2))
corr_types = np.random.choice(types, 10000)
points = [(x, y) for (x, y) in points]

polygons = []
for _ in range(10000):
    ixs = np.random.choice(range(len(points)), 5)
    polygon = []
    for x in ixs:
        polygon.append(points[x])

start = time()
bbox = get_bbox(polygons)
end = time()
print(f"Get bbox used {(end-start):.2f}s")

start = time()
neighbors = get_bbox_neighbors(bbox, 2)
end = time()
print(f"search bbox neighbors used {(end-start):.2f}s")


start = time()
neighbors = get_point_neighbors(points, 10.0)
end = time()
print(f"search point neighbors used {(end-start):.2f}s")

start = time()

cc = CellCombs(types, False)
results = cc.bootstrap(corr_types, neighbors, ignore_self=True)
print(results)

end = time()
print(f"used {(end-start):.2f}s")

s1 = time()
X = [bool(i) for i in np.random.choice([True, False], 10000)]
Y = [bool(i) for i in np.random.choice([True, False], 10000)]
v = comb_bootstrap(X, Y, neighbors, ignore_self=True)
s2 = time()
print(f"used {(s2-s1):.2f}s")
