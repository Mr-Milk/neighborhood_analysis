import numpy as np
from neighborhood_analysis import (CellCombs,
                                   get_bbox,
                                   get_point_neighbors,
                                   get_bbox_neighbors,
                                   comb_bootstrap,
                                   neighbor_components)
from time import time

types = [str(i) for i in range(30)]
points = np.random.randint(0, 1000, (10000, 2))
corr_types = np.random.choice(types, 10000)
points = [(x, y) for (x, y) in points]

polygons = []
for _ in range(100):
    ixs = np.random.choice(range(len(points)), 5)
    polygon = []
    for x in ixs:
        polygon.append(points[x])

start = time()
bbox = get_bbox(polygons)
end = time()
print(f"Get bbox used {(end - start):.5f}s")

start = time()
neighbors = get_bbox_neighbors(bbox, 2, labels=[i for i in range(0, len(bbox))])
end = time()
print(f"search bbox neighbors used {(end - start):.5f}s")

start = time()
neighbors = get_point_neighbors(points, 10.0)
end = time()
print(f"search point neighbors used {(end - start):.5f}s")

start = time()

cc = CellCombs(types, False)
results = cc.bootstrap(corr_types, neighbors, times=1000, ignore_self=True)

end = time()
print(f"neighborhood used {(end - start):.5f}s")

s1 = time()
X = [bool(i) for i in np.random.choice([True, False], 10000)]
Y = [bool(i) for i in np.random.choice([True, False], 10000)]
v = comb_bootstrap(X, Y, neighbors, times=1000, ignore_self=True)
s2 = time()
print(f"marker co-exp used {(s2 - s1):.5f}s")

cents = [i for i in range(len(points))]
s1 = time()
ix, col, data = neighbor_components(dict(zip(cents, neighbors)), dict(zip(cents, corr_types)))
s2 = time()
print(ix, col, data)
print(f"neighbor components used {(s2 - s1):.5f}s")
