import numpy as np
import neighborhood_analysis as na
from neighborhood_analysis import CellCombs, get_neighbors

from time import time

types = ['a', 'b', 'c', 'd', 'e', 'f', 'aa', 'cc', 'dd']
points = np.random.randint(0, 1000, (10000, 2))
corr_types = np.random.choice(types, 10000)
points = [(x, y) for (x, y) in points]
neighbors = get_neighbors(points, 10.0)

start = time()
cc = CellCombs(types)
cc.bootstrap(corr_types, neighbors)
end = time()
print(f"used {(end-start):.2f}s")


