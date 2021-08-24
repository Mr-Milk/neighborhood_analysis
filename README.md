# Neighborhood analysis

> This package is no longer maintained, the functionalities can be found in [spatialtis_core](https://github.com/Mr-Milk/spatialtis-core)

<p align="center">
	<img src="https://img.shields.io/github/workflow/status/Mr-Milk/neighborhood_analysis/Build?style=for-the-badge&logo=github" alt="Build Status"/>
  <img src="img/platform.svg" alt="Platform"/>
	<br/><br/>
	<img src="https://img.shields.io/pypi/v/neighborhood_analysis?style=for-the-badge" alt="PYPI"/>
  <img src="img/python_version.svg" alt="Python version"/>
</p>


A python implementation of neighborhood analysis to profile cell-cell interaction using permutation test.

The parallelism has been implemented, so no need to do multi-processing yourself.

## Installation

Prebuild wheels for: Windwos, MacOS and Linux in 64bit, on Python version 3.6, 3.7, 3.8, 3.9

Requirements: Python >= 3.6

### From pypi

Normally, there should be a compatible wheel for your system. Just run:

```shell script
pip install neighborhood_analysis
```

If not, it will try to compile from a source, you need to install dependencies.

```shell script
# for windows
choco install rustup.install
# for Unix
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

pip install maturin

pip install neighborhood_analysis
```


### From source

Assume you have all the above dependencies, clone the repo, and then run:

```shell script
pip install .
```

## Usage

```python
import numpy as np
from neighborhood_analysis import CellCombs, get_point_neighbors, comb_bootstrap

# Get 10000 points, represent cell location
points = np.random.randint(0, 1000, (10000, 2))

# Assign each cell with a type
types_pool = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
types = np.random.choice(types_pool, 10000)

# Find cell neighbors at radius set at 10
# The points must be a list of tuple, otherwise it will raise TypeError
points = [(x, y) for (x, y) in points]
neighbors = get_point_neighbors(points, 10.0)

cc = CellCombs(types)
result = cc.bootstrap(types, neighbors)
# On my dual-core mac, this step takes 6~7 seconds.

X = [bool(i) for i in np.random.choice([True, False], 10000)]
Y = [bool(i) for i in np.random.choice([True, False], 10000)]
# The types must be a list of bool
z_score = comb_bootstrap(X, Y, neighbors, ignore_self=True)

```

## Documentation

<details>
<summary>Neighborhood analysis</summary>

```python
class CellCombs:

    def __init__(self, types, order=False):
      """
      Constructor function

      Args:
          types: List[str]; All the type of cells in your research
          order: bool (False); If False, A->B and A<-B is the same.

      Return:
          self
      """
    
    def bootstrap(self, types, neighbors, times=500, pval=0.05, method="pval", ignore_self=False):
      """
      Bootstrap functions

      If method is 'pval', 1.0 means association, -1.0 means avoidance, 0.0 means insignificance.
      If method is 'zscore', results is the exact z-score value.

      Args:
          types: List[str]; The type of all the cells
          neighbors: List[List[int]]; eg. {1:[4,5], 2:[6,7]}, cell at index 1 has neighbor cells from index 4 and 5
          times: int (500); How many times to perform bootstrap
          pval: float (0.05); The threshold of p-value
          method: str ('pval'); 'pval' or 'zscore'
          ignore_self: bool (False); Whether to consider self as a neighbor

      Return:
          List of tuples, eg.(('a', 'b'), 1.0), the type a and type b has a relationship as association
      """
```
</details>


<details>
<summary>Neighbors search utility functions</summary>

```python
def get_bbox(points_collections):
  """
  A utility function to return minimum bounding box list of polygons

  Args:
      points_collections: List[List[(float, float)]]; List of 2d points collections

  Return:
      A list of bounding box
  """
```

```python
def get_point_neighbors(points, r):
  """
  A utility function to search for point neighbors using kd-tree

  Args:
      points: List[tuple(float, float)]; Two dimension points
      r: float; The search radius

  Return:
      A list of neighbors' index, return as the order of the input
  """
```

```python
def get_bbox_neighbors(bbox_list, expand=1.0, scale=1.0):
  """
  A utility function to search for bbox neighbors using r-tree

  Args:
      bbox_list: List[tuple(float, float, float, float)]; The minimum bounding box of any polygon
              (minx, miny, maxx, maxy)
      expand: float; The expand unit
      scale: float; The scale fold number

  Return:
      A list of neighbors' index, return as the order of the input
    """
```

</details>

<br><br>
