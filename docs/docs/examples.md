# Examples

Distopia aims to be a small and simple library for calculating kernels common in analysis of molecular simulations. 

Simple examples of use are given below. 

## Distances

First we construct some random coordinates using numpy and calculate distances between them without periodic boundary conditions using `distopia`. You can use `np.float32` for single precision and `np.float64` for double precision. 

```python
import numpy as np

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random(3 * N, dtype=np.float32).reshape(N, 3)
coordinates1 = np.random(3 * N, dtype=np.float32).reshape(N, 3)
result = distopia.calc_bonds_no_box(coordinates0, coordinates1)

# alternatively we can pass in a buffer to use for the results.
result = distopia.calc_bonds_no_box(coordinates0, coordinates1, results=np.empty(N, dtype=np.float32))
```

Using periodic boundary conditions is very similar. For orthorhombic boxes you only need to provide the three box vector lengths `[l, l, l]`

```python
import numpy as np

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random(3 * N, dtype=np.float32).reshape(N, 3)
coordinates1 = np.random(3 * N, dtype=np.float32).reshape(N, 3)
box = np.ararray([10, 10, 10])
result = distopia.calc_bonds_ortho(coordinates0, coordinates1, box)
```

For triclinic boxes the box matrix must be provided in 3x3 matrix form. 

```python
box = np.asarray([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
result = distopia.calc_bonds_ortho(coordinates0, coordinates1, box)
```

## Angles 



## Dihedrals 


## Distance arrays
