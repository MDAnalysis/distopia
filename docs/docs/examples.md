# Examples

Distopia aims to be a small and simple library for calculating kernels common in analysis of molecular simulations. 

Simple examples of use are given below. 

## Distances

First we construct some random coordinates using numpy and calculate distances between them without periodic boundary conditions using `distopia`. You can use `np.float32` for single precision and `np.float64` for double precision. 

### No periodic boundary conditions

```python
import numpy as np
import distopia

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
result = distopia.calc_bonds_no_box(coordinates0, coordinates1)

# alternatively we can pass in a buffer to use for the results.
buffer = np.empty(N, dtype=np.float32)
result = distopia.calc_bonds_no_box(coordinates0, coordinates1, results=buffer)
```

### Orthorhombic periodic boundary conditions

Using periodic boundary conditions is very similar. For orthorhombic boxes you only need to provide the three box vector lengths `[l, l, l]`

```python
import numpy as np
import distopia

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
box = np.asarray([10, 10, 10]).astype(np.float32)
result = distopia.calc_bonds_ortho(coordinates0, coordinates1, box)
```

### Triclinic periodic boundary conditions

For triclinic boxes the box matrix must be provided in 3x3 matrix form. 

```python
import numpy as np
import distopia

N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
box = np.asarray([[10, 0, 0], [0, 10, 0], [0, 0, 10]]).astype(np.float32)
result = distopia.calc_bonds_triclinic(coordinates0, coordinates1, box)
```

### Note

All of the below functions also support orthorhombic and triclinic systems within the same function naming convention but the no-pbc versions are used for demonstration purposes. See the [API](api.md) docs for more details.

## Angles 

Angles function in a similar way to distances, but requiring one more coordinate. 

```python
import numpy as np
import distopia

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates2 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
result = distopia.calc_angles_no_box(coordinates0, coordinates1, coordinates2)
```


## Dihedrals 

Dihedrals require 4 coordinates 

```python
import numpy as np
import distopia

# make N x 3 coordinate arrays
N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates2 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates3 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
result = distopia.calc_dihedrals_no_box(coordinates0, coordinates1, coordinates2, coordinates3)
```

## Distance arrays

You can also do pairwise and self-pairwise distance arrays with `distopia`. 

### Pairwise distances

```python
import numpy as np
import distopia

# make N x 3 and M x 3 coordinate arrays
N = 10000
M = 1000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
coordinates1 = np.random.rand(3 * M).reshape(M, 3).astype(np.float32)
result = distopia.calc_distance_array_no_box(coordinates0, coordinates1)
result # -> will be NxM

# passing in a result buffer is also possible for distance arrays 
buffer = np.empty((N,M), dtype=np.float32)
result = distopia.calc_distance_array_no_box(coordinates0, coordinates1, results=buffer)
```

### Self-pairwise distances

Self distance arrays are similar but use only a single coordinate producing an N*(N-1)/2 array of distances, a flattened upper triangle of the full NxN matrix with the diagonal removed.

```python
import numpy as np
import distopia

# make N x 3 coordinate array
N = 10000
coordinates0 = np.random.rand(3 * N).reshape(N, 3).astype(np.float32)
result = distopia.calc_self_distance_array_no_box(coordinates0)
result # -> will be NxN with result
```

To recover the full NxN matrix use the following

```python
distance_matrix = np.zeros((N,N))
k = 0
for i in range(N):
    for j in range(i + 1, N):
        distance_matrix[i, j] = result[k]
        k += 1
```



## Questions

Please raise any questions or issues on the [issue tracker]((https://github.com/MDAnalysis/distopia/issues)). 