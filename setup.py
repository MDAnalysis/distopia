from setuptools import Extension, setup
import numpy as np
from Cython.Build import cythonize


extensions = [
    Extension('distopia._distopia', ["python/distopia/_distopia.pyx"],
              include_dirs=['/home/richard/code/distopia2_the_highway_warrior/include',
                            np.get_include(),
                            ],
              libraries=['distopia', 'hwy'],
              library_dirs=['/home/richard/code/distopia2_the_highway_warrior/build',
                            '/home/richard/code/distopia2_the_highway_warrior/build/highway',
                            ],
              ),
]


setup(
    name="distopia",
    ext_modules=cythonize(extensions)
)
