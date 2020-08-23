from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize


libdistopia = Extension(
    name='distopia.libdistopia',
    sources=['distopia/libdistopia.pyx'],
    include_dirs=['../'],
    library_dirs=['../build'],
    runtime_library_dirs=['../build/'],
    libraries=['distopia'],
    language='c++',
)


setup(
    name='distopia',
    version='0.0.1',
    packages=find_packages(),
    ext_modules=cythonize([libdistopia],
                          compiler_directives={'language_level': 3}),
)
