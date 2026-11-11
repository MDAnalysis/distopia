# Building and testing

You can build `distopia` from source.

## Checking out the sources
Get the latest sources from the repository https://github.com/MDAnalysis/distopia
```bash
git clone --recurse-submodules https://github.com/MDAnalysis/distopia.git
cd distopia
```

Note that we are using git submodules for the vendored [highway](https://github.com/google/highway)
library (in `libdistopia/highway`, release 1.4.0) and for Google tests and Google benchmarks (in 
`libdistopia/googletest` and `libdistopia/googlebench`).

If the submodules are updated to a more recent version and your local tree has fallen out of sync
then you can run
```bash
git submodule update --init --force --recursive
```


## Building 

Distopia uses the scikit-build build system. You can build the library easily by first installing the prerequisites.

```bash
mamba env create --file devtools/conda_envs/distopia_<platform>-latest.yaml
mamba activate distopia
pip install . 
```

If you want greater control over build and install options you can do a more manual installation specifying CMake args like so:

```bash
python setup.py install -- -DCMAKE_BUILD_TYPE=Release <etc etc>
```

Special CMake flags
- `DISTOPIA_USE_EXTERNAL_HWY`:
  - Off (default) — build with the vendored highway library
  - On — link to an installed libhwy (install with `mamba install libhwy`)
- `DISTOPIA_BUILD_TESTS`:
  - On (default) — build the tests for libdistopia and the benchmarks
  - Off — do not build tests and benchmarks

## Testing

Testing for distopia is done on two levels. The `distopia` python layer has `pytest` tests (which are installed
with `mamba` from conda-forge if you followed the steps above) you can run after installing the package.

```bash
cd distopia/tests
pytest -vvv 
```

Testing for the `libdistopia` C++ layer is done with `googletest` executables. First do a build of the package, then use `ctest` to execute the relevant binaries in `_skbuild/<platform>/cmake-build/libdistopia`. Alternatively you can execute the binaries directly.

```bash
python setup.py build
ctest --test-dir _skbuild/*/cmake-build/libdistopia
```

## Benchmarking

Automated benchmarks for `libdistopia` are available, written using `googlebench`. Similarly to testing first do a build of the library, and then execute the benchmarking binary.

```bash
python setup.py build
cd _skbuild/<platform>/cmake-build/libdistopia
./bench
```
