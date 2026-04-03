# Building and testing

You can build `distopia` from source.

## Checking out the sources
Get the latest sources from the repository https://github.com/MDAnalysis/distopia
```bash
git clone https://github.com/MDAnalysis/distopia.git
cd distopia
```

We also need `highway`, which vendored as a git submodule:
```bash
git submodule update --init libdistopia/highway
```

For testing and benchmarking get `googlebench` and `googletest`:
```bash
git submodule update --init libdistopia/googlebench
git submodule update --init libdistopia/googletest
```


## Building 

Distopia uses the scikit-build build system. You can build the library easily by first installing the prerequisites.

```bash
mamba env create --file devtools/conda_envs/distopia_<platform>-latest.yaml
mamba activate distopia
pip install . 
```

If you want greater control over build and install options you can do a more manual install specifying CMake args like so:

```bash
python setup.py install -- -DCMAKE_BUILD_TYPE=Release <etc etc>
```

## Testing

Testing for distopia is done on two levels. The `distopia` python layer has `pytest` tests you can run after installing the package.

```bash
pip install pytest
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
