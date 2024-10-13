# Building and testing

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
