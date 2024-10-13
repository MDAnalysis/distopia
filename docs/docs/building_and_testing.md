# Building and testing


Distopia uses the scikit-build build system. You can build the library easily by first installing the prerequisites.

```bash
mamba env create --file devtools/conda_envs/distopia_<platform>-latest.yaml
mamba activate distopia
python