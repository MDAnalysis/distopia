from skbuild import setup
import os

description = "Fast distance calculations using explicitly vectorised SIMD"
try:
    readme_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md")
    with open(readme_file) as f:
        long_description = f.read()
except ImportError:
    long_description = description

setup(
    name="distopia",
    packages=['distopia'],
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.20.0"
    ],
)