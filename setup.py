
from skbuild import setup
import versioneer
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
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    long_description_content_type = "text/markdown",
    description=description,
    long_description=long_description,
    author=['Hugo MacDermott-Opeskin', "Richard Gowers"],
    license="MIT",
    packages=['distopia'],
    python_requires=">=3.8",
    keywords=(
        "molecular dynamics distances simulation SIMD"
    ),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    url="https://github.com/MDAnalysis/distopia",
    download_url="https://pypi.org/project/distopia/",
    project_urls={
        "Homepage": "https://github.com/MDAnalysis/distopia",
        "Documentation": "https://www.mdanalysis.org/distopia/",
        "Source Code": "https://github.com/MDAnalysis/distopia",
        "Issue Tracker": "https://github.com/MDAnalysis/distopia/issues",
    },
    install_requires=[
        "numpy>=1.20.0",
        "cython>=0.28.0,<3.0.0",
        "scikit-build",
        "cmake"
    ],
)
