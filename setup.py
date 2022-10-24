
from skbuild import setup
import versioneer
import configparser
from distutils.dist import Distribution

# Get our own instance of Distribution
dist = Distribution()
dist.parse_config_files()
dist.parse_command_line()

def get_dispatch():
    # need to get config value for dispatch before running setup.py
    config = configparser.ConfigParser()
    config.read('./setup.cfg')
    dispatch = config['options']['dispatch'].upper()
    if dispatch == "TRUE":
        dispatch = True
    elif dispatch == "FALSE":
        dispatch = False
    else:
        raise Exception("Invalid option for dispatch in setup.cfg, must be"
                        "True or False")
    if dispatch:
        cmake_dispatch_args = ['-DDISTOPIA_DISPATCH=ON']
    else:
        cmake_dispatch_args = []
    
    return cmake_dispatch_args

cmake_dispatch_args = get_dispatch()

setup(
    name="distopia",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Fast distance calculations using explicitly vectorised SIMD",
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
    cmake_args=cmake_dispatch_args,
    install_requires=[
        "numpy>=1.20.0",
        "cython>=0.28.0,<3.0.0",
        "scikit-build",
        "cmake"
    ],
)
