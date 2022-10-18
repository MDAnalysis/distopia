
from skbuild import setup
import versioneer

setup(
    name="distopia",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Fast distance calculations using explicitly vectorised SIMD",
    author=['Hugo MacDermott-Opeskin', "Richard Gowers"],
    license="MIT",
    packages=['distopia'],
    python_requires=">=3.7",
)