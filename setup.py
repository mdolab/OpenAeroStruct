from setuptools import setup

import re
import os

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("openaerostruct/__init__.py").read(),
)[0]

optional_dependencies = {
    "docs": ["sphinx_mdolab_theme"],
    "test": ["testflo"],
    "ffd": ["pygeo>=1.6.0"],
    "mphys": ["mphys>=1.0.0"],
}

# Add an optional dependency that concatenates all others
optional_dependencies["all"] = sorted(
    [dependency for dependencies in optional_dependencies.values() for dependency in dependencies]
)

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="openaerostruct",
    version=__version__,
    description="OpenAeroStruct",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mdolab/OpenAeroStruct",
    keywords="",
    author="",
    author_email="",
    license="BSD-3",
    packages=[
        "openaerostruct",
        "openaerostruct/docs",
        "openaerostruct/docs/_utils",
        "openaerostruct/geometry",
        "openaerostruct/structures",
        "openaerostruct/aerodynamics",
        "openaerostruct/transfer",
        "openaerostruct/functionals",
        "openaerostruct/integration",
        "openaerostruct/common",
        "openaerostruct/utils",
        "openaerostruct/mphys",
    ],
    # Test files
    package_data={"openaerostruct": ["tests/*.py", "*/tests/*.py", "*/*/tests/*.py"]},
    install_requires=[
        # Remember to update the oldest versions in the GitHub Actions build, the readme, and in docs/installation.rst
        "openmdao>=3.15",
        "numpy>=1.20",
        "scipy>=1.6.0",
        "matplotlib",
    ],
    extras_require=optional_dependencies,
    zip_safe=False,
    # ext_modules=ext,
    entry_points="""
    [console_scripts]
    plot_wing=openaerostruct.utils.plot_wing:disp_plot
    plot_wingbox=openaerostruct.utils.plot_wingbox:disp_plot
    """,
)
