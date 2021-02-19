from setuptools import setup

import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('openaerostruct/__init__.py').read(),
)[0]
setup(name='openaerostruct',
    version=__version__,
    description='OpenAeroStruct',
    author='John Jasa',
    author_email='johnjasa@umich.edu',
    license='BSD-3',
    packages=[
        'openaerostruct',
        'openaerostruct/docs',
        'openaerostruct/docs/_utils',
        'openaerostruct/geometry',
        'openaerostruct/structures',
        'openaerostruct/aerodynamics',
        'openaerostruct/transfer',
        'openaerostruct/functionals',
        'openaerostruct/integration',
        'openaerostruct/common',
        'openaerostruct/utils',
    ],
    # Test files
    package_data={
        'openaerostruct': ['tests/*.py', '*/tests/*.py', '*/*/tests/*.py']
    },
    # TODO: add versions?
    install_requires=[
        'openmdao[docs]>=3.2, <=3.6.0',
        'numpy',
        'scipy',
        'matplotlib',
    ],
    zip_safe=False,
    # ext_modules=ext,
    entry_points="""
    [console_scripts]
    plot_wing=openaerostruct.utils.plot_wing:disp_plot
    plot_wingbox=openaerostruct.utils.plot_wingbox:disp_plot
    """
)
