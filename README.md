OpenAeroStruct
==============

[![Build Status](https://dev.azure.com/mdolab/Public/_apis/build/status%2Fmdolab.OpenAeroStruct?repoName=mdolab%2FOpenAeroStruct&branchName=main)](https://dev.azure.com/mdolab/Public/_build/latest?definitionId=49&repoName=mdolab%2FOpenAeroStruct&branchName=main)
[![codecov](https://codecov.io/gh/mdolab/OpenAeroStruct/branch/main/graph/badge.svg?token=yOxeH7rT2H)](https://codecov.io/gh/mdolab/OpenAeroStruct)
[![Documentation Status](https://readthedocs.com/projects/mdolab-openaerostruct/badge/?version=latest)](https://mdolab-openaerostruct.readthedocs-hosted.com/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/openaerostruct)](https://pypi.org/project/openaerostruct/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/openaerostruct)](https://pypi.org/project/openaerostruct/)

OpenAeroStruct is a lightweight tool that performs aerostructural optimization using OpenMDAO.
It couples a vortex-lattice method (VLM) and a 6 degrees of freedom 3-dimensional spatial beam model to simulate aerodynamic and structural analyses using lifting surfaces.
These simulations are wrapped with an optimizer using NASA's OpenMDAO framework.
The analysis and optimization results can be visualized using included tools, producing figures such as these:

*With a tubular structure*
![Example](https://raw.githubusercontent.com/mdolab/OpenAeroStruct/main/openaerostruct/docs/figures/example.png)

*With a wingbox structure*
![Example2](https://raw.githubusercontent.com/mdolab/OpenAeroStruct/main/openaerostruct/docs/figures/wingbox_fine.png)

Please note that this repository is provided as is without any guaranteed support.
If you would like to highlight issues, ask questions, or make changes, please do so using GitHub Issues and Pull Requests.
The developers will address them at their discretion.

The easiest way to get started is by installing OpenAeroStruct via the Python Package Index with

`pip install openaerostruct`

If you'd like easier access to the examples and source code, install OpenAeroStruct by cloning this repository and entering the folder it generates.
Then do:

`pip install -e .`

Documentation
-------------

Please see the [documentation](https://mdolab-openaerostruct.readthedocs-hosted.com/en/latest/) for more installation details, walkthroughs, and examples.

Citation
--------

For more background, theory, and figures, see the [OpenAeroStruct journal article](https://mdolab.engin.umich.edu/bibliography/Jasa2018a.html).
Please cite this article when using OpenAeroStruct in your research or curricula.

John P. Jasa, John T. Hwang, and Joaquim R. R. A. Martins. "Open-source coupled aerostructural optimization using Python." Structural and Multidisciplinary Optimization 57.4 (2018): 1815-1827. DOI: 10.1007/s00158-018-1912-8

```
@article{Jasa2018a,
	Author = {John P. Jasa and John T. Hwang and Joaquim R. R. A. Martins},
	Doi = {10.1007/s00158-018-1912-8},
	Journal = {Structural and Multidisciplinary Optimization},
	Month = {April},
	Number = {4},
	Pages = {1815--1827},
	Publisher = {Springer},
	Title = {Open-source coupled aerostructural optimization using {Python}},
	Volume = {57},
	Year = {2018}}
```

If using the wingbox model, fuel-weight inertial loads, or structural-weight inertial loads, please cite the following [conference paper](https://www.researchgate.net/publication/327654423_Low-Fidelity_Aerostructural_Optimization_of_Aircraft_Wings_with_a_Simplified_Wingbox_Model_Using_OpenAeroStruct).

Shamsheer S. Chauhan and Joaquim R. R. A. Martins, “Low-Fidelity Aerostructural Optimization of Aircraft Wings with a Simplified Wingbox Model Using OpenAeroStruct,” Proceedings of the 6th International Conference on Engineering Optimization, EngOpt 2018, Springer, Lisbon, Portugal, September 2018, pp. 418–431. doi:10.1007/978-3-319-97773-7 38

```
@inproceedings{Chauhan2018b,
	Author = {Shamsheer S. Chauhan and Joaquim R. R. A. Martins},
	Address = {Lisbon, Portugal},
	Booktitle = {Proceedings of the 6th International Conference on Engineering Optimization, EngOpt 2018},
	Doi = {10.1007/978-3-319-97773-7_38},
	Pages = {418-431},
	Publisher = {Springer},
	Title = {Low-Fidelity Aerostructural Optimization of Aircraft Wings with a Simplified Wingbox Model Using {OpenAeroStruct}},
	Year = {2018}}
```

If using point-mass loads or thrust loads, please cite the following [conference paper](https://www.researchgate.net/publication/333806174_How_Certain_Physical_Considerations_Impact_Aerostructural_Wing_Optimization).

John P. Jasa, Shamsheer S. Chauhan, Justin S. Gray, and Joaquim R. R. A. Martins, “How Certain Physical Considerations Impact Aerostructural Wing Optimization,” AIAA/ISSMO Multidisciplinary Analysis and Optimization Conference, Dallas, TX, 2019. doi:10.2514/6.2019-3242

```
@inproceedings{Jasa2019c,
	Author = {John P. Jasa and Shamsheer S. Chauhan and Justin S. Gray and Joaquim R. R. A. Martins},
	Address = {Dallas, TX},
	Booktitle = {AIAA/ISSMO Multidisciplinary Analysis and Optimization Conference},
	Doi = {10.2514/6.2019-3242},
	Title = {How Certain Physical Considerations Impact Aerostructural Wing Optimization},
	Month = {June},
	Year = {2019}}
```

Version Information
-------------------
The oldest and latest versions of the dependencies that we test regularly are the following (other versions may work, but no guarantees):

| Dependency         | oldest | latest |
|--------------------|--------| ------ |
| Python             | 3.8    | 3.11   |
| NumPy              | 1.20   | latest |
| SciPy              | 1.6.0  | latest |
| OpenMDAO           | 3.35   | latest |
| Matplotlib         | latest | latest |
| MPhys (optional)   | 2.0.0  | latest |
| pyGeo (optional)   | 1.6.0  | latest |
| OpenVSP (optional) | 3.27.1 | 3.27.1 |

If you are looking to use the previous version of OpenAeroStruct which uses OpenMDAO 1.7.4, use OpenAeroStruct 1.0 from [here](https://github.com/mdolab/OpenAeroStruct/releases).

License
-------
Copyright 2018-2023 MDO Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
