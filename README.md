OpenAeroStruct
==============

[![Build Status](https://travis-ci.org/mdolab/OpenAeroStruct.svg?branch=master)](https://travis-ci.org/mdolab/OpenAeroStruct)
[![Coverage Status](https://coveralls.io/repos/github/mdolab/OpenAeroStruct/badge.svg?branch=master)](https://coveralls.io/github/mdolab/OpenAeroStruct?branch=master)

OpenAeroStruct is a lightweight tool that performs aerostructural optimization using OpenMDAO.
It couples a vortex-lattice method (VLM) and a 6 degrees of freedom 3-dimensional spatial beam model to simulate aerodynamic and structural analyses using lifting surfaces.
These simulations are wrapped with an optimizer using NASA's OpenMDAO framework.
The analysis and optimization results can be visualized using included tools, producing figures such as this:

![Example](openaerostruct/docs/example.png)

Install OpenAeroStruct by cloning this repository and entering the folder it generates.
Then do:

`pip install -e .`

Version Information
===================
This version of OpenAeroStruct requires [OpenMDAO](https://github.com/OpenMDAO/openmdao) v2.3+.
If you are looking to use the previous version of OpenAeroStruct which uses OpenMDAO 1.7.4, use OpenAeroStruct v1.0 from [here](https://github.com/mdolab/OpenAeroStruct/releases).

Documentation
=============

Please see the [documentation](https://mdolab.github.io/OpenAeroStruct/) for more installation details, walkthroughs, and examples.

Citing OpenAeroStruct
=====================

For more background, theory, and figures, see the [OpenAeroStruct journal article](http://mdolab.engin.umich.edu/sites/default/files/OAS_SMO_preprint_0.pdf).
Please cite this article when using OpenAeroStruct in your research or curricula.

License
=======
Copyright 2018 MDO Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
