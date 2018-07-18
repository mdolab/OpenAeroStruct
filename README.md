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

Please see the [documentation](https://mdolab.github.io/OpenAeroStruct/) for more installation details, walkthroughs, and examples.
