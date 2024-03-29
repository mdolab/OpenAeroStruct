.. _How_to_Contribute:

How to Contribute to OpenAeroStruct
===================================

.. note::
  This section will be expanded in the future.

OpenAeroStruct is an open-source tool, thus we welcome users to submit additions or fixes to the code to make it better for everybody.

Issues
------
If you have an issue with OpenAeroStruct, a bug to report, or a feature to request, submit an issue on the GitHub repository.
This lets other users know about the issue.
If you are comfortable fixing the issue, please do so and submit a pull request.

Documentation
-------------
When you add or modify code, make sure to provide relevant documentation that explains the new code.
This should be done in code via comments, but also in the Sphinx documentation as well if you add a new feature or capability.
Look at the .rst files in the `docs` section of the repo or click on `view source` on any of the doc pages to see some examples.


Building Docs Locally
---------------------
To build the OpenAeroStruct documentation locally, first install our sphinx theme by ``pip install sphinx_mdolab_theme``.
On Linux or Mac, use the Makefile in the docs folder: ``cd PATH-TO-openaerostruct/docs && make html``.
It is possible to make the documents on Windows without the make utility.
Navigate to the docs folder and invoke :code:`sphinx-build -b html .\ .\_build`

Testing
-------
When you add code or functionality, add tests that cover the new or modified code.
These may be units tests for individual components or regression tests for entire models that use the new functionality.

Each discipline sub-directory has unit tests in the `tests` folder.
For example, `openaerostruct/structures/tests` contains the unit tests for the structural components.
Look at `test_materials_tube.py` within that folder for a simple unit test that you can mimic when you add new components.

Regression tests live in the base `openaerostruct/tests` folder.
If you introduce a new capability or feature, add regression tests to cover this functionality.
These are generally more expensive than a unit test.
They might perform analysis or optimization on a constructed model with certain options enabled, for example.

Pull requests
-------------
Once you have added or modified code, submit a pull request via the GitHub interface.
This will automatically go through all of the tests in the repo to make sure everything is functioning properly.
This also automatically does a coverage test to ensure that any added code is covered in a test.
The main developers of OpenAeroStruct will then merge in the request or provide feedback on how to improve the contribution.

Any code in OpenAeroStruct should adhere to the PEP-8 style.
Upon PRs, we run ``flake8`` linter and ``black`` formatter to check the code style.
To run ``flake8`` and ``black`` locally,

.. code-block:: bash

  $ cd PATH-TO-OpenAeroStruct
  $ pip install flake8==3.9.2
  $ wget https://raw.githubusercontent.com/mdolab/.github/main/.flake8 -O .flake8_mdolab  # download flake8 configuration for OAS
  $ python -m flake8 openaerostruct --append-config .flake8_mdolab --append-config .github/.flake8_oas_specific --count --show-source --statistics
  $ pip install black==22.3.0
  $ black openaerostruct -l 120 --target-version py38
