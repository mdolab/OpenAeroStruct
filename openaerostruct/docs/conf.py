# -*- coding: utf-8 -*-
# This file is execfile()d with the current directory set to its
# containing dir.
import os
import sys
import importlib
from unittest.mock import Mock
from openaerostruct.docs._exts import embed_n2
from openaerostruct.docs._utils.generate_sourcedocs import generate_docs
from sphinx_mdolab_theme.config import *

sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.join("./_exts"))

# Only mock the ones that don't import.
MOCK_MODULES = ["h5py", "petsc4py", "pyoptsparse", "pyDOE2"]
for mod_name in MOCK_MODULES:
    try:
        importlib.import_module(mod_name)
    except ImportError:
        sys.modules[mod_name] = Mock()

# --- General configuration ---

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "numpydoc",
    "sphinx_copybutton",
    "embed_n2",
]

# directories for which to generate sourcedocs
packages = [
    "aerodynamics",
    "functionals",
    "geometry",
    "integration",
    "structures",
    "transfer",
]

generate_docs("..", "../..", packages, project_name="openaerostruct")

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "1.6.2"

numpydoc_show_class_members = False

# General information about the project.
project = "OpenAeroStruct"
copyright = "2018, John Jasa, Dr. John Hwang, Justin S. Gray"
author = "John Jasa, Dr. John Hwang, Justin S. Gray"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open("../__init__.py").read(),
)[0]

# The short X.Y version.
version = __version__

# The full version, including alpha/beta/rc tags.
release = __version__

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%b %d, %Y"

# Output file base name for HTML help builder.
htmlhelp_basename = "OpenAeroStructdoc"

html_extra_path = ["_n2html"]

# The master toctree document.
master_doc = "index"

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "openaerostruct", "OpenAeroStruct Documentation", [author], 1)]
