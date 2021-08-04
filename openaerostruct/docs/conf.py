# -*- coding: utf-8 -*-
# This file is execfile()d with the current directory set to its
# containing dir.
import os
import sys
import openmdao
import importlib
from unittest.mock import Mock
from openmdao.docs.config_params import MOCK_MODULES
from openmdao.docutils import do_monkeypatch, embed_code, embed_options, embed_compare
from openmdao.docs._exts import embed_n2
from openaerostruct.docs._utils.generate_sourcedocs import generate_docs

openmdao_path = os.path.split(os.path.abspath(openmdao.__file__))[0]
sys.path.insert(0, os.path.join(openmdao_path, 'docs', '_exts'))
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('.'))

# Only mock the ones that don't import.
for mod_name in MOCK_MODULES:
    try:
        importlib.import_module(mod_name)
    except ImportError:
        sys.modules[mod_name] = Mock()

do_monkeypatch()

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.coverage',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.githubpages',
              'numpydoc',
              'embed_code',
              'embed_compare',
              'embed_options',
              'embed_n2']

# directories for which to generate sourcedocs
packages = [
    'aerodynamics',
    'functionals',
    'geometry',
    'integration',
    'structures',
    'transfer',
]

generate_docs("..", "../..", packages, project_name='openaerostruct')

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.6.2'

numpydoc_show_class_members = False

# The master toctree document.
#master_doc = 'index'

# General information about the project.
project = u'OpenAeroStruct'
copyright = u'2018, John Jasa, Dr. John Hwang, Justin S. Gray'
author = u'John Jasa, Dr. John Hwang, Justin S. Gray'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('../__init__.py').read(),
)[0]
# The short X.Y version.
version = __version__
# The full version, including alpha/beta/rc tags.
release = __version__

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = '_theme'

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = ['.']
html_static_path = ['_static']
html_context = {
    'css_files': ['_static/style.css', ],
}

# # The name of an image file (relative to this directory) to place at the top
# # of the sidebar.
# html_logo = '_static/OpenMDAO_Logo.png'
#
# # The name of an image file (within the static path) to use as favicon of the
# # docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# # pixels large.
# html_favicon = '_static/OpenMDAO_Favicon.ico'

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# Output file base name for HTML help builder.
htmlhelp_basename = 'OpenAeroStructdoc'

#Customize sidebar
html_sidebars = {
   '**': ['globaltoc.html', 'searchbox.html']
}

html_extra_path = ['_n2html']

# The master toctree document.
master_doc = 'index'
# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'openaerostruct', u'OpenAeroStruct Documentation',
     [author], 1)
]
