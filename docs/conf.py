# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------
import sphinx_rtd_theme

project = 'distopia'
copyright = '2021, Hugo MacDermott-Opeskin Jakub Nabaglo Richard Gowers'
author = 'Hugo MacDermott-Opeskin Jakub Nabaglo Richard Gowers'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['breathe', 'sphinx.ext.autosectionlabel',
              'sphinx_rtd_theme', 'sphinx.ext.napoleon',
              'sphinx.ext.githubpages',
              'sphinx.ext.autosummary',
              'sphinx_sitemap']
breathe_projects = { "distopia": "./doxygen_build/xml" }
breathe_default_project = "distopia"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

site_url = "https://www.mdanalysis.org/distopia/"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'default'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_theme_path = [
    sphinx_rtd_theme.get_html_theme_path()
]

html_theme_options = {
    'canonical_url': '',
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# options common to RTD and MDAnalysis theme

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# For RTD theme: custom.css to override theme defaults.
html_static_path = ['_static']
html_css_files = ['custom.css']
html_logo = "_static/logos/mdanalysis-logo-thin.png"