"""Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

import os
import sys

# import versioningit

sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Tavi"
copyright = "2024, ORNL"  # noqa A001
author = "ORNL"
# version = versioningit.get_version("../../")
# # The full version, including alpha/beta/rc tags
# release = ".".join(version.split(".")[:-1])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.coverage",
    "sphinxcontrib.mermaid",
]

master_doc = "index"
source_suffix = [".rst", ".md"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]
exclude_patterns = ["_build"]
pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"  # pylint: disable=C0103

html_theme_options = {"style_nav_header_background": "#472375"}

html_static_path = ["_static"]

html_css_files = [
    "css/custom.css",
]

epub_show_urls = "footnote"  # pylint: disable=C0103
