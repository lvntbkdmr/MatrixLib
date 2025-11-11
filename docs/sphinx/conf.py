# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Matrix Library'
copyright = '2024, Matrix Library'
author = 'Matrix Library Team'
release = '1.0.0'
version = '1.0.0'

# -- General configuration -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'breathe',  # Integrates Doxygen XML with Sphinx
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

# Breathe configuration
# This tells Breathe where to find the Doxygen XML output
breathe_projects = {
    'MatrixLib': '../doxygen_xml'
}
breathe_default_project = 'MatrixLib'

# Breathe domain-by-extension configuration
# This maps file extensions to C++ domain
breathe_domain_by_extension = {
    'h': 'cpp',
    'hpp': 'cpp',
    'cpp': 'cpp',
}

# Breathe default members
breathe_default_members = ('members', 'undoc-members')

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'  # Read the Docs theme (clean and modern)
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
    'logo_only': False,
    # 'display_version': True,  # Removed - not supported in this theme version
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
}

# HTML title
html_title = 'Matrix Library Documentation'
html_short_title = 'Matrix Library'

# -- Options for LaTeX output --------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-latex-output

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'figure_align': 'htbp',
}

latex_documents = [
    ('index', 'MatrixLibrary.tex', 'Matrix Library Documentation',
     'Matrix Library Team', 'manual'),
]

# -- Options for manual page output --------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-manual-page-output

man_pages = [
    ('index', 'matrixlib', 'Matrix Library Documentation',
     [author], 1)
]

# -- Options for Texinfo output ------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-texinfo-output

texinfo_documents = [
    ('index', 'MatrixLibrary', 'Matrix Library Documentation',
     author, 'MatrixLibrary', 'Avionics-safe matrix and vector library',
     'Miscellaneous'),
]

# -- Extension configuration ---------------------------------------------------

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'cpp': ('https://cppreference.com/', None),
}

