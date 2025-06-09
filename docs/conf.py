# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from importlib.metadata import version

project = "GPSW"
copyright = "2025, Niek Wit"
author = "Niek Wit"
version = version("gpsw")

extensions = [
    "sphinx_design",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = ["force-light-mode.js"]

# -- Theme options -----------------------------------------------------------
html_theme_options = {
    "navbar_end": ["navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/niekwit/gps-orfeome",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
        {
            "name": "Bioconda",
            "url": "https://anaconda.org/bioconda/gpsw",
            "type": "local",
            "icon": "_static/anaconda.svg",
        },
        {
            "name": "Docker Hub",
            "url": "https://hub.docker.com/r/niekwit/gps-orfeome/tags",
            "icon": "fa-brands fa-docker",
            "type": "fontawesome",
        },
    ],
    "logo": {
        "text": "GPSW",
        "image_light": "_static/gpsw_logo.svg",
        "image_dark": "_static/gpsw_logo.svg",
    },
}

html_context = {"default_mode": "light"}

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
