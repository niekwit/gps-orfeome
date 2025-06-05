# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "GPSW"
copyright = "2025, Niek Wit"
author = "Niek Wit"
version = "0.7.0"

extensions = [
    "sphinx_design",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
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
}

html_context = {"default_mode": "light"}

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
