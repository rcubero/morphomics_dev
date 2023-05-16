""" Distribution configuration for Moprhomics
"""
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


config = {
    "description": "morphOMICs: a python package for the topological and statistical  analysis of microglia morphology",
    "author": "Ryan Cubero",
    "url": "https://github.com/rcubero/MorphOMICs",
    "author_email": "ryanjohn.cubero@ist.ac.at",
    "long_description": long_description,
    "long_description_content_type": "text/markdown",
    "install_requires": [
        "TMD>=2.3.0",
        "scipy>=0.13.3",
        "numpy>=1.8.0",
        "scikit-learn>=0.19.1",
        "matplotlib>=3.2.0",
        "umap-learn>=0.3.10",
        "tomli>=2.0.1",
        "ipyvolume>=0.6.1",
        "ipython_genutils>=0.2.0",
    ],
    "setup_requires": ["setuptools_scm"],
    "packages": find_packages(),
    "license": "MIT",
    "scripts": [],
    "name": "morphomics",
    "include_package_data": True,
    "use_scm_version": True,
    "python_requires": ">=3.7",
}


setup(**config)
