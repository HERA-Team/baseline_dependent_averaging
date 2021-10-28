# -*- mode: python; coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License
"""Package for applying baseline-dependent averaging to radio astronomy datasets."""

from setuptools import setup
import glob
import io

with io.open("README.md", "r", encoding="utf-8") as readme_file:
    readme = readme_file.read()

setup_args = {
    "name": "bda",
    "author": "Paul La Plante",
    "url": "https://github.com/HERA-Team/baseline_dependent_averaging",
    "license": "BSD",
    "description": (
        "a tool for applying baseline-dependent averaging to a radio "
        "interferometer dataset"
    ),
    "long_description": readme,
    "long_description_content_type": "text/markdown",
    "package_dir": {"bda": "bda"},
    "packages": ["bda"],
    "scripts": glob.glob("scripts/*"),
    "include_package_data": True,
    "use_scm_version": True,
    "install_requires": ["pyuvdata>=2.2.0", "astropy>=3.0", "setuptools_scm"],
    "extras_require": {"testing": ["pytest>=6.0", "pytest-cov", "pre-commit"]},
    "classifiers": [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    "keywords": "baseline dependent averaging",
}

if __name__ == "__main__":
    setup(**setup_args)
