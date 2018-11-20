# -*- mode: python; coding: utf-8 -*-
# Copyright (c) 2018 Paul La Plante
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import

from setuptools import setup
import glob
import os
import io
import json

from bda import version

data = [version.git_origin, version.git_hash, version.git_description, version.git_branch]
with open(os.path.join('bda', 'GIT_INFO'), 'w') as outfile:
    json.dump(data, outfile)

with io.open('README.md', 'r', encoding='utf-8') as readme_file:
    readme = readme_file.read()

setup_args = {
    'name': 'bda',
    'author': 'Paul La Plante',
    'url': 'https://github.com/plaplant/baseline_dependent_averaging',
    'license': 'BSD',
    'description': 'a tool for applying baseline-dependent averaging to a radio interferometer dataset',
    'long_description': readme,
    'long_description_content_type': 'text/markdown',
    'package_dir': {'bda': 'bda'},
    'packages': ['bda'],
    'scripts': glob.glob('scripts/*'),
    'version': version.version,
    'include_package_data': True,
    'setup_requires': ['pyuvdata>=1.3.3', 'astropy>=2.0', 'six>=1.10'],
    'install_requires': ['pyuvdata>=1.3.3', 'astropy>=2.0', 'six>=1.10'],
    'classifiers': ['Development Status :: 3 - Alpha',
                    'Intended Audience :: Science/Research',
                    'License :: OSI Approved :: BSD License',
                    'Programming Language :: Python :: 2.7',
                    'Programming Language :: Python :: 3.6',
                    'Topic :: Scientific/Engineering :: Astronomy'],
    'keywords': 'baseline dependent averaging'
}

if __name__ == '__main__':
    setup(**setup_args)
