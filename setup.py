#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'numpy',
    'Biopython',
]

# Below needed for Python 2.6
if sys.version_info < (2,7,):
    requirements.append('argparse')

test_requirements = [
    'testfixtures',
]

setup(
    name='mutator',
    version='0.1.0',
    description="Generate mutated sequence files from a reference genome.",
    long_description=readme + '\n\n' + history,
    author="Hugh Rand",
    author_email='Hugh.Rand@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/mutator',
    packages=[
        'mutator',
    ],
    package_dir={'mutator':
                 'mutator'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['bioinformatics', 'NGS', 'mutator'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    scripts=[
        'mutator/mutator.py',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
