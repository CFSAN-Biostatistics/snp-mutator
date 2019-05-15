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
if sys.version_info < (2, 7, ):
    requirements.append('argparse')

test_requirements = [
    'testfixtures',
]

setup(
    name='snp-mutator',
    version='1.2.0',
    description="Generate mutated sequence files from a reference genome.",
    long_description=readme + '\n\n' + history,
    author="Hugh Rand",
    author_email='Hugh.Rand@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/snp-mutator',
    packages=[
        'snpmutator',
    ],
    package_dir={'snpmutator':
                 'snpmutator'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['bioinformatics', 'NGS', 'snp-mutator'],
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
    entry_points={'console_scripts': ['snpmutator = snpmutator.script:main']},
    test_suite='tests',
    tests_require=test_requirements
)
