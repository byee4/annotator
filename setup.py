#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='annotator',
    version='0.0.5',
    url='',
    license='',
    author='brianyee',
    author_email='',
    description='Annotate BED files with genic information given a priority list',
    packages=['annotator'],
    package_dir={
        'annotator': 'annotator',
    },
    entry_points = {
        'console_scripts': [
            'annotate-bed = annotator.annotate:main',
        ]
    }
)
