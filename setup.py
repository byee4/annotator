#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='annotator',
    version='0.0.2',
    url='',
    license='',
    author='brianyee',
    author_email='',
    description='Annotate stuff',
    packages=['annotator'],
    package_dir={
        'annotator': 'annotator',
    },
    entry_points = {
        'console_scripts': [
            'anno = annotator.annotate:main',
        ]
    }
)