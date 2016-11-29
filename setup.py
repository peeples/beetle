#!/usr/bin/env python

from setuptools import setup

setup(
    name='beetle',
    author='Molly Peeples',
    author_email='molly@stsci.edu',
    version='0.1',
    py_modules=['beetle',
                'calculate_rates',
                'create_galaxies',
                'constants',
                'make_metals',
                'make_plots',
                'read_in_galaxies',
                'peeples2014'],  # Python modules to install (without the .py in the filename)
    scripts=['beetle']  # This is the full name of the script "chester"; this will be installed to a bin/ directory
)
