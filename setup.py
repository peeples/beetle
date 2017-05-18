from setuptools import setup

setup(
    name='beetle',
    author='Molly Peeples',
    author_email='molly@stsci.edu',
    url="https://github.com/peeples/beetle",
    version='0.1',
    py_modules=['beetle',
                'calculate_rates',
                'create_galaxies',
                'constants',
                'make_metals',
                'make_plots',
                'read_in_galaxies',
                'peeples2014'],  # Python modules to install (without the .py in the filename)
    scripts=['beetle.py']  # This is the full name of the script "beetle"; this will be installed to a bin/ directory
)
