"""
Installation script for mapalign
"""

from setuptools import setup, find_packages

from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

with open(path.join(here, 'requirements.txt')) as f:
    requirements = [val.strip() for val in f.readlines()]

setup(
    name='mapalign',

    version='0.2.1',
    description='Mapalign',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/satra/mapalign',

    # Author details
    author='satra',

    # Choose your license
    license='Apache License, 2.0',
    provides=['mapalign'],
    packages=find_packages(),
    py_modules=["align", "dist", "embed"],
    
    install_requires=requirements[:-1],
    tests_require=requirements[-1:]
)
