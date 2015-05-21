#! /usr/bin/python
"""Setuptools-based setup script for MDGenesis.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup

setup(name='mdgenesis',
      version='0.1.0-dev',
      author='Chris Ing', 
      author_email='ing.chris@gmail.com',
      packages=['mdgenesis','mdgenesis.modules'],
      license='GPL 2',
      long_description=open('README.rst').read(),
      requires=['MDSynthesis', 'MDAnalysis']
     )
