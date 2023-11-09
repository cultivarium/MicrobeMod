#!/usr/bin/env python

from setuptools import setup, find_packages

from MicrobeMod._version import __version__

setup(name='MicrobeMod',
      version=__version__,
      description='Calculation of strain-level metrics',
      url='https://github.com/cultivarium/MicrobeMod',
      author='Alex Crits-Christoph',
      author_email='alex@cultivarium.org',
      license='MIT License',
      package_data={'MicrobeMod': ['./db/restriction_metadata.csv','db/HMMs/*','db/rebase_blast/*']},
      packages=find_packages(exclude=["tests"]),
      scripts=['bin/MicrobeMod'],
      python_requires='>=3.4.0',
      install_requires=[
          'pandas>=1.5.0',
          'biopython>=1.81'
      ],
      zip_safe=False)