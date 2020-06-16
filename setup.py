#!/usr/bin/env python3
"""
ScopeSim: A python package to simulate telescope observations

How to compile and put these on pip::

    $ python setup.py sdist
    $ twine upload dist/*

"""
from setuptools import setup, find_packages


with open('README.md') as f:
    __readme__ = f.read()

with open('LICENCE') as f:
    __license__ = f.read()

with open('scopesim/version.py') as f:
    __version__ = f.readline().split("'")[1]


def setup_package():
    setup(name='ScopeSim',
          version=__version__,
          description="Generalised telescope observation simulator",
          long_description=__readme__,
          long_description_content_type='text/markdown',
          author="Kieran Leschinski",
          author_email="kieran.leschinski@unive.ac.at",
          url="https://github.com/astronomyk/ScopeSim",
          license="GNU General Public License",
          package_dir={'scopesim': 'scopesim'},
          include_package_data=True,
          packages=find_packages(exclude=('tests', 'data', 'docs_to_be_sorted',
                                          'misc', 'OLD_code', )),
          install_requires=["numpy>=1.16",
                            "scipy>=1.0.0",
                            "astropy>=2.0",
                            "matplotlib>=1.5",
                            "pyyaml>5.1",
                            "requests>=2.20",
                            "beautifulsoup4>=4.4",
                            "synphot>0.1.3",
                            "skycalc_ipy",
                            "anisocado",
                            ],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: MIT License",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
