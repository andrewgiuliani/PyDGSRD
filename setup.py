from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import setuptools

__version__ = '0.0.1'

setup(
    name='PyDGSRD1D',
    long_description='',
    install_requires=['numpy', 'scipy', 'argparse', 'matplotlib', 'pytest'],
    packages = ["pydgsrd1d"],
    package_dir = {"pydgsrd1d": "pydgsrd1d"},
    zip_safe=False,
)
