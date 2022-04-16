#!/usr/bin/env python
# -*- encoding: utf8 -*-
import glob
import inspect
import io
import os

from setuptools import find_packages
from setuptools import setup


long_description = """
Source code: https://github.com/chenyk1990/pyortho""".strip() 

def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")).read()

setup(
    name="pyortho",
    version="0.0.3",
    license='GNU General Public License, Version 3 (GPLv3)',
    description="A python package for local signal-and-noise orthogonalization and local similarity calculation",
    long_description=long_description,
    author="pyortho developing team",
    author_email="chenyk2016@gmail.com",
    url="https://github.com/chenyk1990/pyortho",
    packages=['pyortho'],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
    ],
    keywords=[
        "seismology", "exploration seismology", "array seismology", "denoising", "science", "orthogonalization", "local orthogonalization", "local similarity", "signal-to-noise ratio", "damped rank reduction method"
    ],
    install_requires=[
        "numpy", "scipy", "matplotlib"
    ],
    extras_require={
        "docs": ["sphinx", "ipython", "runipy"]
    }
)
