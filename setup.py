#!/usr/bin/env python

from distutils.core import setup

setup(name='gandalf',
      version='0.4.0',
      package_dir={'gandalf': ''},
      packages=['gandalf','gandalf.analysis','gandalf.analysis.swig_generated'],
      package_data={'gandalf.analysis.swig_generated': ['_SphSim.so']}
      )