#!/usr/bin/env python

from distutils.core import setup, Extension


FactorGraph_module = Extension('_FactorGraph',
                           sources=['FactorGraph.cpp', 'FactorCreator.cpp', 'Factor.cpp', 'FactorGraph_wrap.cxx'],
                           extra_compile_args=['-std=c++11', '-I/usr/local/include','-DPERFECT_MATCHING_DOUBLE'],
                           )

setup (name = 'FactorGraph',
       version = '0.1',
       author      = "Zhen Zhang",
       description = """Python module for Factor Graph Inference""",
       ext_modules = [FactorGraph_module],
       py_modules = ["FactorGraph"],
       )
