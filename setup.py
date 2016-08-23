#!/usr/bin/env python

from distutils.core import setup, Extension


FactorGraph_module = Extension('_FactorGraph',
                           sources=['FactorGraph.cpp', 'FactorCreator.cpp', 'Factor.cpp', 'FactorGraph_wrap.cxx', 'Auction.cpp'],
                           extra_compile_args=['-O3', '-I/usr/local/include','-DPERFECT_MATCHING_DOUBLE', '-ffast-math', '-pipe', '-fomit-frame-pointer', '-std=c++11'],
                           #extra_compile_args=['-g', '-I/usr/local/include','-DPERFECT_MATCHING_DOUBLE', '-std=c++11'],
                           )

setup (name = 'FactorGraph',
       version = '0.1',
       author      = "Zhen Zhang",
       description = """Python module for Factor Graph Inference""",
       ext_modules = [FactorGraph_module],
       py_modules = ["FactorGraph"],
       )
