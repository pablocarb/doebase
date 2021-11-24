# @package sbol
# A Python wrapper for libSBOLc, a module for reading, writing, and constructing
# genetic designs according to the standardized specifications of the Synthetic Biology Open Language
#
# @file sbol.py
# Implements a high-level, Pythonic interface for the SWIG-Python classes in libsbol
#
# @namespace sbol.sbol
# High level wrappers for libSBOLc
#
# @namespace sbol.libsbol
# Low level SWIG-Python wrappers for libSBOLc
#
# @namespace sbol.sbol_test
# Unit tests

from platform import system
if system() == 'Darwin':
    platform_system = 'macOS'
elif system() == 'Windows':
    platform_system = 'Windows'
elif system() == 'Linux':
    platform_system = 'Linux'


from sys import path as sys_path
from os import path as os_path
here = os_path.abspath(os_path.dirname(__file__))
# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(
    1,
    os_path.join(
        here,
        'lib',
        platform_system
    )
)


from libsbol import *
