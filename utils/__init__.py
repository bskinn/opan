#-------------------------------------------------------------------------------
# Name:        __init__
# Purpose:     opan.utils.__init__ file
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     5 Oct 2015
# Copyright:   (c) Brian Skinn 2015
# License:     The MIT License; see "license.txt" for full license terms
#                   and contributor agreement.
#
#       This file is part of opan (OpenAnharmonic), a system for automated
#       computation of anharmonic properties of molecular systems via wrapper
#       calls to computational/quantum chemical software packages.
#
#       http://www.github.com/bskinn/opan
#
#-------------------------------------------------------------------------------

"""Utility functions for OpenAnharmonic, including execution automation.


Functions
---------
check_geom       -- Confirm two OpenBabel geometries (atom types and
                        coordinates) match to within a specified tolerance
delta_fxn        -- Generalized Kronecker delta function
make_timestamp   -- Construct a string time-elapsed timestamp in h/m/s format
pack_tups        -- Pack an arbitrary combination of iterables and non-
                        iterables into a list of tuples
safe_cast        -- Robustified casting with a post-check to confirm the cast
                        actually resulted in the proper type
template_subst   -- Perform a field-based substitution into a string template

Sub-Modules
-----------
execute          -- Functions for execution of external computational
                        software packages
inertia          -- Inertia-related tools (c.o.m., rotational constants, etc.)
symm             -- Molecular symmetry utility functions
vector           -- Vector utility functions

"""

from __future__ import absolute_import

__all__ = ['decorate', 'execute', 'inertia', 'symm', 'vector']

from . import *
from .base import check_geom, delta_fxn, make_timestamp, pack_tups
from .base import safe_cast, template_subst



