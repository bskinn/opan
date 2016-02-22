#-------------------------------------------------------------------------------
# Name:        __init__
# Purpose:     opan.utils.__init__ file
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     5 Oct 2015
# Copyright:   (c) Brian Skinn 2016
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


**Sub-Modules**

:mod:`~opan.utils.decorate` -- Custom OpenAnharmonic decorators

:mod:`~opan.utils.execute` -- Functions for execution of external computational
software packages

:mod:`~opan.utils.inertia` -- Inertia-related tools (center of mass,
rotational constants, principal moments/axes, etc.)

:mod:`~opan.utils.symm` -- Molecular symmetry utility functions
(\ **INCOMPLETE**\ )

:mod:`~opan.utils.vector` -- Vector utility functions

**Functions**

.. autofunction:: opan.utils.base.assert_npfloatarray

.. autofunction:: opan.utils.base.check_geom(c1,a1,c2,a2[,tol])

.. autofunction:: opan.utils.base.delta_fxn(a,b)

.. autofunction:: opan.utils.base.iterable

.. autofunction:: opan.utils.base.make_timestamp

.. autofunction:: opan.utils.base.pack_tups

.. autofunction:: opan.utils.base.safe_cast

.. autofunction::
    opan.utils.base.template_subst(template, subs[, delims=['<', '>']])


"""

from __future__ import absolute_import

__all__ = ['decorate', 'execute', 'inertia', 'symm', 'vector']

from . import *
from .base import check_geom, delta_fxn, make_timestamp, pack_tups
from .base import safe_cast, template_subst, iterable
from .base import assert_npfloatarray



