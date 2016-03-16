#-------------------------------------------------------------------------------
# Name:        __init__
# Purpose:     opan.vpt2.__init__ file
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

"""Submodule implementing VPT2 anharmonic computations.

.. warning::

    This module is under active development. API &c. may change
    with little notice.

**Sub-Modules**

:mod:`~opan.vpt2.repo` -- HDF5 repository for :class:`OpanVPT2`


**Classes**

:class:`~opan.vpt2.base.OpanVPT2` -- Core driver class for VPT2
anharmonic calculations.


**API**

.. autoclass:: opan.vpt2.base.OpanVPT2


"""

from __future__ import absolute_import

__all__ = ['repo']

from . import *
from .base import OpanVPT2




