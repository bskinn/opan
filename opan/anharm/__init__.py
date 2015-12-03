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

"""**NEEDS DOCSTRING**, add submodules and base members to this docstring

Sub-Modules
-----------

:mod:`~opan.anharm.repo` -- HDF5 repository for :class:`OPAN_ANHARM`

"""

from __future__ import absolute_import

__all__ = ['repo']

from . import *
from .base import OpanAnharm




