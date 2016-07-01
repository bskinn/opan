#-------------------------------------------------------------------------------
# Name:        utils.freq
# Purpose:     Submodule containing utility functions dealing with calculations
#               of inertia tensor-related properties
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Jun 2016
# Copyright:   (c) Brian Skinn 2016
# License:     The MIT License; see "license.txt" for full license terms
#                   and contributor agreement.
#
#       This file is part of opan (Open Anharmonic), a system for automated
#       computation of anharmonic properties of molecular systems via wrapper
#       calls to computational/quantum chemical software packages.
#
#       http://www.github.com/bskinn/opan
#
#-------------------------------------------------------------------------------

"""Utilities for calculation of harmonic-oscillator vibrational frequencies,
and quantities related to / derived from them.

These functions are housed separately from the :mod:`opan.vpt2` VPT2 module
since they may have broader applicability to other envisioned capabilities of
Open Anharmonic, and because they will likely be independently valuable
to users.

**Functions**

.. autofunction:: mw_hess(hess, masses)

"""

from opan.utils.decorate import arraysqueeze as _arraysqueeze

@_arraysqueeze(1)
def mw_hess(hess, masses):
    pass





if __name__ == '__main__':  # pragma: no cover
    print("Module not executable.")

