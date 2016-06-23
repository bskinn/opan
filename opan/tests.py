#-------------------------------------------------------------------------------
# Name:        tests
# Purpose:     Definitions of all unit tests for the Open Anharmonic package.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Jul 2015
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

# Module-level imports
import unittest, opan.test


# ==========================  Helper Functions  ============================= #

def any_params(d, k):
    """ Shortcut function for finding truthy values in a |dict|

    If the |dict| `d` contains truthy values at one or more of the keys
    in `k`, returns |True|. Otherwise, |False|.

    Parameters
    ----------
    d
        |dict| -- Dictionary of values to evaluate.

    k
        iterable -- Keys to search for truthy values.

    Returns
    -------
    any_found
        |bool| -- |True| if any truthy values were found at the indicated
        keys; |False| otherwise.

    Raises
    ------
    ~exceptions.KeyError
        If any of the values held in `k` is not a valid key in `d`.

    """

    return any(d[x] for x in k)

## end def any_params


if __name__ == '__main__':

    import os, sys, argparse as ap

    # Arguments for selecting test suites
    ALL = 'all'
    BASE = 'base'
    CONST = 'const'
    ERROR = 'error'
    SUPERS = 'supers'
    UTILS = 'utils'
    UTILS_BASE = 'utils_base'
    UTILS_DECORATE = 'utils_decorate'
    UTILS_INERTIA = 'utils_inertia'
    UTILS_VECTOR = 'utils_vector'
    XYZ = 'xyz'
    XYZ_FILEDATA = 'xyz_filedata'
    XYZ_DIRECTDATA = 'xyz_directdata'
    XYZ_USAGE = 'xyz_usage'

    ORCA = 'orca'
    ORCA_ENGRAD = 'orca_engrad'
    ORCA_HESS = 'orca_hess'

    # Prefix string for long arguments
    PFX = "--{0}"

    # Create the parser
    prs = ap.ArgumentParser(description="Run tests for Open Anharmonic.")

    # Create the top-level test groups
    gp_global = prs.add_argument_group(title="Global Test Options")
    gp_const = prs.add_argument_group(title="opan.const Tests")
    gp_error = prs.add_argument_group(title="opan.error Tests")
    gp_utils = prs.add_argument_group(title="opan.utils Tests")
    gp_xyz = prs.add_argument_group(title="opan.xyz Tests")

    gp_orca = prs.add_argument_group(title="ORCA Object Tests")

    # ====  OPTIONS  ==== #
    prs.add_argument('-v', action='store_true',
            help="Show verbose output")

    # ====  GLOBAL  ==== #
    # Add the all-tests argument
    gp_global.add_argument(PFX.format(ALL),
            action='store_true',
            help="Run all tests (overrides any other selections)")
    gp_global.add_argument(PFX.format(SUPERS),
            action='store_true',
            help="Run tests for superclasses")
    gp_global.add_argument(PFX.format(BASE),
            action='store_true',
            help="Run opan base-package tests")

    # ====  CONST  ==== #
    gp_const.add_argument(PFX.format(CONST),
            action='store_true', help="Run all opan.const tests")

    # ====  ERROR  ==== #
    gp_error.add_argument(PFX.format(ERROR),
            action='store_true', help="Run all opan.error tests")

    # ====  UTILS  ==== #
    # Add the arguments for the various suite cases
    gp_utils.add_argument(PFX.format(UTILS),
            action='store_true', help="Run all opan.utils tests")
    gp_utils.add_argument(PFX.format(UTILS_BASE),
            action='store_true', help="Run opan.utils.base tests")
    gp_utils.add_argument(PFX.format(UTILS_INERTIA),
            action='store_true', help="Run opan.utils.inertia tests")
    gp_utils.add_argument(PFX.format(UTILS_DECORATE),
            action='store_true', help="Run opan.utils.decorate tests")
    gp_utils.add_argument(PFX.format(UTILS_VECTOR),
            action='store_true', help="Run opan.utils.vector tests")

    # ====  XYZ  ==== #
    gp_xyz.add_argument(PFX.format(XYZ),
            action='store_true', help="Run all opan.xyz tests")
    gp_xyz.add_argument(PFX.format(XYZ_FILEDATA),
            action='store_true', help="Run opan.xyz tests on loading data from disk")
    gp_xyz.add_argument(PFX.format(XYZ_DIRECTDATA),
            action='store_true', help="Run opan.xyz tests on loading data dynamically")
    gp_xyz.add_argument(PFX.format(XYZ_USAGE),
            action='store_true', help="Run opan.xyz usage tests")

    # ====  ORCA Objects  ==== #
    gp_orca.add_argument(PFX.format(ORCA),
            action='store_true', help="Run all ORCA object tests")
    gp_orca.add_argument(PFX.format(ORCA_ENGRAD),
            action='store_true', help="Run OrcaEngrad tests")
    gp_orca.add_argument(PFX.format(ORCA_HESS),
            action='store_true', help="Run OrcaHess tests")

    # Pull the dictionary of the stored flags, with the unused args,
    #  and update sys.argv
    ns, args_left = prs.parse_known_args()
    params = vars(ns)
    sys.argv = sys.argv[:1] + args_left

    # Create the test suite to be compiled for running
    TestMasterSuite = unittest.TestSuite()

    # Base tests
    if any_params(params, [ALL, BASE]):
        TestMasterSuite.addTest(opan.test.opan_base.suite())

    # Superclasses
    if any_params(params, [ALL, SUPERS]):
        TestMasterSuite.addTest(opan.test.opan_supers.suite())

    # opan.const
    if any_params(params, [ALL, CONST]):
        TestMasterSuite.addTest(opan.test.opan_const.suite())

    # opan.error
    if any_params(params, [ALL, ERROR]):
        TestMasterSuite.addTest(opan.test.opan_error.suite())

    # opan.utils.base
    if any_params(params, [ALL, UTILS, UTILS_BASE]):
        TestMasterSuite.addTest(opan.test.opan_utils_base.suite())

    # opan.utils.inertia
    if any_params(params, [ALL, UTILS, UTILS_INERTIA]):
        TestMasterSuite.addTest(opan.test.opan_utils_inertia.suite())

    # opan.utils.decorate
    if any_params(params, [ALL, UTILS, UTILS_DECORATE]):
        TestMasterSuite.addTest(opan.test.opan_utils_decorate.suite())

    # opan.utils.vector
    if any_params(params, [ALL, UTILS, UTILS_VECTOR]):
        TestMasterSuite.addTest(opan.test.opan_utils_vector.suite())

    # opan.xyz (file data)
    if any_params(params, [ALL, XYZ, XYZ_FILEDATA]):
        TestMasterSuite.addTest(opan.test.opan_xyz.suite_FileData())

    # opan.xyz (direct data)
    if any_params(params, [ALL, XYZ, XYZ_DIRECTDATA]):
        TestMasterSuite.addTest(opan.test.opan_xyz.suite_DirectData())

    # opan.xyz (usage)
    if any_params(params, [ALL, XYZ, XYZ_USAGE]):
        TestMasterSuite.addTest(opan.test.opan_xyz.suite_Usage())


    # OrcaEngrad
    if any_params(params, [ALL, ORCA, ORCA_ENGRAD]):
        TestMasterSuite.addTest(opan.test.orca_engrad.suite())

    # OrcaHess
    if any_params(params, [ALL, ORCA, ORCA_HESS]):
        TestMasterSuite.addTest(opan.test.orca_hess.suite())


    # Create the text test runner and execute
    ttr = unittest.TextTestRunner(buffer=True,
                verbosity=(2 if params['v'] else 1))
    success = ttr.run(TestMasterSuite).wasSuccessful()

    # Store success/fail code and return it
    sys.exit(0 if success else 1)

## end main block

