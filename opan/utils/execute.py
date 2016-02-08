#-------------------------------------------------------------------------------
# Name:        utils.execute
# Purpose:     Submodule containing utility functions for execution of
#               computational chemistry software packages for OpenAnharmonic
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

""" Functions to enable OpenAnharmonic to execute external software packages.

**Functions**

"""



# Module-level imports
from ..const import DEF as _DEF, EnumSoftware as _E_SW, EnumFileType as _E_FT


# Functions

def execute_orca(inp_tp, work_dir, exec_cmd, subs=None, subs_delims=("<",">"), \
            sim_name="orcarun", \
            inp_ext=_DEF.FILE_EXTS[_E_SW.ORCA][_E_FT.INPUTFILE], \
            out_ext=_DEF.FILE_EXTS[_E_SW.ORCA][_E_FT.OUTPUT], \
            wait_to_complete=True, \
            bohrs=False):
    """Executes |orca| on a dynamically constructed input file.

    .. warning:: Function is still under active development! Execution with
        `wait_to_complete` == |True| should be robust, however.

    **Execution**

    Generates an |orca| input file dynamically from information passed into the
    various arguments, performs the run, and returns with exit info and
    computation results *(in some fashion; still under development)*.
    Any required resources (GBW,
    XYZ, etc.) MUST already be present in `work_dir`. No check for pre-existing
    files of the same base name is made; any such will be overwritten.

    |orca| MUST be called using a wrapper script; this function does not
    implement the redirection necessary to send output from a direct |orca| call
    to a file on disk.

    If `wait_to_complete` is |True|, the :func:`subprocess.call` syntax will
    be used and the function will not return until execution of the
    wrapper script completes.
    If False, *[indicate what will be returned if not waiting]*.

    *#!DOC: execute_orca: The different output modes, depending on waiting
    or not.*

    The command to call |orca| must be specified in the parameter list syntax of
    the `args` argument to the :class:`subprocess.Popen` constructor.
    The implementation is flexible and general, to allow interface with local
    scripts for, e.g., submission to a job queue in a shared-resource
    environment.

    Valid |orca| input syntax of the resulting text is NOT checked
    before calling |orca|.

    No mechanism is implemented to detect hangs of |orca|. Periodic manual
    oversight is recommended.

    **Template Substitution**

    See :func:`utils.template_subst <opan.utils.base.template_subst>` for
    implementation details of the tag substitution mechanism.

    Here, in addition to performing any substitutions on the input file
    template as indicated by `subs`,
    the special tags **INP** and **OUT**, enclosed with the
    `subs_delims` delimiters, will be replaced with ``sim_name + '.' + inp_ext``
    and ``sim_name + '.' + out_ext``, respectively, in the input template and
    in all elements of
    `exec_cmd` before executing the call. In the special case of
    `inp_ext` == |None|, the **INP** tag will be replaced with just
    `sim_name` (no extension), and similarly for **OUT** if
    `out_ext` == |None|. The tag **NAME** will be replaced just with
    `sim_name` in all cases.

    `inp_ext` and `out_ext` must be different, to avoid collisions.

    **Return Values**

    The information returned depends on the value of `wait_to_complete`:

    If `wait_to_complete` == |True|:

        A `tuple` of objects is returned, with elements of type ::

            (OrcaOutput, OpanXYZ, OrcaEngrad, OrcaHess)

        These objects contain the corresponding results from the computation,
        if the latter exist, or |None| if they are missing.

    If `wait_to_complete` == |False|:

        **TBD**, but current intention is to return the PID of the spawned
        subprocess.


    **Signature**


    Parameters
    ----------
    inp_tp  : str
        Template text for the input file to be generated.

    work_dir : str
        Path to base working directory. Must already exist and contain any
        resource files (.gbw, .xyz, etc.) required for the calculation.

    exec_cmd : list of str
        Sequence of strings defining the |orca| execution call in the syntax of
        the :class:`~subprocess.Popen` constructor. This call must
        be to a local script; stream redirection of the forked process
        is not supported in this function.

    subs    : `dict` of |str|, optional
        Substitutions to be performed in the template (see *Template
        Substitution*, above).

    subs_delims : 2-tuple of str, optional
        Tag delimiters passed directly to
        :func:`~opan.utils.base.template_subst`.  Defaults to `('<','>')`.

    sim_name : str, optional
        Basename to use for the input/output/working files.
        If omitted, "orcarun" will be used.

    inp_ext : str
        Extension to be used for the input file generated.

    out_ext : str
        Extension to be used for the output file generated.


    Returns
    -------
    `tuple` of objects or |int|_ PID
        Varies depending on `wait_to_complete`; see *Return Values* above



    Raises
    ------
    ~exceptions.ValueError
        If `inp_ext` and `out_ext` are identical.

    ~exceptions.KeyError
        If special tag names **INP**, **OUT**, or **NAME** are defined in `subs`

    ~exceptions.TypeError
        If any elements in `subs` are not tuples

    """

    # Imports
    import os, subprocess as sp
    from ..output import OrcaOutput
    from ..xyz import OpanXYZ
    from ..grad import OrcaEngrad
    from ..hess import OrcaHess
    from ..utils import template_subst

    # Special key constants
    INPKEY = "INP"
    OUTKEY = "OUT"
    NAMEKEY = "NAME"

    # Store old dir; switch to new; default exception fine for
    #  handling case of invalid dir.
    olddir = os.getcwd()
    os.chdir(work_dir)

    # Check for inp_ext identical to out_ext
    if inp_ext == out_ext:
        raise(ValueError("'inp_ext' and 'out_ext' cannot be identical."))
    ##end if

    # Build the input and output file names and store the substitution keys
    if inp_ext:
        inp_fname = sim_name + '.' + inp_ext
    else:
        inp_fname = sim_name
    ##end if
    if out_ext:
        out_fname = sim_name + '.' + out_ext
    else:
        out_fname = sim_name
    ##end if

    SPECIALSUBS = {INPKEY: inp_fname, OUTKEY: out_fname, NAMEKEY: sim_name}

    # Complain if special tags used in subs
    if not set(SPECIALSUBS.keys()).isdisjoint(subs):
        raise(KeyError("Redefinition of special tag(s) is forbidden: {0}"
                .format(list(set(SPECIALSUBS.keys()).intersection(subs)))))
    ## end if

    # Perform the replacement into the exec_cmd of the input and output
    #  filenames.
    exec_cmd_subs = [
                s.replace(subs_delims[0] + INPKEY + subs_delims[1], inp_fname)
            for s in exec_cmd]
    exec_cmd_subs = [
                s.replace(subs_delims[0] + OUTKEY + subs_delims[1], out_fname)
            for s in exec_cmd_subs]
    exec_cmd_subs = [
                s.replace(subs_delims[0] + NAMEKEY + subs_delims[1], sim_name)
            for s in exec_cmd_subs]

    # Perform the content substitutions into the template string
    augsubs = subs.copy().update(SPECIALSUBS)
    input_text = template_subst(inp_tp, augsubs, delims=subs_delims)

    # Create and write the input file
    with open(inp_fname, 'w') as input_file:
        input_file.write(input_text)
    ##end with

    # Perform the ORCA call; collect return values as appropriate
    #!TODO: execute_orca: Implement non-waiting return
    if wait_to_complete:
        # Run ORCA
        sp.call(exec_cmd_subs, cwd=os.getcwd())

        # Bind ORCA_XXXXX objects and return. Have to address possibility of
        #  any or all of these not existing.
        try:
            o_out = OrcaOutput(out_fname,'file')
        except IOError:
            o_out = None
        ## end try

        try:
            o_xyz = OpanXYZ(path=os.path.join(work_dir, sim_name + ".xyz"), \
                                                                bohrs=bohrs)
        except IOError:
            o_xyz = None
        ## end try

        try:
            o_trj = OpanXYZ(path=os.path.join(work_dir, sim_name + ".trj"), \
                                                                bohrs=bohrs)
        except IOError:
            o_trj = None
        ## end try

        try:
            o_engrad = OrcaEngrad(os.path.join(work_dir, sim_name + ".engrad"))
        except IOError:
            o_engrad = None
        ## end try

        try:
            o_hess = OrcaHess(os.path.join(work_dir, sim_name + ".hess"))
        except IOError:
            o_hess = None
        ## end try

    else:
        raise(NotImplementedError("Background execution not yet implemented."))
    ## end if

    # Return to prior working directory
    os.chdir(olddir)

    # Return something appropriate, either computation results or information
    #  on the queued computation.
    #TODO: execute_orca: Must refine this, esp for the different exec modes
    return o_out, o_xyz, o_engrad, o_hess

## end def execute_orca


if __name__ == '__main__':      # pragma: no cover
    print("Module not executable.")
