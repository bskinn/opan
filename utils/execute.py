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


# Module-level imports
from ..const import DEF, E_Software as E_SW, E_FileType as E_FT
from ..const import E_DispDirection


# Functions

def execute_orca(inp_tp, work_dir, exec_cmd, subs=None, subs_delims=("<",">"), \
            sim_name="orcarun", \
            inp_ext=DEF.File_Extensions[E_SW.ORCA][E_FT.inputfile], \
            out_ext=DEF.File_Extensions[E_SW.ORCA][E_FT.output], \
            wait_to_complete=True, \
            bohrs=False):
    """Executes ORCA on a dynamically constructed input file.

    Generates an ORCA input file dynamically from information passed into the
    various arguments, performs the run, and returns with exit info and
    computation results (in some fashion...).  Any required resources (GBW,
    XYZ, etc.) MUST already be present in 'work_dir'. No check for pre-existing
    files of the same base name is made; any such will be overwritten.

    ORCA MUST be called using a wrapper script; this method does not
    implement the redirection necessary to send output from a direct ORCA call
    to a file on disk.

    If 'wait_to_complete' is True, the subprocess.call() syntax will be used
    and the function will not return until ORCA execution completes. If False,
    [indicate what will be returned if not waiting.]

    #!DOC: execute_orca: The different output modes, depending on waiting or not.

    The command to call ORCA must be specified in the parameter list syntax of
    the subprocess.Popen constructor (https://docs.python.org/2/library/
    subprocess.html#popen-constructor).  The execution is flexible and general,
    to allow interface with local scripts for, e.g., submission to a job queue
    in a shared-resource environment.  The special tags "INP" and "OUT",
    enclosed with the 'subs_delims' delimiters, will be replaced with
        sim_name+'.'+inp_ext
    and
        sim_name+'.'+out_ext
    respectively in all elements of 'exec_cmd' before executing the call.
    In the special cases of inp_ext=None, the INP tag will be replaced with
    just sim_name (no extension), and similarly if out_ext=None. Similarly, the
    tag NAME will be replaced just with 'sim_name' in all cases.

    inp_ext and out_ext must be different, to avoid collisions.

    Valid ORCA input syntax of the resulting text is NOT checked before calling
    ORCA.

    No mechanism is implemented to detect hangs of ORCA. Periodic manual
    oversight is recommended.

    Parameters
    ----------
    inp_tp  : str
        Template text for the input file to be generated.
    work_dir : str
        Path to base working directory. Must already exist and contain any
        resource files (GBW, XYZ, etc.) required for the calculation.
    exec_cmd : list of str
        Sequence of strings defining the ORCA execution call in the syntax of
        the subprocess.Popen constructor. This call must be to a local script;
        stream redirection of the forked process is not supported in this
        method.
    subs    : iterable of (str, str) tuples/lists, optional
        Substitutions to be performed in template. The first element of each
        tuple corresponds to the occurrence of one or more delimited
        substitution markers in the template to be substituted for the entire
        text of the second element. Start and end delimiters are set using
        'subs_delims'.
        For example, the tuple ("ABC", "text") would convert:
            The <ABC> is working
        to:
            The text is working
        Substitutions are performed in the order they appear in the 'subs'
        tuple. Recursive substitution as the tag parsing proceeds is thus
        feasible.
    subs_delims : 2-tuple of str, optional
        Indicates starting and ending strings to be used as delimiters for the
        substitution operations into the input template; the 0th and 1st
        elements are, respectively, the starting and ending strings.  Defaults
        are "<" and ">", respectively.  Strings of any length are acceptable.
        For example, to substitute a tag of the form {|TAG|}, the tuple
        ("{|","|}") should be passed to subs_delims.  Any elements past the
        second are ignored. No checking is performed for whether the delimiters
        are "sensible" or not.
    sim_name : str, optional
        Basename to use for the input/output/working files.
        If omitted, "orcarun" will be used.
    inp_ext : str
        Extension to be used for the input file generated.
    out_ext : str
        Extension to be used for the output file generated.

    Returns
    -------
    if wait_to_complete = True:
        oo : ORCA_OUTPUT
            Collected object containing output results

    #!DOC: Return values for execute_orca, once they exist

    Raises
    ------
    ValueError : If inp_ext and out_ext are identical.

    """

    # Imports
    import os, subprocess as sp
    from ..output import ORCA_OUTPUT
    from ..xyz import OPAN_XYZ
    from ..grad import ORCA_ENGRAD
    from ..hess import ORCA_HESS
    from ..utils import template_subst

    # Switch to dir; default exception fine for handling case of invalid dir.
    os.chdir(work_dir)

    # Check for inp_ext identical to out_ext
    if inp_ext == out_ext:
        raise(ValueError("'inp_ext' and 'out_ext' cannot be identical."))
    ##end if

    # Build the input and output file names
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

    # Perform the replacement into the exec_cmd of the input and output
    #  filenames.
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "INP" + subs_delims[1], inp_fname) \
            for s in exec_cmd]
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "OUT" + subs_delims[1], out_fname) \
            for s in exec_cmd_subs]
    exec_cmd_subs = [ \
                s.replace(subs_delims[0] + "NAME" + subs_delims[1], sim_name) \
            for s in exec_cmd_subs]

    # Perform the content substitutions into the template string
    input_text = template_subst(inp_tp, subs, subs_delims=subs_delims)

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
        #  things other than the output not existing
        o_out = ORCA_OUTPUT(out_fname,'file')
        try:
            o_xyz = OPAN_XYZ(path=os.path.join(work_dir, sim_name + ".xyz"), \
                                                                bohrs=bohrs)
        except IOError:
            o_xyz = None
        ## end try
        try:
            o_trj = OPAN_XYZ(path=os.path.join(work_dir, sim_name + ".trj"), \
                                                                bohrs=bohrs)
        except IOError:
            o_trj = None
        ## end try
        try:
            o_engrad = ORCA_ENGRAD(os.path.join(work_dir, sim_name + ".engrad"))
        except IOError:
            o_engrad = None
        ## end try
        try:
            o_hess = ORCA_HESS(os.path.join(work_dir, sim_name + ".hess"))
        except IOError:
            o_hess = None
        ## end try

    else:
        raise(NotImplementedError("Background execution not yet implemented."))
    ## end if

    # Return something appropriate
    #TODO: execute_orca: Must refine this, esp for the different exec modes
    return o_out, o_xyz, o_engrad, o_hess

## end def execute_orca


if __name__ == '__main__':
    print("Module not executable.")
