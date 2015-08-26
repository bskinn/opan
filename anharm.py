#-------------------------------------------------------------------------------
# Name:        anharm
# Purpose:     Central module for VPT2 anharmonic wrapper.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     19 Jul 2015
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

# Imports


class ORCA_ANHARM(object):
    """ #DOC: ORCA_ANHARM class-level docstring
    """

    # Imports
    from .const import E_Software

    def __init__(self):
        """ #DOC: ORCA_ANHARM.__init__ docstring (if needed)?
        """

        # Declare the holding variables for the XYZ, GRAD, HESS
        self.w_xyz = None
        self.w_grad = None
        self.w_hess = None

        # Declare the holding variable for the repo
        self.repo = None

    ## end def __init__


    def new_from_files(self, basepath, basename, repo, \
                    bohrs=False, \
                    software=E_Software.ORCA, \
                    repo_clobber=False, **kwargs):
        """ #DOC: new_from_files docstring
                    xyz_ext='xyz', \
                    grad_ext='engrad', \
                    hess_ext='hess', \
                    kwargs
        """

        # Imports
        import os
        from os import path as osp
        from .xyz import OPAN_XYZ as OX
        from .grad import ORCA_ENGRAD as OE
        from .hess import ORCA_HESS as OH
        from .repo import OPAN_REPO as OR
        from .const import E_DispDirection as E_DDir, E_FileType as E_FT
        from .const import E_Software as E_SW
        from .const import DEF
        from .error import ANHARMError as ANHErr

##        # Store working directory for restore?
##        prev_dir = os.getcwd()

        # Complain if anything is already bound
        if not self.w_xyz == None:
            raise(ANHErr(ANHErr.status, \
                    "XYZ object is already bound", \
                    ""))
        ## end if
        if not self.w_grad == None:
            raise(ANHErr(ANHErr.status, \
                    "GRAD object is already bound", \
                    ""))
        ## end if
        if not self.w_hess == None:
            raise(ANHErr(ANHErr.status, \
                    "HESS object is already bound", \
                    ""))
        ## end if
        if not self.repo == None:
            raise(ANHErr(ANHErr.status, \
                    "Repository object is already bound", \
                    ""))
        ## end if

        # RESUME: anharm--factor for loading from different software pkgs

        # Load the three data files
        self.w_xyz = OX( osp.join(basepath, \
                basename + osp.extsep + xyz_ext) )
        self.w_grad = OE( osp.join(basepath, \
                basename + osp.extsep + engrad_ext), \
                0, E_DDir.NoDisp, 0.0 )
        self.w_hess = OH( osp.join(basepath, \
                basename + osp.extsep + hess_ext), \
                0, E_DDir.NoDisp, 0.0 )

        # Only accept new repos for now
        if not isinstance(repo, str):
            raise(TypeError("Must create new repository when loading " + \
                    "a new dataset."))
        ## end if

        # Repo is string, treat as filename and try to load
        # Check if it's a complete path
        # If it's a relative path, prepend the basepath
        if osp.split(repo[0]) > 0 and not osp.isabs(repo):
            repo = osp.join(basepath, repo)
        ## end if

        # Complain if it's a directory
        if osp.isdir(repo):
            raise(IOError("Cannot bind repository -- specified " + \
                    "location is a directory"))
        ## end if

        # If file exists ...
        if osp.isfile(repo):
                # Depending on clobber, either delete existing or raise error
                if repo_clobber:
                    # Clobber old repo
                    os.remove(repo)
                else:
                    # Raise error
                    raise(IOError("Target repository file exists and " + \
                            "clobber is disabled."))
                ## end if
            ## end if

        # Should be good to create the repo
        self.repo = OR(repo)

        # Store the geometry info, grad, and Hessian to the repo

        # Think that's all here?

    ## end def new_from_files

## end class ORCA_ANHARM



if __name__ == '__main__':
    print('Module not executable.')
