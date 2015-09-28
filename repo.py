#-------------------------------------------------------------------------------
# Name:        repo
# Purpose:     Encapsulation of input/output from HDF5 repository.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Jul 2015
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


class OPAN_REPO(object):
    """ #DOC: OPAN_REPO docstring

    Operations here **DO NOT** ensure consistency with the surrounding data.
    Such consistency must be checked/handled at a higher level. This is just
    a wrapper class to facilitate I/O interactions!

    Currently no chunking or any filters are used. May or may not be worth
    robustification. Things stored are not likely to be huge....

    --- to be filled ---

    Instantiation
    -------------
    __init__(...)
        --- to be completed ---

    Class Variables
    ---------------


    Instance Variables
    ------------------


    Methods
    -------


    Generators
    ----------
    (none)

    """

    # Imports
    from .const import E_DispDirection as E_DispDir


    # Class variables
    # Data group names
    G_param = 'param'
    G_geom_ref = 'geom_ref'
    G = frozenset([
                G_geom_ref,
                G_param
                    ])

    # Dynamic formatting names
    F_mode_fmt = 'm%05d%c'

    # String for 'dirty' flag
    N_dirty = 'dirty'

    # dict for translating DispDir enum to direction code
    dircode = {
                E_DispDir.Negative: 'n',
                E_DispDir.NoDisp: 'o',
                E_DispDir.Positive: 'p'
                }


    # Instance methods
    def __init__(self, fname=None):
        """ #DOC: REPO.__init__ docstring

        """

        # Imports
        import h5py as h5

        # If string passed, try opening h5.File; otherwise init with repo
        #  link as None
        if isinstance(fname, str):
            self.fname = fname
            self._repo = h5.File(fname)
        elif fname == None:
            self._repo = None
            self.fname = None
        else:
            raise(TypeError("Invalid filename type: " + str(type(fname))))
        ## end if

    ## end def __init__


    def close(self):
        """ #DOC: REPO.close docstring
        """

        # Imports
        import h5py as h5
        from .error import REPOError

        # Close the repo if it's bound, and wipe link; else complain
        if self._repo != None:
            self._repo.close()
            self._repo = None
            # Leave self.fname defined for potential easy re-opening later on.
        else:
            raise(REPOError(REPOError.status, \
                    "Cannot close; no repository open", \
                    "Last repo: " + str(self.fname) ))
        ## end if

    ## end def close


    def load(self, fname):
        """ #DOC: REPO.load docstring
        """

        # Imports
        import h5py as h5
        from .error import REPOError

        # If repo not None, complain
        if not self._repo == None:
            raise(REPOError(REPOError.status, \
                    "Repository already open", \
                    "File: " + str(self.fname) ))
        ## end if

        # If string passed, try opening h5.File; otherwise complain
        if isinstance(fname, str):
            self.fname = fname
            self._repo = h5.File(fname)
        else:
            raise(TypeError("Invalid filename type: " + str(type(fname))))
        ## end if

        # Dirty/clean status will be as it was when repo was last closed

    ## end def load


    def is_open(self):
        """ #DOC: is_open docstring
        """

        # Imports
        import h5py as h5

        # Return whether repo is bound; depends on the 'None' assignment in
        #  'close'
        retval = isinstance(self._repo, h5.File)
        return retval

    ## end def isopen



    def store_data(self, data, datatype, mode, disp, clobber=False):
        """ #DOC: store_data docstring
        """

        # Imports
        import h5py as h5, numpy as np
        from .error import REPOError as RErr
        from .const import E_DispDirection as E_DispDir
        from .const import ERA_Data

        # Must be valid mode
        if not (mode >=0 and isinstance(mode, int)):
            raise(ValueError("Mode must be a non-negative integer"))
        ## end if

        # Must be a valid disp direction
        if not disp in E_DispDir:
            raise(ValueError("'" + str(disp) + "' is not a valid " + \
                    "displacement enum value"))
        ## end if

        # Must be a valid repo data type
        if not datatype in ERA_Data:
            raise(ValueError("'" + str(datatype) + "' is not a valid " + \
                    "data type enum value"))
        ## end if

        # Get the appropriate geom group name
        if disp == E_DispDir.NoDisp:
            grpname = self.G_geom_ref
        else:
            grpname = self.F_mode_fmt % \
                                    (mode, self.dircode[disp])
        ## end if

        # Get the group, creating if absent
        try:
            grp = self._repo.require_group(grpname)
        except AttributeError:
            # Presume repo not open/attached
            raise(RErr(RErr.status, \
                        "Cannot store; no repository open", ""))
        ## end try

        # If dataset exists in repo group, obliterate or complain. Can't use
        #  'require_dataset' since the object could be changing dimension(s)
        #  and h5py can't do that
        if datatype in grp.keys():
            if clobber:
                grp.pop(datatype)
            else:
                raise(RErr(RErr.data, \
                        "Dataset to be stored exists and clobber == False", \
                        self._repo.filename))
            ## end if
        ## end if

        # Store the new data. DOES NOT ENSURE CONSISTENCY with any
        #  other data in the repository.
        grp.create_dataset(datatype, data=data)

        # Set as dirty and flush the repo
        self.set_dirty(True)
        self._repo.flush()

    ## end def store_geom


    def get_data(self, datatype, mode, disp):
        """ #DOC: docstring for get_data
        """

        # Imports
        import os, numpy as np
        from .const import E_DispDirection as E_DispDir
        from .const import ERA_Data
        from .error import REPOError as RErr

        # Must be valid mode
        if not (mode >=0 and isinstance(mode, int)):
            raise(ValueError("Mode must be a non-negative integer"))
        ## end if

        # Must be a valid disp direction
        if not disp in E_DispDir:
            raise(ValueError("'" + str(disp) + "' is not a valid " + \
                    "displacement enum value"))
        ## end if

        # Must be a valid data type
        if not datatype in ERA_Data:
            raise(ValueError("'" + str(datatype) + "' is not a valid " + \
                    "repository data type enum value"))
        ## end if

        # Get the appropriate geom group name
        if disp == E_DispDir.NoDisp:
            grpname = self.G_geom_ref
        else:
            grpname = self.F_mode_fmt % \
                                    (mode, self.dircode[disp])
        ## end if

        # Get the group, complaining if repo not bound
        try:
            grp = self._repo.get(grpname)
        except AttributeError:
            raise(RErr(RErr.status, \
                        "Cannot load; no repository open", ""))
        ## end try

        # If succeeded, check if group not found
        if grp == None:
            raise(RErr(RErr.group, \
                    "Group '" + grpname + "' not found", \
                    self.fname))
        ## end if

        # Group found, try loading the data object; complain if not found
        try:
            out_data = grp.get(datatype).value
        except AttributeError:
            raise(RErr(RErr.data, \
                    "Dataset '" + datatype + "' not found", \
                    self.fname))
        ## end try

        # Matrixify and return
        out_data = np.matrix(out_data)
        return out_data

    ## end def get_data


    def has_data(self, datatype, mode, disp):
        """ #DOC: has_data docstring
        """

        # Imports
        from .error import REPOError

        # Just try to get the data. Simply pass up all exceptions except for
        #  those if the group or data doesn't exist. DO re-raise a non-bound
        #  repo. Initialize value to DATA PRESENT (True) for case where data
        #  is present and the 'try' succeeds.
        retval = True
        try:
            self.get_data(datatype, mode, disp)
        except REPOError as RErr:
            if RErr.tc in (RErr.data, RErr.group):
                retval = False
            ## end if
        ## end try

        # Should be good to return
        return retval

    ## end def has_data


    def store_param(self, value, param, clobber=False):
        """ #DOC: store_param docstring
        """

        # Imports
        from .const import ERA_Param
        from .error import REPOError as RErr

        # Must be a valid parameter name
        if not param in ERA_Param:
            raise(ValueError("'" + str(param) + "' is not a valid " + \
                    "parameter enum value"))
        ## end if

        # Get the params group, complaining if repo not bound
        try:
            grp = self._repo.require_group(self.G_param)
        except AttributeError:
            raise(RErr(RErr.status, \
                        "Cannot store; no repository open", ""))
        ## end try

        # If succeeded, check if group not found
        if grp == None:
            raise(RErr(RErr.group, \
                    "Parameters group not found", \
                    self.fname))
        ## end if

        # If dataset exists in repo group, obliterate or complain. Can't use
        #  'require_dataset' since the object could be changing dimension(s)
        #  and h5py can't do that
        if param in grp.keys():
            if clobber:
                grp.pop(param)
            else:
                raise(RErr(RErr.data, \
                        "Parameter to be stored exists and clobber == False", \
                        self._repo.filename))
            ## end if
        ## end if

        # Store the new data. DOES NOT ENSURE CONSISTENCY with any
        #  other data in the repository.
        grp.create_dataset(param, data=value)

        # Set as dirty and flush the repo
        self.set_dirty(True)
        self._repo.flush()

    ## end def store_param


    def get_param(self, param):
        """ #DOC: docstring for get_param
        """

        # Imports
        import os, numpy as np
        from .const import ERA_Param
        from .error import REPOError as RErr

        # Must be a valid parameter name
        if not param in ERA_Param:
            raise(ValueError("'" + str(param) + "' is not a valid " + \
                    "parameter enum value"))
        ## end if

        # Get the params group, complaining if repo not bound
        try:
            grp = self._repo.require_group(self.G_param)
        except AttributeError:
            raise(RErr(RErr.status, \
                        "Cannot load; no repository open", ""))
        ## end try

        # Group should be guaranteed present with 'require_group', try
        #  loading the parameter; complain if not found.
        try:
            out_param = grp.get(param).value
        except AttributeError:
            raise(RErr(RErr.data, \
                    "Parameter '" + param + "' not found", \
                    self.fname))
        ## end try

        # If param is an array, matrixify.
        if len(grp.get(param).shape) != 0:
            out_param = np.matrix(out_param)
        ## end if

        # Return
        self._repo.flush()
        return out_param

    ## end def get_param


    def has_param(self, param):
        """ #DOC: has_param docstring
        """

        # Imports
        from .error import REPOError

        # Try to get the param; pass along all errors, except 'data' error
        #  from REPOError
        retval = True
        try:
            self.get_param(param)
        except REPOError as RErr:
            if RErr.tc == RErr.data:
                retval = False
            ## end if
        ## end try

        # Should be good to return
        return retval

    ## end def has_param


    def is_dirty(self):
        """ #DOC: docstring for is_dirty
        """

        # Imports
        from .error import REPOError as RErr

        # Get the return value from the dataset, complaining if repo not
        #  bound. Using 'require_dataset' since any repo w/o a defined
        #  'dirty' value is just going to be assumed to be dirty.
        try:
            retval = self._repo.require_dataset(self.N_dirty, \
                                shape=(), dtype=bool, data=True).value
        except AttributeError:
            raise(RErr(RErr.status, \
                        "Cannot report dirty status; no repository open", ""))
        ## end try

        # Either way it evaluated, should be good to return. Flush first.
        self._repo.flush()
        return retval

    ## end def is_dirty


    def set_dirty(self, dirty):
        """ #DOC: set_clean docstring
        """

        # Complain if 'dirty' isn't boolean
        if not isinstance(dirty, bool):
            raise(ValueError("'dirty' must be Boolean"))
        ## end if

        # Try to retrieve the dataset; complain if repo not bound.
        try:
            dset = self._repo.require_dataset(self.N_dirty, \
                                shape=(), dtype=bool)
        except AttributeError:
            raise(RErr(RErr.status, \
                        "Cannot set dirty status; no repository open", ""))
        ## end try

        # Change the value to the indicated value
        dset[()] = dirty

        # Done!
    ## end def set_dirty


    def get_XYZ(self, mode, disp):
        """ #DOC: docstring for get_xyz
        """

        # Imports
        from .xyz import OPAN_XYZ as XYZ
        from .const import ERA_Param, ERA_Data

        # Generate XYZ and return
        out_XYZ = XYZ(atom_syms=self.get_param(ERA_Param.atoms), \
                        coords=self.get_data(ERA_Data.geom, mode, disp))
        return out_XYZ

    ## end def get_XYZ






if __name__ == '__main__':
    print("Module not executable.")
