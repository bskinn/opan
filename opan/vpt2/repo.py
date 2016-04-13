#-------------------------------------------------------------------------------
# Name:        repo
# Purpose:     Encapsulation of input/output from HDF5 repository.
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     12 Jul 2015
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

""" Sub-module for HDF5 repo interactions for VPT2 anharmonic calculations.

.. warning::

    Module is under active development. API &c. may change
    with little notice.

**Classes**

.. autoclass:: OpanAnharmRepo


"""

# Imports


class OpanAnharmRepo(object):
    """ HDF5 repo interface for VPT2 anharmonic calculations.

    Operations here **DO NOT** ensure consistency with the surrounding data.
    Such consistency must be checked/handled at a higher level. This is just
    a wrapper class to facilitate I/O interactions!

    Currently no chunking or any filters are used. May or may not be worth
    robustification. Things stored are not likely to be huge....

    *#DOC: --- to be filled ---*

    **Instantiation**

    .. automethod:: __init__(...)


    **Class Variables**

    *To be documented*

    **Instance Variables**

    *To be documented*

    **Methods**

    *To be documented*

    """

    # Imports
    from ..const import EnumDispDirection as _E_DD


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
                _E_DD.NEGATIVE: 'n',
                _E_DD.NO_DISP: 'o',
                _E_DD.POSITIVE: 'p'
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
        elif fname is None:
            self._repo = None
            self.fname = None
        else:
            raise TypeError("Invalid filename type: {0}".format(type(fname)))
        ## end if

    ## end def __init__


    def close(self):
        """ #DOC: REPO.close docstring
        """

        # Imports
        import h5py as h5
        from ..error import RepoError

        # Close the repo if it's bound, and wipe link; else complain
        if self._repo != None:
            self._repo.close()
            self._repo = None
            # Leave self.fname defined for potential easy re-opening later on.
        else:
            raise RepoError(RepoError.STATUS,
                    "Cannot close; no repository open",
                    "Last repo: {0}".format(self.fname))
        ## end if

    ## end def close


    def load(self, fname):
        """ #DOC: REPO.load docstring
        """

        # Imports
        import h5py as h5
        from ..error import RepoError

        # If repo not None, complain
        if not self._repo == None:
            raise RepoError(RepoError.STATUS,
                    "Repository already open",
                    "File: {0}".format(self.fname))
        ## end if

        # If string passed, try opening h5.File; otherwise complain
        if isinstance(fname, str):
            self.fname = fname
            self._repo = h5.File(fname)
        else:
            raise TypeError("Invalid filename type: {0}".format(type(fname)))
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
        import h5py as h5
        from ..error import RepoError as RErr
        from ..const import EnumDispDirection as _E_DD
        from ..const import EnumAnharmRepoData

        # Must be valid mode
        if not (mode >=0 and isinstance(mode, int)):
            raise ValueError("Mode must be a non-negative integer")
        ## end if

        # Must be a valid disp direction
        if not disp in _E_DD:
            raise ValueError("'{0}' is not a valid " +
                    "displacement enum value".format(disp))
        ## end if

        # Must be a valid repo data type
        if not datatype in EnumAnharmRepoData:
            raise ValueError("'{0}' is not a valid " +
                    "data type enum value".format(datatype))
        ## end if

        # Get the appropriate geom group name
        if disp == _E_DD.NO_DISP:
            grpname = self.G_geom_ref
        else:
            grpname = self.F_mode_fmt % (mode, self.dircode[disp])
        ## end if

        # Get the group, creating if absent
        try:
            grp = self._repo.require_group(grpname)
        except AttributeError:
            # Presume repo not open/attached
            raise RErr(RErr.STATUS,
                        "Cannot store; no repository open", "")
        ## end try

        # If dataset exists in repo group, obliterate or complain. Can't use
        #  'require_dataset' since the object could be changing dimension(s)
        #  and h5py can't do that
        if datatype in grp.keys():
            if clobber:
                grp.pop(datatype)
            else:
                raise RErr(RErr.DATA,
                        "Dataset to be stored exists and clobber == False",
                        self._repo.filename)
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
        import os
        from ..const import EnumDispDirection as _E_DD
        from ..const import EnumAnharmRepoData
        from ..error import RepoError as RErr

        # Must be valid mode
        if not (mode >=0 and isinstance(mode, int)):
            raise ValueError("Mode must be a non-negative integer")
        ## end if

        # Must be a valid disp direction
        if not disp in _E_DD:
            raise ValueError("'{0}' is not a valid " +
                    "displacement enum value".format(disp))
        ## end if

        # Must be a valid data type
        if not datatype in EnumAnharmRepoData:
            raise ValueError("'{0}' is not a valid " +
                    "repository data type enum value".format(datatype))
        ## end if

        # Get the appropriate geom group name
        if disp == _E_DD.NO_DISP:
            grpname = self.G_geom_ref
        else:
            grpname = self.F_mode_fmt % \
                                    (mode, self.dircode[disp])
        ## end if

        # Get the group, complaining if repo not bound
        try:
            grp = self._repo.get(grpname)
        except AttributeError:
            raise RErr(RErr.STATUS,
                        "Cannot load; no repository open", "")
        ## end try

        # If succeeded, check if group not found
        if grp is None:
            raise RErr(RErr.GROUP,
                    "Group '" + grpname + "' not found",
                    self.fname)
        ## end if

        # Group found, try loading the data object; complain if not found
        try:
            out_data = grp.get(datatype).value
        except AttributeError:
            raise RErr(RErr.DATA,
                    "Dataset '" + datatype + "' not found",
                    self.fname)
        ## end try

        # Return as-is (h5py creates NumPy arrays of appropriate dimensions)
        return out_data

    ## end def get_data


    def has_data(self, datatype, mode, disp):
        """ #DOC: has_data docstring
        """

        # Imports
        from ..error import RepoError

        # Just try to get the data. Simply pass up all exceptions except for
        #  those if the group or data doesn't exist. DO re-raise a non-bound
        #  repo. Initialize value to DATA PRESENT (True) for case where data
        #  is present and the 'try' succeeds.
        retval = True
        try:
            self.get_data(datatype, mode, disp)
        except RepoError as RErr:
            if RErr.tc in (RErr.DATA, RErr.GROUP):
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
        from ..const import EnumAnharmRepoParam
        from ..error import RepoError as RErr

        # Must be a valid parameter name
        if not param in EnumAnharmRepoParam:
            raise ValueError("'{0}' is not a valid " +
                    "parameter enum value".format(param))
        ## end if

        # Get the params group, complaining if repo not bound
        try:
            grp = self._repo.require_group(self.G_param)
        except AttributeError:
            raise RErr(RErr.STATUS,
                        "Cannot store; no repository open", "")
        ## end try

        # If succeeded, check if group not found
        if grp is None:
            raise RErr(RErr.GROUP,
                    "Parameters group not found",
                    self.fname)
        ## end if

        # If dataset exists in repo group, obliterate or complain. Can't use
        #  'require_dataset' since the object could be changing dimension(s)
        #  and h5py can't do that
        if param in grp.keys():
            if clobber:
                grp.pop(param)
            else:
                raise RErr(RErr.DATA,
                        "Parameter to be stored exists and clobber == False",
                        self._repo.filename)
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
        import os
        from ..const import EnumAnharmRepoParam
        from ..error import RepoError as RErr

        # Must be a valid parameter name
        if not param in EnumAnharmRepoParam:
            raise ValueError("'{0}' is not a valid " +
                    "parameter enum value".format(param))
        ## end if

        # Get the params group, complaining if repo not bound
        try:
            grp = self._repo.require_group(self.G_param)
        except AttributeError:
            raise RErr(RErr.STATUS,
                        "Cannot load; no repository open", "")
        ## end try

        # Group should be guaranteed present with 'require_group', try
        #  loading the parameter; complain if not found.
        try:
            out_param = grp.get(param).value
        except AttributeError:
            raise RErr(RErr.DATA,
                    "Parameter '" + param + "' not found",
                    self.fname)
        ## end try

        # Return
        self._repo.flush()
        return out_param

    ## end def get_param


    def has_param(self, param):
        """ #DOC: has_param docstring
        """

        # Imports
        from ..error import RepoError

        # Try to get the param; pass along all errors, except 'data' error
        #  from RepoError
        retval = True
        try:
            self.get_param(param)
        except RepoError as RErr:
            if RErr.tc == RErr.DATA:
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
        from ..error import RepoError as RErr

        # Get the return value from the dataset, complaining if repo not
        #  bound. Using 'require_dataset' since any repo w/o a defined
        #  'dirty' value is just going to be assumed to be dirty.
        try:
            retval = self._repo.require_dataset(self.N_dirty, \
                                shape=(), dtype=bool, data=True).value
        except AttributeError:
            raise RErr(RErr.STATUS,
                        "Cannot report dirty status; no repository open", "")
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
            raise ValueError("'dirty' must be Boolean")
        ## end if

        # Try to retrieve the dataset; complain if repo not bound.
        try:
            dset = self._repo.require_dataset(self.N_dirty, \
                                shape=(), dtype=bool)
        except AttributeError:
            raise RErr(RErr.STATUS,
                        "Cannot set dirty status; no repository open", "")
        ## end try

        # Change the value to the indicated value
        dset[()] = dirty

        # Done!
    ## end def set_dirty


    def get_XYZ(self, mode, disp):
        """ #DOC: docstring for get_xyz
        """

        # Imports
        from ..xyz import OpanXYZ as XYZ
        from ..const import EnumAnharmRepoParam, EnumAnharmRepoData

        # Generate XYZ and return
        out_XYZ = XYZ(atom_syms=self.get_param(EnumAnharmRepoParam.atoms), \
                        coords=self.get_data(EnumAnharmRepoData.geom, mode, disp))
        return out_XYZ

    ## end def get_XYZ






if __name__ == '__main__':
    print("Module not executable.")
