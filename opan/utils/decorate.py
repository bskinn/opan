#-------------------------------------------------------------------------------
# Name:        utils.decorate
# Purpose:     Submodule containing custom decorators for Open Anharmonic
#
# Author:      Brian Skinn
#                bskinn@alum.mit.edu
#
# Created:     30 Oct 2015
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

""" Custom decorators defined for Open Anharmonic.

All of these decorators were built using the class form per the exposition
`here <http://www.artima.com/weblogs/viewpost.jsp?thread=240845>`__
|extlink|.

Decorators
----------

.. autoclass:: arraysqueeze

.. autoclass:: kwarg_fetch

"""

# Imports
from functools import wraps as _wraps


# Decorators
class arraysqueeze(object):
    """ Converts selected arguments to squeezed np.arrays

    Pre-applies an ``np.asarray(...).squeeze()`` conversion to all positional
    arguments according to integer indices passed, and to any keyword arguments
    according to any strings passed.

    Each |int| argument passed instructs the decorator to convert the
    corresponding positional argument in the function definition.

    Each |str| argument passed instructs the decorator to convert the
    corresponding keyword argument.

    |str| parameters corresponding to keyword arguments absent in a
    particular function call and positional/optional argument indices
    beyond the range of the actual `\*args` of the decorated
    function are ignored.

    .. warning:: Likely fragile with optional arguments; needs to be tested.

    Arguments
    ---------
    \*args : |int| or |str|
        Arguments to convert to squeezed |nparray|.

    """

    def __init__(self, *args):
        """ Pass the positions of arguments to be arraysqueezed as integers
        """

        # Check for integers and strings
        for arg in args:
            if not (isinstance(arg, int) or isinstance(arg, str)):
                raise ValueError("Invalid decorator argument: {0}".format(arg))

        # If all ok, store
        self.arglist = args

    def __call__(self, f):
        """ Call the wrapped function after arraysqueezing selected arguments.

        Absent keyword arguments and positional/optional argument indices
        beyond the range of the actual `*args` of `f` are ignored.

        """
        @_wraps(f)
        def wrapped_f(*args, **kwargs):

            # Must import numpy
            import numpy as np

            # Working list of args values, since args is a tuple
            w_args = list(args)

            # Parse the arguments and arraysqueeze them. If squeezed to a
            #  singleton array, rewrap as a dimension-one array.
            for mod_arg in self.arglist:
                if isinstance(mod_arg, int) and mod_arg < len(w_args):
                    w_args[mod_arg] = np.asarray(w_args[mod_arg]).squeeze()
                    if not w_args[mod_arg].shape:
                        w_args[mod_arg] = w_args[mod_arg][np.newaxis]

                elif isinstance(mod_arg, str) and mod_arg in kwargs:
                    kwargs[mod_arg] = np.asarray(kwargs[mod_arg]).squeeze()
                    if not kwargs[mod_arg].shape:
                        kwargs[mod_arg] = kwargs[mod_arg][np.newaxis]

                # no 'else:' since type checked in __init__

            # Execute the function and return its result
            return f(*w_args, **kwargs)

        # end def wrapped_f

        # Return the decorated function
        return wrapped_f

    # end def __call__

# end class arraysqueeze


class kwarg_fetch(object):
    """ Fetch a function call to a missing or |None| keyword argument

    Arguments
    ---------
    (fill in here)

    *A 'kwarg' argument

    .. warning::

        Fragile to needing to pass optional arguments to the callable?

    .. Decorator built using the class form per the exposition
        `here <http://www.artima.com/weblogs/viewpost.jsp?thread=240845>`__
        |extlink|.

    """

    import keyword as _kw

    @classmethod
    def ok_kwarg(cls, val):
        """ DOCSTRING """
        try:
            return str.isidentifier(val) and not cls._kw.iskeyword(val)
        except TypeError:
            # Catch integer values in particular; fail if found
            return False

    @classmethod
    def ok_tuplearg(cls, val):
        """ DOCSTRING """

        if not (isinstance(val, tuple) and len(val) == 2 and
                isinstance(val[0], int) and cls.ok_kwarg(val[1])):
            return False
        return True

    @classmethod
    def ok_argarg(cls, val):
        """ DOCSTRING """
        try:
            ok_kw = str.isidentifier(val) and not cls._kw.iskeyword(val)
        except TypeError:
            # Again, catch integer values in particular and fail if found
            ok_kw = False

        ok_pos = isinstance(val, int)

        ok_tup = cls.ok_tuplearg(val)

        return ok_kw or ok_pos or ok_tup

    def __init__(self, *args, **kwargs):
        """Initialize with the keyword, callable, and relevant arguments

        """

        # Don't want named arguments anywhere in this initializer, since
        #  that would constrain the keywords allowable for calls to the
        #  fetching callable. 'kw' and 'c' probably aren't all that common,
        #  but better to make it more robust/flexible, especially since
        #  the fix is pretty simple.

        # Convert args to a list
        args = list(args)

        # Retrieve the keyword id and the callable
        kw = args.pop(0)
        c = args.pop(0)

        # Proof and store all of the things
        if not self.ok_kwarg(kw):
            raise ValueError("'kw' argument must be a valid non-keyword "
                             "identifier")
        self.kw = kw

        if not callable(c):
            raise ValueError("'c' argument must be callable")
        self.c = c

        if not all(map(self.ok_argarg, args)):
            raise ValueError("All 'args' must be valid non-keyword "
                             "identifier strings, integers, or (int, str) "
                             "tuples (str as valid identifiers)")
        if not all(map(self.ok_argarg, kwargs.values())):
            raise ValueError("All 'kwargs' values must be valid non-keyword "
                             "identifier strings, integers, or (int, str) "
                             "tuples (str as valid identifiers)")
        self.arglist = args
        self.kwarglist = kwargs

    def __call__(self, f):
        """Call the wrapped function after any needed fetch."""

        @_wraps(f)
        def wrapped_f(*args, **kwargs):

            # Check for if the target kwarg is missing
            if self.kw not in kwargs:
                # Missing. Must fetch.
                # Initialize as empty the arguments list to pass to the
                # callable
                fetch_args = []

                # Populate fetch_args sequentially in the order specified by
                # the particular decorator constructor
                for a in self.arglist:
                    # Simple type checks on the argument specifiers
                    # should suffice since they were checked at
                    # construction.
                    if isinstance(a, str):
                        # Keyword argument; must handle possible absence
                        if a in kwargs:
                            fetch_args.append(kwargs[a])
                        else:
                            fetch_args.append(None)
                    elif isinstance(a, tuple):
                        # Tuple argument for optional-positional args.
                        # Could be present as positional or as keyword,
                        # or could be absent.
                        if len(args) > a[0]:
                            # Sufficient positional arguments; assume
                            # present and passed as optional-positional
                            fetch_args.append(args[a[0]])
                        elif a[1] in kwargs:
                            # Present in the passed-in kwargs
                            fetch_args.append(kwargs[a[1]])
                        else:
                            # Not found
                            fetch_args.append(None)
                    else:
                        # Positional argument; integer value
                        fetch_args.append(args[a])

                # Populate fetch_kwargs according to what was specified
                # at decorator construction
                fetch_kwargs = {}
                for item in self.kwarglist.items():
                    # Same as above -- simple type checks should suffice
                    if isinstance(item[1], str):
                        # Keyword argument; must handle possible absence
                        if item[1] in kwargs:
                            fetch_kwargs.update({item[0]: kwargs[item[1]]})
                        else:
                            fetch_kwargs.update({item[0]: None})
                    elif isinstance(item[1], tuple):
                        # Optional-positional
                        if len(args) > item[1][0]:
                            # Sufficient positional args
                            fetch_kwargs.update({item[0]: args[item[1][0]]})
                        elif item[1][1] in kwargs:
                            # Insufficient positional; add from kwargs
                            # if present
                            fetch_kwargs.update({item[0]: kwargs[item[1][1]]})
                        else:
                            # Not found; store None
                            fetch_kwargs.update({item[0]: None})
                    else:
                        # Assume positional, integer value
                        fetch_kwargs.update({item[0]: args[item[1]]})

                # Call the callable and store the result into the target
                # keyword
                c_result = self.c(*fetch_args, **fetch_kwargs)
                kwargs.update({self.kw: c_result})

            # Whether the target kwarg was present or generated/injected,
            # call the wrapped function
            return f(*args, **kwargs)

        # end def wrapped_f

        # Return the wrapped function
        return wrapped_f

    # end def __call__

# end class kwarg_fetch


if __name__ == '__main__':
    print("Module not executable.")


