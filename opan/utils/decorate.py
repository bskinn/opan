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

.. autoclass:: kwargfetch

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


class kwargfetch(object):
    """Fetch a missing keyword argument with a custom callable & arguments

    This decorator implements a form of non-persistent memoization for
    use in networks of inter-related and/or nested functions, where:

    * External users may have reason to call any of the functions directly
    * Most or all of the functions call one or more of the same specific
      "supporting" functions that potentially represent significant
      computational overhead
    * Calls with identical function arguments are not likely to recur
      in typical use by external users, and thus fully persistent memoization
      would in general be a waste of memory

    The memoization is implemented via injection of a specific keyword
    argument into a call to the wrapped function, where the inserted value
    is obtained from a call in turn to a specified callable using
    arguments drawn from the wrapped call.  If the target keyword argument
    is already present in the wrapped call, no action is taken.

    .. note::

        The API description below is wholly non-intuitive and likely
        impossible to follow. The examples provided in the
        :doc:`User's Guide </userguide/usage/utils/decorate>` will probably
        be much more illuminating.

    Arguments
    ---------
    args[0]
        |str| --
        Name of the keyword argument to be injected into the call to the
        wrapped function

    args[1]
        |callable| --
        Object to call to generate value to be injected into the target
        keyword argument (`args[0]`)

    args[2..n]
        |int| or |str| --
        Indicate which positional (|int|) and keyword (|str|) parameters
        of the wrapped function call are to be passed to the
        |callable| of `args[1]` as positional parameters, in the
        order provided within `args[2..n]`

    kwargs
        |int| or |str| --
        Indicate which positional (|int|) and keyword (|str|) parameters
        of the wrapped function call are to be passed to the
        |callable| of `args[1]` as keyword parameters, where the keys
        indicated here in `kwargs` are those used in the call to
        `args[1]`


    .. Decorator built using the class form per the exposition
       `here <http://www.artima.com/weblogs/viewpost.jsp?thread=240845>`__
       |extlink|.

    """

    @staticmethod
    def ok_kwarg(val):
        """Helper method for screening keyword arguments"""

        import keyword

        try:
            return str.isidentifier(val) and not keyword.iskeyword(val)
        except TypeError:
            # Non-string values are never a valid keyword arg
            return False

    @classmethod
    def ok_argarg(cls, val):
        """Helper method for screening valid arguments of any type"""

        return cls.ok_kwarg(val) or isinstance(val, int)

    def __init__(self, *args, **kwargs):
        """Initialize with the keyword, callable, and relevant arguments"""

        # Don't want named arguments anywhere in this initializer, since
        #  that would constrain the keywords allowable for calls to the
        #  fetching callable. 'kw' and 'c' probably aren't all that common,
        #  but better to make it more robust/flexible, especially since
        #  the fix is pretty simple.

        # Convert args to a list
        args = list(args)

        # Retrieve the keyword id and the callable
        try:
            kw = args.pop(0)
        except IndexError as e:
            raise TypeError("'Target keyword' argument absent") from e

        try:
            c = args.pop(0)
        except IndexError as e:
            raise TypeError("'Callable' argument absent") from e

        # Proof and store all of the things
        if not self.ok_kwarg(kw):
            raise ValueError("'kw' argument must be a valid non-keyword "
                             "identifier")
        self.kw = kw

        if not callable(c):
            raise TypeError("'c' argument must be callable")
        self.c = c

        if not all(map(self.ok_argarg, args)):
            raise ValueError("All 'args' must be valid non-keyword "
                             "identifier strings or integers")
        if not all(map(self.ok_argarg, kwargs.values())):
            raise ValueError("All 'kwargs' values must be valid non-keyword "
                             "identifier strings or integers")
        self.arglist = args
        self.kwarglist = kwargs

    def __call__(self, f):
        """Call the wrapped function after any needed fetch"""

        from inspect import signature as sig

        @_wraps(f)
        def wrapped_f(*args, **kwargs):

            # Check for if the target kwarg is missing
            if self.kw not in kwargs:
                # Missing. Must fetch.

                # Retrieve and materialize the enumerated arguments list.
                # Depends on the parameters being contained in an OrderedDict
                # so that signature order is retained.
                params = list(sig(f).parameters)

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
                        # Keyword argument; handle possible absence with get()
                        fetch_args.append(kwargs.get(a))
                    else: 
                        # Integer argument for (optional-)positional args.
                        # Could be present as positional or as keyword,
                        # or could be absent.
                        if len(args) > a:
                            # Sufficient positional arguments; assume
                            # present and passed as optional-positional
                            fetch_args.append(args[a])
                        else:
                            # The **kwargs is not valid for this, so exclude
                            pname = params[:-1][a]
                            if pname in kwargs:
                                # Present in the passed-in kwargs
                                fetch_args.append(kwargs[pname])
                            else:
                                # Not found; pass the function default
                                fetch_args.append(sig(f).parameters[pname]
                                                  .default)

                # Populate fetch_kwargs according to what was specified
                # at decorator construction
                fetch_kwargs = {}
                for item in self.kwarglist.items():
                    # Same as above -- simple type checks should suffice
                    if isinstance(item[1], str):
                        # Keyword argument; handle possible absence with get()
                        fetch_kwargs.update({item[0]: kwargs.get(item[1])})
                    else:
                        # Optional-positional
                        if len(args) > item[1]:
                            # Sufficient positional args
                            fetch_kwargs.update({item[0]: args[item[1]]})
                        else:
                            # The **kwargs is not valid for this, so exclude
                            pname = params[:-1][item[1]]
                            if pname in kwargs:
                                # Insufficient positional; add from kwargs
                                # if present
                                fetch_kwargs.update({item[0]: kwargs[pname]})
                            else:
                                # Not found; store the function default
                                fetch_kwargs.update({item[0]:
                                                    sig(f)
                                                    .parameters[pname].default})

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


