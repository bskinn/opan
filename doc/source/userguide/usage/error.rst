.. Usage for error

error
=====

.. note:: Most interactive use of ``opan`` will not require detailed
          knowledge of the custom errors in this module.

The custom exceptions in this module are all subclassed from
:class:`opan.error.OpanError`, which itself is a subclass of
:class:`python:Exception`.  In addition to the
typical error message included as part of initializing
an :class:`python:Exception`, the custom error subclasses of
:class:`~opan.error.OpanError` also define a typecode and a source attribute
(typically a filename or other data source) to allow more finely-grained
definition of error origins and types.  In the below example, attempting
to import the source file for this usage page as an OpenBabel XYZ file
quite sensibly results in an error:

.. testsetup:: errordemo

    import os
    os.chdir('source\\userguide\\usage')

.. testcleanup:: errordemo

    os.chdir('..\\..\\..')

.. doctest:: errordemo

    >>> x = opan.xyz.OpanXYZ(path='error.rst')
    Traceback (most recent call last):
    ...
    XYZError: (XYZFILE) No valid geometry found: XYZ file: error.rst

The custom exception :exc:`~opan.error.XYZError` is raised
with typecode :attr:`~opan.error.XYZError.XYZFILE`, indicating a problem
with the indicated input file.  The external data source causing the
exception is included after the final colon (``error.rst``, this file).
If no data source is relevant to a given exception, it is omitted.

The subclasses of :class:`opan.const.OpanEnum` are equipped with
membership testing of and iteration over their respective typecodes:

.. testsetup:: error

    from opan.error import XYZError

.. doctest:: error

    >>> 'XYZFILE' in XYZError
    True
    >>> [tc for tc in sorted(XYZError)]
    ['DIHED', 'NONPRL', 'OVERWRITE', 'XYZFILE']

Raising these exceptions follows standard syntax, save for the extra
'typecode' and 'source' parameters:

.. doctest:: error

    >>> raise XYZError(XYZError.OVERWRITE, "Spurious overwrite", "Console")
    Traceback (most recent call last):
    ...
    XYZError: (OVERWRITE) Spurious overwrite: Console


