.. Discussion of notational and other conventions in the opan documentation

Notational and Algebraic Conventions
====================================

In order to aid comprehension, this documentation strives to obey the following
formatting conventions.

**In text:**

 * Function/method parameters are set in *italics*.

 * When not used for emphasis, **bold** text is used to set apart
   words/phrases with specific contextual meaning.

 * Code snippets are set in fixed font, colored red, and
   boxed: ``[v for v in range(6)]``.

 * References to Python objects are set in fixed font, colored black
   (in most cases), and boxed: :class:`~opan.xyz.OpanXYZ`. Where practical,
   they will be linked to the relevant documentation (via
   :mod:`~sphinx.ext.intersphinx`).

 * Engineering units are set in upright Roman (serif) equation type:
   :math:`\left(\mathrm{E_h\over B^2}\right)`

 * Code examples are set in fixed font, formatted like a Python console
   session, and boxed the full width of the page (as per
   :mod:`~sphinx.ext.doctest`):

   >>> 1 + 2 + 3
   6



**In equations:**

 * x

|

 * *Upright sans-serif* in equations when referring to code variable names

 * Roman (serif) text in equations when referring to code functions/methods

 * Italics in text when referring to mathematical symbols

 * Italics in text when referring to code functions or names

 * N always refers to the number of atoms in the system of interest

 * Where relevant, G refers to the number of geometries in an OpenBabel XYZ file

 * Arbitrary integers represented as R, S, ...

 * Where included, capitalized instructions [add list] are to be interpreted in
   accordance with `RFC 2119 <http://tools.ietf.org/html/rfc2119>`__

 * [Others?]


.. toctree:
    :maxdepth: 2

