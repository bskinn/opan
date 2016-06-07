.. Discussion of notational and other conventions in the opan documentation

Notational and Algebraic Conventions
====================================

In order to aid comprehension, this documentation strives to obey the following
formatting/notation conventions.

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

 * Code examples are set in fixed font and boxed the full width of the
   page::

       from opan.utils.vector import proj
       pvec = proj(vec1, vec2)

   Interactive usage examples are formatted like a Python console
   session, as per :mod:`~sphinx.ext.doctest`:

   >>> 1 + 2 + 3
   6

 * Where included, the key words "MUST", "MUST NOT", "REQUIRED", "SHALL",
   "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED",  "MAY", and
   "OPTIONAL" in this documentation are to be interpreted as described in
   :rfc:`2119`.

 * Mathematical symbols and engineering units are set as in equations, below.



**In equations:**

 * Roman variables are set in serif italics:
   :math:`x + C`

 * Lowercase Greek variables are set in italics:
   :math:`\theta + \pi \over \gamma`

 * Uppercase Greek variables are set in upright serif:
   :math:`\Phi + \Theta`

 * Vectors are set in bold, upright, serif Roman symbols:
   :math:`\mathbf{r}_2 - \mathbf{r}_1`

 * Matrices are set as bold, uppercase, serif Roman or Greek symbols:
   :math:`\mathbf{M}\!^{^{^{-1}\!/_2}} \mathbf{F}
   \mathbf{M}\!^{^{^{-1}\!/_2}}`

 * Engineering units are set in upright Roman (serif) equation type:
   :math:`\left(\mathrm{E_h\over B^2}\right)`


**Common symbols:**

 * :math:`N` -- the number of atoms in the system of interest

 * :math:`G` -- the number of geometries in an OpenBabel XYZ file

 * :math:`R, S, ...` -- Arbitrary integers 



.. toctree:
    :maxdepth: 2

