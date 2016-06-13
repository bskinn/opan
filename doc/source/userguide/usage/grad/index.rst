.. Usage for grad core/base/superclass

grad
====

*Draft scratch content...*

Gradient objects are a thing that have to be properly done in order for
various of the automation functiony applications of the ``opan`` to work right,
since different softwares make their things in different ways, but
the core automating thingy has to not care what software made the data but
the data has to be presented uniformly.  Most console interactions with
``opan`` won't care about this much, but it's worth noting here the things
that can be expected from ALL thingies, even though many other thingies
will likely be available for any given other various software.

Firstly, the instance members specified as having to be there by
:class:`~opan.grad.SuperOpanGrad`:

 * :attr:`gradient` -- 1-D |nparray| of |npfloat| -- Gradient data in
   :math:`\left(\mathrm{E_h\over B}\right)` |units|

 * :attr:`geom` -- 1-D |nparray| of |npfloat| -- Geometry data in
   :math:`\mathrm B` |units|

 * :attr:`atom_syms` -- |list| of |str| -- Atomic symbols in **ALL CAPS**

There will need to be a private ``_load`` method, but that shouldn't ever be
useful or usable, interactively.

Other than that, the below subpages describe the software-specific data
available in the specific subclass objects of
:class:`~opan.grad.SuperOpanGrad`.

**Implemented Subclasses**

.. toctree::
    :maxdepth: 1

    orca_engrad

