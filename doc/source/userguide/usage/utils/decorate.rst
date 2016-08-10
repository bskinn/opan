.. Usage for utils.decorate

utils.decorate
==============

**arraysqueeze**

*Stuff*


**kwargfetch**

This decorator is intended for use within sets of nested functions,
where a call to an "outer" function results in a chain of calls into
the nesting structure.  If each of these nested calls requires

Dummy test block:

.. testsetup:: kwargfetch

    from opan.utils.decorate import kwargfetch

.. doctest:: kwargfetch

    >>> def f_2p(a, b):
    ...     return 2 * a + 3 * b

    >>> @kwargfetch('kw', f_2p, 1, 'm')
    ... def testfxn(x, y, **kwargs):
    ...     return (y - x) * kwargs['kw']

    >>> testfxn(3, 7, m=4)
    ... # kw=f_2p(7, 4) = 26
    ... # testfxn then returns 4*26 = 104
    104

    >>> testfxn(y=7, m=4, x=3)
    ... # kwargfetch is robust against positional arguments
    ... # passed as keyword arguments
    104







