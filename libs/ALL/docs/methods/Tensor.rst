.. _tensorDetails:

Tensor method
=============

In order to equalize the load of individual domains, the assumption is
made that this can be achieved by equalizing the work in each cartesian
direction, i.e. the work of all domains having the same coordinate in a
cartesian direction is collected and the width of all these domains in
this direction is adjusted by comparing this collective work with the
collective work of the neighboring domains. This is done independently
for each cartesian direction in the system.

Required number of vertices
 - two, one describing the lower left front point and one describing the
   upper right back point of the domain

Advantages
 - topology of the system is maintained (orthogonal domains, neighbor
   relations)
 - no update for the neighbor relations is required in the calling
   code
 - if a code was able to deal with orthogonal domains, only small
   changes are expected to include this strategy

Disadvantages
 - due to the comparison of collective work loads in cartesian layers
   and restrictions resulting from this construction, the final result
   might lead to a sub-optimal domain distribution

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
