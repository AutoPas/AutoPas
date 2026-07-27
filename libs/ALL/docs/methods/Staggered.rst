.. _staggeredDetails:

Staggered method
================

The staggered grid approach is a hierarchical one. In a first step the
work of all domains sharing the same cartesian coordinate with respect
to the highest dimension (z in three dimensions) is collected. Then,
like in the TENSOR strategy the layer width in this dimension is
adjusted based on comparison of the collected work with the collected
work of the neighboring domains in the same dimension. As a second step
each of these planes is divided into a set of columns, where all domains
share the same cartesian coordinate in the next lower dimension (y in
three dimensions). For each of these columns the before described
procedure is repeated, i.e. work collected and the width of the columns
adjusted accordingly. Finally, in the last step the work of individual
domains is compared to direct neighbors and the width in the last
dimension (x in three dimensions) adjusted. This leads to a staggered
grid, that is much better suited to describe inhomogeneous work
distributions than the TENSOR strategy.

Required number of vertices
 - two, one describing the lower left front point and one describing the
   upper right back point of the domain

Advantages
 - very good equalization results for the work load
 - maintains orthogonal domains

Disadvantages
 - changes topology of the domains and requires adjustment of neighbor
   relations
 - communication pattern in the calling code might require adjustment to
   deal with changing neighbors

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
