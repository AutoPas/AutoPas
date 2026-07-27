.. _histogramDetails:

Histogram method
================

The histogram-based method works differently than the first two methods
in such a way, that it is a global method that requires three distinct
steps for a single adjustment. In each of these steps the following
takes place: a partial histogram needs to be created over the workload,
e.g. number of particles, along one direction on each domain, then these
are supplied to the method, where a global histogram is computed. With
this global histogram a distribution function is created. This is used
to compute the optimal (possible) width of domains in that direction.
For the second and third steps the computation of the global histograms
and distribution functions take place in subsystems, being the results
of the previous step.  The result is the most optimal distribution of
domains according to the Staggered-grid method, at the cost of global
exchange of work, due to the global adjustment, which makes the method
not well suited to highly dynamic systems, due to the need of frequent
updates. On the other hand the method is well suited for static
problems, e.g. grid-based simulations.

Note: Currently the order of dimensions is: z-y-x.

Required number of vertices
 - two, one describing the lower left front point and one describing the
   upper right back point of the domain

Additional requirements
 - partial histogram created over the workload on the local domain in
   the direction of the current correction step

Advantages
 - supplies an optimal distribution of domains (restricted by width of
   bins used for the histogram)
 - only three steps needed to acquire result

Disadvantages
 - expensive in cost of communication, due to global shifts
 - requires preparation of histogram and shift of work between each of
   the three correction steps

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
