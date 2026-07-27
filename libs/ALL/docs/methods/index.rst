The Loadbalancing Methods
=========================

.. toctree::
  :hidden:
  
  Tensor.rst
  Staggered.rst
  Histogram.rst

A short overview of the methods is given below and more details can be
found on their respective page.

:ref:`Tensor product<tensorDetails>`
  The work on all processes is reduced over the cartesian planes in the
  systems. This work is then equalized by adjusting the borders of the
  cartesian planes.

:ref:`Staggered grid<staggeredDetails>`
  A 3-step hierarchical approach is applied, where:

  1. work over the cartesian planes is reduced, before the borders of
     these planes are adjusted
  2. in each of the cartesian planes the work is reduced for each
     cartesian column. These columns are then adjusted to each other to
     homogenize the work in each column
  3. the work between neighboring domains in each column is adjusted.
     Each adjustment is done locally with the neighboring planes,
     columns or domains by adjusting the adjacent boundaries.

:ref:`Topological mesh<topologicalDetails>` (WIP)
  In contrast to the previous methods this method adjusts domains not by
  moving boundaries but vertices, i.e.  corner points, of domains. For
  each vertex a force, based on the differences in work of the neighboring
  domains, is computed and the vertex is shifted in a way to equalize the
  work between these neighboring domains.

:ref:`Voronoi mesh<voronoiDetails>` (WIP)
  Similar to the topological mesh method, this method computes a force,
  based on work differences. In contrast to the topological mesh method,
  the force acts on a Voronoi point rather than a vertex, i.e. a point
  defining a Voronoi cell, which describes the domain. Consequently, the
  number of neighbors is not a conserved quantity, i.e. the topology may
  change over time. ALL uses the Voro++ library published by the Lawrence
  Berkeley National Laboratory for the generation of the Voronoi mesh.

:ref:`Histogram based staggered grid<histogramDetails>`
  Resulting in the same grid, as the staggered grid scheme, this scheme
  uses the cumulative work function in each of the three cartesian
  directions in order to generate this grid. Using histograms and the
  previously defined distribution of process domains in a cartesian grid,
  this scheme generates in three steps a staggered-grid result, in which
  the work is distributed as evenly as the resolution of the underlying
  histogram allows. In contrast to the previously mentioned schemes this
  scheme depends on a global exchange of work between processes.

:ref:`Orthogonal recursive bisection<orthogonalDetails>` (Planned)
  Comparable to the histogram based staggered grid scheme, this scheme
  uses cumulative work functions evaluated by histograms in order to find
  a new distribution of workload in a hierarchical manner. In contrast to
  the histogram based staggered grid scheme, the subdivision of domains is
  not based on a cartesian division of processes, but on the prime
  factorization of the total number of processes. During each step, the
  current slice of the system is distributed into smaller sub-slices along
  the longest edge, that roughly contain the same workload. The number of
  sub-slices in then determined by the corresponding prime number of that
  step, i.e. when trying to establish a distribution for 24 (3 * 2 * 2 *
  2) processes, the bisection would use four steps. In the first step the
  system would be subdivided into 3 sub domains along the largest edge,
  each containing roughly the same workload.  Each of these sub domains
  then is divided into two sub sub domains independently. This procedure
  is then repeated for all the prime numbers in the prime factorization.

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
