..  In ReStructured Text (ReST) indentation and spacing are very important (it is how ReST knows what to do with your
    document). For ReST to understand what you intend and to render it correctly please to keep the structure of this
    template. Make sure that any time you use ReST syntax (such as for ".. sidebar::" below), it needs to be preceded
    and followed by white space (if you see warnings when this file is built they this is a common origin for problems).

..  We allow the template to be standalone, so that the library maintainers add it in the right place

..  Firstly, let's add technical info as a sidebar and allow text below to wrap around it. This list is a work in
    progress, please help us improve it. We use *definition lists* of ReST_ to make this readable.

..  sidebar:: Software Technical Information

  Name
    A Load Balancing Library (ALL)/GMPM-PoC

  Language
    Fortran/C/C++

  Licence
    `BSD 3-Clause <https://choosealicense.com/licenses/bsd-3-clause/>`_   

  Documentation Tool
    In source documentation using Doxygen, additional man pages and plain
    text

  Application Documentation
    Non public/in repo

  Relevant Training Material
    `Webinar (YT) <https://www.youtube.com/watch?v=2K2YFdzIJF4&list=PLmhmpa4C4MzY02eaacXImTts2aGJHrdwQ&index=3>`_

  Software Module Developed by
    Stephan Schulz


..  In the next line you have the name of how this module will be referenced in the main documentation (which you  can
    reference, in this case, as ":ref:`example`"). You *MUST* change the reference below from "example" to something
    unique otherwise you will cause cross-referencing errors. The reference must come right before the heading for the
    reference to work (so don't insert a comment between).

.. _all_mpm_integration:

###############
MPM Integration
###############

..  Let's add a local table of contents to help people navigate the page

..  contents:: :local:

..  Add an abstract for a *general* audience here. Write a few lines that explains the "helicopter view" of why you are
    creating this module. For example, you might say that "This module is a stepping stone to incorporating XXXX effects
    into YYYY process, which in turn should allow ZZZZ to be simulated. If successful, this could make it possible to
    produce compound AAAA while avoiding expensive process BBBB and CCCC."


The material point method (MPM) is used to simulate continuous matter and
is especially suited for the simulation of large deformations. Once large
deformation are present, a dynamic load balancing solution is sensible to
efficiently simulate large systems. Even if the initial work distribution
is good, it is very often the case, that this distribution is much less so
during the simulation run itself. The load balancing library ALL provides
an easy plug and play solution to this problem and this module describes
the details in how the library is integrated. Thanks to the good load
balancing provided by the library larger systems can be simulated with
less computational cost.

Purpose of Module
_________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

This module shows the straight forwardness of including the load balancing
library into already existing code. Depending on the simulation code
additional precautions must be taken and those needed for the MPM
simulation code are presented here. The prerequisites for the simulation
code are also shown. Looking at these will help determine whether a
simulation code is particularly suited for integrating ALL or if some
further work is needed when integrating.

.. This module also shows a real world application of the :ref:`Fortran
   interface<all_fortran_interface>` provided with
   :ref:`ALL<ALL_background>`.

This module also shows a real world application of the Fortran interface
provided with :ref:`ALL<ALL_background>` (documented in :ref:`all_fortran_interface`).

The MPM simulation code with integrated ALL is used by Stephan Schulz in
his thesis.

.. If needed you can include latex mathematics like
   :math:`\frac{ \sum_{t=0}^{N}f(t,k) }{N}`
   which won't show up on GitLab/GitHub but will in final online documentation.


Background Information
______________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

.. If the modifications are to an existing code base (which is typical)
   then this would be the place to name that application. List any
   relevant urls and explain how to get access to that code. There needs
   to be enough information here so that the person reading knows where to
   get the source code for the application, what version this information
   is relevant for, whether this requires any additional patches/plugins,
   etc.

.. Overall, this module is supposed to be self-contained, but linking to
   specific URLs with more detailed information is encouraged. In other
   words, the reader should not need to do a websearch to understand the
   context of this module, all the links they need should be already in
   this module.

The load balancing library ALL is integrated into the material point
method simulation code GMPM-PoC, which is written by Stephan Schulz during
his thesis. The simulation code will be released to the public in the
future.

Certain aspects of the simulation code require additional treatment of the
library, or additional features of the library. First, the open boundaries
of the simulation code require continuous updates of the outer domain
boundaries of the boundary domains. The system extent is adapted to the
particle's bounding box each time step. This also means, the geometry
generated in the last balance step by the library cannot be used directly.
It is therefore saved by the simulation code, adapted to the new system
extent and then given to the library as the basis for the new geometry.

The communication is based on grid halos and only accommodates nearest
neighbor communication. This causes the minimum domain size to be the
width of exactly this halo. The library supports this feature and only the
aforementioned outer domain bounds must be checked for compliance with the
minimum size. The other domain boundaries are automatically sufficiently
large due to the library.

And lastly, the particle communication at the end of each time step also
only accounts for nearest neighbor communication. This means, that a
domain's boundary must not change so much, that it needs to retrieve
particles from a process that is not its nearest neighbor. Due to the way
the library moves boundaries in the staggered grid and tensor approaches,
this is also already guaranteed to be true. There is always an overlap
between the old domain's volume and the new domain's.

However, the library also has a few requirements of the simulation code.
Due to changing domains, particles must be able to be communicated across
processes, which is implemented for all particle codes with moving
particles. Depending on the load balancing method this communication may
be more involved. In the case of the tensor method the topology does not
change and every process has the same 26 neighbors during the entire
simulation. If, however, the staggered grid approach is used, the
communication must also handle changing number of neighbors and determine
where they are and what regions they belong to. For example it is common,
that one half of a boundary must be communicated to one process and the
other to a different one. So the fixed relationship between boundaries and
neighbors is broken up. The GMPM-PoC code was already designed with such a
communication scheme in mind and provided the necessary flexibility to
simply enable the staggered grid method after fixing a few communication
bugs.


Building and Testing
____________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

.. Provide the build information for the module here and explain how tests are run. This needs to be adequately detailed,
   explaining if necessary any deviations from the normal build procedure of the application (and links to information
   about the normal build process needs to be provided).

To build the code just run ``make LB=ALL`` and everything should be build
automatically including dependencies. Make sure the correct compiler are
found in the path and if you want to use Intel compilers you need to set
``COMPILER=INTEL`` as well. The normal caveats and required modules for
some HPC systems are the described in the main code's ``README``.

Source Code
___________

.. Notice the syntax of a URL reference below `Text <URL>`_ the backticks matter!

.. Here link the source code *that was created for the module*. If you are using Github or GitLab and the `Gitflow Workflow
  <https://www.atlassian.com/git/tutorials/comparing-workflows#gitflow-workflow>`_ you can point to your feature branch.
  Linking to your pull/merge requests is even better. Otherwise you can link to the explicit commits.

.. * `Link to a merge request containing my source code changes
     <https://github.com/easybuilders/easybuild-easyblocks/pull/1106>`_

.. There may be a situation where you cannot do such linking. In this case, I'll go through an example that uses a patch
  file to highlight my source code changes, for that reason I would need to explain what code (including exact version
  information), the source code is for.

.. You can create a similar patch file by (for example if you are using git for your version control) making your changes
  for the module in a feature branch and then doing something like the following:

The main changes are the replacement of the original domain decomposition
function which used to equi partition the system extent. Now, ALL is
called to update the domain geometry.


.. literalinclude:: ./MPMIntegration.F90
   :linenos:
   :language: Fortran

To include the library and its VTK dependency into the existing make build
system, the following snippets were used. This builds a 'minimal' VTK and
links ALL against this. During the linking of the main simulation code VTK
is linked using ``$(VTK_LIB)`` where the order is very important. The
calling of make in this Makefile is deprecated and should be replaced by
appropriate calls to ``cmake --build`` and ``cmake --install``.


.. literalinclude:: ./MPMIntegration.makefile
   :linenos:
   :language: Makefile

..    :emphasize-lines: 2,9-11
.. If the patch is very long you will probably want to add it as a subpage which can be done as follows

.. .. toctree::
      :glob:
      :maxdepth: 1

..    patch

..  Remember to change the reference "patch" for something unique in your patch file subpage or you will have
    cross-referencing problems

.. you can reference it with :ref:`patch`

.. Here are the URL references used (which is alternative method to the one described above)

.. _ReST: http://www.sphinx-doc.org/en/stable/rest.html
.. _Sphinx: http://www.sphinx-doc.org/en/stable/markup/index.html

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
