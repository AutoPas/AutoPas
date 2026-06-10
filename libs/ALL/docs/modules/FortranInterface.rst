..  In ReStructured Text (ReST) indentation and spacing are very important (it is how ReST knows what to do with your
    document). For ReST to understand what you intend and to render it correctly please to keep the structure of this
    template. Make sure that any time you use ReST syntax (such as for ".. sidebar::" below), it needs to be preceded
    and followed by white space (if you see warnings when this file is built they this is a common origin for problems).

..  We allow the template to be standalone, so that the library maintainers add it in the right place

..  Firstly, let's add technical info as a sidebar and allow text below to wrap around it. This list is a work in
    progress, please help us improve it. We use *definition lists* of ReST_ to make this readable.

..  sidebar:: Software Technical Information

  Name
    A Load Balancing Library (ALL)

  Language
    C++, Fortran interfaces available

  Licence
    `BSD 3-Clause <https://choosealicense.com/licenses/bsd-3-clause/>`_

  Documentation Tool
    In source provided by Doxygen, additional using Sphinx

  Application Documentation
    http://slms.pages.jsc.fz-juelich.de/websites/all-website/sphinx/api/ALL_module.html

  Relevant Training Material
    `Webinar (YT) <https://www.youtube.com/watch?v=2K2YFdzIJF4&list=PLmhmpa4C4MzY02eaacXImTts2aGJHrdwQ&index=3>`_

  Software Module Developed by
    Stephan Schulz


..  In the next line you have the name of how this module will be referenced in the main documentation (which you  can
    reference, in this case, as ":ref:`example`"). You *MUST* change the reference below from "example" to something
    unique otherwise you will cause cross-referencing errors. The reference must come right before the heading for the
    reference to work (so don't insert a comment between).

.. _all_fortran_interface:

#####################
ALL Fortran interface
#####################

..  Let's add a local table of contents to help people navigate the page

..  contents:: :local:

..  Add an abstract for a *general* audience here. Write a few lines that explains the "helicopter view" of why you are
    creating this module. For example, you might say that "This module is a stepping stone to incorporating XXXX effects
    into YYYY process, which in turn should allow ZZZZ to be simulated. If successful, this could make it possible to
    produce compound AAAA while avoiding expensive process BBBB and CCCC."

There is still a lot of Fortran code in use and the low entry barrier to
the language makes it an easy choice for beginning scientific programmers.
To be able to utilize the features of the library in Fortran code this
interface is provided. It provides basically the same functionality as the
C++ interface.

Purpose of Module
_________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

This module is necessary for any Fortran developers trying to use this
library.

.. It is currently in use by the Fortran Multi Particle Method written for
   the thesis of Stephan Schulz. This application of the interface is
   documented in the according :ref:`module<all_mpm_integration>`.

It is currently in use by the Fortran Multi Particle Method written for
the thesis of Stephan Schulz. This application of the interface is
documented in the module :ref:`all_mpm_integration`.

.. TODO:

.. * If there are published results obtained using this code, describe them briefly in terms readable for non-expert users.
  If you have few pictures/graphs illustrating the power or utility of the module, please include them with
  corresponding explanatory captions.

.. If you want to add a citation, such as [CIT2009]_, please check the source code to see how this is done. Note that
.. citations may get rearranged, e.g., to the bottom of the "page".

.. .. [CIT2009] This is a citation (as often used in journals).

Background Information
______________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

The interface is part of ALL which can be found at
https://gitlab.version.fz-juelich.de/SLMS/loadbalancing in the ``src``
subdirectory. It is called ``ALL_module.F90``. Additionally, an internal
C wrapper is used for the C++ class, since Fortran can only interoperate
with C.


Building and Testing
____________________

.. Keep the helper text below around in your module by just adding "..  " in front of it, which turns it into a comment

The Fortran module must be explicitly enabled when building the library.
This is done by setting the CMake variable ``CM_ALL_FORTRAN`` to ``ON``.
The use of the ``mpi_f08`` module can also be enabled with
``CM_ALL_USE_F08``. Then the new MPI derived types can be used directly.
The requisite compilers and MPI implementations must be present. Also
note, that the ALL module must be compiled by the same Fortran compiler as
the application (and MPI implementation if any MPI module is used). More
information is available in the libraries documentation in the code's
repository in ``docs/Install.rst``.


Source Code
___________

.. Notice the syntax of a URL reference below `Text <URL>`_ the backticks matter!

The source code for this interface consists of the C wrapper
`src/ALL_fortran.cpp <https://gitlab.version.fz-juelich.de/SLMS/loadbalancing/-/blob/master/src/ALL_fortran.cpp>`_
and the Fortran module `ALL`
`src/ALL_module.F90 <https://gitlab.version.fz-juelich.de/SLMS/loadbalancing/-/blob/master/src/ALL_module.F90>`_.

.. vim: et sw=2 ts=2 tw=74 spell spelllang=en_us:
