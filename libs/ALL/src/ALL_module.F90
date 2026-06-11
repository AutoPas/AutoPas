!Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
!Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany
!Copyright 2019-2020 Stephan Schulz, Ruhr-Universität Bochum, Germany
!
!Redistribution and use in source and binary forms, with or without modification, are
!permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list
!   of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice, this
!   list of conditions and the following disclaimer in the documentation and/or
!   other materials provided with the distribution.
!
!3. Neither the name of the copyright holder nor the names of its contributors
!   may be used to endorse or promote products derived from this software without
!   specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
!EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
!SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
!BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
!SUCH DAMAGE.
! -DCMAKE_Fortran_FLAGS='-std=f2003 -Wall -Wextra -fbacktrace' -DCM_ALL_FORTRAN=ON -DCMAKE_BUILD_TYPE=Debug -DCM_ALL_DEBUG=ON

!> The function API for ALL
!!
!! The ``*_c`` functions should never be called directly, but their non suffixed
!! counterparts should. Safety checks and conveniences provided by Fortran are
!! available there. For a simple example on how to use this see ``examples/``
!! directory.
module ALL_module
   use iso_c_binding
   use iso_fortran_env
   implicit none
   private
   ! The different loadbalancing methods. Must be kept in sync with C side.
   integer(c_int), public, parameter :: ALL_STAGGERED = 0
   integer(c_int), public, parameter :: ALL_TENSOR = 1
   integer(c_int), public, parameter :: ALL_FORCEBASED = 2
#ifdef ALL_VORONOI_ACTIVE
   integer(c_int), public, parameter :: ALL_VORONOI = 3
#endif
   integer(c_int), public, parameter :: ALL_HISTOGRAM = 4
   integer(c_int), public, parameter :: ALL_TENSOR_MAX = 5
   integer(c_int), public, parameter :: ALL_UNIMPLEMENTED = 128
   ! The different error IDs. Must be kept in sync with CustomException class
   integer(c_int), public, parameter :: ALL_ERROR_GENERIC = 1
   integer(c_int), public, parameter :: ALL_ERROR_POINTDIMENSIONMISSMATCH = 2
   integer(c_int), public, parameter :: ALL_ERROR_INVALIDCOMMTYPE = 3
   integer(c_int), public, parameter :: ALL_ERROR_INVALIDARGUMENT = 4
   integer(c_int), public, parameter :: ALL_ERROR_OUTOFBOUNDS = 5
   integer(c_int), public, parameter :: ALL_ERROR_INTERNAL = 6
   integer(c_int), public, parameter :: ALL_ERROR_FILESYSTEM = 6
   ! Maximum character length of error descriptions
   integer, public, parameter :: ALL_ERROR_LENGTH = 1024
   ! Direct interface with C wrapper
   ! TODO add intent() to dummy arguments
   interface
      function all_init_c(method, dim, gamma) result(this) bind(C)
         use iso_c_binding
         integer(c_int), value :: method
         integer(c_int), value :: dim
         real(c_double), value :: gamma
         type(c_ptr) :: this
      end function
      subroutine all_finalize_c(obj) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
      end subroutine
      subroutine all_set_gamma_c(obj, gamma) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         real(c_double), value :: gamma
      end subroutine
      subroutine all_set_proc_grid_params_c(obj, nloc, loc, nsize, size) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: nloc
         integer(c_int), dimension(nloc) :: loc
         integer(c_int), value :: nsize
         integer(c_int), dimension(nsize) :: size
      end subroutine
      subroutine all_set_proc_tag_c(obj, tag) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: tag
      end subroutine
      subroutine all_set_min_domain_size_c(obj, dim, domain_size) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: dim
         real(c_double), dimension(dim) :: domain_size
      end subroutine
      subroutine all_set_work_c(obj, work) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         real(c_double), value :: work
      end subroutine
      subroutine all_set_work_multi_c(obj, work, dim) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: dim
         real(c_double), dimension(dim) :: work
      end subroutine
      subroutine all_set_vertices_c(obj, n, dim, vertices) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: n
         integer(c_int), value :: dim
         real(c_double), dimension(n*dim) :: vertices
      end subroutine
      subroutine all_set_communicator_c(obj, comm) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer, value :: comm
      end subroutine
      subroutine all_set_sys_size_c(obj, size, dim) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: dim
         real(c_double), dimension(dim) :: size
      end subroutine
      subroutine all_set_method_data_histogram_c(obj, nbins) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), dimension(3) :: nbins
      end subroutine
      subroutine all_setup_c(obj) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
      end subroutine
      subroutine all_balance_c(obj) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
      end subroutine
      subroutine all_get_gamma_c(obj, gamma) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         real(c_double) :: gamma
      end subroutine
      subroutine all_get_number_of_vertices_c(obj, n) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int) :: n
      end subroutine
      subroutine all_get_vertices_c(obj, n, vertices) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: n
         real(c_double), dimension(*) :: vertices
      end subroutine
      subroutine all_get_prev_vertices_c(obj, n, vertices) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: n
         real(c_double), dimension(*) :: vertices
      end subroutine
      subroutine all_get_dimension_c(obj, dim) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int) :: dim
      end subroutine
      subroutine all_get_length_of_work_c(obj, length) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int) :: length
      end subroutine
      subroutine all_get_work_c(obj, work) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         real(c_double) :: work
      end subroutine
      subroutine all_get_work_array_c(obj, work, length) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: length
         real(c_double), dimension(length) :: work
      end subroutine
      subroutine all_get_number_of_neighbors_c(obj, count) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int) :: count
      end subroutine
      subroutine all_get_neighbors_c(obj, neighbors, count) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: count
         integer(c_int), dimension(count) :: neighbors
      end subroutine
#ifdef ALL_VTK_OUTPUT
      subroutine all_print_vtk_outlines_c(obj, step) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: step
      end subroutine
      subroutine all_print_vtk_vertices_c(obj, step) bind(C)
         use iso_c_binding
         type(c_ptr), value :: obj
         integer(c_int), value :: step
      end subroutine
#endif
      function all_errno_c() result(errno) bind(C)
         use iso_c_binding
         integer(c_int) :: errno
      end function
      subroutine all_reset_errno_c() bind(C)
      end subroutine
      subroutine all_errdesc_c(text, length) bind(C)
         use iso_c_binding
         character(len=1, kind=c_char) :: text(*)
         integer(c_size_t), value :: length
      end subroutine
   end interface

   !> The object oriented API is this object. It contains all relevant functions
   type, public :: ALL_t
      type(c_ptr), private :: object = C_NULL_PTR !< pointer to C++ object used on C side
      integer, private :: dim = 0 !< dimensionality of system, set during init
   contains
      procedure :: init => ALL_init
      procedure :: finalize => ALL_finalize
      procedure :: set_gamma => ALL_set_gamma
      procedure :: set_proc_grid_params => ALL_set_proc_grid_params
      procedure :: set_proc_tag => ALL_set_proc_tag
      procedure :: set_min_domain_size => ALL_set_min_domain_size
      procedure :: set_work => ALL_set_work
      procedure :: set_work_multi => ALL_set_work_multi
      procedure :: set_vertices => ALL_set_vertices
#ifdef ALL_USE_F08
      procedure :: set_communicator => ALL_set_communicator_f08
#else
      procedure :: set_communicator => ALL_set_communicator_int
#endif
      procedure :: set_sys_size => ALL_set_sys_size
      procedure :: set_method_data_histgram => ALL_set_method_data_histgram
      procedure :: setup => ALL_setup
      procedure :: balance => ALL_balance
      procedure :: get_gamma => ALL_get_gamma
      procedure :: get_number_of_vertices => ALL_get_number_of_vertices
      procedure :: get_vertices => ALL_get_vertices
      procedure :: get_vertices_alloc => ALL_get_vertices_alloc
      procedure :: get_prev_vertices => ALL_get_prev_vertices
      procedure :: get_dimension => ALL_get_dimension
      procedure :: get_length_of_work => ALL_get_length_of_work
      procedure :: get_work => ALL_get_work
      procedure :: get_work_array => ALL_get_work_array
      procedure :: get_number_of_neighbors => ALL_get_number_of_neighbors
      procedure :: get_neighbors => ALL_get_neighbors
#ifdef ALL_VTK_OUTPUT
      procedure :: print_vtk_outlines => ALL_print_vtk_outlines
      procedure :: print_vtk_vertices => ALL_print_vtk_vertices
#endif
   end type

   ! This is the old, but still supported API of separate functions
   public :: ALL_init
   public :: ALL_finalize
   public :: ALL_set_gamma
   public :: ALL_set_proc_grid_params
   public :: ALL_set_proc_tag
   public :: ALL_set_min_domain_size
   public :: ALL_set_work
   public :: ALL_set_work_multi
   public :: ALL_set_vertices
   public :: ALL_set_communicator
   public :: ALL_set_sys_size
   public :: ALL_set_method_data_histgram
   public :: ALL_setup
   public :: ALL_balance
   public :: ALL_get_gamma
   public :: ALL_get_number_of_vertices
   public :: ALL_get_vertices
   public :: ALL_get_vertices_alloc
   public :: ALL_get_prev_vertices
   public :: ALL_get_dimension
   public :: ALL_get_length_of_work
   public :: ALL_get_work
   public :: ALL_get_work_array
   public :: ALL_get_number_of_neighbors
   public :: ALL_get_neighbors
#ifdef ALL_VTK_OUTPUT
   public :: ALL_print_vtk_outlines
   public :: ALL_print_vtk_vertices
#endif
   public :: ALL_error
   public :: ALL_reset_error
   public :: ALL_error_description

   interface ALL_set_communicator
      module procedure ALL_set_communicator_int
#ifdef ALL_USE_F08
      module procedure ALL_set_communicator_f08
#endif
   end interface
contains
   !> Initialises the loadbalancer
   subroutine ALL_init(this, method, dim, gamma)
      class(ALL_t), intent(out) :: this !< teh ALL object is returned
      integer, intent(in) :: method !< Must be one of the ALL_* method values
      integer, intent(in) :: dim !< dimensionality of system
      real(c_double), intent(in) :: gamma !< gamma value for load balancer (ignored for TENSOR and STAGGERED)
      this%object = all_init_c(method, dim, gamma)
      this%dim = dim
   end subroutine
   !> Delete the loadbalance object
   subroutine ALL_finalize(this)
      class(ALL_t), intent(in) :: this
      call all_finalize_c(this%object)
   end subroutine
   !> Set gamma value for load balancer
   subroutine ALL_set_gamma(this, gamma)
      class(ALL_t), intent(in) :: this
      real(c_double), intent(in) :: gamma
      call all_set_gamma_c(this%object, gamma)
   end subroutine
   !> Set the grid parameters for this process
   subroutine ALL_set_proc_grid_params(this, loc, ranks)
      class(ALL_t), intent(in) :: this
      integer, dimension(this%dim), intent(in) :: loc !< index of domain in `dim` directions (0-indexed!)
      integer, dimension(this%dim), intent(in) :: ranks !< total number of domains in `dim` directions
      call all_set_proc_grid_params_c(this%object, this%dim, loc, this%dim, ranks)
   end subroutine
   !> Set process identifying tag for output
   subroutine ALL_set_proc_tag(this, tag)
      class(ALL_t), intent(in) :: this
      integer, intent(in) :: tag !< tag of local process, only output in VTK outlines
      call all_set_proc_tag_c(this%object, tag)
   end subroutine
   !> Set a minimum domain size
   subroutine ALL_set_min_domain_size(this, domain_size)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(this%dim), intent(in) :: domain_size !< minimum domain size
      call all_set_min_domain_size_c(this%object, this%dim, domain_size)
   end subroutine
   !> Set the work of this process
   subroutine ALL_set_work(this, work)
      class(ALL_t), intent(in) :: this
      real(c_double), intent(in) :: work !< work of this domain
      call all_set_work_c(this%object, work)
   end subroutine
   !> Set multi dimensional work of this process
   subroutine ALL_set_work_multi(this, work)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:), intent(in) :: work !< multi dimensional work of this domain
#ifndef NDEBUG
      if (size(work) /= this%dim) then
         write (error_unit, '(a)') "ALL_set_work_multi: Wrong dimension for work"
         stop
      end if
#endif
      call all_set_work_multi_c(this%object, work, size(work))
   end subroutine
   !> Set new vertices
    !!
    !! @todo write interface for single rank array of vertices
   subroutine ALL_set_vertices(this, vertices)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:, :), intent(in) :: vertices !< vertices of domain, for `n` domains must have shape: (dim,n)
#ifndef NDEBUG
      if (size(vertices, 1) /= this%dim) then
         write (error_unit, '(a)') "ALL_set_vertices: Wrong dimension for vertices"
         stop
      end if
#endif
      call all_set_vertices_c(this%object, size(vertices, 2), this%dim, vertices)
   end subroutine
   !> Set the MPI communicator with an mpi_f08 object
#ifdef ALL_USE_F08
   subroutine ALL_set_communicator_f08(this, comm)
      use mpi_f08
      class(ALL_t), intent(in) :: this
      type(MPI_Comm), intent(in) :: comm !< MPI Communicator
      call all_set_communicator_c(this%object, comm%mpi_val)
   end subroutine
#endif
   !> Set the MPI communicator with an ``mpi`` oder ``mpif.h`` communicator
   subroutine ALL_set_communicator_int(this, comm)
      class(ALL_t), intent(in) :: this
      integer, intent(in) :: comm !< MPI Communicator, not type checked!
      call all_set_communicator_c(this%object, comm)
   end subroutine
   !> Set size of system, which is required for the histogram method
   subroutine ALL_set_sys_size(this, syssize)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:), intent(in) :: syssize !< system size
#ifndef NDEBUG
      if (size(syssize) /= this%dim) then
         write (error_unit, '(a)') "ALL_set_size: Wrong dimension for size"
         stop
      end if
#endif
      call all_set_sys_size_c(this%object, syssize, size(syssize))
   end subroutine
   !> Set number of bins for histogram method
   subroutine ALL_set_method_data_histgram(this, nbins)
      class(ALL_t), intent(in) :: this
      integer(c_int), dimension(3), intent(in) :: nbins !< Number of bins per dimension
      call all_set_method_data_histogram_c(this%object, nbins)
   end subroutine
   !> Set up the loadbalancer
   subroutine ALL_setup(this)
      class(ALL_t), intent(in) :: this
      call all_setup_c(this%object)
   end subroutine
   !> Run loadbalancer and calculate new vertices
   subroutine ALL_balance(this)
      class(ALL_t), intent(in) :: this
      call all_balance_c(this%object)
   end subroutine
   !> Retrieve currently set gamma value of balancer
   subroutine ALL_get_gamma(this, gamma)
      class(ALL_t), intent(in) :: this
      real(c_double), intent(out) :: gamma
      call all_get_gamma_c(this%object, gamma)
   end subroutine
   !> Retrieve number of vertices for current vertices
   subroutine ALL_get_number_of_vertices(this, n)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(out) :: n !< set to number of new vertices
      call all_get_number_of_vertices_c(this%object, n)
   end subroutine
   !> Retrieve current vertices
   subroutine ALL_get_vertices(this, vertices)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:, :), intent(out) :: vertices !< set to new vertices, must be exact size (dim,n)
      integer(c_int) :: n_verts
      call this%get_number_of_vertices(n_verts)
#ifndef NDEBUG
      if (size(vertices, 1) /= this%dim .or. size(vertices, 2) < n_verts) then
         write (error_unit, '(a)') "ALL_get_vertices: vertices array not large enough!"
         stop
      end if
#endif
      call all_get_vertices_c(this%object, n_verts, vertices)
   end subroutine
   !> Same function as get_vertices, but takes an allocatable array, and resizes automatically
   subroutine ALL_get_vertices_alloc(this, vertices)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:, :), allocatable, intent(inout) :: vertices !< set to new vertices, may be reallocated to fit
      integer(c_int) :: n_verts
      call this%get_number_of_vertices(n_verts)
      if (allocated(vertices)) then
         if (size(vertices, 1) /= this%dim .or. size(vertices, 2) < n_verts) then
            deallocate (vertices)
         end if
      end if
      if (.not. allocated(vertices)) allocate (vertices(this%dim, n_verts))
      call all_get_vertices_c(this%object, n_verts, vertices)
   end subroutine
   !> Retrieve vertices from before loadbalancing
   subroutine ALL_get_prev_vertices(this, vertices)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:, :), intent(out) :: vertices !< set to prev vertices, must be exact size (dim,n), unchecked!
#ifndef NDEBUG
      ! TODO Check against number of vertices not only dimensionality
      if (size(vertices, 1) /= this%dim) then
         write (error_unit, '(a)') "ALL_get_prev_vertices: vertices array not large enough!"
         stop
      end if
#endif
      call all_get_prev_vertices_c(this%object, size(vertices, 2), vertices)
   end subroutine
   !> Get current dimension from loadbalancer
   subroutine ALL_get_dimension(this, dim)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(out) :: dim
      call all_get_dimension_c(this%object, dim)
   end subroutine
   !> Retrieve length of work array
   subroutine ALL_get_length_of_work(this, length)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(out) :: length
      call all_get_length_of_work_c(this%object, length)
   end subroutine
   !> Retrieve first element of work array
   subroutine ALL_get_work(this, work)
      class(ALL_t), intent(in) :: this
      real(c_double), intent(out) :: work
      call all_get_work_c(this%object, work)
   end subroutine
   !> Retrieve work array, which must already be the correct size
   subroutine ALL_get_work_array(this, work)
      class(ALL_t), intent(in) :: this
      real(c_double), dimension(:), intent(out) :: work
      integer :: length
      call ALL_get_length_of_work(this, length)
      if (size(work) /= length) then
         write (error_unit, '(a)') "ALL_get_work_array: work has wrong length!"
         stop
      end if
      call all_get_work_array_c(this%object, work, size(work))
   end subroutine
   !> Retrieve number of neighbors (i.e. length of neighbors list)
   subroutine ALL_get_number_of_neighbors(this, count)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(out) :: count
      call all_get_number_of_vertices_c(this%object, count)
   end subroutine
   !> Retrieve list of neighboring ranks (must be correct size already)
   subroutine ALL_get_neighbors(this, neighbors)
      class(ALL_t), intent(in) :: this
      integer(c_int), dimension(:), intent(out) :: neighbors
      integer :: count
      call ALL_get_number_of_neighbors(this, count)
      if (size(neighbors) /= count) then
         write (error_unit, '(a)') "ALL_get_neighbors: neighbors has wrong length!"
         stop
      end if
      call all_get_neighbors_c(this%object, neighbors, size(neighbors))
   end subroutine

#ifdef ALL_VTK_OUTPUT
   !> Print VTK outlines (must be enabled in build step)
   subroutine ALL_print_vtk_outlines(this, step)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(in) :: step !< current step, used for filename
      call all_print_vtk_outlines_c(this%object, step)
   end subroutine
   !> Print VTK domain vertices (must be enabled in build step)
   subroutine ALL_print_vtk_vertices(this, step)
      class(ALL_t), intent(in) :: this
      integer(c_int), intent(in) :: step !< current step, used for filename
      call all_print_vtk_vertices_c(this%object, step)
   end subroutine
#endif
   function ALL_error() result(error)
      integer :: error
      error = all_errno_c()
   end function
   subroutine ALL_reset_error()
      call all_reset_errno_c()
   end subroutine
   function ALL_error_description() result(desc)
      character(len=ALL_ERROR_LENGTH) :: desc
      call all_errdesc_c(desc, int(ALL_ERROR_LENGTH, c_size_t))
   end function
end module

! VIM modeline for indentation
! vim: sw=4 et
