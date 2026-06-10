!Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
!Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany
!Copyright 2020-2020 Stephan Schulz, Forschungszentrum Juelich GmbH, Germany
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

program ALL_Staggered_f
   use ALL_module
   use iso_c_binding
   use, intrinsic :: iso_fortran_env, only: stdin=>input_unit&
      , stdout=>output_unit&
      , stderr=>error_unit&
      , file_storage_size
#ifndef ALL_USE_F08
   use mpi
#else
   use mpi_f08
#endif
   implicit none
   integer, parameter :: real64 = selected_real_kind(15)

   call main()
contains
   subroutine print_loc(my_rank, maximum_rank, my_location, number_of_processes)
      integer, intent(in) :: my_rank, maximum_rank
      integer, dimension(3), intent(in) :: my_location, number_of_processes
      integer :: rank, error
      do rank=0, maximum_rank-1
         if(rank == my_rank) then
            write(stdout,'(a,i3.3,a,i3,",",i3,",",i3,a,i3,",",i3,",",i3,a)')&
               "[", my_rank, "] Location: (", my_location, ")/(", number_of_processes, ")"
            flush(stdout)
            call MPI_Barrier(MPI_COMM_WORLD, error)
         else
            call MPI_Barrier(MPI_COMM_WORLD, error)
         endif
      enddo
   end subroutine

   subroutine print_domain(my_rank, maximum_rank, domain_vertices)
      integer, intent(in) :: my_rank, maximum_rank
      real(real64), dimension(3,2), intent(in) :: domain_vertices
      integer :: rank, error
      do rank=0, maximum_rank-1
         if(rank == my_rank) then
            write(stdout,'("[",i3.3,"]",a,es12.4,a,es12.4,a,es12.4)')&
               my_rank, " Lower: ",&
               domain_vertices(1,1), achar(9),&
               domain_vertices(2,1), achar(9),&
               domain_vertices(3,1)
            write(stdout,'("[",i3.3,"]",a,es12.4,a,es12.4,a,es12.4)')&
               my_rank, " Upper: ",&
               domain_vertices(1,2), achar(9),&
               domain_vertices(2,2), achar(9),&
               domain_vertices(3,2)
            flush(stdout)
            call MPI_Barrier(MPI_COMM_WORLD,error)
         else
            call MPI_Barrier(MPI_COMM_WORLD,error)
         endif
      enddo
   end subroutine

   subroutine print_work(my_rank, maximum_rank, my_work)
      integer, intent(in) :: my_rank, maximum_rank
      real(real64), intent(in) :: my_work
      integer :: rank, error
      do rank=0, maximum_rank-1
         if(rank == my_rank) then
            write(stdout,'(a,i3.3,a,es12.4)')&
               "[", my_rank, "] Work: ", my_work
            flush(stdout)
            call MPI_Barrier(MPI_COMM_WORLD, error)
         else
            call MPI_Barrier(MPI_COMM_WORLD, error)
         endif
      enddo
   end subroutine

   subroutine print_testing_output(my_rank, maximum_rank, new_vertices, timestep)
      integer, intent(in) :: my_rank, maximum_rank, timestep
      real(real64), dimension(3,2), intent(in) :: new_vertices
      integer :: rank, error
      do rank=0, maximum_rank-1
         if(rank == my_rank) then
            !write(stdout,'(a,i4,",",i3.3,a,3f11.6)')&
            !    "[", timestep, my_rank, "] Result Width: ",&
            !    new_vertices(:,2)-new_vertices(:,1)
            !flush(stdout)
            write(stdout,'(a,i4,",",i3.3,",",i2.2,a,3f11.6)')&
               "[", timestep, my_rank, 0, "] Result Vertex: ",&
               new_vertices(:,1)
            flush(stdout)
            write(stdout,'(a,i4,",",i3.3,",",i2.2,a,3f11.6)')&
               "[", timestep, my_rank, 1, "] Result Vertex: ",&
               new_vertices(:,2)
            flush(stdout)
            call MPI_Barrier(MPI_COMM_WORLD, error)
         else
            call MPI_Barrier(MPI_COMM_WORLD, error)
         endif
      enddo
   end subroutine

   subroutine print_binary(my_rank, my_work, new_vertices, my_location, number_of_processes, fh, timestep)
      use iso_fortran_env, only: ERROR_UNIT

#ifndef ALL_USE_F08
      integer, intent(in) :: fh
      integer, dimension(MPI_STATUS_SIZE) :: state
#else
      type(MPI_File), intent(in) :: fh
      type(MPI_Status) :: state
#endif
      integer, intent(in) :: my_rank, timestep
      integer, dimension(3), intent(in) :: my_location, number_of_processes
      real(real64), intent(in) :: my_work
      real(real64), dimension(3,2), intent(in) :: new_vertices
      integer(MPI_OFFSET_KIND) :: offset
      real(real64), dimension(11) :: buffer
      integer :: error

      offset = ((timestep-1) * number_of_processes(3) * 11 + my_rank * 11) * 8

      buffer(1) = real(my_rank,real64)
      buffer(2) = my_work
      buffer(3) = new_vertices(1,1)
      buffer(4) = new_vertices(2,1)
      buffer(5) = new_vertices(3,1)
      buffer(6) = new_vertices(1,2)
      buffer(7) = new_vertices(2,2)
      buffer(8) = new_vertices(3,2)
      buffer(9) = my_location(1)
      buffer(10) = my_location(2)
      buffer(11) = my_location(3)

      call MPI_FILE_WRITE_AT(fh, offset, buffer, 11, MPI_DOUBLE_PRECISION, state, error);
      !write(ERROR_UNIT,"(11f10.3,a,i0,a,i0)") buffer," ",offset, " ", timestep

   end subroutine

   subroutine main()
      integer :: error
      integer :: current_step
      integer, parameter :: number_of_steps = 50
      integer, parameter :: dimensions = 3
      real(real64), parameter :: loadbalancer_gamma = 0. ! ignored for staggered method
      integer, dimension(dimensions) :: my_location, number_of_processes
      real(real64), dimension(dimensions) :: minimum_domain_size
      real(real64), dimension(dimensions,2) :: domain_vertices, new_vertices
      real(real64), parameter :: domain_size = 1.0
      real(real64) :: my_work
      integer :: my_rank, maximum_rank
      integer :: i,j
      type(ALL_t) :: jall

      ! number of command line arguments
      integer :: n_clargs
      ! command line arguments
      character(256), dimension(:), allocatable :: clargs
      character(256) :: filename, filename2
      integer(MPI_OFFSET_KIND) :: offset
      integer :: test_file
      real(real64) :: d
#ifndef ALL_USE_F08
      integer :: fh
      integer, dimension(MPI_STATUS_SIZE) :: state
#else
      type(MPI_File) :: fh
      type(MPI_Status) :: state
#endif

#ifdef ALL_VTK_OUTPUT_EXAMPLE
      character (len=ALL_ERROR_LENGTH) :: error_msg
#endif

      call MPI_Init(error)

      ! The ALL_TENSOR method can be used just like the staggered method.
      ! Just exchange ALL_STAGGERED with ALL_TENSOR in the init call.
      call jall%init(ALL_STAGGERED, dimensions, loadbalancer_gamma)
      my_location(:) = 0
      ! All domains are placed along a line in z direction, evne though they are three dimensional
      call MPI_Comm_rank(MPI_COMM_WORLD, my_location(3), error)
      my_rank = my_location(3)

      number_of_processes(:) = 1
      call MPI_Comm_size(MPI_COMM_WORLD, number_of_processes(3), error)
      maximum_rank = number_of_processes(3)

      n_clargs = command_argument_count()
      if (n_clargs > 1) then
         call get_command_argument(2,clargs(2))
         write(filename,"(a,i0,a)") trim(clargs(2)),number_of_processes,trim(".bin")
      else
         write(filename,"(a,i0,a)") trim("./ALL_Staggered_f_"),number_of_processes(3),trim(".bin")
      endif
      call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,IOR(MPI_MODE_CREATE, MPI_MODE_WRONLY), MPI_INFO_NULL, fh, error)

      if(my_rank == 0) then
         write(stdout,'(a,i0)') "Ranks: ", maximum_rank
         write(stdout,'(a,i0)') "Number of Steps: ", number_of_steps
         flush(stdout)
      endif

      call MPI_Barrier(MPI_COMM_WORLD, error)
      call print_loc(my_rank, maximum_rank, my_location, number_of_processes)

      ! For a cartesian communicator this is not required, but we are using
      ! MPI_COMM_WORLD here.
      call jall%set_proc_grid_params(my_location, number_of_processes)

      ! If we want ot set a minimum domain size:
      minimum_domain_size(:) = 0.1
      call jall%set_min_domain_size(minimum_domain_size)

      call jall%set_communicator(MPI_COMM_WORLD)

      ! We also set the optional process tag for the output.
      ! This can be useful if we want to know which of 'our'
      ! ranks is which in the output produces by the library.
      ! The ranks used inside the library do not necessarily
      ! match our own numbering.
      call jall%set_proc_tag(my_rank)

      call jall%setup()

      ! A first domain distribution must be given to the balancer.
      ! We use the provided ALL::Point class to define the vertices,
      ! but a simple double array can also be used. We need 2 vertices
      ! which correspond to the two opposing corners.
      ! We create a cubic domain initially
      do i=1, 2
         do j=1,dimensions
            domain_vertices(j,i) = (my_location(j)+i-1) * domain_size
         enddo
      enddo
      call print_domain(my_rank, maximum_rank, domain_vertices)
      call jall%set_vertices(domain_vertices)

      ! Calculate the work fo our domain. Here we just use
      my_work = my_rank + 1.
      call jall%set_work(my_work)
      call print_work(my_rank, maximum_rank, my_work)
      do current_step=1, number_of_steps
         ! In a real code we need to set the updated work in each
         ! iteration before calling balance()
         if(my_rank == 0) then
            write(stdout,'(a,i0,"/",i0)') "Starting step: ", current_step, number_of_steps
            flush(stdout)
         endif
#ifdef ALL_VTK_OUTPUT_EXAMPLE
         call ALL_reset_error()
         call jall%print_vtk_outlines(current_step)
         if(ALL_error() /= 0) then
            error_msg = ALL_error_description()
            write(stdout, '(a)') error_msg
            ! Maybe also abort if the output is necesssary or handle this
            ! some other way.
            call MPI_Abort(MPI_COMM_WORLD, 2, error)
            stop
         endif
#endif
         call jall%balance()

         call jall%get_vertices(new_vertices)

         !call print_testing_output(my_rank, maximum_rank, new_vertices, current_step)
         call print_binary(my_rank, my_work, new_vertices, my_location, number_of_processes, fh, current_step)

         call MPI_Barrier(MPI_COMM_WORLD, error)
      enddo
#ifdef ALL_VTK_OUTPUT_EXAMPLE
      call ALL_reset_error()
      call jall%print_vtk_outlines(current_step)
      if(ALL_error() /= 0) then
         error_msg = ALL_error_description()
         write(stdout, '(a)') error_msg
         ! Maybe also abort if the output is necesssary or handle this
         ! some other way.
         call MPI_Abort(MPI_COMM_WORLD, 2, error)
         stop
      endif
#endif

      call jall%finalize()
      call MPI_FILE_CLOSE(fh, error)

      call MPI_Finalize(error)
   end subroutine
end program
