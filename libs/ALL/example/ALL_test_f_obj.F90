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

! module for ALL access from Fortran
program ALL_test_f
    use ALL_module
    use iso_c_binding
#ifndef ALL_USE_F08
    use mpi
#else
    use mpi_f08
#endif
    implicit none

#ifndef ALL_USE_F08
    integer                             ::  cart_comm
#else
    type(MPI_Comm)                      ::  cart_comm
#endif
    integer                             ::  dims(3)
    integer                             ::  coords(3)
    logical                             ::  period(3)
    real(8)                             ::  length(3)
    real(8)                             ::  vertices(3,2)
    integer                             ::  error
    integer                             ::  rank, n_ranks
    real(8),dimension(:,:),allocatable  ::  new_vertices
    integer                             ::  i

    type(ALL_t)     ::  lb

    call MPI_INIT(error)

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,n_ranks,error)

    ! create cartesian communicator
    dims = 0
    period = .true.
    call MPI_DIMS_CREATE(n_ranks,3,dims, error)
    call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, period, .true., cart_comm, error)

    ! compute cell length
    length = 1.0d0 / real(dims,8)

    ! compute vertices of local domain
    call MPI_CART_COORDS(cart_comm, rank, 3, coords, error)
    vertices(:,1) = real(coords,8) * length
    vertices(:,2) = real(coords+1,8) * length

    do i = 0, n_ranks-1
        if (i == rank) then
            write(*,"(a,i7,2(a,3es14.7))") "rank: ", rank, " old vertices: ", vertices(:,1), ", ", vertices(:,2)
            call MPI_BARRIER(cart_comm, error)
        else
            call MPI_BARRIER(cart_comm, error)
        end if
    end do

    call lb%init(ALL_STAGGERED, 3, 4.0d0)
    call lb%set_proc_tag(rank)
    call lb%set_work(real( product(coords,1)*64,8))
    call lb%set_vertices(vertices)
    ! if using a non cartesian communicator, this call would be required
    ! call lb%set_proc_grid_params(coords,dims)
    call lb%set_communicator(cart_comm)

    call lb%setup()
    call lb%balance()

    ! If not allocatable as target, must use
    ! lb%get_number_of_vertices and lb%get_vertices instead!
    call lb%get_vertices_alloc(new_vertices)

    do i = 0, n_ranks-1
        if (i == rank) then
            write(*,"(2(a,i7),2(a,3es14.7))") "rank: ", rank, " ", size(new_vertices,2), " new vertices: ", new_vertices(:,1), &
                ", ", new_vertices(:,2)
            call MPI_BARRIER(cart_comm, error)
        else
            call MPI_BARRIER(cart_comm, error)
        end if
    end do
    deallocate(new_vertices)
    call lb%finalize()
    call MPI_FINALIZE(error);

END PROGRAM

