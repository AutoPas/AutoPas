!Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
!Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany
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
PROGRAM ALL_test_f
    USE ALL_module
    USE ISO_C_BINDING
    use mpi
    IMPLICIT NONE

    INTEGER                             ::  cart_comm
    INTEGER                             ::  dims(3)
    INTEGER                             ::  coords(3)
    LOGICAL                             ::  period(3)
    REAL(8)                             ::  length(3)
    REAL(8)                             ::  vertices(3,2)
    INTEGER                             ::  error 
    INTEGER                             ::  rank, n_ranks
    INTEGER                             ::  n_vertices
    REAL(8),DIMENSION(:,:),ALLOCATABLE  ::  new_vertices
    INTEGER                             ::  i 
    CHARACTER(LEN=ALL_ERROR_LENGTH)     ::  all_error_desc

    TYPE(ALL_t)     ::  obj

    CALL MPI_INIT(error)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,n_ranks,error)

    ! create cartesian communicator
    dims = 0
    period = .true.
    CALL MPI_DIMS_CREATE(n_ranks,3,dims, error) 
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, period, .true., cart_comm, error)

    ! compute cell length
    length = 1.0d0 / real(dims,8)

    ! compute vertices of local domain
    CALL MPI_CART_COORDS(cart_comm, rank, 3, coords, error)
    vertices(:,1) = real(coords,8) * length
    vertices(:,2) = real(coords+1,8) * length

    DO i = 0, n_ranks-1
        IF (i == rank) THEN
            WRITE(*,"(a,i7,2(a,3es14.7))") "rank: ", rank, " old vertices: ", vertices(:,1), ", ", vertices(:,2)
            call MPI_BARRIER(cart_comm, error)
        ELSE
            call MPI_BARRIER(cart_comm, error)
        END IF
    END DO

    CALL ALL_reset_error()
    CALL ALL_init(obj,ALL_STAGGERED,3,4.0d0)
    IF (ALL_error() /= 0) THEN
        print*, "Error: ", ALL_error()
        all_error_desc= ALL_error_description()
        print*, "Message: ", all_error_desc
        call MPI_Abort(MPI_COMM_WORLD, 1, error)
        stop
    END IF
    CALL ALL_set_proc_tag(obj, rank)
    CALL ALL_set_work(obj, real( product(coords,1)*64,8) )
    CALL ALL_set_vertices(obj,vertices)
    ! if using a non cartesian communicator, this call would be required
    !CALL ALL_set_proc_grid_params(obj,coords,dims)
    CALL ALL_set_communicator(obj,cart_comm)

    CALL ALL_setup(obj)
    CALL ALL_balance(obj)

    CALL ALL_get_number_of_vertices(obj,n_vertices)
    ALLOCATE(new_vertices(3,n_vertices))
    CALL ALL_get_vertices(obj,new_vertices)

    DO i = 0, n_ranks-1
        IF (i == rank) THEN
            WRITE(*,"(2(a,i7),2(a,3es14.7))") "rank: ", rank, " ", n_vertices, " new vertices: ", new_vertices(:,1), &
                ", ", new_vertices(:,2)
            call MPI_BARRIER(cart_comm, error)
        ELSE
            call MPI_BARRIER(cart_comm, error)
        END IF
    END DO
    DEALLOCATE(new_vertices)
    CALL ALL_finalize(obj);
    CALL MPI_FINALIZE(error);

END PROGRAM

