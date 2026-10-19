module fortran_interface
    use iso_c_binding
    implicit none

    type LoopData
        integer :: n                    ! #iterations (for reinforcement learning agent)
        integer :: p                    ! #threads (for reinforcement learning agent)
        integer :: autoSearch           ! flag: continue automatic search or not
        integer :: cDLS                 ! current DLS
        integer :: bestDLS              ! Best DLS
        integer :: searchTrials         ! number of trials to find the best DLS
        integer :: cChunk               ! current chunk size
        integer :: bestChunk            ! chunk size of the best DLS technique
        double precision :: cTime       ! current DLS time
        double precision :: bestTime    ! loop time of the best DLS
        double precision :: cLB         ! load imbalance of the current DLS
        double precision :: bestLB      ! load imbalance of the best DLS
        double precision :: cRobust     ! current Robustness
        double precision :: cStdDev     ! standard deviation of the current DLS
        double precision :: bestStdDev  ! standard deviation of the best DLS
        double precision :: cSkew       ! skewness of the current DLS
        double precision :: bestSkew    ! skewness of the best DLS
        double precision :: cKurt       ! kurtosis of the current DLS
        double precision :: bestKurt    ! kurtosis of the best DLS
        double precision :: cCOV        ! coefficient of variation of the current DLS
        double precision :: bestCOV     ! coefficient of variation of the best DLS
    end type LoopData

    interface
        function my_cpp_function(arg) bind(C, name = "my_cpp_function")
            use iso_c_binding
            integer(c_int), value :: arg
            integer(c_int) :: my_cpp_function
        end function my_cpp_function

        integer(c_int) function search(loop_id, agent_type, info, dimension) bind(c, name = "search")
            use iso_c_binding, only : c_int
            import :: LoopData
            character(len = 1), intent(in) :: loop_id
            integer, value :: agent_type
            type (LoopData), intent (inout) :: info
            integer, value :: dimension
        end function search

    end interface
end module fortran_interface