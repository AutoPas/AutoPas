! The following additional functions were used:
!
! Additional information:
! - LB_METHOD_PRE is set by the preprocessor to ALL_STAGGERED or ALL_TENSOR.
! - The domain_bounds_old will only be used during initialisation for the
!   initial domain configuration.
! - ``this_image()`` returns the current image index, i.e. current MPI rank+1.
! - The work ist estimated using ``lb_estimate_work`` which takes the current
!   domain size and number of particles as arguments.

    function domain_decomposition_jall(bounds, dh, num_images3, domain_bounds_old, work, output) result(domain_bounds)
        use ISO_C_BINDING
        type(boundingbox), intent(in) :: bounds !< simulation bounds
        real (kind = real_kind), intent(in) :: dh !< grid width
        integer, dimension(3), intent(in) :: num_images3 !< the 1 indexed number of images in 3D
        type(boundingbox_aligned), intent(in) :: domain_bounds_old !< current domain bounds
        real(real_kind), intent(in) :: work !< work of this domain
        logical, intent(in) :: output !< output domain bounds to `vtk_outline` directory
        type(boundingbox_aligned) :: domain_bounds
        type(boundingbox), save :: domain_old
        type(ALL_t), save :: jall ! ALL object which is initialized once
        real (kind = real_kind), dimension(3,2) :: verts
        integer, dimension(3), save :: this_image3 ! the 1 indexed image number in 3D
        logical, save :: first_run = .true.
        integer, save :: step = 1
        logical, dimension(2,3), save :: domain_at_sim_bound ! true if domain is at the lower/upper simulation boundary
        real (kind = real_kind), dimension(3), save :: min_size
        integer(c_int), parameter :: LB_METHOD = LB_METHOD_PRE
        character (len=ALL_ERROR_LENGTH) :: error_msg
        if(first_run) then
            ! calculate this_image3
            block
                integer :: x,y,z, cnt
                cnt = 1
                do z=1,num_images3(3)
                    do y=1,num_images3(2)
                        do x=1,num_images3(1)
                            if(this_image()==cnt) this_image3 = (/ x,y,z /)
                            cnt = cnt + 1
                        enddo
                    enddo
                enddo
            end block
            call jall%init(LB_METHOD,3,4.0d0)
            call jall%set_proc_grid_params(this_image3-1, num_images3)
            call jall%set_proc_tag(this_image())
            call jall%set_communicator(MPI_COMM_WORLD)
            min_size(:) = (abs(Rcont_min)+abs(Rcont_max))*dh
            call jall%set_min_domain_size(min_size)
            domain_old%bounds_unaligned = domain_bounds_old%bounds_aligned
            domain_at_sim_bound(1,:) = this_image3==1 ! first image in a direction is automatically at sim bound
            domain_at_sim_bound(2,:) = this_image3==num_images3 ! last image likewise at sim bound
            call jall%setup()
        endif
        call jall%set_work(real(work,real_kind))
        !! The `domain_old` bounds are not the actual domain bounds, which
        !! are aligned to grid widths, but what we got from the previous
        !! iteration of load balancing. However, the simulation boundaries are
        !! unchanged by the load balancing.
        block
            type(boundingbox_aligned) :: aligned_bnds
            real (kind = real_kind), dimension(3) :: lowest_upper_bound, highest_lower_bound
            !> Align the simulation boundaries to the grid and add an additional
            !! grid width on the top. These may be used instead of our current
            !! bounds, so they should align properly on the upper bound, if we
            !! are a simulation boundary. If the simulation bounds have not
            !! changed they should still coincide with the domain bounds.
            aligned_bnds%bounds_aligned = floor(bounds%bounds_unaligned/dh)*dh
            aligned_bnds%bounds_aligned(2,:) = aligned_bnds%bounds_aligned(2,:) + dh
            !> To make sure, the shrinking domain is still always large enough
            !! and in particular is not shrunk into the neighbouring domain.
            !! This can happen if the bounding box is not present in the current
            !! domain, so the outer bound is moved across the inner bound. This
            !! must be avoided at all cost. Additionally, we also need to ensure
            !! the minimum domain width. Also, the outer bound of all boundary
            !! domains, must be the same. To achieve this, the outermost inner
            !! bound is calculated in each direction. This then allows us to
            !! compute the innermost position any outer bound may have to still
            !! be the required distance from every next inner bound.
            ! For the lowest domains:
            lowest_upper_bound = comm_co_min_f(domain_old%bounds_unaligned(2,:))
            aligned_bnds%bounds_aligned(1,:) = min(lowest_upper_bound-min_size, aligned_bnds%bounds_aligned(1,:))
            ! For the highest domains:
            highest_lower_bound = comm_co_max_f(domain_old%bounds_unaligned(1,:))
            aligned_bnds%bounds_aligned(2,:) = max(highest_lower_bound+min_size, aligned_bnds%bounds_aligned(2,:))
            ! And now set the boundary domains outer bounds to the new, fixed bounds
            where(domain_at_sim_bound)
                domain_old%bounds_unaligned = aligned_bnds%bounds_aligned
            end where
        end block
        !> Make sure that the old domain bounds are sensible. we are only
        !! updating them, based in the previous value. This also means
        !! the first call must already contain a sensible approximation
        !! (the equidistant (geometric) distribution suffices for that).
        verts = transpose(domain_old%bounds_unaligned)
        call jall%set_vertices(verts)
        call jall%balance()
        call jall%get_result_vertices(verts)
        domain_bounds%bounds_aligned = transpose(verts)
        domain_old%bounds_unaligned = domain_bounds%bounds_aligned
        domain_bounds%bounds_aligned = nint(domain_bounds%bounds_aligned/dh)*dh
        if(output) then
            call ALL_reset_error()
            call jall%print_vtk_outlines(step)
            if(ALL_error() /= 0) then
                error_msg = ALL_error_description()
                print*, "Error in ALL detected:"
                print*, error_msg
            endif
        endif
        first_run = .false.
        step = step + 1
        call assert_domain_width(domain_bounds, dh)
    end function

    !> Estimate local work
    function lb_estimate_work(n_part, domain_bounds_old) result(work)
        integer, intent(in) :: n_part !< number of particles of this domain
        type(boundingbox_aligned), intent(in) :: domain_bounds_old !< domain bounds
        real(real_kind) :: work
        real(real_kind), parameter :: beta = 0.128 ! empirically determined
        work = n_part + beta*product( domain_bounds_old%bounds_aligned(2,:)-domain_bounds_old%bounds_aligned(1,:) )/grid%dh**3
    end function
