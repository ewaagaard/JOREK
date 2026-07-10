!*******************************************************************************
!* Module: mod_boundary_ndotB                                                  *
!*                                                                             *
!* Store and retrieve n.B/|B| values per boundary node for stellarator SBC.    *
!* Uses MPI_Allreduce to share values across all ranks.                        *
!* Supports toroidal variation: stores n.B per (node, plane) and provides      *
!* Fourier-transformed coefficients for each toroidal harmonic.                *
!*******************************************************************************

module mod_boundary_ndotB

  implicit none
  private
  
  public :: init_boundary_ndotB
  public :: accumulate_ndotB_at_node
  public :: accumulate_ndotB_at_node_plane
  public :: finalize_boundary_ndotB
  public :: get_ndotB_at_node
  public :: get_ndotB_fourier_at_node
  public :: get_vpar_target_for_column
  
  ! Per-node storage (toroidally averaged, for backward compatibility)
  real*8, allocatable, save :: ndotB_per_node(:)
  real*8, allocatable, save :: ndotB_count_per_node(:)
  
  ! Per-node-plane storage (for toroidal variation)
  real*8, allocatable, save :: ndotB_per_node_plane(:,:)
  real*8, allocatable, save :: ndotB_count_per_node_plane(:,:)
  
  ! Fourier coefficients (computed in finalize)
  real*8, allocatable, save :: ndotB_fourier_cos(:,:)  ! (n_nodes, n_tor)
  real*8, allocatable, save :: ndotB_fourier_sin(:,:)  ! (n_nodes, n_tor)
  
  ! Fourier coefficients of vpar_target (with tanh smoothing applied in physical space)
  real*8, allocatable, save :: vpar_target_fourier_cos(:,:)  ! (n_nodes, n_tor)
  real*8, allocatable, save :: vpar_target_fourier_sin(:,:)  ! (n_nodes, n_tor)
  
  logical, save :: ndotB_initialized = .false.
  logical, save :: ndotB_finalized = .false.
  logical, save :: ndotB_computed_once = .false.  ! Track if ndotB has been computed at least once
  
  integer, save :: n_nodes_stored = 0
  integer, save :: n_plane_stored = 0
  integer, save :: n_tor_stored = 0

  public :: is_ndotB_frozen  ! Query function for frozen state

contains

!---------------------------------------------------------------------------
! Check if ndotB is frozen (already computed and ndotB_evolving = .false.)
!---------------------------------------------------------------------------
function is_ndotB_frozen() result(frozen)
  use phys_module, only: ndotB_evolving
  implicit none
  logical :: frozen
  
  frozen = ndotB_computed_once .and. (.not. ndotB_evolving)
end function is_ndotB_frozen

subroutine init_boundary_ndotB(n_nodes, n_plane_in, n_tor_in)
  use phys_module, only: ndotB_evolving
  implicit none
  integer, intent(in) :: n_nodes
  integer, intent(in), optional :: n_plane_in, n_tor_in
  integer :: n_plane_local, n_tor_local
  
  ! If ndotB is frozen (already computed once and not evolving), skip reinitialization
  ! This preserves the stored values from the first timestep
  if (ndotB_computed_once .and. (.not. ndotB_evolving)) then
    ! Just reset the finalized flag to allow boundary_conditions to use stored values
    ndotB_finalized = .true.  ! Keep finalized state (values already valid)
    return
  endif
  
  ! Handle optional arguments
  n_plane_local = 1
  n_tor_local = 1
  if (present(n_plane_in)) n_plane_local = n_plane_in
  if (present(n_tor_in)) n_tor_local = n_tor_in
  
  ! Deallocate if already allocated
  if (allocated(ndotB_per_node)) deallocate(ndotB_per_node)
  if (allocated(ndotB_count_per_node)) deallocate(ndotB_count_per_node)
  if (allocated(ndotB_per_node_plane)) deallocate(ndotB_per_node_plane)
  if (allocated(ndotB_count_per_node_plane)) deallocate(ndotB_count_per_node_plane)
  if (allocated(ndotB_fourier_cos)) deallocate(ndotB_fourier_cos)
  if (allocated(ndotB_fourier_sin)) deallocate(ndotB_fourier_sin)
  if (allocated(vpar_target_fourier_cos)) deallocate(vpar_target_fourier_cos)
  if (allocated(vpar_target_fourier_sin)) deallocate(vpar_target_fourier_sin)
  
  ! Allocate per-node (backward compatible)
  allocate(ndotB_per_node(n_nodes))
  allocate(ndotB_count_per_node(n_nodes))
  
  ! Allocate per-node-plane
  allocate(ndotB_per_node_plane(n_nodes, n_plane_local))
  allocate(ndotB_count_per_node_plane(n_nodes, n_plane_local))
  
  ! Allocate Fourier coefficients
  allocate(ndotB_fourier_cos(n_nodes, n_tor_local))
  allocate(ndotB_fourier_sin(n_nodes, n_tor_local))
  allocate(vpar_target_fourier_cos(n_nodes, n_tor_local))
  allocate(vpar_target_fourier_sin(n_nodes, n_tor_local))
  
  ndotB_per_node = 0.d0
  ndotB_count_per_node = 0.d0
  ndotB_per_node_plane = 0.d0
  ndotB_count_per_node_plane = 0.d0
  ndotB_fourier_cos = 0.d0
  ndotB_fourier_sin = 0.d0
  vpar_target_fourier_cos = 0.d0
  vpar_target_fourier_sin = 0.d0
  
  n_nodes_stored = n_nodes
  n_plane_stored = n_plane_local
  n_tor_stored = n_tor_local
  
  ndotB_initialized = .true.
  ndotB_finalized = .false.
  
end subroutine init_boundary_ndotB

subroutine accumulate_ndotB_at_node(inode, ndotB_value)
  implicit none
  integer, intent(in) :: inode
  real*8, intent(in) :: ndotB_value
  
  if (.not. ndotB_initialized) return
  if (inode < 1 .or. inode > n_nodes_stored) return
  
  ndotB_per_node(inode) = ndotB_per_node(inode) + ndotB_value
  ndotB_count_per_node(inode) = ndotB_count_per_node(inode) + 1.d0
  
end subroutine accumulate_ndotB_at_node

subroutine accumulate_ndotB_at_node_plane(inode, mp, ndotB_value)
  !---------------------------------------------------------------------------
  ! Accumulate ndotB at a specific node and toroidal plane
  ! This is called for each plane during element assembly
  ! Skip if ndotB is frozen (already computed and not evolving)
  !---------------------------------------------------------------------------
  use phys_module, only: ndotB_evolving
  implicit none
  integer, intent(in) :: inode, mp
  real*8, intent(in) :: ndotB_value
  
  ! If frozen, skip accumulation (use stored values from first timestep)
  if (ndotB_computed_once .and. (.not. ndotB_evolving)) return
  
  if (.not. ndotB_initialized) return
  if (inode < 1 .or. inode > n_nodes_stored) return
  if (mp < 1 .or. mp > n_plane_stored) return
  
  ndotB_per_node_plane(inode, mp) = ndotB_per_node_plane(inode, mp) + ndotB_value
  ndotB_count_per_node_plane(inode, mp) = ndotB_count_per_node_plane(inode, mp) + 1.d0
  
  ! Also accumulate into the toroidally-averaged array for backward compatibility
  ndotB_per_node(inode) = ndotB_per_node(inode) + ndotB_value
  ndotB_count_per_node(inode) = ndotB_count_per_node(inode) + 1.d0
  
end subroutine accumulate_ndotB_at_node_plane

subroutine finalize_boundary_ndotB()
  use mpi
  use data_structure
  use mod_parameters, only: n_period
  use corr_neg, only: corr_neg_temp
  use phys_module, only: vpar_sbc_alpha0, vpar_sbc_strength, vpar_sbc_smooth_sign, &
                         vpar_sbc_angle_scale, T_0, T_1, GAMMA, ndotB_evolving, loop_voltage, &
                         sbc_use_local_T
  use mod_model_settings, only: var_T
  implicit none
  real*8, parameter :: pi = 3.14159265358979d0
  integer :: i, mp, in, n_with_data, ierr, my_id
  real*8 :: ndotB_sum, ndotB_avg, ndotB_min, ndotB_max
  real*8 :: phi, ndotB_val, ndotB_scaled, cos_n, sin_n
  real*8 :: alpha0_rad, alpha_rad, factor_sbc, vpar_target_val, cs, T_local
  real*8, allocatable :: ndotB_global(:), count_global(:)
  real*8, allocatable :: ndotB_plane_global(:,:), count_plane_global(:,:)
  
  if (.not. ndotB_initialized) return
  if (ndotB_finalized) return
  
  ! If frozen (already computed once and not evolving), skip re-finalization
  if (ndotB_computed_once .and. (.not. ndotB_evolving)) then
    ndotB_finalized = .true.
    return
  endif
  
  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
  
  ! Allocate temporary arrays for global reduction
  allocate(ndotB_global(n_nodes_stored))
  allocate(count_global(n_nodes_stored))
  allocate(ndotB_plane_global(n_nodes_stored, n_plane_stored))
  allocate(count_plane_global(n_nodes_stored, n_plane_stored))
  
  ! Sum up toroidally-averaged values from all ranks
  call MPI_Allreduce(ndotB_per_node, ndotB_global, n_nodes_stored, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(ndotB_count_per_node, count_global, n_nodes_stored, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Sum up per-plane values from all ranks
  call MPI_Allreduce(ndotB_per_node_plane, ndotB_plane_global, &
                     n_nodes_stored * n_plane_stored, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(ndotB_count_per_node_plane, count_plane_global, &
                     n_nodes_stored * n_plane_stored, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  ! Compute averages for toroidally-averaged values
  n_with_data = 0
  ndotB_sum = 0.d0
  ndotB_min = 1.d10
  ndotB_max = -1.d10
  
  do i = 1, n_nodes_stored
    if (count_global(i) > 0.5d0) then
      ndotB_per_node(i) = ndotB_global(i) / count_global(i)
      n_with_data = n_with_data + 1
      ndotB_sum = ndotB_sum + ndotB_per_node(i)
      if (ndotB_per_node(i) < ndotB_min) ndotB_min = ndotB_per_node(i)
      if (ndotB_per_node(i) > ndotB_max) ndotB_max = ndotB_per_node(i)
    else
      ndotB_per_node(i) = 0.d0
    endif
    
    ! Compute averages for per-plane values
    do mp = 1, n_plane_stored
      if (count_plane_global(i, mp) > 0.5d0) then
        ndotB_per_node_plane(i, mp) = ndotB_plane_global(i, mp) / count_plane_global(i, mp)
      else
        ndotB_per_node_plane(i, mp) = 0.d0
      endif
    enddo
  enddo
  
  deallocate(ndotB_global, count_global)
  deallocate(ndotB_plane_global, count_plane_global)
  
  ! Compute Fourier coefficients from per-plane data
  ! f(phi) = a_0 + sum_n [a_n * cos(n*n_period*phi) + b_n * sin(n*n_period*phi)]
  ndotB_fourier_cos = 0.d0
  ndotB_fourier_sin = 0.d0
  vpar_target_fourier_cos = 0.d0
  vpar_target_fourier_sin = 0.d0
  
  ! Use node-local temperature if needed
  if (sbc_use_local_T .and. var_T .gt. 0) then
    ! Note: 2T mode (var_T=0) always falls back to T_1.
    ! For physical 2T SBC, this should use Ti+Te at the boundary.
    ! This requires separate implementation when 2T SBC is needed.
    T_local = corr_neg_temp(node_list%node(i)%values(1,1,var_T))
  else
    T_local = corr_neg_temp(T_1)  ! Use normalized SOL temperature
  endif

  ! Compute SBC parameters once
  alpha0_rad = vpar_sbc_alpha0 * pi / 180.d0
  cs = sqrt(GAMMA * T_local)  ! Use reference sound speed
  
  do i = 1, n_nodes_stored
    do in = 1, n_tor_stored
      do mp = 1, n_plane_stored
        ! Mode in=2 has n_period oscillations over 2pi (one per field period) -- lowest stellarator mode.
        phi = 2.d0 * pi * dble(mp-1) / dble(n_plane_stored * n_period)
        ndotB_val = ndotB_per_node_plane(i, mp)
        
        ! Compute vpar_target in physical space
        if (vpar_sbc_smooth_sign) then
          ! Smooth sign formulation: vpar = cs * tanh(ndotB_scaled / sin(alpha0))
          ! This avoids sign() discontinuity that causes Gibbs phenomenon
          ! ndotB = sin(alpha), so sin(alpha0) normalizes to make transition at alpha0
          ndotB_scaled = ndotB_val * vpar_sbc_angle_scale
          vpar_target_val = cs * tanh(ndotB_scaled / sin(alpha0_rad)) * vpar_sbc_strength
        else
          ! Original formulation: vpar = sign(ndotB) * cs * tanh(|alpha|/alpha0)
          ! Has sign() discontinuity causing Gibbs overshoot at ndotB sign changes
          ndotB_scaled = ndotB_val * vpar_sbc_angle_scale
          alpha_rad = asin(min(1.d0, max(-1.d0, abs(ndotB_scaled))))
          factor_sbc = tanh(alpha_rad / alpha0_rad)
          vpar_target_val = sign(1.d0, ndotB_scaled) * cs * factor_sbc * vpar_sbc_strength
        endif
        
        ! Fourier coefficient for mode (in-1) [0-indexed internally]
        ! in=1 is n=0 mode (constant), in=2 is n=1 mode, etc.
        cos_n = cos(dble(in-1) * dble(n_period) * phi)
        sin_n = sin(dble(in-1) * dble(n_period) * phi)
        
        ! Fourier transform of ndotB (for diagnostics)
        ndotB_fourier_cos(i, in) = ndotB_fourier_cos(i, in) + ndotB_val * cos_n
        ndotB_fourier_sin(i, in) = ndotB_fourier_sin(i, in) + ndotB_val * sin_n
        
        ! Fourier transform of vpar_target (for BC application - nonlinear!)
        vpar_target_fourier_cos(i, in) = vpar_target_fourier_cos(i, in) + vpar_target_val * cos_n
        vpar_target_fourier_sin(i, in) = vpar_target_fourier_sin(i, in) + vpar_target_val * sin_n
      enddo
      
      ! Normalize by number of planes (DFT normalization)
      ndotB_fourier_cos(i, in) = ndotB_fourier_cos(i, in) / dble(n_plane_stored)
      ndotB_fourier_sin(i, in) = ndotB_fourier_sin(i, in) / dble(n_plane_stored)
      vpar_target_fourier_cos(i, in) = vpar_target_fourier_cos(i, in) / dble(n_plane_stored)
      vpar_target_fourier_sin(i, in) = vpar_target_fourier_sin(i, in) / dble(n_plane_stored)
      
      ! For n=0 mode (in=1), the normalization is just the average
      ! For n>0 modes, multiply by 2 (standard DFT convention for real signals)
      if (in > 1) then
        ndotB_fourier_cos(i, in) = ndotB_fourier_cos(i, in) * 2.d0
        ndotB_fourier_sin(i, in) = ndotB_fourier_sin(i, in) * 2.d0
        vpar_target_fourier_cos(i, in) = vpar_target_fourier_cos(i, in) * 2.d0
        vpar_target_fourier_sin(i, in) = vpar_target_fourier_sin(i, in) * 2.d0
      endif
    enddo
  enddo
  
  if (n_with_data > 0) then
    ndotB_avg = ndotB_sum / dble(n_with_data)
  else
    ndotB_avg = 0.d0
  endif
  
  if (my_id == 0) then
    write(*,'(A,I6,A,E12.4,A,E12.4,A,E12.4)') &
      " ndotB stats (MPI reduced): n_nodes=", n_with_data, " min=", ndotB_min, &
      " max=", ndotB_max, " avg=", ndotB_avg
    if (.not. ndotB_evolving) then
      write(*,'(A)') "   ndotB_evolving=.false.: values frozen for subsequent timesteps"
    endif

    if (loop_voltage .ne. 0.d0 .and. .not. ndotB_evolving) then
      write(*,'(A)') "WARNING: loop_voltage != 0 but ndotB_evolving=.false."
      write(*,'(A)') "  SBC will use frozen t=0 field geometry. Set ndotB_evolving=.true. for consistency."
    endif
  endif
  
  ndotB_finalized = .true.
  ndotB_computed_once = .true.  ! Mark that ndotB has been computed at least once
  
end subroutine finalize_boundary_ndotB

function get_ndotB_at_node(inode) result(ndotB_value)
  implicit none
  integer, intent(in) :: inode
  real*8 :: ndotB_value
  
  if (.not. ndotB_initialized) then
    ndotB_value = 0.d0
    return
  endif
  
  if (inode < 1 .or. inode > n_nodes_stored) then
    ndotB_value = 0.d0
    return
  endif
  
  ndotB_value = ndotB_per_node(inode)
  
end function get_ndotB_at_node

subroutine get_ndotB_fourier_at_node(inode, in, ndotB_cos, ndotB_sin)
  !---------------------------------------------------------------------------
  ! Get Fourier coefficients of n.B/|B| at a boundary node for harmonic 'in'
  ! Returns both cosine and sine components for the toroidal mode
  ! in=1 corresponds to n=0 (constant), in=2 to n=1*n_period, etc.
  !---------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: inode, in
  real*8, intent(out) :: ndotB_cos, ndotB_sin
  
  if (.not. ndotB_initialized) then
    ndotB_cos = 0.d0
    ndotB_sin = 0.d0
    return
  endif
  
  if (inode < 1 .or. inode > n_nodes_stored) then
    ndotB_cos = 0.d0
    ndotB_sin = 0.d0
    return
  endif
  
  if (in < 1 .or. in > n_tor_stored) then
    ndotB_cos = 0.d0
    ndotB_sin = 0.d0
    return
  endif
  
  ndotB_cos = ndotB_fourier_cos(inode, in)
  ndotB_sin = ndotB_fourier_sin(inode, in)
  
end subroutine get_ndotB_fourier_at_node

function get_vpar_target_for_column(inode, in) result(vpar_target)
  !---------------------------------------------------------------------------
  ! Return vpar_target for JOREK column index 'in' at boundary node.
  ! Handles column-to-harmonic mapping and cos/sin selection internally.
  ! JOREK basis: in=1 DC, (2,3) first harmonic cos/-sin, (4,5) second, etc.
  !---------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: inode, in
  real*8 :: vpar_target
  integer :: k_fourier

  if (.not. ndotB_initialized .or. inode < 1 .or. inode > n_nodes_stored) then
    vpar_target = 0.d0
    return
  endif

  if (in .eq. 1) then
    vpar_target = vpar_target_fourier_cos(inode, 1)   ! DC; sin(0)=0 identically
    return
  endif

  k_fourier = in / 2 + 1
  if (k_fourier > n_tor_stored) then
    vpar_target = 0.d0
    return
  endif

  if (mod(in, 2) .eq. 0) then
    ! Even in: cosine component
    vpar_target = vpar_target_fourier_cos(inode, k_fourier)
  else
    ! Odd in: sine component; JOREK basis is -sin
    vpar_target = -vpar_target_fourier_sin(inode, k_fourier)
  endif

end function get_vpar_target_for_column

end module mod_boundary_ndotB
