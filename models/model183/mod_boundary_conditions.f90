!*******************************************************************************
!* Subroutine: boundary_condition                                              *
!*******************************************************************************
!*                                                                             *
!* Add boundary condition on the matrix.                                       *
!*                                                                             *
!* Parameters:                                                                 *
!*   my_id        - Identifier of the node in MPI_COMM_WORLD                   *
!*   node_list    - List of nodes                                              *
!*   element_list - List of all elements                                       *
!*   local_elms   - List of local elements                                     *
!*   n_local_elms - Number of local elements                                   *
!*   index_min    - Minimal index of local elements                            *
!*   index_max    - Maximal index of local elements                            *
!*   xpoint2      -                                                            *
!*   xcase2       -                                                            *
!*   psi_axis     -                                                            *
!*   psi_bnd      -                                                            *
!*   Z_xpoint     -                                                            *
!*   gmres        - boolean indicating if we are using GMRES method            *
!*   solve_only   - Indicate if we want to perform only solve                  *
!*                                                                             *
!*******************************************************************************
module mod_boundary_conditions
implicit none
contains
  subroutine boundary_conditions( my_id, node_list, element_list, bnd_node_list, local_elms,& 
                                n_local_elms, index_min, index_max, rhs_loc, xpoint2,     &
                                xcase2, R_axis, Z_axis, psi_axis, psi_bnd,                &
                                R_xpoint, Z_xpoint, psi_xpoint, a_mat)

    use mod_assembly, only : boundary_conditions_add_one_entry, boundary_conditions_add_RHS

    use phys_module, only: F0, bc_natural_open, GAMMA, T_0, &
                           vpar_sbc_enable, vpar_sbc_alpha0, vpar_sbc_strength, &
                           particle_flux_sbc_enable, heat_flux_sbc_enable, &
                           loop_voltage, tstep, central_density, central_mass, &
                           sbc_use_local_T
    use mod_boundary_ndotB, only: get_ndotB_at_node, get_ndotB_fourier_at_node, &
                                  get_vpar_target_fourier_at_node
    use mod_model_settings, only: var_Psi, var_Phi, var_zj, var_w, var_rho, var_T, &
                                  var_Vpar, var_Ti, var_Te, n_var
    use vacuum, only: is_freebound
    use constants, only: MU_ZERO, ATOMIC_MASS_UNIT
    use mpi_mod
    use mod_locate_irn_jcn
    use mod_integer_types
    use data_structure

    implicit none

    ! --- Routine parameters
    integer,                            intent(in)    :: my_id
    type (type_node_list),              intent(in)    :: node_list
    type (type_element_list),           intent(in)    :: element_list
    type (type_bnd_node_list),          intent(in)    :: bnd_node_list
    integer,                            intent(in)    :: local_elms(*)
    integer,                            intent(in)    :: n_local_elms
    integer,                            intent(in)    :: index_min
    integer,                            intent(in)    :: index_max
    logical,                            intent(in)    :: xpoint2
    integer,                            intent(in)    :: xcase2
    real*8,                             intent(in)    :: R_axis
    real*8,                             intent(in)    :: Z_axis
    real*8,                             intent(in)    :: psi_axis
    real*8,                             intent(in)    :: psi_bnd
    real*8,                             intent(in)    :: R_xpoint(2)
    real*8,                             intent(in)    :: Z_xpoint(2)
    real*8,                             intent(in)    :: psi_xpoint(2)
    real*8,                             intent(inout) :: rhs_loc(*)
    type(type_SP_MATRIX)                              :: a_mat

    ! Internal parameters
    real*8                :: zbig
    integer               :: i, in, iv, inode, k
    integer               :: ielm
    integer               :: index_node
    
    ! v_par SBC variables
    integer               :: k_fourier
    real*8                :: vpar_target, vpar_current, delta_vpar
    real*8                :: ndotB_norm, ndotB_cos, ndotB_sin, T_local
    real*8                :: cs, alpha_rad, alpha0_rad, factor_sbc
    real*8                :: vpar_target_cos, vpar_target_sin
    real*8, parameter     :: pi = 3.14159265358979d0

    zbig = 1.d12
       do i=1, n_local_elms

          ielm = local_elms(i)

          do iv=1, n_vertex_max

             inode = element_list%element(ielm)%vertex(iv)

             if (node_list%node(inode)%boundary .ne. 0) then

                do in=a_mat%i_tor_min, a_mat%i_tor_max

                   do k=1, n_var
                     if (bc_natural_open .and. k .eq. var_zj) cycle
                     ! Skip Dirichlet BC for density when particle flux SBC is enabled
                     ! (weak-form natural BC handled in mod_boundary_matrix_open.f90)
                     if (particle_flux_sbc_enable .and. k .eq. var_rho) cycle
                     ! Skip Dirichlet BC for temperature when heat flux SBC is enabled
                     if (heat_flux_sbc_enable .and. k .eq. var_T) cycle

                      !------------------------------------ boundary nodes (types 1, 2, 3)
                      ! Type 1: Open field lines (tokamak divertor targets)
                      ! Type 2: Wall/limiter (stellarator with divertor geometry)
                      ! Type 3: Corner nodes (both target and wall)
                      if ((node_list%node(inode)%boundary .eq. 1) .or. &
                          (node_list%node(inode)%boundary .eq. 2) .or. &
                          (node_list%node(inode)%boundary .eq. 3)) then

                         if ((k .eq. var_Psi) .or. (k .eq. var_Phi) .or. (k .eq. var_zj) .or. &
                              (k .eq. var_w) .or. (k .eq. var_rho) .or. (k .eq. var_T) .or. (k .eq. var_Vpar) .or. &
                              (k .eq. var_Ti) .or. (k .eq. var_Te)) then
 
                          if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)
                            
                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)
                            
                            ! v_par Dirichlet BC: add RHS term for non-zero target
                            ! For n=0 mode, use toroidally-averaged ndotB
                            ! For n>0 modes, use Fourier coefficients of ndotB
                            if (k .eq. var_Vpar) then
                              ! Compute angle-dependent target if SBC enabled
                              if (vpar_sbc_enable) then
                                ! Temperature for sound speed: local edge T or core T_0
                                if (sbc_use_local_T) then
                                  T_local = node_list%node(inode)%values(1,1,var_T)
                                  if (T_local < 1.d-8) T_local = T_0  ! Fallback for numerical zeros
                                else
                                  T_local = T_0  ! Use core temperature (default)
                                endif
                                cs = sqrt(GAMMA * T_local)
                                alpha0_rad = vpar_sbc_alpha0 * pi / 180.d0
                                
                                if (in .eq. 1) then
                                  ! n=0 mode: use toroidally-averaged ndotB
                                  ndotB_norm = get_ndotB_at_node(inode)
                                  alpha_rad = asin(min(1.d0, max(-1.d0, abs(ndotB_norm))))
                                  factor_sbc = tanh(alpha_rad / alpha0_rad)
                                  vpar_target = sign(1.d0, ndotB_norm) * cs * factor_sbc * vpar_sbc_strength
                                  
                                  ! Get current vpar value at this node (n=0 mode)
                                  vpar_current = node_list%node(inode)%values(1,1,var_Vpar)
                                  delta_vpar = vpar_target - vpar_current
                                  
                                  call boundary_conditions_add_RHS(                     &
                                         index_node, k, in, index_min, index_max,       &
                                         rhs_loc, zbig * delta_vpar,                    &
                                         a_mat%i_tor_min, a_mat%i_tor_max)
                                else
                                  ! n>0 modes: use precomputed vpar_target Fourier coefficients
                                  ! Map JOREK harmonic index (in) to Fourier mode index:
                                  ! JOREK layout: in=1: n=0, in=2: +cos(nfp*phi),
                                  !   in=3: -sin(nfp*phi), in=4: +cos(2*nfp*phi), ...
                                  ! Fourier storage: k=1: n=0, k=2: n=nfp, k=3: n=2*nfp, ...
                                  k_fourier = in / 2 + 1
                                  call get_vpar_target_fourier_at_node(inode, k_fourier, vpar_target_cos, vpar_target_sin)
                                  
                                  ! Select cos or sin and apply JOREK -sin convention
                                  if (mod(in, 2) .eq. 0) then
                                    ! Even in: cosine component
                                    delta_vpar = vpar_target_cos - node_list%node(inode)%values(in,1,var_Vpar)
                                  else
                                    ! Odd in: sine component (JOREK uses -sin basis)
                                    delta_vpar = -vpar_target_sin - node_list%node(inode)%values(in,1,var_Vpar)
                                  endif
                                  
                                  call boundary_conditions_add_RHS(                     &
                                         index_node, k, in, index_min, index_max,       &
                                         rhs_loc, zbig * delta_vpar,                    &
                                         a_mat%i_tor_min, a_mat%i_tor_max)
                                endif
                              else
                                ! SBC disabled: zero BC for all modes (no RHS addition needed)
                              endif
                            endif

                            index_node = node_list%node(inode)%index(2)

                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)
                            
                            ! Derivative BC: zero for constant v_par (no RHS addition)
                            
                          endif
                        endif
                      endif

                      !------------------------------------ wall aligned with fluxsurface (in case of x-point grid)
                      if ((node_list%node(inode)%boundary .eq. 2) .or. (node_list%node(inode)%boundary .eq. 3)) then

                         if ( (.not. is_freebound(in,k)) ) then ! apply fixed boundary conditions where necessary

                            index_node = node_list%node(inode)%index(1)

                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)
                            
                            ! v_par Dirichlet BC: add RHS term for non-zero target
                            ! For n=0 mode, use toroidally-averaged ndotB
                            ! For n>0 modes, use Fourier coefficients of ndotB
                            if (k .eq. var_Vpar) then
                              ! Compute angle-dependent target if SBC enabled
                              if (vpar_sbc_enable) then
                                ! Temperature for sound speed: local edge T or core T_0
                                if (sbc_use_local_T) then
                                  T_local = node_list%node(inode)%values(1,1,var_T)
                                  if (T_local < 1.d-8) T_local = T_0  ! Fallback for numerical zeros
                                else
                                  T_local = T_0  ! Use core temperature (default)
                                endif
                                cs = sqrt(GAMMA * T_local)
                                alpha0_rad = vpar_sbc_alpha0 * pi / 180.d0
                                
                                if (in .eq. 1) then
                                  ! n=0 mode: use toroidally-averaged ndotB
                                  ndotB_norm = get_ndotB_at_node(inode)
                                  alpha_rad = asin(min(1.d0, max(-1.d0, abs(ndotB_norm))))
                                  factor_sbc = tanh(alpha_rad / alpha0_rad)
                                  vpar_target = sign(1.d0, ndotB_norm) * cs * factor_sbc * vpar_sbc_strength
                                  
                                  ! Get current vpar value at this node (n=0 mode)
                                  vpar_current = node_list%node(inode)%values(1,1,var_Vpar)
                                  delta_vpar = vpar_target - vpar_current
                                  
                                  call boundary_conditions_add_RHS(                     &
                                         index_node, k, in, index_min, index_max,       &
                                         rhs_loc, zbig * delta_vpar,                    &
                                         a_mat%i_tor_min, a_mat%i_tor_max)
                                else
                                  ! n>0 modes: use precomputed vpar_target Fourier coefficients
                                  k_fourier = in / 2 + 1
                                  call get_vpar_target_fourier_at_node(inode, k_fourier, vpar_target_cos, vpar_target_sin)
                                  if (mod(in, 2) .eq. 0) then
                                    delta_vpar = vpar_target_cos - node_list%node(inode)%values(in,1,var_Vpar)
                                  else
                                    delta_vpar = -vpar_target_sin - node_list%node(inode)%values(in,1,var_Vpar)
                                  endif
                                  
                                  call boundary_conditions_add_RHS(                     &
                                         index_node, k, in, index_min, index_max,       &
                                         rhs_loc, zbig * delta_vpar,                    &
                                         a_mat%i_tor_min, a_mat%i_tor_max)
                                endif
                              endif
                            endif

                            index_node = node_list%node(inode)%index(3)

                            call boundary_conditions_add_one_entry(                 &
                                   index_node, k, in, index_node, k, in,            &
                                   zbig, index_min, index_max, a_mat)

                         endif

                      endif
                      
                      !------------------------------------ divertor/frozen nodes (type 4)
                      ! Type 4: Divertor region nodes (frozen, no evolution)
                      ! All variables frozen at current values: delta = 0
                      ! Used for extended grids where divertor region should not evolve
                      if (node_list%node(inode)%boundary .eq. 4) then
                         ! Freeze ALL variables for this node
                         if ((k .eq. var_Psi) .or. (k .eq. var_Phi) .or. (k .eq. var_zj) .or. &
                              (k .eq. var_w) .or. (k .eq. var_rho) .or. (k .eq. var_T) .or. (k .eq. var_Vpar) .or. &
                              (k .eq. var_Ti) .or. (k .eq. var_Te)) then
                           
                           ! Apply Dirichlet BC: delta = 0 (no change from current value)
                           ! For all 4 degrees of freedom (value and 3 derivatives)
                           index_node = node_list%node(inode)%index(1)
                           call boundary_conditions_add_one_entry(                 &
                                  index_node, k, in, index_node, k, in,            &
                                  zbig, index_min, index_max, a_mat)
                           ! RHS = 0 (default) → delta = 0 → frozen at current value
                           
                           index_node = node_list%node(inode)%index(2)
                           call boundary_conditions_add_one_entry(                 &
                                  index_node, k, in, index_node, k, in,            &
                                  zbig, index_min, index_max, a_mat)
                           
                           index_node = node_list%node(inode)%index(3)
                           call boundary_conditions_add_one_entry(                 &
                                  index_node, k, in, index_node, k, in,            &
                                  zbig, index_min, index_max, a_mat)
                           
                           index_node = node_list%node(inode)%index(4)
                           call boundary_conditions_add_one_entry(                 &
                                  index_node, k, in, index_node, k, in,            &
                                  zbig, index_min, index_max, a_mat)
                         endif
                      endif

                   enddo  ! k=1,n_var (variables loop)
                   
                   ! Apply loop voltage to drive psi evolution (n=0 mode at boundary)
                   ! Physics: loop_voltage drives Ohmic current via Faraday's law:
                   !   d(psi)/dt = -V_loop => psi(t) = psi(0) - V_loop * t
                   ! This causes magnetic flux to evolve, changing topology over time.
                   ! Used for: tearing mode studies, current ramp scenarios, stellarator
                   ! transport with evolving fields.
                   if ( loop_voltage .ne. 0.d0 ) then
                      if ( in == 1 ) then  ! n=0 mode only
                         if ( (.not. is_freebound(in, var_psi)) ) then
                            index_node = node_list%node(inode)%index(1)
                            call boundary_conditions_add_RHS(       &
                                      index_node, var_psi, in,      &
                                      index_min, index_max,         &
                                      RHS_loc, zbig*loop_voltage*sqrt(MU_ZERO*central_density*central_mass*ATOMIC_MASS_UNIT*1.d20)*tstep, &
                                      a_mat%i_tor_min, a_mat%i_tor_max)
                         endif
                      endif
                   endif

                enddo  ! in (toroidal modes)
             endif
          enddo
       enddo

    return
  end subroutine boundary_conditions
end module mod_boundary_conditions
