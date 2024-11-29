!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! Author: Joseph.Mouallem @FV3 GFDL/Princeton
! Initial release: ../../....
module duogrid_mod

  use mpp_mod, only: FATAL, mpp_error
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
  use mpp_mod, only: mpp_pe, mpp_npes, mpp_root_pe, mpp_sync, mpp_send, mpp_recv
  use mpp_domains_mod, only: mpp_global_max, mpp_global_min, mpp_get_tile_id, mpp_global_sum
  use constants_mod, only: OMEGA, grav, rdgas, rvgas, pi => pi_8, kappa
  use mpp_domains_mod, only: BITWISE_EXACT_SUM, domain2d, BITWISE_EFP_SUM

  use mpp_domains_mod, only: mpp_update_domains

  use lib_grid_mod, only: R_GRID, RADIUS
  use lib_grid_mod, only: lib_4pt_area

  use global_grid_mod, only: global_grid_type
  use global_grid_mod, only: global_grid_init, global_grid_end

  use fv_arrays_mod, only: duogrid_type, fv_grid_bounds_type, fv_flags_type, fv_atmos_type, fv_grid_type
  !use fv_grid_utils_mod, only: c2l_ord2, great_circle_dist
  use fv_timing_mod, only: timing_on, timing_off

  use mpp_parameter_mod, only: DGRID_NE, BGRID_NE, BGRID_SW, AGRID, SCALAR_PAIR, NORTH, EAST, WEST, SOUTH, CGRID_NE

  !use fv_grid_utils_mod, only: mid_pt_cart, latlon2xyz, vect_cross, normalize_vect
  implicit none

  private

  public :: duogrid_init
  public :: duogrid_end
  public :: ext_scalar
  public :: ext_vector
  public :: fill_corner_region, lagrange_poly_interp

  interface fill_corner_region
    module procedure fill_corner_region_2d
    module procedure fill_corner_region_3d
    module procedure fill_corner_region_4d
  end interface fill_corner_region
  interface lagrange_poly_interp
    module procedure lagrange_poly_interp_2d
    module procedure lagrange_poly_interp_3d
    module procedure lagrange_poly_interp_4d
  end interface lagrange_poly_interp

  interface fill_corners_domain_decomp
    module procedure fill_corners_domain_decomp_2d
    module procedure fill_corners_domain_decomp_3d
  end interface fill_corners_domain_decomp

  interface ext_scalar
    module procedure ext_scalar_2d
    module procedure ext_scalar_3d
    module procedure ext_scalar_4d
  end interface ext_scalar

  integer :: interporder = 3 !order is interporder+1
  logical :: laginter = .true. ! use the laginterp function

contains

  subroutine duogrid_init(Atm)

    type(fv_atmos_type), intent(INOUT) :: Atm
    ! !--- local
    type(global_grid_type) :: gg
    integer :: ng, res

    !--- check status
    if (atm%gridstruct%dg%is_initialized) return

    !     if (.not.mp%is_initialized) &
    !         call mpp_error(FATAL, 'duogrid_init: mp needs to be intialized')

    !set indices for extended bd of duo
    call set_bd_ext_duo(atm%gridstruct%dg)

    !--- alloc grid
    call duogrid_alloc(atm%gridstruct%dg, atm%bd, atm%flagstruct)

    !--- create global_grid and initialize with it
    res = atm%flagstruct%npx - 1
    ng = atm%bd%ng + 2 !need more haloes for staggered points
    ng = atm%gridstruct%dg%bd%ng + 2 !need more haloes for staggered points

    call global_grid_init(gg, res, ng, atm%gridstruct%grid_type, atm%flagstruct%do_schmidt, &
                          atm%flagstruct%target_lon, atm%flagstruct%target_lat)

    call duogrid_init_with_global_grid(atm%gridstruct%dg, gg, atm%gridstruct%dg%bd, atm%flagstruct, atm%global_tile)

    !--- coriolis forces
    call duogrid_init_coriolis(atm%gridstruct%dg, atm%gridstruct%dg%bd)

    !--- A 2 CD grid metrics
    call a2stag_metrics(atm%gridstruct, atm%gridstruct%dg%bd, atm%gridstruct%dg)

    !--- lagrange corners coeff
    !call compute_lagrange_coeff(atm%gridstruct%dg%ss_store, atm%bd, atm%gridstruct%dg)
    call compute_lagrange_coeff(atm%gridstruct%dg%xp, atm%gridstruct%dg%xm, &
                                atm%gridstruct%dg%yp, atm%gridstruct%dg%ym, atm%bd, atm%gridstruct%dg)

    !--- end global_grid
    call global_grid_end(gg)

    !--- set status
    atm%gridstruct%dg%is_initialized = .true.

  end subroutine duogrid_init
  !===========================================================================
  subroutine set_bd_ext_duo(dg)
    type(duogrid_type), intent(inout) :: dg
    integer :: is, ie, js, je, isd, ied, jsd, jed
    call mpp_get_compute_domain(dg%domain_for_duo, is, ie, js, je)
    call mpp_get_data_domain(dg%domain_for_duo, isd, ied, jsd, jed)
    dg%bd%is = is
    dg%bd%ie = ie
    dg%bd%js = js
    dg%bd%je = je
    dg%bd%isd = isd
    dg%bd%ied = ied
    dg%bd%jsd = jsd
    dg%bd%jed = jed
    dg%bd%ng = 4

  end subroutine set_bd_ext_duo
  !> @brief end method for duogrid_type
  subroutine duogrid_end(dg)
    type(duogrid_type), intent(inout) :: dg

    !--- check status
    if (.not. dg%is_initialized) return

    !--- dealloc grid
    call duogrid_dealloc(dg)

    !--- set status
    dg%is_initialized = .false.

  end subroutine duogrid_end
  !===========================================================================
  subroutine duogrid_init_coriolis(dg, bd)
    type(duogrid_type), intent(inout) :: dg
    type(fv_grid_bounds_type), intent(in) :: bd
    !--- local
    real(kind=R_GRID) :: lon, lat
    real :: alpha
    integer :: i, j
    integer :: isd, ied, jsd, jed

    !--- assign dims
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    !--- zero out f
    dg%a_f(:, :) = 0.
    dg%a_kik_f(:, :) = 0.

    !alpha = dg%fg%alpha
    alpha = 0

    !--- initialize f
    do j = jsd, jed
      do i = isd, ied
        ! ext
        lon = dg%a_pt(1, i, j)
        lat = dg%a_pt(2, i, j)
        dg%a_f(i, j) = 2.*OMEGA*(-1.*cos(lon)*cos(lat)*sin(alpha) + &
                                 sin(lat)*cos(alpha))
      end do
    end do

  end subroutine duogrid_init_coriolis

  subroutine duogrid_init_with_global_grid(dg, gg, bd, flagstruct, tile)
    type(duogrid_type), intent(inout) :: dg
    type(global_grid_type), intent(inout) :: gg
    type(fv_grid_bounds_type), intent(in) :: bd
    type(fv_flags_type), intent(inout) :: flagstruct
    integer, intent(in) :: tile
    !--- local
    integer :: i, j, n, ii, jj
    integer :: is, ie, js, je, isd, ied, jsd, jed, npx, npy

    !--- index rules
    !    a_pt ii = i*2   jj = j*2
    !    b_pt ii = i*2-1 jj = j*2-1
    !    c_pt ii = i*2-1 jj = j*2
    !    d_pt ii = i*2   jj = j*2-1

    !--- assign dims
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    n = tile

    npx = flagstruct%npx
    npy = flagstruct%npy

    if (.not. dg%k2e_nord == gg%k2e_nord) then
      call mpp_error(FATAL, 'k2e_nord should be concsistent between dg and gg')
    end if

    !--- assign a-grid
    do j = jsd, jed
      do i = isd, ied
        ! point
        ii = i*2; jj = j*2
        dg%a_pt(:, i, j) = gg%pt_ext(:, ii, jj, n)
        dg%a_x(i, j) = gg%ext_x(ii, jj, n)
        dg%a_y(i, j) = gg%ext_y(ii, jj, n)

        dg%a_kik_x(i, j) = gg%kik_x(ii, jj, n)
        dg%a_kik_y(i, j) = gg%kik_y(ii, jj, n)

        dg%a_gco(:, :, i, j) = gg%ext_g_co(:, :, ii, jj, n)
        dg%a_gct(:, :, i, j) = gg%ext_g_ctr(:, :, ii, jj, n)

        dg%a_c2l(:, :, i, j) = gg%ext_c2l(:, :, ii, jj, n)
        dg%a_l2c(:, :, i, j) = gg%ext_l2c(:, :, ii, jj, n)

        dg%a_sina(i, j) = gg%ext_sina(ii, jj, n)
        dg%a_cosa(i, j) = gg%ext_cosa(ii, jj, n)

        !--- length and area
        ii = i*2 - 1; jj = j*2
        dg%a_dx(i, j) = gg%ext_dx(ii, jj, n) + gg%ext_dx(ii + 1, jj, n)
        ii = i*2; jj = j*2 - 1
        dg%a_dy(i, j) = gg%ext_dy(ii, jj, n) + gg%ext_dy(ii, jj + 1, n)
        ii = i*2 - 1; jj = j*2 - 1
        dg%a_da(i, j) = gg%ext_da(ii, jj, n) + gg%ext_da(ii, jj + 1, n) + &
                        gg%ext_da(ii + 1, jj, n) + gg%ext_da(ii + 1, jj + 1, n)
        dg%rda(i, j) = 1.0/dg%a_da(i, j)
        dg%rdx(i, j) = 1.0/dg%a_dx(i, j)
        dg%rdy(i, j) = 1.0/dg%a_dy(i, j)

        !--- k2e
        dg%k2e_loc(i, j) = gg%k2e_loc(i, j, n)
        dg%k2e_coef(:, i, j) = gg%k2e_coef(:, i, j, n)
      end do
    end do

    !--- assign b-grid
    do j = jsd, jed + 1
      do i = isd, ied + 1
        ii = i*2 - 1; jj = j*2 - 1
        dg%b_pt(:, i, j) = gg%pt_ext(:, ii, jj, n)
        dg%b_x(i, j) = gg%ext_x(ii, jj, n)
        dg%b_y(i, j) = gg%ext_y(ii, jj, n)

        dg%b_kik_x(i, j) = gg%kik_x(ii, jj, n)
        dg%b_kik_y(i, j) = gg%kik_y(ii, jj, n)

        dg%b_gco(:, :, i, j) = gg%ext_g_co(:, :, ii, jj, n)
        dg%b_gct(:, :, i, j) = gg%ext_g_ctr(:, :, ii, jj, n)

        dg%b_c2l(:, :, i, j) = gg%ext_c2l(:, :, ii, jj, n)
        dg%b_l2c(:, :, i, j) = gg%ext_l2c(:, :, ii, jj, n)

        dg%b_sina(i, j) = gg%ext_sina(ii, jj, n)
        dg%b_cosa(i, j) = gg%ext_cosa(ii, jj, n)

        !--- length and area
        ii = i*2 - 1 - 1; jj = j*2 - 1
        dg%b_dx(i, j) = gg%ext_dx(ii, jj, n) + gg%ext_dx(ii + 1, jj, n)
        ii = i*2 - 1; jj = j*2 - 1 - 1
        dg%b_dy(i, j) = gg%ext_dy(ii, jj, n) + gg%ext_dy(ii, jj + 1, n)
        ii = i*2 - 1 - 1; jj = j*2 - 1 - 1
        dg%b_da(i, j) = gg%ext_da(ii, jj, n) + gg%ext_da(ii, jj + 1, n) + &
                        gg%ext_da(ii + 1, jj, n) + gg%ext_da(ii + 1, jj + 1, n)
        dg%b_rda(i, j) = 1.0/dg%b_da(i, j)
        dg%b_rdx(i, j) = 1.0/dg%b_dx(i, j)
        dg%b_rdy(i, j) = 1.0/dg%b_dy(i, j)

        !--- k2e
        dg%k2e_loc_b(i, j) = gg%k2e_loc_b(i, j, n)
        dg%k2e_coef_b(:, i, j) = gg%k2e_coef_b(:, i, j, n)
      end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! WRITE GRIDS FOR JUPYTER NOTEBOOK !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do j = jsd, jed + 1
!    do i = isd, ied + 1
!      write (300 + mpp_pe(), *), 'ijdpt1', i, j, dg%b_pt(1, i, j)
!      write (400 + mpp_pe(), *), 'ijdpt2', i, j, dg%b_pt(2, i, j)
!    end do
!    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define a_grid_recalculate
#ifdef a_grid_recalculate
    do j = jsd, jed
      do i = isd, ied
        dg%a_da(i, j) = lib_4pt_area( &
                        dg%b_pt(:, i, j), dg%b_pt(:, i + 1, j), dg%b_pt(:, i + 1, j + 1), dg%b_pt(:, i, j + 1), &
                        RADIUS)
        dg%rda(i, j) = 1.0/dg%a_da(i, j)
      end do
    end do
#endif

!TO CHECK: revisit the c/d grid assignement for some variables
    !--- assign c-grid
    do j = jsd, jed + 1
      do i = isd, ied + 1

        !--- point
        ii = i*2 - 1; jj = j*2
        dg%c_pt(:, i, j) = gg%pt_ext(:, ii, jj, n)
        dg%c_x(i, j) = gg%ext_x(ii, jj, n)
        dg%c_y(i, j) = gg%ext_y(ii, jj, n)

        dg%c_kik_x(i, j) = gg%kik_x(ii, jj, n)
        dg%c_kik_y(i, j) = gg%kik_y(ii, jj, n)

        dg%c_gco(:, :, i, j) = gg%ext_g_co(:, :, ii, jj, n)
        dg%c_gct(:, :, i, j) = gg%ext_g_ctr(:, :, ii, jj, n)

        dg%c_c2l(:, :, i, j) = gg%ext_c2l(:, :, ii, jj, n)
        dg%c_l2c(:, :, i, j) = gg%ext_l2c(:, :, ii, jj, n)

        dg%c_sina(i, j) = gg%ext_sina(ii, jj, n)
        dg%c_cosa(i, j) = gg%ext_cosa(ii, jj, n)

        !--- length and area
        ii = i*2 - 1 - 1; jj = j*2
        dg%c_dx(i, j) = gg%ext_dx(ii, jj, n) + gg%ext_dx(ii + 1, jj, n)
        ii = i*2 - 1; jj = j*2 - 1
        dg%c_dy(i, j) = gg%ext_dy(ii, jj, n) + gg%ext_dy(ii, jj + 1, n)
        ii = i*2 - 1 - 1; jj = j*2 - 1
        dg%c_da(i, j) = gg%ext_da(ii, jj, n) + gg%ext_da(ii, jj + 1, n) + &
                        gg%ext_da(ii + 1, jj, n) + gg%ext_da(ii + 1, jj + 1, n)
        dg%c_rda(i, j) = 1.0/dg%c_da(i, j)
        dg%c_rdx(i, j) = 1.0/dg%c_dx(i, j)
        dg%c_rdy(i, j) = 1.0/dg%c_dy(i, j)
      end do
    end do
    !--- k2e
    do j = jsd, jed
      do i = isd, ied + 1
        dg%k2e_loc_c_x(i, j) = gg%k2e_loc_c_x(i, j, n)
        dg%k2e_coef_c_x(:, i, j) = gg%k2e_coef_c_x(:, i, j, n)
      end do
    end do
    do j = jsd, jed + 1
      do i = isd, ied
        dg%k2e_loc_c_y(i, j) = gg%k2e_loc_c_y(i, j, n)
        dg%k2e_coef_c_y(:, i, j) = gg%k2e_coef_c_y(:, i, j, n)
      end do
    end do

    !--- assign d-grid
    do j = jsd, jed + 1
      do i = isd, ied + 1

        !--- point
        ii = i*2; jj = j*2 - 1
        dg%d_pt(:, i, j) = gg%pt_ext(:, ii, jj, n)
        dg%d_x(i, j) = gg%ext_x(ii, jj, n)
        dg%d_y(i, j) = gg%ext_y(ii, jj, n)

        dg%d_kik_x(i, j) = gg%kik_x(ii, jj, n)
        dg%d_kik_y(i, j) = gg%kik_y(ii, jj, n)

        dg%d_gco(:, :, i, j) = gg%ext_g_co(:, :, ii, jj, n)
        dg%d_gct(:, :, i, j) = gg%ext_g_ctr(:, :, ii, jj, n)

        dg%d_c2l(:, :, i, j) = gg%ext_c2l(:, :, ii, jj, n)
        dg%d_l2c(:, :, i, j) = gg%ext_l2c(:, :, ii, jj, n)

        dg%d_sina(i, j) = gg%ext_sina(ii, jj, n)
        dg%d_cosa(i, j) = gg%ext_cosa(ii, jj, n)

        !--- length
        ii = i*2 - 1; jj = j*2 - 1
        dg%d_dx(i, j) = gg%ext_dx(ii, jj, n) + gg%ext_dx(ii + 1, jj, n)
        ii = i*2; jj = j*2 - 1 - 1
        dg%d_dy(i, j) = gg%ext_dy(ii, jj, n) + gg%ext_dy(ii, jj + 1, n)
        ii = i*2 - 1; jj = j*2 - 1 - 1
        dg%d_da(i, j) = gg%ext_da(ii, jj, n) + gg%ext_da(ii, jj + 1, n) + &
                        gg%ext_da(ii + 1, jj, n) + gg%ext_da(ii + 1, jj + 1, n)
        dg%d_rda(i, j) = 1.0/dg%c_da(i, j)
        dg%d_rdx(i, j) = 1.0/dg%c_dx(i, j)

      end do
    end do

    !--- k2e
    do j = jsd, jed
      do i = isd, ied + 1
        dg%k2e_loc_d_x(i, j) = gg%k2e_loc_d_x(i, j, n)
        dg%k2e_coef_d_x(:, i, j) = gg%k2e_coef_d_x(:, i, j, n)
      end do
    end do
    do j = jsd, jed + 1
      do i = isd, ied
        dg%k2e_loc_d_y(i, j) = gg%k2e_loc_d_y(i, j, n)
        dg%k2e_coef_d_y(:, i, j) = gg%k2e_coef_d_y(:, i, j, n)
      end do
    end do

    ! Set edge pe logicals
    if (is == 1) dg%rmp_w = .true.
    if (ie == npx - 1) dg%rmp_e = .true.
    if (js == 1) dg%rmp_s = .true.
    if (je == npy - 1) dg%rmp_n = .true.
    ! Set corner pe logicals
    dg%rmp_sw = (dg%rmp_s .and. dg%rmp_w)
    dg%rmp_se = (dg%rmp_s .and. dg%rmp_e)
    dg%rmp_nw = (dg%rmp_n .and. dg%rmp_w)
    dg%rmp_ne = (dg%rmp_n .and. dg%rmp_e)

    if (mpp_pe() == mpp_root_pe()) print *, "DUO GRID variables: DONE!"
  end subroutine duogrid_init_with_global_grid
  !===========================================================================

  !===========================================================================
  !> @brief update scalar block halo and remap to ext grid when needed
  subroutine ext_scalar_2d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(INOUT) :: dg
    type(domain2d), intent(INOUT)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: istag, jstag
    real, dimension(:, :) :: var

    integer, dimension(:, :), allocatable :: k2e_loc_u
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_u
    integer :: isd, ied, jsd, jed, i, j

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    !--- update domains
   ! if (istag == 1 .and. jstag == 1) then
   !   call mpp_update_domains(var, domain, complete=.true., position=NORTH + EAST)
   ! end if
   ! if (istag == 0 .and. jstag == 1) then
   !   call mpp_update_domains(var, domain, complete=.true., position=NORTH)
   ! end if
   ! if (istag == 1 .and. jstag == 0) then
   !   call mpp_update_domains(var, domain, complete=.true., position=EAST)
   ! end if
    if (istag == 0 .and. jstag == 0) then
      call mpp_update_domains(var, domain, complete=.true.)
    else
      call mpp_error(FATAL, 'ext_scalar_2d for istag and jstag /=0 not implemented')
    end if

    allocate (k2e_coef_u(dg%k2e_nord, isd:ied, jsd:jed))
    allocate (k2e_loc_u(isd:ied, jsd:jed))
    do i = isd, ied
    do j = jsd, jed
      k2e_loc_u(i, j) = dg%k2e_loc(i, j)
      k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
    end do
    end do
    call cube_rmp(var, dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
    call fill_corners_domain_decomp(var, dg, bd, domain, istag, jstag)
    call fill_corner_region(var, bd, dg, istag, jstag)

    deallocate (k2e_coef_u)
    deallocate (k2e_loc_u)
  end subroutine ext_scalar_2d
  !===========================================================================
  !> @brief update scalar block halo and remap to ext grid when needed
  subroutine ext_scalar_3d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(INOUT) :: dg
    type(domain2d), intent(INOUT)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: istag, jstag
    real, dimension(:, :, :) :: var
    !--- local
    integer :: i, j, isd, ied, jsd, jed, k, npz
    integer :: dims(3)

    integer, dimension(:, :), allocatable :: k2e_loc_u
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_u

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    if ((istag == 0 .and. jstag == 1) .or. (istag == 1 .and. jstag == 0)) then
      call mpp_error(FATAL, 'ext_scalar_3d for stag fields not implemented')
    end if

    !--- update domains
    if (istag == 1 .and. jstag == 1) then
      allocate (k2e_coef_u(dg%k2e_nord, isd:ied + 1, jsd:jed + 1))
      allocate (k2e_loc_u(isd:ied + 1, jsd:jed + 1))
!      k2e_loc_u = dg%k2e_loc_b
!      k2e_coef_u = dg%k2e_coef_b

      do i = isd, ied + 1
      do j = jsd, jed + 1
        k2e_loc_u(i, j) = dg%k2e_loc_b(i, j)
        k2e_coef_u(:, i, j) = dg%k2e_coef_b(:, i, j)
      end do
      end do
      call mpp_update_domains(var, domain, complete=.true., position=NORTH + EAST)
    end if

    if (istag == 0 .and. jstag == 0) then
      allocate (k2e_coef_u(dg%k2e_nord, isd:ied, jsd:jed))
      allocate (k2e_loc_u(isd:ied, jsd:jed))
      do i = isd, ied
      do j = jsd, jed
        k2e_loc_u(i, j) = dg%k2e_loc(i, j)
        k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
      end do
      end do
      !k2e_loc_u = dg%k2e_loc
      !k2e_coef_u = dg%k2e_coef
      call mpp_update_domains(var, domain, complete=.true.)
    end if

    !--- remap to ext
    dims = shape(var)
    npz = dims(3)

    do k = 1, npz
      call cube_rmp(var(:, :, k), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
    end do
    call fill_corners_domain_decomp(var, dg, bd, domain, istag, jstag)
    call fill_corner_region(var, bd, dg, istag, jstag)

    deallocate (k2e_coef_u)
    deallocate (k2e_loc_u)
  end subroutine ext_scalar_3d
  !===========================================================================
  !> @brief update scalar block halo and remap to ext grid when needed
  subroutine ext_scalar_4d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(INOUT) :: dg
    type(domain2d), intent(INOUT)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: istag, jstag
    real, dimension(:, :, :, :) :: var
    !--- local
    integer :: k, n, npz, nq
    integer :: i, j, isd, ied, jsd, jed
    integer :: dims(4)
    integer, dimension(:, :), allocatable :: k2e_loc_u
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_u

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    !--- update domains
    if (istag == 0 .and. jstag == 0) then
      allocate (k2e_coef_u(dg%k2e_nord, isd:ied, jsd:jed))
      allocate (k2e_loc_u(isd:ied, jsd:jed))
      do i = isd, ied
      do j = jsd, jed
        k2e_loc_u(i, j) = dg%k2e_loc(i, j)
        k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
      end do
      end do

      call mpp_update_domains(var, domain, complete=.true.)
    else
      call mpp_error(FATAL, 'ext_scalar_4d for istag and jstag /=0 not implemented')
    end if

    !--- remap to ext
    dims = shape(var)
    npz = dims(3)
    nq = dims(4)

    do n = 1, nq
      do k = 1, npz
        call cube_rmp(var(:, :, k, n), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
      end do
      call fill_corners_domain_decomp(var(:, :, :, n), dg, bd, domain, istag, jstag)
    end do
    call fill_corner_region(var, bd, dg, istag, jstag)

    if (istag == 0 .and. jstag == 0) then
      deallocate (k2e_coef_u)
      deallocate (k2e_loc_u)
    end if
  end subroutine ext_scalar_4d
  !===========================================================================
  !> @brief update vector block halo and remap to ext grid when needed
  subroutine ext_vector(u_in, v_in, dg, bd, domain, gridstruct, flagstruct, ieu_stag, jeu_stag, iev_stag, jev_stag)
    type(duogrid_type), intent(INOUT) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: ieu_stag, jeu_stag, iev_stag, jev_stag
    type(domain2d), intent(INOUT)         :: domain
    type(fv_flags_type), intent(IN) :: flagstruct
    type(fv_grid_type), intent(INOUT) :: gridstruct
    real, dimension(bd%isd:bd%ied + ieu_stag, bd%jsd:bd%jed + jeu_stag, flagstruct%npz) :: u_in
    real, dimension(bd%isd:bd%ied + iev_stag, bd%jsd:bd%jed + jev_stag, flagstruct%npz) :: v_in
    !--- local
    real, dimension(bd%isd:bd%ied + ieu_stag, bd%jsd:bd%jed + jeu_stag, flagstruct%npz) :: u
    real, dimension(bd%isd:bd%ied + iev_stag, bd%jsd:bd%jed + jev_stag, flagstruct%npz) :: v
    real, dimension(bd%isd - 1:bd%ied + 1 + ieu_stag, bd%jsd - 1:bd%jed + 1 + jeu_stag, flagstruct%npz) :: up1
    real, dimension(bd%isd - 1:bd%ied + 1 + iev_stag, bd%jsd - 1:bd%jed + 1 + jev_stag, flagstruct%npz) :: vp1
    integer :: is, ie, js, je, isd, ied, jsd, jed, isdp1, iedp1, jsdp1, jedp1, ng, npz
    integer :: i, j, k
    real, dimension(:, :, :), allocatable :: ull, vll, ullp1, vllp1

    real, dimension(:, :, :, :), allocatable :: c2l, l2c
    integer, dimension(:, :), allocatable :: k2e_loc_u
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_u
    integer, dimension(:, :), allocatable :: k2e_loc_v
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_v

    !--- assign parameters
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    isdp1 = bd%isd - 1
    iedp1 = bd%ied + 1
    jsdp1 = bd%jsd - 1
    jedp1 = bd%jed + 1
    ng = bd%ng

    npz = flagstruct%npz

    allocate (ull(isd:ied, jsd:jed, npz))
    allocate (vll(isd:ied, jsd:jed, npz))
    allocate (ullp1(isd - 1:ied + 1, jsd - 1:jed + 1, npz))
    allocate (vllp1(isd - 1:ied + 1, jsd - 1:jed + 1, npz))
    ull = -99999.
    vll = -99999.
    allocate (c2l(2, 2, isd:ied + iev_stag, jsd:jed + jev_stag))
    allocate (l2c(2, 2, isd:ied + iev_stag, jsd:jed + jev_stag))
    c2l = -99999.
    l2c = -99999.

    ! Since we are transforming C/D variables to A using c2l, we only use the A-grid
    ! remapping coeff. C/D extensions coeff are available, using them however with the current
    ! logic produced more noised compared to the first method.
    ! To CHECK:  do C/D on the fly without passing by A.

!Agrid
    if ((ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0)) then  !AGRID CGRID
      allocate (k2e_coef_u(dg%k2e_nord, isd:ied, jsd:jed))
      allocate (k2e_loc_u(isd:ied, jsd:jed))
      allocate (k2e_coef_v(dg%k2e_nord, isd:ied, jsd:jed))
      allocate (k2e_loc_v(isd:ied, jsd:jed))
      k2e_loc_u = -99999
      k2e_coef_u = -99999
      k2e_loc_v = -99999
      k2e_coef_v = -99999
      do j = jsd, jed
      do i = isd, ied
        c2l(:, :, i, j) = dg%a_c2l(:, :, i, j)
        l2c(:, :, i, j) = dg%a_l2c(:, :, i, j)
        k2e_loc_u(i, j) = dg%k2e_loc(i, j)
        k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
        k2e_loc_v(i, j) = dg%k2e_loc(i, j)
        k2e_coef_v(:, i, j) = dg%k2e_coef(:, i, j)
      end do
      end do

    else

      allocate (k2e_coef_u(dg%k2e_nord, isd - 1:ied + 1, jsd - 1:jed + 1))
      allocate (k2e_loc_u(isd - 1:ied + 1, jsd - 1:jed + 1))
      allocate (k2e_coef_v(dg%k2e_nord, isd - 1:ied + 1, jsd - 1:jed + 1))
      allocate (k2e_loc_v(isd - 1:ied + 1, jsd - 1:jed + 1))
      k2e_loc_u = -99999
      k2e_coef_u = -99999
      k2e_loc_v = -99999
      k2e_coef_v = -99999
      do j = jsd - 1, jed + 1
      do i = isd - 1, ied + 1
        !c2l (:,:,i,j) = dg%a_c2l(:,:,i,j)
        !l2c(:,:,i,j) = dg%a_l2c(:,:,i,j)
        k2e_loc_u(i, j) = dg%k2e_loc(i, j)
        k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
        k2e_loc_v(i, j) = dg%k2e_loc(i, j)
        k2e_coef_v(:, i, j) = dg%k2e_coef(:, i, j)
      end do
      end do

    end if
    if (ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0) then  !AGRID
      do k = 1, npz
        do j = js, je
          do i = is, ie
            ull(i, j, k) = c2l(1, 1, i, j)*u_in(i, j, k) + c2l(1, 2, i, j)*v_in(i, j, k)
          end do
        end do
        do j = js, je
          do i = is, ie
            vll(i, j, k) = c2l(2, 1, i, j)*u_in(i, j, k) + c2l(2, 2, i, j)*v_in(i, j, k)
          end do
        end do
      end do
    end if

    if (ieu_stag == 0 .and. jeu_stag == 1 .and. iev_stag == 1 .and. jev_stag == 0) then  !DGRID
      call mpp_update_domains(u_in, v_in, domain, gridtype=DGRID_NE) ! update interior pe first
      !call c2l_ord2(u_in, v_in, ull, vll, gridstruct, npz, 0, bd, .false.) !zonal-merdional vel
      call c2l_ord2(u_in, v_in, ull, vll, gridstruct, npz, 0, bd, .true.) !zonal-merdional vel

      do k = 1, npz
        do j = js, je + 1
          do i = is, ie + 1
            ullp1(i, j, k) = ull(i, j, k)
          end do
        end do
        do j = js, je + 1
          do i = is, ie + 1
            vllp1(i, j, k) = vll(i, j, k)
          end do
        end do
      end do

    end if

    if (ieu_stag == 1 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 1) then  !CGRID
      call mpp_update_domains(u_in, v_in, domain, gridtype=CGRID_NE) !update interior pe first
      call c2l_ord2_cgrid(u_in, v_in, ull, vll, gridstruct, npz, 0, bd, .true.) !zonal-merdional vel
      do k = 1, npz
        do j = js, je + 1
          do i = is, ie + 1
            ullp1(i, j, k) = ull(i, j, k)
          end do
        end do
        do j = js, je + 1
          do i = is, ie + 1
            vllp1(i, j, k) = vll(i, j, k)
          end do
        end do
      end do

    end if  !C-grid

    !Update latlonwinds haloes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0)) then  ! AGRID
      call mpp_update_domains(ull, domain, complete=.false.)
      call mpp_update_domains(vll, domain, complete=.true.)
    else !CDGRID
      call mpp_update_domains(ullp1, dg%domain_for_duo, complete=.false.)
      call mpp_update_domains(vllp1, dg%domain_for_duo, complete=.true.)
    end if

    !--- remap to ext
      !!!!!!!!!!!!!!!!!!
    if ((ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0)) then   !AGRID
      do k = 1, npz
        call cube_rmp(ull(:, :, k), dg, bd, domain, 0, 0, k2e_loc_u, k2e_coef_u)
        call cube_rmp(vll(:, :, k), dg, bd, domain, 0, 0, k2e_loc_v, k2e_coef_v)
      end do
      call fill_corners_domain_decomp(ull, dg, bd, domain, 0, 0)
      call fill_corners_domain_decomp(vll, dg, bd, domain, 0, 0)

    else !CDGRID

      do k = 1, npz
        call cube_rmp(ullp1(:, :, k), dg, dg%bd, dg%domain_for_duo, 0, 0, k2e_loc_u, k2e_coef_u)
        call cube_rmp(vllp1(:, :, k), dg, dg%bd, dg%domain_for_duo, 0, 0, k2e_loc_v, k2e_coef_v)
      end do
      call fill_corners_domain_decomp(ullp1, dg, dg%bd, dg%domain_for_duo, 0, 0)
      call fill_corners_domain_decomp(vllp1, dg, dg%bd, dg%domain_for_duo, 0, 0)

      if (ieu_stag == 1 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 1) then  !CGRID
        call cubed_a2c_halo(npz, ullp1, vllp1, up1, vp1, dg)
      else
        call cubed_a2d_halo(npz, ullp1, vllp1, up1, vp1, dg)
      end if

      do k = 1, npz
      do j = jsd, jed + jeu_stag
      do i = isd, ied + ieu_stag
        u(i, j, k) = up1(i, j, k)
      end do
      end do
      do j = jsd, jed + jev_stag
      do i = isd, ied + iev_stag
        v(i, j, k) = vp1(i, j, k)
      end do
      end do
      end do
    end if

    if (ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0) then  !AGRID
      !--- convert to wind
      ! need to update the interior pes as well, so no switch
      ! south
      do k = 1, npz
        do j = jsd, js - 1
          do i = isd, ied
            u_in(i, j, k) = l2c(1, 1, i, j)*ull(i, j, k) + l2c(1, 2, i, j)*vll(i, j, k)
          end do
        end do
        do j = jsd, js - 1
          do i = isd, ied
            v_in(i, j, k) = l2c(2, 1, i, j)*ull(i, j, k) + l2c(2, 2, i, j)*vll(i, j, k)
          end do
        end do
      end do

      ! north
      do k = 1, npz
        do j = je + 1, jed
          do i = isd, ied
            u_in(i, j, k) = l2c(1, 1, i, j)*ull(i, j, k) + l2c(1, 2, i, j)*vll(i, j, k)
          end do
        end do
        do j = je + 1, jed
          do i = isd, ied
            v_in(i, j, k) = l2c(2, 1, i, j)*ull(i, j, k) + l2c(2, 2, i, j)*vll(i, j, k)
          end do
        end do
      end do

      ! west
      do k = 1, npz
        do j = js, je
          do i = isd, is - 1
            u_in(i, j, k) = l2c(1, 1, i, j)*ull(i, j, k) + l2c(1, 2, i, j)*vll(i, j, k)
          end do
        end do
        do j = js, je
          do i = isd, is - 1
            v_in(i, j, k) = l2c(2, 1, i, j)*ull(i, j, k) + l2c(2, 2, i, j)*vll(i, j, k)
          end do
        end do
      end do

      ! east
      do k = 1, npz
        do j = js, je
          do i = ie + 1, ied
            u_in(i, j, k) = l2c(1, 1, i, j)*ull(i, j, k) + l2c(1, 2, i, j)*vll(i, j, k)
          end do
        end do
        do j = js, je
          !do i = ie+ieu_stag+1,ied+iev_stag
          do i = ie + 1, ied
            v_in(i, j, k) = l2c(2, 1, i, j)*ull(i, j, k) + l2c(2, 2, i, j)*vll(i, j, k)
          end do
        end do
      end do

    else ! CD

      if (dg%rmp_s) then
        ! south
        do k = 1, npz
          do j = jsd, js - 1
            do i = isd, ied + ieu_stag
              u_in(i, j, k) = u(i, j, k)
            end do
          end do
          do j = jsd, js - 1
            do i = isd, ied + iev_stag
              v_in(i, j, k) = v(i, j, k)
            end do
          end do
        end do
      end if

      if (dg%rmp_n) then
        ! north
        do k = 1, npz
          do j = je + jeu_stag + 1, jed + jeu_stag
            do i = isd, ied + ieu_stag
              u_in(i, j, k) = u(i, j, k)
            end do
          end do
          !do j = je+jeu_stag+1,jed+jev_stag
          do j = je + jev_stag + 1, jed + jev_stag
            do i = isd, ied + iev_stag
              v_in(i, j, k) = v(i, j, k)
            end do
          end do
        end do
      end if

      if (dg%rmp_w) then
        ! west
        do k = 1, npz
          do j = jsd, jed + jeu_stag
            do i = isd, is - 1
              u_in(i, j, k) = u(i, j, k)
            end do
          end do
          do j = jsd, jed + jev_stag
            do i = isd, is - 1
              v_in(i, j, k) = v(i, j, k)
            end do
          end do
        end do
      end if

      if (dg%rmp_e) then
        ! east
        do k = 1, npz
          do j = jsd, jed + jeu_stag
            do i = ie + ieu_stag + 1, ied + ieu_stag
              u_in(i, j, k) = u(i, j, k)
            end do
          end do
          do j = jsd, jed + jev_stag
            !do i = ie+ieu_stag+1,ied+iev_stag
            do i = ie + iev_stag + 1, ied + iev_stag
              v_in(i, j, k) = v(i, j, k)
            end do
          end do
        end do
      end if

    end if !if C-grid

!Fill corner region
!!!!!!!!!!!!!!!!!!!
    call fill_corner_region(u_in, bd, dg, ieu_stag, jeu_stag)
    call fill_corner_region(v_in, bd, dg, iev_stag, jev_stag)

    !--- deallocate wk arrays
    deallocate (ull)
    deallocate (vll)
    deallocate (ullp1)
    deallocate (vllp1)
    deallocate (k2e_coef_u)
    deallocate (k2e_coef_v)
    deallocate (k2e_loc_u)
    deallocate (k2e_loc_v)
    deallocate (c2l)
    deallocate (l2c)

  end subroutine ext_vector
  !===========================================================================
  subroutine cube_rmp(var, dg, bd, domain, istag, jstag, k2e_loc, k2e_coef)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(domain2d), intent(INOUT) :: domain
    integer, intent(in) :: istag, jstag
    !integer, intent(in) :: ext_x, ext_y
    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: var
    ! local
    real, dimension(:, :), allocatable :: var_kik
    logical :: rmp_w, rmp_e, rmp_s, rmp_n, rmp_sw, rmp_se, rmp_ne, rmp_nw
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: loc, lo, offset, jmin, jmax
    !   real, dimension(:,:), allocatable :: k2e_loc
    !   real, dimension(:,:,:), allocatable :: k2e_coef
    !  integer, dimension(bd%isd:bd%ied+ext_x, bd%jsd:bd%jed+ext_y) :: k2e_loc
    !  real, dimension(dg%k2e_nord,bd%isd:bd%ied+ext_x, bd%jsd:bd%jed+ext_y) :: k2e_coef
    integer, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: k2e_loc
    real(kind=R_GRID), dimension(dg%k2e_nord, bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: k2e_coef

    !--- assign parameters
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    ng = bd%ng

    rmp_w = dg%rmp_w
    rmp_e = dg%rmp_e
    rmp_s = dg%rmp_s
    rmp_n = dg%rmp_n

    rmp_sw = dg%rmp_sw
    rmp_se = dg%rmp_se
    rmp_ne = dg%rmp_ne
    rmp_nw = dg%rmp_nw

    offset = dg%k2e_nord/2

    !--- allocata wk array
    allocate (var_kik(isd:ied + istag, jsd:jed + jstag))

    !--- copy neighbor to var_kik
    do ii = 1, ng
      do i = isd, ied + istag
        j = js - ii
        var_kik(i, j) = var(i, j)
        j = je + jstag + ii
        var_kik(i, j) = var(i, j)
      end do
      do j = js, je + jstag
        i = is - ii
        var_kik(i, j) = var(i, j)
        i = ie + istag + ii
        var_kik(i, j) = var(i, j)
      end do
    end do

    !--- south
    if (rmp_s) then
      do ii = 1, ng

        j = js - ii

        do i = is, ie + istag

          loc = k2e_loc(i, j)
          lo = loc - offset

          var(i, j) = 0.
          do n = 1, dg%k2e_nord
            var(i, j) = var(i, j) + var_kik(lo + n, j)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !--- north
    if (rmp_n) then
      do ii = 1, ng

        j = je + jstag + ii

        do i = is, ie + istag

          loc = k2e_loc(i, j)
          lo = loc - offset

          var(i, j) = 0.
          do n = 1, dg%k2e_nord
            var(i, j) = var(i, j) + var_kik(lo + n, j)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !--- west
    if (rmp_w) then
      jmin = js
      jmax = je + jstag
! For few pe layouts we can extend the rmp process
! to fill the pe corners at the edges; but for most of them
! the lo+n goes beyond the var_kik dimensions, meaning the remapped
! value at the corners needs to grab data outside of what is available on its pe
! so we end up filling the pe corner region from the neighbooring pe.

      do ii = 1, ng

        i = is - ii

        !do j = js,je+jstag
        do j = jmin, jmax

          loc = k2e_loc(i, j)
          lo = loc - offset

          var(i, j) = 0.
          do n = 1, dg%k2e_nord
            var(i, j) = var(i, j) + var_kik(i, lo + n)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !--- east
    if (rmp_e) then
      do ii = 1, ng

        i = ie + istag + ii

        do j = js, je + jstag

          loc = k2e_loc(i, j)
          lo = loc - offset

          var(i, j) = 0.
          do n = 1, dg%k2e_nord
            var(i, j) = var(i, j) + var_kik(i, lo + n)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !move this outside of cube_rmp to do the comms 3D
    !call fill_corners_domain_decomp(var, dg, bd, domain, istag, jstag)

    deallocate (var_kik)

  end subroutine cube_rmp
  !===========================================================================

! The subroutine fill_corners_domain_decomp below fills the corner regions
! of the pes at the edges from neighboring pes.
! Below is a schematic of how pes X and Y corners are filled at the west cubed sphere edge
! first, in the first code section of rmp_w, we fill the corner region y of peX from peY data (region y).
! then we fill the corner region x of pe Y from pe X region X.
!
! we loop over pes from lowest up to higher ranks, simultaneously over the six tiles.
! if there is a pe Z above pe X, we first do peY/X comm then pe X/Z comm
! Same logic applies for east, north, south. Works for all grid staggering :-)
!
!        cubed edge
!           |
!           |
!     x x x |
!     x x x |    pe X
!     x x x |
!    -------------- pe X lower edge - pe Y upper edge
!     y y y |
!     y y y |
!     y y y |
!           |
!
!
!           |
!           |
!     x x x |
!     x x x |
!     x x x |
!    -------------- pe Y upper edge - pe X lower edge
!     y y y |
!     y y y |    pe Y
!     y y y |
!           |
!
!
!

  subroutine fill_corners_domain_decomp_2d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(domain2d), intent(IN)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:) :: var
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, layoutx, layouty, tope, frompe, gid, kk
    integer :: tile(1), npes_per_tile
    real, dimension(1:bd%ng, 1:bd%ng) :: corner

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    ng = bd%ng

    tile = mpp_get_tile_id(domain)
    npes_per_tile = mpp_npes()/6.
    gid = mpp_pe()
    layoutx = dg%layout(1)
    layouty = dg%layout(2)

    if (layoutx + layouty > 2) then

!!!!!!!!!!!!!!!!!!
!!!!!! WEST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_w .and. layouty > 1) then

        do kk = 1, layouty - 1
          tope = (kk)*layoutx + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(isd + i - 1, je - ng + j)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(isd + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          frompe = (kk)*layoutx + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(isd + i - 1, js + jstag + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(isd + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

        end do !enddo kk

      end if !endif rmp_w

!!!!!!!!!!!!!!!!!!
!!!!!! EAST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_e .and. layouty > 1) then

        do kk = 1, layouty - 1

          tope = (kk)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1
          frompe = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1

          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(ie + istag + 1 + i - 1, je - ng + j)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(ie + istag + 1 + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          frompe = (kk)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1
          tope = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1

          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(ie + istag + 1 + i - 1, js + jstag + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

        end do
      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!SOUTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_s .and. layoutx > 1) then

        do kk = 1, layoutx - 1

          tope = (kk) + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
          do i = 1, ng
            do j = 1, ng
              corner(i, j) = var(ie + istag - ng + i, jsd + j - 1)
              corner(i, j) = var(ie - ng + i, jsd + j - 1)
            end do
          end do
          call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(isd + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          frompe = (kk) + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
          do i = 1, ng
            do j = 1, ng
              corner(i, j) = var(is + istag + i - 1, jsd + j - 1)
            end do
          end do
          call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(ie + istag + 1 + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

        end do
      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!NORTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_n .and. layoutx > 1) then

        do kk = 1, layoutx - 1

          tope = (kk) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(ie - ng + i, je + jstag + 1 + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(isd + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

          frompe = (kk) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do i = 1, ng
              do j = 1, ng
                corner(i, j) = var(is + istag + i - 1, je + jstag + 1 + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do i = 1, ng
              do j = 1, ng
                var(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

        end do ! kk
      end if ! rmpn

    end if ! layoux+layouty>2

  end subroutine fill_corners_domain_decomp_2d

  subroutine fill_corners_domain_decomp_3d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(domain2d), intent(IN)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, 1:) :: var
    real, allocatable, dimension(:, :, :) :: corner
    ! real, dimension(1:bd%ng, 1:bd%ng) :: corner

    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, layoutx, layouty, tope, frompe, gid, k, kk
    integer :: tile(1), npes_per_tile, npz, dims(3)

    dims = shape(var)
    npz = dims(3)

    allocate (corner(1:bd%ng, 1:bd%ng, 1:npz))

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    ng = bd%ng

    tile = mpp_get_tile_id(domain)
    npes_per_tile = mpp_npes()/6.
    gid = mpp_pe()
    layoutx = dg%layout(1)
    layouty = dg%layout(2)

    if (layoutx + layouty > 2) then

!!!!!!!!!!!!!!!!!!
!!!!!! WEST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_w .and. layouty > 1) then

        do kk = 1, layouty - 1
          tope = (kk)*layoutx + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(isd + i - 1, je - ng + j, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(isd + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          frompe = (kk)*layoutx + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(isd + i - 1, js + jstag + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(isd + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

        end do !enddo kk

      end if !endif rmp_w

!!!!!!!!!!!!!!!!!!
!!!!!! EAST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_e .and. layouty > 1) then

        do kk = 1, layouty - 1

          tope = (kk)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1
          frompe = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(ie + istag + 1 + i - 1, je - ng + j, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(ie + istag + 1 + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          frompe = (kk)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1
          tope = (kk - 1)*layoutx + (tile(1) - 1)*npes_per_tile + layoutx - 1

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(ie + istag + 1 + i - 1, js + jstag + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

        end do
      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!SOUTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_s .and. layoutx > 1) then

        do kk = 1, layoutx - 1

          tope = (kk) + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(ie - ng + i, jsd + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(isd + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          frompe = (kk) + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(is + istag + i - 1, jsd + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(ie + istag + 1 + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

        end do
      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!NORTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_n .and. layoutx > 1) then

        do kk = 1, layoutx - 1

          tope = (kk) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          frompe = (kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(ie - ng + i, je + jstag + 1 + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(isd + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          frompe = (kk) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile
          tope = (kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1)*npes_per_tile

          if (gid == frompe) then
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  corner(i, j, k) = var(is + istag + i - 1, je + jstag + 1 + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=tope)
          end if

          if (gid == tope) then
            call mpp_recv(corner, size(corner), from_pe=frompe)
            do k = 1, npz
              do i = 1, ng
                do j = 1, ng
                  var(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

        end do ! kk
      end if ! rmpn

    end if ! layoux+layouty>2

    deallocate (corner)

  end subroutine fill_corners_domain_decomp_3d

! subroutine fill_corner_region will fill the corner regions of pes lying on
! the corners of the cubes sphere nw, ne, se, sw.
! the data is extrapolated from the haloe regions present on the pe following
! a logic similar to ZA22's using a lagrangian interpolation function.

! Ideally, values should be interoplated (no extrapolation needed) using the exact
! ZA22 formulation; however, this is not done here yet and requires more complicated work.

! In addition, using the lagrangian interpolation function for extrapolation
! purposes is not highly recommended in general; however, it is tolerable in this case given that
! the extrapolation is performed for couple points just outside of the domain on which the data reside.

  subroutine fill_corner_region_2d(vel, bd, dg, istag, jstag)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(in) :: dg
    integer, intent(IN) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:) :: vel

    integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
    !integer :: interporder = 3 !order is interporder+1

    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: veltemp
    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: veltempp
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    !lastpoint in compute domain
    ie = ie + istag
    je = je + jstag

    ! NE corner
!     |   1--2--3
!     |   |     |
!     |   4  5  6
!     |   |     |
!     |   7--8--9
!     -------------

    if (dg%rmp_ne) then
      call lagrange_poly_interp(vel, ie + 1, je + 2, bd, dg, 'X+', istag, jstag, interporder) !4
      call lagrange_poly_interp(vel, ie + 1, je + 3, bd, dg, 'X+', istag, jstag, interporder) !1
      call lagrange_poly_interp(vel, ie + 2, je + 3, bd, dg, 'X+', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, ie + 2, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, ie + 3, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !9
      call lagrange_poly_interp(vel, ie + 3, je + 2, bd, dg, 'Y+', istag, jstag, interporder) !6

      veltemp = vel
      veltempp = vel

      !diag
      i = ie + 1
      j = je + 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))

      i = ie + 3
      j = je + 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))

      i = ie + 2
      j = je + 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! NW
!        1--2--3 |
!        |     | |
!        4  5  6 |
!        |     | |
!        7--8--9 |
!        ---------

    if (dg%rmp_nw) then
      call lagrange_poly_interp(vel, is - 1, je + 2, bd, dg, 'X-', istag, jstag, interporder) !6
      call lagrange_poly_interp(vel, is - 1, je + 3, bd, dg, 'X-', istag, jstag, interporder) !3
      call lagrange_poly_interp(vel, is - 2, je + 3, bd, dg, 'X-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, is - 2, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, is - 3, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !7
      call lagrange_poly_interp(vel, is - 3, je + 2, bd, dg, 'Y+', istag, jstag, interporder) !4

      veltemp = vel
      veltempp = vel

      !diag
      i = is - 1
      j = je + 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))

      i = is - 3
      j = je + 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))

      i = is - 2
      j = je + 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SE corner
!     --------------
!     |   1--2--3
!     |   |     |
!     |   4  5  6
!     |   |     |
!     |   7--8--9
!
    if (dg%rmp_se) then
      call lagrange_poly_interp(vel, ie + 1, js - 2, bd, dg, 'X+', istag, jstag, interporder) !4
      call lagrange_poly_interp(vel, ie + 1, js - 3, bd, dg, 'X+', istag, jstag, interporder) !7
      call lagrange_poly_interp(vel, ie + 2, js - 3, bd, dg, 'X+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, ie + 2, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, ie + 3, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !3
      call lagrange_poly_interp(vel, ie + 3, js - 2, bd, dg, 'Y-', istag, jstag, interporder) !6

      veltemp = vel
      veltempp = vel

      !diag
      i = ie + 1
      j = js - 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
      i = ie + 3
      j = js - 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
      i = ie + 2
      j = js - 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SW corner
!     -------------
!        1--2--3  |
!        |     |  |
!        4  5  6  |
!        |     |  |
!        7--8--9  |
!

    if (dg%rmp_sw) then
      call lagrange_poly_interp(vel, is - 1, js - 2, bd, dg, 'X-', istag, jstag, interporder) !6
      call lagrange_poly_interp(vel, is - 1, js - 3, bd, dg, 'X-', istag, jstag, interporder) !9
      call lagrange_poly_interp(vel, is - 2, js - 3, bd, dg, 'X-', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, is - 2, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, is - 3, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !1
      call lagrange_poly_interp(vel, is - 3, js - 2, bd, dg, 'Y-', istag, jstag, interporder) !4

      veltemp = vel
      veltempp = vel

      !diag
      i = is - 1
      j = js - 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
      i = is - 3
      j = js - 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
      i = is - 2
      j = js - 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j) = 0.5*(veltemp(i, j) + veltempp(i, j))
    end if

  end subroutine fill_corner_region_2d

  subroutine fill_corner_region_3d(vel, bd, dg, istag, jstag)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(in) :: dg
    integer, intent(IN) :: istag, jstag
    ! real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: vel
    real, dimension(bd%isd:, bd%jsd:, 1:) :: vel

    integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, npz, dims(3)
    !integer :: interporder = 3 !order is interporder+1

    real, allocatable, dimension(:, :, :) :: veltemp
    real, allocatable, dimension(:, :, :) :: veltempp

    dims = shape(vel)
    npz = dims(3)

! any other solutions to get the dimensions right?
    allocate (veltemp(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag, 1:npz))
    allocate (veltempp(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag, 1:npz))

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    !lastpoint in compute domain
    ie = ie + istag
    je = je + jstag

    ! NE corner
!     |   1--2--3
!     |   |     |
!     |   4  5  6
!     |   |     |
!     |   7--8--9
!     -------------

    if (dg%rmp_ne) then
      call lagrange_poly_interp(vel, ie + 1, je + 2, bd, dg, 'X+', istag, jstag, interporder) !4
      call lagrange_poly_interp(vel, ie + 1, je + 3, bd, dg, 'X+', istag, jstag, interporder) !1
      call lagrange_poly_interp(vel, ie + 2, je + 3, bd, dg, 'X+', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, ie + 2, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, ie + 3, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !9
      call lagrange_poly_interp(vel, ie + 3, je + 2, bd, dg, 'Y+', istag, jstag, interporder) !6

      veltemp = vel
      veltempp = vel

      !diag
      i = ie + 1
      j = je + 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))

      i = ie + 3
      j = je + 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))

      i = ie + 2
      j = je + 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! NW
!        1--2--3 |
!        |     | |
!        4  5  6 |
!        |     | |
!        7--8--9 |
!        ---------

    if (dg%rmp_nw) then
      call lagrange_poly_interp(vel, is - 1, je + 2, bd, dg, 'X-', istag, jstag, interporder) !6
      call lagrange_poly_interp(vel, is - 1, je + 3, bd, dg, 'X-', istag, jstag, interporder) !3
      call lagrange_poly_interp(vel, is - 2, je + 3, bd, dg, 'X-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, is - 2, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, is - 3, je + 1, bd, dg, 'Y+', istag, jstag, interporder) !7
      call lagrange_poly_interp(vel, is - 3, je + 2, bd, dg, 'Y+', istag, jstag, interporder) !4

      veltemp = vel
      veltempp = vel

      !diag
      i = is - 1
      j = je + 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))

      i = is - 3
      j = je + 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))

      i = is - 2
      j = je + 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y+', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SE corner
!     --------------
!     |   1--2--3
!     |   |     |
!     |   4  5  6
!     |   |     |
!     |   7--8--9
!
    if (dg%rmp_se) then
      call lagrange_poly_interp(vel, ie + 1, js - 2, bd, dg, 'X+', istag, jstag, interporder) !4
      call lagrange_poly_interp(vel, ie + 1, js - 3, bd, dg, 'X+', istag, jstag, interporder) !7
      call lagrange_poly_interp(vel, ie + 2, js - 3, bd, dg, 'X+', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, ie + 2, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, ie + 3, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !3
      call lagrange_poly_interp(vel, ie + 3, js - 2, bd, dg, 'Y-', istag, jstag, interporder) !6

      veltemp = vel
      veltempp = vel

      !diag
      i = ie + 1
      j = js - 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
      i = ie + 3
      j = js - 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
      i = ie + 2
      j = js - 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X+', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SW corner
!     -------------
!        1--2--3  |
!        |     |  |
!        4  5  6  |
!        |     |  |
!        7--8--9  |
!

    if (dg%rmp_sw) then
      call lagrange_poly_interp(vel, is - 1, js - 2, bd, dg, 'X-', istag, jstag, interporder) !6
      call lagrange_poly_interp(vel, is - 1, js - 3, bd, dg, 'X-', istag, jstag, interporder) !9
      call lagrange_poly_interp(vel, is - 2, js - 3, bd, dg, 'X-', istag, jstag, interporder) !8
      call lagrange_poly_interp(vel, is - 2, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !2
      call lagrange_poly_interp(vel, is - 3, js - 1, bd, dg, 'Y-', istag, jstag, interporder) !1
      call lagrange_poly_interp(vel, is - 3, js - 2, bd, dg, 'Y-', istag, jstag, interporder) !4

      veltemp = vel
      veltempp = vel

      !diag
      i = is - 1
      j = js - 1
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
      i = is - 3
      j = js - 3
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
      i = is - 2
      j = js - 2
      call lagrange_poly_interp(veltemp, i, j, bd, dg, 'X-', istag, jstag, interporder)
      call lagrange_poly_interp(veltempp, i, j, bd, dg, 'Y-', istag, jstag, interporder)
      vel(i, j, :) = 0.5*(veltemp(i, j, :) + veltempp(i, j, :))
    end if

    deallocate (veltemp, veltempp)
  end subroutine fill_corner_region_3d

!!!!!!!!!!
!!slower!!
!!!!!!!!!!
!  subroutine fill_corner_region_3d(vel, bd, dg, istag, jstag)
!    type(fv_grid_bounds_type), intent(IN) :: bd
!    type(duogrid_type), intent(in) :: dg
!    integer, intent(IN) :: istag, jstag
!    real, dimension(bd%isd:, bd%jsd:, 1:) :: vel
!
!    integer :: k, npz, dims(3)
!
!    dims = shape(vel)
!    npz = dims(3)
!
!    do k = 1, npz
!      call fill_corner_region_2d(vel(:, :, k), bd, dg, istag, jstag)
!    end do
!  end subroutine fill_corner_region_3d

  subroutine fill_corner_region_4d(vel, bd, dg, istag, jstag)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(in) :: dg
    integer, intent(IN) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, 1:, 1:) :: vel

    integer :: k, npz, q, nq, dims(4)

    dims = shape(vel)
    npz = dims(3)
    nq = dims(4)

    do q = 1, nq
      call fill_corner_region_3d(vel(:, :, :, q), bd, dg, istag, jstag)
    end do

  end subroutine fill_corner_region_4d

! Below is the algorithm of the lagrangian interpolation function used to compute
! the cubed sphere corner region and the farthest haloe lines in cubed_a2d_halo and
! cubed_a2c_halo

! This will take any field of any staggering and computed the value at (i,j)
! along the given direction.

! along X+: uses data along X in the positive direction: so (i,j) is at the right side
!                ---x----x---x---(i,j)
! along X-: uses data along X in the negative direction: so (i,j) is at the left side
!           (i,j)---x----x---x---

! In general, this is made for the corner region, where data, for example at ie+2 is computed
! from data edged at bd%ie
! For the a2d/a2c halo data, an optional argument corner=.false. could be used to
! interpolate data from the closest neighbooring point
! for instance for ie+3 along X+, the last used point is ie+2.
! the interoplation order depicts how many points are used. if order=1, uses 2 points.

! compute coeff a priori
! TODO rearrange ss_store
  subroutine compute_lagrange_coeff(xp, xm, yp, ym, bd, dg)

    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(inout) :: dg
    real(kind=R_GRID), intent(out) :: xp(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: xm(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: yp(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: ym(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4)

    !local
    integer ii, iii, it, jstag, istag
    integer is, js, ie, je, i, j
    real(kind=R_GRID) :: ss, S

    do istag = 0, 1
    do jstag = 0, 1
    do i = bd%ie - interporder, bd%ied + istag
    do j = bd%jsd, bd%jed + jstag
      is = bd%is
      js = bd%js
      !  ie = bd%ie + istag
      !  je = bd%je + jstag
      do ii = bd%ie - interporder, bd%ie
        ss = 1.0
        do iii = bd%ie - interporder, bd%ie
          if (ii < iii) it = -1
          if (ii > iii) it = 1
          if (ii .eq. iii) cycle
          if (j > bd%je - 1) then
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, iii, j))/ &
                        great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j)))
          else
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
                        great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
          end if
        end do
        xp(bd%ie - ii + 1, i, j, 2*istag + jstag + 1) = ss
      end do

    end do
    end do

    do i = bd%isd, bd%is + interporder
    do j = bd%jsd, bd%jed + jstag

      do ii = is + interporder, is, -1
        ss = 1.0
        do iii = is + interporder, is, -1
          if (ii < iii) it = 1
          if (ii > iii) it = -1
          if (ii .eq. iii) cycle
          if (j > bd%je - 1) then
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, iii, j))/ &
                        great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j)))
          else
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
                        great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
          end if

        end do
        xm(is + interporder - ii + 1, i, j, 2*istag + jstag + 1) = ss
      end do
    end do
    end do

    do i = bd%isd, bd%ied + istag
    do j = bd%je - interporder, bd%jed + jstag
      do ii = bd%je - interporder, bd%je
        ss = 1.0
        do iii = bd%je - interporder, bd%je
          if (ii < iii) it = -1
          if (ii > iii) it = 1
          if (ii .eq. iii) cycle
          if (i > bd%ie - 1) then
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, i, iii))/ &
                        great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
          else
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, i - istag, iii))/ &
                        great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))

          end if

        end do
        yp(bd%je - ii + 1, i, j, 2*istag + jstag + 1) = ss
      end do
    end do
    end do

    do i = bd%isd, bd%ied + istag
    do j = bd%jsd, bd%js + interporder
      do ii = js + interporder, js, -1
        ss = 1.0
        do iii = js + interporder, js, -1
          if (ii < iii) it = 1
          if (ii > iii) it = -1
          if (ii .eq. iii) cycle
          if (i > bd%ie - 1) then
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, i, iii))/ &
                        great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
          else
            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, i - istag, iii))/ &
                        great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))
          end if

        end do
        ym(js + interporder - ii + 1, i, j, 2*istag + jstag + 1) = ss
      end do
    end do
    end do
    end do
    end do

    if (mpp_pe() == mpp_root_pe()) print *, "DUO lag coeff: DONE!"
  end subroutine compute_lagrange_coeff

!  subroutine compute_lagrange_poly_interp_2d(ss_store, i, j, bd, dg, direction, istag, jstag, order, corner)
!
!    integer, intent(in) :: i, j, istag, jstag, order
!    type(fv_grid_bounds_type), intent(IN) :: bd
!    !real, dimension(bd%isd:, bd%jsd:) :: field
!    type(duogrid_type), intent(in) :: dg
!    character(len=*), intent(in) :: direction
!    logical, optional, intent(in) :: corner
!    real(kind=R_GRID), intent(out) :: ss_store(5, 1:interporder+1, bd%isd:bd%ied+1, bd%jsd:bd%jed+1, 0:1,0:1)
!
!    !local
!    integer ii, iii, it
!    integer is, js, ie, je
!    real(kind=R_GRID) :: ss, S
!
!    is = bd%is
!    js = bd%js
!    ie = bd%ie + istag
!    je = bd%je + jstag
!    it = 1
!
!!TO check: compute weight coefficients a priori?!
!
!        do ii = bd%ie - interporder, bd%ie
!          ss = 1.0
!          do iii = bd%ie - interporder, bd%ie
!            if (ii < iii) it = -1
!            if (ii > iii) it = 1
!            if (ii .eq. iii) cycle
!            if (j > bd%je - 1) then
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, iii, j))/great_circle_dist(dg%a_pt(:, ii, j), &
!                                                                                                            dg%a_pt(:, iii, j)))
!            else
!            ss = ss*it*(great_circle_dist(dg%a_pt(:, i-istag, j-jstag), dg%a_pt(:, iii, j-jstag))/great_circle_dist(dg%a_pt(:, ii, j-jstag),   dg%a_pt(:, iii, j-jstag)))
!            end if
!          end do
!          ss_store(1,bd%ie-ii+1,i,j,istag,jstag)=ss
!        end do
!
!
!        do ii = is + order, is, -1
!          ss = 1.0
!          do iii = is + order, is, -1
!            if (ii < iii) it = 1
!            if (ii > iii) it = -1
!            if (ii .eq. iii) cycle
!            if (j > bd%je - 1) then
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, iii, j))/great_circle_dist(dg%a_pt(:, ii, j), &
!                                                                                                    dg%a_pt(:, iii, j)))
!            else
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, iii, j - jstag))/great_circle_dist(dg%a_pt(:, ii, j - jstag), &
!                                                                                                        dg%a_pt(:, iii, j - jstag)))
!            end if
!          end do
!          ss_store(2,is+order-ii+1,i,j,istag,jstag)=ss
!        end do
!
!
!        do ii = bd%je - order, bd%je
!          ss = 1.0
!          do iii = bd%je - order, bd%je
!            if (ii < iii) it = -1
!            if (ii > iii) it = 1
!            if (ii .eq. iii) cycle
!            if (i > bd%ie - 1) then
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, i, iii))/great_circle_dist(dg%a_pt(:, i, ii), &
!                                                                                                            dg%a_pt(:, i, iii)))
!            else
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i-istag, j-jstag), dg%a_pt(:, i-istag, iii))/great_circle_dist(dg%a_pt(:, i-istag, ii), &
!                                                                                                        dg%a_pt(:, i - istag, iii)))
!
!            end if
!          end do
!          ss_store(3,bd%je-ii+1,i,j,istag,jstag)=ss
!        end do
!
!        do ii = js + order, js, -1
!          ss = 1.0
!          do iii = js + order, js, -1
!            if (ii < iii) it = 1
!            if (ii > iii) it = -1
!            if (ii .eq. iii) cycle
!            if (i > bd%ie - 1) then
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, i, iii))/great_circle_dist(dg%a_pt(:, i, ii), &
!                                                                                                    dg%a_pt(:, i, iii)))
!            else
!              ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, i - istag, iii))/great_circle_dist(dg%a_pt(:, i - istag, ii), &
!                                                                                                        dg%a_pt(:, i - istag, iii)))
!            end if
!          end do
!          ss_store(4,js+order-ii+1,i,j,istag,jstag)=ss
!        end do
!
!  end subroutine compute_lagrange_poly_interp_2d

  subroutine lagrange_poly_interp_2d(field, i, j, bd, dg, direction, istag, jstag, order, corner)

    integer, intent(in) :: i, j, istag, jstag, order
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, dimension(bd%isd:, bd%jsd:) :: field
    type(duogrid_type), target, intent(in) :: dg
    character(len=*), intent(in) :: direction
    logical, optional, intent(in) :: corner

    !local
    integer ii, iii, it
    integer is, js, ie, je
    real(kind=R_GRID) :: ss, S
    !real(kind=R_GRID) :: ss_store(5, 1:interporder+1, bd%isd:bd%ied+1, bd%jsd:bd%jed+1, 0:1,0:1)

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xm
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: yp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ym

    xp => dg%xp
    xm => dg%xm
    yp => dg%yp
    ym => dg%ym

    S = 0. ! here we sum all the polynomial*solution

    is = bd%is
    js = bd%js
    ie = bd%ie + istag
    je = bd%je + jstag
    it = 1

    call timing_on('lag')
    !call compute_lagrange_poly_interp_2d(ss_store, i, j, bd, dg, direction, istag, jstag, order, corner)

!TO check: compute weight coefficients a priori?!

    if (direction == 'X+') then

      do ii = bd%ie - interporder, bd%ie
        S = S + xp(bd%ie - ii + 1, i, j, 2*istag + jstag + 1)*field(ii + istag, j)
      end do

      field(i, j) = S

    end if

    if (direction == 'X-') then

      do ii = is + order, is, -1
        S = S + xm(is + order - ii + 1, i, j, 2*istag + jstag + 1)*field(ii, j)
      end do

      field(i, j) = S

    end if

    if (direction == 'Y+') then

      do ii = bd%je - order, bd%je
        S = S + yp(bd%je - ii + 1, i, j, 2*istag + jstag + 1)*field(i, ii + jstag)
      end do

      field(i, j) = S

    end if

    if (direction == 'Y-') then

      do ii = js + order, js, -1
        S = S + ym(js + order - ii + 1, i, j, 2*istag + jstag + 1)*field(i, ii)
      end do

      field(i, j) = S

    end if
    nullify (xp, xm, yp, ym)

    call timing_off('lag')
  end subroutine lagrange_poly_interp_2d

  subroutine lagrange_poly_interp_3d(field, i, j, bd, dg, direction, istag, jstag, order, corner)
    integer, intent(in) :: i, j, istag, jstag, order
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, dimension(bd%isd:, bd%jsd:, 1:) :: field
    type(duogrid_type), target, intent(in) :: dg
    character(len=*), intent(in) :: direction
    logical, optional, intent(in) :: corner

    !local
    integer ii, iii, it
    integer is, js, ie, je
    integer dims(3), k, npz
    real(kind=R_GRID) :: ss
    real, allocatable, dimension(:) :: S
    !real(kind=R_GRID) :: ss_store(5, 1:interporder+1, bd%isd:bd%ied+1, bd%jsd:bd%jed+1, 0:1,0:1)

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xm
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: yp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ym

    xp => dg%xp
    xm => dg%xm
    yp => dg%yp
    ym => dg%ym

    dims = shape(field)
    npz = dims(3)
    allocate (S(1:npz))

    S = 0. ! here we sum all the polynomial*solution

    is = bd%is
    js = bd%js
    ie = bd%ie + istag
    je = bd%je + jstag
    it = 1

    call timing_on('lag')
    !call compute_lagrange_poly_interp_2d(ss_store, i, j, bd, dg, direction, istag, jstag, order, corner)

    do k = 1, npz
      if (direction == 'X+') then

        do ii = bd%ie - interporder, bd%ie
          S(k) = S(k) + xp(bd%ie - ii + 1, i, j, 2*istag + jstag + 1)*field(ii + istag, j, k)
        end do

        field(i, j, k) = S(k)

      end if

      if (direction == 'X-') then

        do ii = is + order, is, -1
          S(k) = S(k) + xm(is + order - ii + 1, i, j, 2*istag + jstag + 1)*field(ii, j, k)
        end do

        field(i, j, k) = S(k)

      end if

      if (direction == 'Y+') then

        do ii = bd%je - order, bd%je
          S(k) = S(k) + yp(bd%je - ii + 1, i, j, 2*istag + jstag + 1)*field(i, ii + jstag, k)
        end do

        field(i, j, k) = S(k)

      end if

      if (direction == 'Y-') then

        do ii = js + order, js, -1
          S(k) = S(k) + ym(js + order - ii + 1, i, j, 2*istag + jstag + 1)*field(i, ii, k)
        end do

        field(i, j, k) = S(k)

      end if
    end do
    nullify (xp, xm, yp, ym)
    deallocate (S)

    call timing_off('lag')

  end subroutine lagrange_poly_interp_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this is much slower!?!?!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine lagrange_poly_interp_3d(field, i, j, bd, dg, direction, istag, jstag, order, corner)
!
!    integer, intent(in) :: i, j, istag, jstag, order
!    type(fv_grid_bounds_type), intent(IN) :: bd
!    real, dimension(bd%isd:, bd%jsd:, 1:) :: field
!    type(duogrid_type), intent(in) :: dg
!    character(len=*), intent(in) :: direction
!    logical, optional, intent(in) :: corner
!
!    !local
!    integer k, npz, dims(3)
!
!    dims = shape(field)
!    npz = dims(3)
!
!    do k = 1, npz
!      call lagrange_poly_interp_2d(field(:, :, k), i, j, bd, dg, direction, istag, jstag, order, corner)
!    end do
!
!  end subroutine lagrange_poly_interp_3d

  subroutine lagrange_poly_interp_4d(field, i, j, bd, dg, direction, istag, jstag, order, corner)

    integer, intent(in) :: i, j, istag, jstag, order
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, dimension(bd%isd:, bd%jsd:, 1:, 1:) :: field
    type(duogrid_type), intent(in) :: dg
    character(len=*), intent(in) :: direction
    logical, optional, intent(in) :: corner

    !local
    integer k, q, npz, nq, dims(4)

    dims = shape(field)
    npz = dims(3)
    nq = dims(4)

    do q = 1, nq
      call lagrange_poly_interp_3d(field(:, :, :, nq), i, j, bd, dg, direction, istag, jstag, order, corner)
    end do

  end subroutine lagrange_poly_interp_4d

  !===========================================================================
  !===========================================================================
  !===========================================================================

  subroutine cubed_a2c_halo(npz, uatemp, vatemp, uc, vc, dg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HERE WE USE THE BD OF ng=4 domain, so isd,ied,jsd,jed are offsetted by one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Purpose; Transform wind on A grid to C grid

    type(duogrid_type), intent(in), target :: dg
    integer, intent(in):: npz
    real, intent(inout), dimension(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed, npz):: uatemp, vatemp
    real, intent(inout):: uc(dg%bd%isd:dg%bd%ied + 1, dg%bd%jsd:dg%bd%jed, npz)
    real, intent(inout):: vc(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed + 1, npz)
! local:
    real v3(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)
    real ue(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)    ! 3D winds at edges
    real ve(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)    ! 3D winds at edges

    integer i, j, k

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ew, es
    real(kind=R_GRID), pointer, dimension(:, :, :)   :: vlon, vlat

    integer :: isd, ied, jsd, jed

    isd = dg%bd%isd
    ied = dg%bd%ied
    jsd = dg%bd%jsd
    jed = dg%bd%jed

    vlon => dg%vlon_ext
    vlat => dg%vlat_ext
    ew => dg%ew_ext
    es => dg%es_ext

! check if this ok
    call fill_corner_region(vatemp, dg%bd, dg, 0, 0)
    call fill_corner_region(uatemp, dg%bd, dg, 0, 0)

    do k = 1, npz
! Compute 3D wind on A grid
      do j = jsd, jed
        do i = isd, ied
          v3(1, i, j) = uatemp(i, j, k)*vlon(i, j, 1) + vatemp(i, j, k)*vlat(i, j, 1)
          v3(2, i, j) = uatemp(i, j, k)*vlon(i, j, 2) + vatemp(i, j, k)*vlat(i, j, 2)
          v3(3, i, j) = uatemp(i, j, k)*vlon(i, j, 3) + vatemp(i, j, k)*vlat(i, j, 3)
        end do
      end do

! A --> C
! Interpolate to cell edges
      do j = jsd, jed
        do i = isd + 1, ied
          ue(1, i, j) = 0.5*(v3(1, i - 1, j) + v3(1, i, j))
          ue(2, i, j) = 0.5*(v3(2, i - 1, j) + v3(2, i, j))
          ue(3, i, j) = 0.5*(v3(3, i - 1, j) + v3(3, i, j))
        end do
      end do

      do j = jsd + 1, jed
        do i = isd, ied
          ve(1, i, j) = 0.5*(v3(1, i, j - 1) + v3(1, i, j))
          ve(2, i, j) = 0.5*(v3(2, i, j - 1) + v3(2, i, j))
          ve(3, i, j) = 0.5*(v3(3, i, j - 1) + v3(3, i, j))
        end do
      end do

      do j = jsd, jed
        do i = isd + 1, ied
          uc(i, j, k) = ue(1, i, j)*ew(1, i, j, 1) + &
                        ue(2, i, j)*ew(2, i, j, 1) + &
                        ue(3, i, j)*ew(3, i, j, 1)
        end do
      end do
      do j = jsd + 1, jed
        do i = isd, ied
          vc(i, j, k) = ve(1, i, j)*es(1, i, j, 2) + &
                        ve(2, i, j)*es(2, i, j, 2) + &
                        ve(3, i, j)*es(3, i, j, 2)
        end do
      end do

    end do         ! k-loop

  end subroutine cubed_a2c_halo

  subroutine cubed_a2d_halo(npz, uatemp, vatemp, ud, vd, dg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HERE WE USE THE BD OF ng=4 domain, so isd,ied,jsd,jed are offsetted by one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Purpose; Transform wind on A grid to D grid

    type(duogrid_type), intent(IN), target :: dg
    integer, intent(in):: npz
    real, intent(inout), dimension(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed, npz):: uatemp, vatemp
    real, intent(inout):: ud(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed + 1, npz)
    real, intent(inout):: vd(dg%bd%isd:dg%bd%ied + 1, dg%bd%jsd:dg%bd%jed, npz)
! local:
    real v3(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)
    real ue(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)    ! 3D winds at edges
    real ve(3, dg%bd%is - dg%bd%ng:dg%bd%ie + 1 + dg%bd%ng, dg%bd%js - dg%bd%ng:dg%bd%je + 1 + dg%bd%ng)    ! 3D winds at edges

    integer i, j, k

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ew, es
    real(kind=R_GRID), pointer, dimension(:, :, :) :: vlon, vlat

    integer :: isd, ied, jsd, jed

    isd = dg%bd%isd
    ied = dg%bd%ied
    jsd = dg%bd%jsd
    jed = dg%bd%jed

    ud = -999.
    vd = -888.

    vlon => dg%vlon_ext
    vlat => dg%vlat_ext
    ew => dg%ew_ext
    es => dg%es_ext

! check if this ok
    call fill_corner_region(vatemp, dg%bd, dg, 0, 0)
    call fill_corner_region(uatemp, dg%bd, dg, 0, 0)

    do k = 1, npz
! Compute 3D wind on A grid
      do j = jsd, jed
        do i = isd, ied
          v3(1, i, j) = uatemp(i, j, k)*vlon(i, j, 1) + vatemp(i, j, k)*vlat(i, j, 1)
          v3(2, i, j) = uatemp(i, j, k)*vlon(i, j, 2) + vatemp(i, j, k)*vlat(i, j, 2)
          v3(3, i, j) = uatemp(i, j, k)*vlon(i, j, 3) + vatemp(i, j, k)*vlat(i, j, 3)
        end do
      end do

! A --> D
! Interpolate to cell edges
      do j = jsd + 1, jed
        do i = isd, ied
          ue(1, i, j) = 0.5*(v3(1, i, j - 1) + v3(1, i, j))
          ue(2, i, j) = 0.5*(v3(2, i, j - 1) + v3(2, i, j))
          ue(3, i, j) = 0.5*(v3(3, i, j - 1) + v3(3, i, j))
        end do
      end do

      do j = jsd, jed
        do i = isd + 1, ied
          ve(1, i, j) = 0.5*(v3(1, i - 1, j) + v3(1, i, j))
          ve(2, i, j) = 0.5*(v3(2, i - 1, j) + v3(2, i, j))
          ve(3, i, j) = 0.5*(v3(3, i - 1, j) + v3(3, i, j))
        end do
      end do

      do j = jsd + 1, jed
        do i = isd, ied
          ud(i, j, k) = ue(1, i, j)*es(1, i, j, 1) + &
                        ue(2, i, j)*es(2, i, j, 1) + &
                        ue(3, i, j)*es(3, i, j, 1)
        end do
      end do
      do j = jsd, jed
        do i = isd + 1, ied
          vd(i, j, k) = ve(1, i, j)*ew(1, i, j, 2) + &
                        ve(2, i, j)*ew(2, i, j, 2) + &
                        ve(3, i, j)*ew(3, i, j, 2)
        end do
      end do

    end do         ! k-loop

  end subroutine cubed_a2d_halo

  subroutine c2l_ord2_cgrid(u, v, ua, va, gridstruct, km, grid_type, bd, do_halo)
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: km, grid_type
    real, intent(in) ::  u(bd%isd:bd%ied + 1, bd%jsd:bd%jed, km)
    real, intent(in) ::  v(bd%isd:bd%ied, bd%jsd:bd%jed + 1, km)
    type(fv_grid_type), intent(IN), target :: gridstruct
    logical, intent(in) :: do_halo
!
    real, intent(out):: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    real, intent(out):: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    real wu(bd%is - 1:bd%ie + 2, bd%js - 1:bd%je + 1)
    real wv(bd%is - 1:bd%ie + 1, bd%js - 1:bd%je + 2)
    real u1(bd%is - 1:bd%ie + 1), v1(bd%is - 1:bd%ie + 1)
    integer i, j, k
    integer :: is, ie, js, je

    real, dimension(:, :), pointer :: a11, a12, a21, a22
    real, dimension(:, :), pointer :: dx, dy, rdxa, rdya

    a11 => gridstruct%a11
    a12 => gridstruct%a12
    a21 => gridstruct%a21
    a22 => gridstruct%a22

    dx => gridstruct%dx
    dy => gridstruct%dy
    rdxa => gridstruct%rdxa
    rdya => gridstruct%rdya

    if (do_halo) then
      is = bd%is - 1
      ie = bd%ie + 1
      js = bd%js - 1
      je = bd%je + 1
    else
      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
    end if

!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
    do k = 1, km
      if (grid_type < 4) then
        do j = js, je
          do i = is, ie + 1
            wu(i, j) = u(i, j, k)*dy(i, j)
          end do
        end do
        do j = js, je + 1
          do i = is, ie
            wv(i, j) = v(i, j, k)*dx(i, j)
          end do
        end do

        do j = js, je
          do i = is, ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = 2.*(wu(i, j) + wu(i + 1, j))/(dy(i, j) + dy(i + 1, j))
            v1(i) = 2.*(wv(i, j) + wv(i, j + 1))/(dx(i, j) + dx(i, j + 1))
! Cubed (cell center co-variant winds) to lat-lon:
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i) !org
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i) !org
          end do
        end do
      else
! 2nd order:
        do j = js, je
          do i = is, ie
            va(i, j, k) = 0.5*(v(i, j, k) + v(i, j + 1, k))
            ua(i, j, k) = 0.5*(u(i, j, k) + u(i + 1, j, k))
          end do
        end do
      end if
    end do

  end subroutine c2l_ord2_cgrid

  subroutine unit_vect_latlon_ext(pp, elon, elat)
    real(kind=R_GRID), intent(IN)  :: pp(2)
    real(kind=R_GRID), intent(OUT) :: elon(3), elat(3)

    real(selected_real_kind(20)):: lon, lat
    real(selected_real_kind(20)):: sin_lon, cos_lon, sin_lat, cos_lat

    lon = pp(1)
    lat = pp(2)

    sin_lon = sin(lon)
    cos_lon = cos(lon)
    sin_lat = sin(lat)
    cos_lat = cos(lat)

    elon(1) = -sin_lon
    elon(2) = cos_lon
    elon(3) = 0.d0

    elat(1) = -sin_lat*cos_lon
    elat(2) = -sin_lat*sin_lon
    elat(3) = cos_lat

  end subroutine unit_vect_latlon_ext

  subroutine a2stag_metrics(gridstruct, bd, dg)

    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(INOUT), target :: dg
    type(fv_grid_type), intent(IN), target :: gridstruct
! local:

    integer i, j, k

    real(kind=R_GRID), dimension(bd%isd:bd%ied, bd%jsd:bd%jed, 3)   :: vlon, vlat
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ew, es

    real(kind=R_GRID), pointer, dimension(:, :, :) :: aagrid_local
    real(kind=R_GRID), pointer, dimension(:, :, :) :: bbgrid_local

    real(kind=R_GRID) p1(3), p2(3), p3(3), pp(3)
    real(kind=R_GRID) grid3(3, bd%isd:bd%ied + 1, bd%jsd:bd%jed + 1)

    integer :: isd, ied, jsd, jed
    integer :: l, ll

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    allocate (ew(3, isd:ied + 1, jsd:jed, 2))
    allocate (es(3, isd:ied, jsd:jed + 1, 2))

    allocate (aagrid_local(isd:ied, jsd:jed, 2))
    allocate (bbgrid_local(isd:ied + 1, jsd:jed + 1, 2))

    do j = jsd, jed
      do i = isd, ied
        do l = 1, 2
          aagrid_local(i, j, l) = dg%a_pt(l, i, j)
        end do
      end do
    end do

    do j = jsd, jed + 1
      do i = isd, ied + 1
        do l = 1, 2
          bbgrid_local(i, j, l) = dg%b_pt(l, i, j)
        end do
      end do
    end do

    do j = jsd, jed + 1
      do i = isd, ied + 1
        call latlon2xyz(bbgrid_local(i, j, 1:2), grid3(1, i, j))
      end do
    end do

    do j = jsd, jed
      do i = isd + 1, ied
        call mid_pt_cart(bbgrid_local(i, j, 1:2), bbgrid_local(i, j + 1, 1:2), pp)
        call latlon2xyz(aagrid_local(i - 1, j, 1:2), p3)
        call latlon2xyz(aagrid_local(i, j, 1:2), p1)
        call vect_cross(p2, p3, p1)
        call vect_cross(ew(1:3, i, j, 1), p2, pp)
        call normalize_vect(ew(1:3, i, j, 1))
!---
        call vect_cross(p1, grid3(1, i, j), grid3(1, i, j + 1))
        call vect_cross(ew(1:3, i, j, 2), p1, pp)
        call normalize_vect(ew(1:3, i, j, 2))
      end do
    end do

    do j = jsd + 1, jed
      do i = isd, ied
        call mid_pt_cart(bbgrid_local(i, j, 1:2), bbgrid_local(i + 1, j, 1:2), pp)
        call latlon2xyz(aagrid_local(i, j, 1:2), p1)
        call latlon2xyz(aagrid_local(i, j - 1, 1:2), p3)
        call vect_cross(p2, p3, p1)
        call vect_cross(es(1:3, i, j, 2), p2, pp)
        call normalize_vect(es(1:3, i, j, 2))
!---
        call vect_cross(p3, grid3(1, i, j), grid3(1, i + 1, j))
        call vect_cross(es(1:3, i, j, 1), p3, pp)
        call normalize_vect(es(1:3, i, j, 1))
      end do
    end do

! Compute 3D wind on A grid
    do j = jsd, jed
      do i = isd, ied
        call unit_vect_latlon_ext(aagrid_local(i, j, 1:2), vlon(i, j, 1:3), vlat(i, j, 1:3))
      end do
    end do

! Populate in dg
    do l = 1, 3
      do j = jsd, jed
        do i = isd + 1, ied
          dg%ew_ext(l, i, j, 1) = ew(l, i, j, 1)
          dg%ew_ext(l, i, j, 2) = ew(l, i, j, 2)
        end do
      end do

      do j = jsd + 1, jed
        do i = isd, ied
          dg%es_ext(l, i, j, 1) = es(l, i, j, 1)
          dg%es_ext(l, i, j, 2) = es(l, i, j, 2)
        end do
      end do

      do j = jsd, jed
        do i = isd, ied
          dg%vlon_ext(i, j, l) = vlon(i, j, l)
          dg%vlat_ext(i, j, l) = vlat(i, j, l)
        end do
      end do
    end do

    deallocate (ew)
    deallocate (es)

    deallocate (aagrid_local)
    deallocate (bbgrid_local)

  end subroutine a2stag_metrics

  !===========================================================================
  !===========================================================================
  !===========================================================================

  subroutine duogrid_alloc(dg, bd, flagstruct)
    type(duogrid_type), intent(INOUT), target :: dg
    type(fv_flags_type), intent(INOUT) :: flagstruct
    type(fv_grid_bounds_type), intent(IN) :: bd
    !--- local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: npx, npy

    !--- assign dims
    npx = flagstruct%npx
    npy = flagstruct%npy

    ng = dg%bd%ng
    is = dg%bd%is
    ie = dg%bd%ie
    js = dg%bd%js
    je = dg%bd%je
    isd = dg%bd%isd
    ied = dg%bd%ied
    jsd = dg%bd%jsd
    jed = dg%bd%jed

    !--- allocation

    ! lon lat
    allocate (dg%a_pt(2, isd:ied, jsd:jed))
    allocate (dg%b_pt(2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_pt(2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_pt(2, isd:ied + 1, jsd:jed + 1))

    ! scalars
    allocate (dg%a_x(isd:ied, jsd:jed))
    allocate (dg%a_y(isd:ied, jsd:jed))
    allocate (dg%a_kik_x(isd:ied, jsd:jed))
    allocate (dg%a_kik_y(isd:ied, jsd:jed))
    allocate (dg%a_dx(isd:ied, jsd:jed))
    allocate (dg%a_dy(isd:ied, jsd:jed))
    allocate (dg%a_da(isd:ied, jsd:jed))
    allocate (dg%rda(isd:ied, jsd:jed))
    allocate (dg%rdx(isd:ied, jsd:jed))
    allocate (dg%rdy(isd:ied, jsd:jed))

    allocate (dg%a_sina(isd:ied, jsd:jed))
    allocate (dg%a_cosa(isd:ied, jsd:jed))

    ! matrices
    allocate (dg%a_gco(2, 2, isd:ied, jsd:jed))
    allocate (dg%a_gct(2, 2, isd:ied, jsd:jed))

    ! allocate( dg%c_ct2ort_x(2,2,isd:ied+1,jsd:jed  ) )
    ! allocate( dg%c_ort2ct_x(2,2,isd:ied+1,jsd:jed  ) )
    ! allocate( dg%d_ct2ort_y(2,2,isd:ied,  jsd:jed+1) )
    ! allocate( dg%d_ort2ct_y(2,2,isd:ied,  jsd:jed+1) )

    allocate (dg%a_c2l(2, 2, isd:ied, jsd:jed))
    allocate (dg%a_l2c(2, 2, isd:ied, jsd:jed))

    allocate (dg%a_kik_gco(2, 2, isd:ied, jsd:jed))
    allocate (dg%a_kik_gct(2, 2, isd:ied, jsd:jed))
    allocate (dg%a_kik_c2l(2, 2, isd:ied, jsd:jed))
    allocate (dg%a_kik_l2c(2, 2, isd:ied, jsd:jed))

    ! B-grid
    allocate (dg%b_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_kik_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_kik_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_gco(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_gct(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_c2l(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_l2c(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_sina(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_cosa(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_dx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_dy(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_da(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_rda(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_rdx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%b_rdy(isd:ied + 1, jsd:jed + 1))

    ! C-grid
    allocate (dg%c_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_kik_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_kik_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_gco(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_gct(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_c2l(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_l2c(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_sina(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_cosa(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_dx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_dy(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_da(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_rda(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_rdx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%c_rdy(isd:ied + 1, jsd:jed + 1))

    ! D-grid
    allocate (dg%d_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_kik_x(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_kik_y(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_gco(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_gct(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_c2l(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_l2c(2, 2, isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_sina(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_cosa(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_dx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_dy(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_da(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_rda(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_rdx(isd:ied + 1, jsd:jed + 1))
    allocate (dg%d_rdy(isd:ied + 1, jsd:jed + 1))

    ! k2e
    allocate (dg%k2e_loc(isd:ied, jsd:jed))
    allocate (dg%k2e_coef(dg%k2e_nord, isd:ied, jsd:jed))
    allocate (dg%k2e_loc_b(isd:ied + 1, jsd:jed + 1))
    allocate (dg%k2e_coef_b(dg%k2e_nord, isd:ied + 1, jsd:jed + 1))
    allocate (dg%k2e_loc_c_x(isd:ied + 1, jsd:jed))
    allocate (dg%k2e_coef_c_x(dg%k2e_nord, isd:ied + 1, jsd:jed))
    allocate (dg%k2e_loc_c_y(isd:ied, jsd:jed + 1))
    allocate (dg%k2e_coef_c_y(dg%k2e_nord, isd:ied, jsd:jed + 1))
    allocate (dg%k2e_loc_d_x(isd:ied + 1, jsd:jed))
    allocate (dg%k2e_coef_d_x(dg%k2e_nord, isd:ied + 1, jsd:jed))
    allocate (dg%k2e_loc_d_y(isd:ied, jsd:jed + 1))
    allocate (dg%k2e_coef_d_y(dg%k2e_nord, isd:ied, jsd:jed + 1))

    ! Coriolis
    allocate (dg%a_f(isd:ied, jsd:jed))
    allocate (dg%a_kik_f(isd:ied, jsd:jed))

    ! extended grid metrics for A 2 C/D transform
    ! these are used with domainforduo ng=4, so offset by one
    allocate (dg%vlon_ext(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed, 3))
    allocate (dg%vlat_ext(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed, 3))
    allocate (dg%ew_ext(3, dg%bd%isd:dg%bd%ied + 1, dg%bd%jsd:dg%bd%jed, 2))
    allocate (dg%es_ext(3, dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%jed + 1, 2))

    ! corner lag coeff
    ! these are used with domain ng=3
    allocate (dg%xp(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4))
    allocate (dg%xm(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4))
    allocate (dg%yp(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4))
    allocate (dg%ym(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4))

  end subroutine duogrid_alloc
  !===========================================================================
  subroutine duogrid_dealloc(dg)
    type(duogrid_type), intent(inout)  :: dg

    !--- deallocation
    deallocate (dg%a_f)
    deallocate (dg%a_kik_f)

    deallocate (dg%k2e_loc)
    deallocate (dg%k2e_coef)
    deallocate (dg%k2e_loc_b)
    deallocate (dg%k2e_coef_b)
    deallocate (dg%k2e_loc_c_x)
    deallocate (dg%k2e_coef_c_x)
    deallocate (dg%k2e_loc_c_y)
    deallocate (dg%k2e_coef_c_y)
    deallocate (dg%k2e_loc_d_x)
    deallocate (dg%k2e_coef_d_x)
    deallocate (dg%k2e_loc_d_y)
    deallocate (dg%k2e_coef_d_y)

    deallocate (dg%a_gco)
    deallocate (dg%a_gct)

    !  deallocate( dg%c_ct2ort_x )
    !  deallocate( dg%c_ort2ct_x )
    !  deallocate( dg%d_ct2ort_y )
    !  deallocate( dg%d_ort2ct_y )

    deallocate (dg%a_c2l)
    deallocate (dg%a_l2c)

    deallocate (dg%a_kik_gco)
    deallocate (dg%a_kik_gct)
    deallocate (dg%a_kik_c2l)
    deallocate (dg%a_kik_l2c)

    deallocate (dg%a_x)
    deallocate (dg%a_y)
    deallocate (dg%a_kik_x)
    deallocate (dg%a_kik_y)
    deallocate (dg%a_dx)
    deallocate (dg%a_dy)
    deallocate (dg%a_da)
    deallocate (dg%rda)
    deallocate (dg%rdx)
    deallocate (dg%rdy)

    deallocate (dg%a_sina)
    deallocate (dg%a_cosa)

    ! lon lat
    deallocate (dg%a_pt)
    deallocate (dg%b_pt)

    ! B-grid
    deallocate (dg%b_x)
    deallocate (dg%b_y)
    deallocate (dg%b_kik_x)
    deallocate (dg%b_kik_y)
    deallocate (dg%b_gco)
    deallocate (dg%b_gct)
    deallocate (dg%b_c2l)
    deallocate (dg%b_l2c)
    deallocate (dg%b_sina)
    deallocate (dg%b_cosa)
    deallocate (dg%b_dx)
    deallocate (dg%b_dy)
    deallocate (dg%b_da)
    deallocate (dg%b_rda)
    deallocate (dg%b_rdx)
    deallocate (dg%b_rdy)

    ! C-grid
    deallocate (dg%c_x)
    deallocate (dg%c_y)
    deallocate (dg%c_kik_x)
    deallocate (dg%c_kik_y)
    deallocate (dg%c_gco)
    deallocate (dg%c_gct)
    deallocate (dg%c_c2l)
    deallocate (dg%c_l2c)
    deallocate (dg%c_sina)
    deallocate (dg%c_cosa)
    deallocate (dg%c_dx)
    deallocate (dg%c_dy)
    deallocate (dg%c_da)
    deallocate (dg%c_rda)
    deallocate (dg%c_rdx)
    deallocate (dg%c_rdy)

    ! D-grid
    deallocate (dg%d_x)
    deallocate (dg%d_y)
    deallocate (dg%d_kik_x)
    deallocate (dg%d_kik_y)
    deallocate (dg%d_gco)
    deallocate (dg%d_gct)
    deallocate (dg%d_c2l)
    deallocate (dg%d_l2c)
    deallocate (dg%d_sina)
    deallocate (dg%d_cosa)
    deallocate (dg%d_dx)
    deallocate (dg%d_dy)
    deallocate (dg%d_da)
    deallocate (dg%d_rda)
    deallocate (dg%d_rdx)
    deallocate (dg%d_rdy)

    ! extended grid metrics for A 2 C/D transform
    deallocate (dg%vlon_ext)
    deallocate (dg%vlat_ext)
    deallocate (dg%ew_ext)
    deallocate (dg%es_ext)
  end subroutine duogrid_dealloc
  !===========================================================================


 subroutine c2l_ord2(u, v, ua, va, gridstruct, km, grid_type, bd, do_halo)
 type(fv_grid_bounds_type), intent(IN) :: bd
  integer, intent(in) :: km, grid_type
  real, intent(in) ::  u(bd%isd:bd%ied,bd%jsd:bd%jed+1,km)
  real, intent(in) ::  v(bd%isd:bd%ied+1,bd%jsd:bd%jed,km)
 type(fv_grid_type), intent(IN), target :: gridstruct
 logical, intent(in) :: do_halo
!
  real, intent(out):: ua(bd%isd:bd%ied, bd%jsd:bd%jed,km)
  real, intent(out):: va(bd%isd:bd%ied, bd%jsd:bd%jed,km)
!--------------------------------------------------------------
! Local
  real wu(bd%is-1:bd%ie+1,  bd%js-1:bd%je+2)
  real wv(bd%is-1:bd%ie+2,  bd%js-1:bd%je+1)
  real u1(bd%is-1:bd%ie+1), v1(bd%is-1:bd%ie+1)
  integer i, j, k
  integer :: is,  ie,  js,  je

  real, dimension(:,:), pointer :: a11, a12, a21, a22
  real, dimension(:,:), pointer :: dx, dy, rdxa, rdya

  a11 => gridstruct%a11
  a12 => gridstruct%a12
  a21 => gridstruct%a21
  a22 => gridstruct%a22

  dx   => gridstruct%dx
  dy   => gridstruct%dy
  rdxa => gridstruct%rdxa
  rdya => gridstruct%rdya

  if (do_halo) then
     is  = bd%is-1
     ie  = bd%ie+1
     js  = bd%js-1
     je  = bd%je+1
  else
     is  = bd%is
     ie  = bd%ie
     js  = bd%js
     je  = bd%je
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,km,grid_type,u,dx,v,dy,ua,va,a11,a12,a21,a22) &
!$OMP                          private(u1, v1, wu, wv)
  do k=1,km
     if ( grid_type < 4 ) then
       do j=js,je+1
          do i=is,ie
             wu(i,j) = u(i,j,k)*dx(i,j)
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             wv(i,j) = v(i,j,k)*dy(i,j)
          enddo
       enddo

       do j=js,je
          do i=is,ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
             u1(i) = 2.*(wu(i,j) + wu(i,j+1)) / (dx(i,j)+dx(i,j+1))
             v1(i) = 2.*(wv(i,j) + wv(i+1,j)) / (dy(i,j)+dy(i+1,j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
             ua(i,j,k) = a11(i,j)*u1(i) + a12(i,j)*v1(i)
             va(i,j,k) = a21(i,j)*u1(i) + a22(i,j)*v1(i)
          enddo
       enddo
     else
! 2nd order:
       do j=js,je
          do i=is,ie
             ua(i,j,k) = 0.5*(u(i,j,k)+u(i,  j+1,k))
             va(i,j,k) = 0.5*(v(i,j,k)+v(i+1,j,  k))
          enddo
       enddo
     endif
  enddo

 end subroutine c2l_ord2


 real function great_circle_dist( q1, q2, radius )
      real(kind=R_GRID), intent(IN)           :: q1(2), q2(2)
      real(kind=R_GRID), intent(IN), optional :: radius

      real (kind=R_GRID):: p1(2), p2(2)
      real (kind=R_GRID):: beta
      integer n

      do n=1,2
         p1(n) = q1(n)
         p2(n) = q2(n)
      enddo

      beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                         sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

      if ( present(radius) ) then
           great_circle_dist = radius * beta
      else
           great_circle_dist = beta   ! Returns the angle
      endif

  end function great_circle_dist

 subroutine mid_pt_cart(p1, p2, e3)
    real(kind=R_GRID), intent(IN)  :: p1(2), p2(2)
    real(kind=R_GRID), intent(OUT) :: e3(3)
!-------------------------------------
    real(kind=R_GRID) e1(3), e2(3)

    call latlon2xyz(p1, e1)
    call latlon2xyz(p2, e2)
    call mid_pt3_cart(e1, e2, e3)

 end subroutine mid_pt_cart

 subroutine mid_pt3_cart(p1, p2, e)
       real(kind=R_GRID), intent(IN)  :: p1(3), p2(3)
       real(kind=R_GRID), intent(OUT) :: e(3)
!
       real (kind=R_GRID):: q1(3), q2(3)
       real (kind=R_GRID):: dd, e1, e2, e3
       integer k

       do k=1,3
          q1(k) = p1(k)
          q2(k) = p2(k)
       enddo

       e1 = q1(1) + q2(1)
       e2 = q1(2) + q2(2)
       e3 = q1(3) + q2(3)

       dd = sqrt( e1**2 + e2**2 + e3**2 )
       e1 = e1 / dd
       e2 = e2 / dd
       e3 = e3 / dd

       e(1) = e1
       e(2) = e2
       e(3) = e3

 end subroutine mid_pt3_cart



 subroutine latlon2xyz(p, e, id)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real(kind=R_GRID), intent(in) :: p(2)
 real(kind=R_GRID), intent(out):: e(3)
 integer, optional, intent(in):: id   ! id=0 do nothing; id=1, right_hand

 integer n
 real (kind=R_GRID):: q(2)
 real (kind=R_GRID):: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz

 subroutine vect_cross(e, p1, p2)
 real(kind=R_GRID), intent(in) :: p1(3), p2(3)
 real(kind=R_GRID), intent(out):: e(3)
!
! Perform cross products of 3D vectors: e = P1 X P2
!
      e(1) = p1(2)*p2(3) - p1(3)*p2(2)
      e(2) = p1(3)*p2(1) - p1(1)*p2(3)
      e(3) = p1(1)*p2(2) - p1(2)*p2(1)

 end subroutine vect_cross


 subroutine normalize_vect(e)
!                              Make e an unit vector
 real(kind=R_GRID), intent(inout):: e(3)
 real(kind=R_GRID):: pdot
 integer k

    pdot = e(1)**2 + e(2)**2 + e(3)**2
    pdot = sqrt( pdot )

    do k=1,3
       e(k) = e(k) / pdot
    enddo

 end subroutine normalize_vect


end module duogrid_mod

