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

  !use lib_grid_mod, only: R_GRID, RADIUS
  !use lib_grid_mod, only: lib_4pt_area

  use fv_arrays_mod, only: R_GRID

  use global_grid_mod, only: global_grid_type
  use global_grid_mod, only: global_grid_init, global_grid_end

  use fv_arrays_mod, only: duogrid_type, fv_grid_bounds_type, fv_flags_type, fv_atmos_type, fv_grid_type
  !use fv_grid_utils_mod, only: c2l_ord2, great_circle_dist
  use fv_timing_mod, only: timing_on, timing_off

  use mpp_parameter_mod, only: DGRID_NE, BGRID_NE, BGRID_SW, AGRID, SCALAR_PAIR, NORTH, EAST, WEST, SOUTH, CGRID_NE

  ! use fv_grid_utils_mod, only: mid_pt_cart, latlon2xyz, vect_cross, normalize_vect
  implicit none

  private

#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  public :: duogrid_init
  public :: duogrid_end
  public :: ext_scalar
  public :: ext_vector
!  public :: fill_corner_region, lagrange_poly_interp
!  public :: cubed_a2d_halo

!  interface fill_corner_region
!    module procedure fill_corner_region_2d
!    module procedure fill_corner_region_3d
!    module procedure fill_corner_region_4d
!  end interface fill_corner_region
!  interface lagrange_poly_interp
!    module procedure lagrange_poly_interp_2d
!    module procedure lagrange_poly_interp_3d
!    module procedure lagrange_poly_interp_4d
!  end interface lagrange_poly_interp
!
!  interface fill_corners_domain_decomp
!    module procedure fill_corners_domain_decomp_2d
!    module procedure fill_corners_domain_decomp_3d
!  end interface fill_corners_domain_decomp
!
  interface fill_buff_corners_domain_decomp
    module procedure fill_buff_corners_domain_decomp_2d
    module procedure fill_buff_corners_domain_decomp_3d
  end interface fill_buff_corners_domain_decomp

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
! extrapolation for divergence
    call compute_lagrange_coeff_extra(atm%gridstruct%dg%xp3, atm%gridstruct%dg%xm3, &
                                      atm%gridstruct%dg%yp3, atm%gridstruct%dg%ym3, atm%bd, atm%gridstruct%dg, gg)
    call compute_lagrange_coeff_extra(atm%gridstruct%dg%xp, atm%gridstruct%dg%xm, &
                                      atm%gridstruct%dg%yp, atm%gridstruct%dg%ym, atm%gridstruct%dg%bd, atm%gridstruct%dg, gg)

    !--- end global_grid
    !call global_grid_end(gg)

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
        dg%a_kik(:, i, j) = gg%pt_kik(:, ii, jj, n)
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
!#ifdef a_grid_recalculate
!    do j = jsd, jed
!      do i = isd, ied
!        dg%a_da(i, j) = lib_4pt_area( &
!                        dg%b_pt(:, i, j), dg%b_pt(:, i + 1, j), dg%b_pt(:, i + 1, j + 1), dg%b_pt(:, i, j + 1), &
!                        RADIUS)
!        dg%rda(i, j) = 1.0/dg%a_da(i, j)
!      end do
!    end do
!#endif

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
    integer :: isd, ied, jsd, jed, i, j
    integer, dimension(bd%isd:bd%ied, bd%jsd:bd%jed) :: k2e_loc_u, k2e_loc_v
    real(kind=R_GRID), dimension(dg%k2e_nord, bd%isd:bd%ied, bd%jsd:bd%jed) :: k2e_coef_u, k2e_coef_v
    real, dimension(bd%isd:bd%ied, bd%je:bd%jed) :: S_N
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%js) :: S_S
    real, dimension(bd%ie:bd%ied, bd%jsd:bd%jed) :: S_E
    real, dimension(bd%isd:bd%is, bd%jsd:bd%jed) :: S_W

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    if (istag == 0 .and. jstag == 0) then
      call mpp_update_domains(var, domain, complete=.true.)
    else
      call mpp_error(FATAL, 'ext_scalar_2d for istag and jstag /=0 not implemented')
    end if

    do i = isd, ied
    do j = jsd, jed
      k2e_loc_u(i, j) = dg%k2e_loc(i, j)
      k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
    end do
    end do

    call create_buff_2d(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)
    call cube_rmp_buff(S_N(:, :), S_E(:, :), S_S(:, :), S_W(:, :), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
    call fill_corner_buffer(S_N(:, :), S_E(:, :), S_S(:, :), S_W(:, :), bd, dg, istag, jstag)
    call fill_buff_corners_domain_decomp(S_N, S_E, S_S, S_W, dg, bd, domain, istag, jstag)
    call apply_buff_2d(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)

  end subroutine ext_scalar_2d
  !===========================================================================
  !> @brief update scalar block halo and remap to ext grid when needed
  subroutine ext_scalar_3d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(INOUT) :: dg
    type(domain2d), intent(INOUT)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, :) :: var
    real, dimension(:, :, :), allocatable :: S_N, S_S, S_E, S_W ! scalar buffers
    real, dimension(:, :, :), allocatable :: varp1
    !--- local
    integer :: i, j, isd, ied, jsd, jed, k, npz
    integer :: dims(3)

    !real, dimension(bd%isd:bd%ied + istag, bd%je + jstag:bd%jed + jstag, dg%npz) :: S_N
    !real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%js, dg%npz) :: S_S
    !real, dimension(bd%ie + istag:bd%ied + istag, bd%jsd:bd%jed + jstag, dg%npz) :: S_E
    !real, dimension(bd%isd:bd%is, bd%jsd:bd%jed + jstag, dg%npz) :: S_W

    integer, dimension(:, :), allocatable :: k2e_loc_u
    real(kind=R_GRID), dimension(:, :, :), allocatable :: k2e_coef_u

    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    if ((istag == 0 .and. jstag == 1) .or. (istag == 1 .and. jstag == 0)) then
      call mpp_error(FATAL, 'ext_scalar_3d for stag fields not implemented')
    end if

      call timing_on('ext_scalar3d')
    !--- update domains
    if (istag == 1 .and. jstag == 1) then
      allocate (k2e_coef_u(dg%k2e_nord, isd:ied + 1, jsd:jed + 1))
      allocate (k2e_loc_u(isd:ied + 1, jsd:jed + 1))
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
      call mpp_update_domains(var, domain, complete=.true.)
    end if

    !--- remap to ext
    dims = shape(var)
    npz = dims(3) !!!! SOME FIELDS ARE NPZ+1

    allocate (S_N(bd%isd:bd%ied + istag, bd%je + jstag:bd%jed + jstag, npz))
    allocate (S_S(bd%isd:bd%ied + istag, bd%jsd:bd%js, npz))
    allocate (S_E(bd%ie + istag:bd%ied + istag, bd%jsd:bd%jed + jstag, npz))
    allocate (S_W(bd%isd:bd%is, bd%jsd:bd%jed + jstag, npz))

    if (istag == 0 .and. jstag == 0) then
      call create_buff(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)

      call timing_on('cube_buffrmp&cornerbuff')
      !$OMP parallel do default(none) shared(npz,S_N,S_E,S_S,S_W,dg,bd,domain,istag,jstag,k2e_loc_u,k2e_coef_u)
      do k = 1, npz
     call cube_rmp_buff(S_N(:, :, k), S_E(:, :, k), S_S(:, :, k), S_W(:, :, k), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
        call fill_corner_buffer(S_N(:, :, k), S_E(:, :, k), S_S(:, :, k), S_W(:, :, k), bd, dg, istag, jstag)
      end do
      call timing_off('cube_buffrmp&cornerbuff')
      call timing_on('fillcornerdomaindecomp')
      call fill_buff_corners_domain_decomp(S_N, S_E, S_S, S_W, dg, bd, domain, istag, jstag)
      call timing_off('fillcornerdomaindecomp')
      call apply_buff(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)
    end if

    if (istag == 1 .and. jstag == 1) then
      call create_buff(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)
      ! FLIP order TO EXTRAPOLATE or INTERPOLATE DIVERGENCE AT CORNERS (also check the compute lagrange coeff xp3)
      ! compute_lagrange_coeff_extra: with cube_rmp then fill_corner => extrapolate div from remaped data
      ! compute_lagrange_coeff: with fill_corner then cube_rmp=> interpolate div from kinked data
      call timing_on('cube_buffrmp&cornerbuff')
      !$OMP parallel do default(none) shared(npz,S_N,S_E,S_S,S_W,dg,bd,domain,istag,jstag,k2e_loc_u,k2e_coef_u)
      do k = 1, npz
     call cube_rmp_buff(S_N(:, :, k), S_E(:, :, k), S_S(:, :, k), S_W(:, :, k), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
        call fill_corner_buffer(S_N(:, :, k), S_E(:, :, k), S_S(:, :, k), S_W(:, :, k), bd, dg, istag, jstag)
      end do
      call timing_off('cube_buffrmp&cornerbuff')
      call timing_on('fillcornerdomaindecomp')
      call fill_buff_corners_domain_decomp(S_N, S_E, S_S, S_W, dg, bd, domain, istag, jstag)
      call timing_off('fillcornerdomaindecomp')
      call apply_buff(var, dg, bd, S_N, S_E, S_S, S_W, istag, jstag)
    end if

    deallocate (k2e_coef_u)
    deallocate (k2e_loc_u)
    deallocate (S_N, S_E, S_S, S_W)
      call timing_off('ext_scalar3d')
  end subroutine ext_scalar_3d
  !===========================================================================
  !> @brief update scalar block halo and remap to ext grid when needed
  subroutine ext_scalar_4d(var, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(INOUT) :: dg
    type(domain2d), intent(INOUT)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: istag, jstag
    real, dimension(:, :, :, :) :: var
    real, dimension(:, :, :, :), allocatable :: S_N, S_S, S_E, S_W ! scalar buffers
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

      call timing_on('ext_scalar4d')
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

    allocate (S_N(bd%isd:bd%ied + istag, bd%je + jstag:bd%jed + jstag, npz, nq))
    allocate (S_S(bd%isd:bd%ied + istag, bd%jsd:bd%js, npz, nq))
    allocate (S_E(bd%ie + istag:bd%ied + istag, bd%jsd:bd%jed + jstag, npz, nq))
    allocate (S_W(bd%isd:bd%is, bd%jsd:bd%jed + jstag, npz, nq))

    do n = 1, nq
    call create_buff(var(:, :, :, nq), dg, bd, S_N(:, :, :, nq), S_E(:, :, :, nq), S_S(:, :, :, nq), S_W(:, :, :, nq), istag, jstag)
    end do !nq

    do n = 1, nq
!$OMP parallel do default(none) shared(npz,nq,S_N,S_E,S_S,S_W,dg,bd,domain,istag,jstag,k2e_loc_u,k2e_coef_u)
      do k = 1, npz
       call cube_rmp_buff(S_N(:, :, k,nq), S_E(:, :, k,nq), S_S(:, :, k,nq), S_W(:, :, k,nq), dg, bd, domain, istag, jstag, k2e_loc_u, k2e_coef_u)
        call fill_corner_buffer(S_N(:, :, k, nq), S_E(:, :, k, nq), S_S(:, :, k, nq), S_W(:, :, k, nq), bd, dg, istag, jstag)
      end do
      call fill_buff_corners_domain_decomp(S_N(:,:,:,nq), S_E(:,:,:,nq), S_S(:,:,:,nq), S_W(:,:,:,nq), dg, bd, domain, istag, jstag)
     call apply_buff(var(:, :, :, nq), dg, bd, S_N(:, :, :, nq), S_E(:, :, :, nq), S_S(:, :, :, nq), S_W(:, :, :, nq), istag, jstag)
    end do !nq

    if (istag == 0 .and. jstag == 0) then
      deallocate (k2e_coef_u)
      deallocate (k2e_loc_u)
      deallocate (S_E, S_N, S_W, S_S)
    end if
      call timing_off('ext_scalar4d')
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
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng, npz
    integer :: i, j, k, k2e_nord

    real, dimension(bd%isd:bd%ied, bd%jsd:bd%jed, flagstruct%npz) :: ull, vll
    real, dimension(bd%isd - 1:bd%ied + 1, bd%jsd - 1:bd%jed + 1, flagstruct%npz) :: ullp1, vllp1
! for Agrid 3 haloes
    real, dimension(bd%isd:bd%ied, bd%je:bd%jed, flagstruct%npz) :: u_N
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%js, flagstruct%npz) :: u_S
    real, dimension(bd%ie:bd%ied, bd%jsd:bd%jed, flagstruct%npz) :: u_E
    real, dimension(bd%isd:bd%is, bd%jsd:bd%jed, flagstruct%npz) :: u_W
    real, dimension(bd%isd:bd%ied, bd%je:bd%jed, flagstruct%npz) :: v_N
    real, dimension(bd%isd:bd%ied, bd%jsd:bd%js, flagstruct%npz) :: v_S
    real, dimension(bd%ie:bd%ied, bd%jsd:bd%jed, flagstruct%npz) :: v_E
    real, dimension(bd%isd:bd%is, bd%jsd:bd%jed, flagstruct%npz) :: v_W

! for CDgrid 4 haloes
    real, dimension(dg%bd%isd:dg%bd%ied, dg%bd%je:dg%bd%jed, flagstruct%npz) :: u_N1
    real, dimension(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%js, flagstruct%npz) :: u_S1
    real, dimension(dg%bd%ie:dg%bd%ied, dg%bd%jsd:dg%bd%jed, flagstruct%npz) :: u_E1
    real, dimension(dg%bd%isd:dg%bd%is, dg%bd%jsd:dg%bd%jed, flagstruct%npz) :: u_W1
    real, dimension(dg%bd%isd:dg%bd%ied, dg%bd%je:dg%bd%jed, flagstruct%npz) :: v_N1
    real, dimension(dg%bd%isd:dg%bd%ied, dg%bd%jsd:dg%bd%js, flagstruct%npz) :: v_S1
    real, dimension(dg%bd%ie:dg%bd%ied, dg%bd%jsd:dg%bd%jed, flagstruct%npz) :: v_E1
    real, dimension(dg%bd%isd:dg%bd%is, dg%bd%jsd:dg%bd%jed, flagstruct%npz) :: v_W1

    real(kind=R_GRID), dimension(2, 2, bd%isd:bd%ied, bd%jsd:bd%jed) :: c2l, l2c

    integer, dimension(bd%isd - 1:bd%ied + 1, bd%jsd - 1:bd%jed + 1) :: k2e_loc_u, k2e_loc_v
    real(kind=R_GRID), dimension(dg%k2e_nord, bd%isd - 1:bd%ied + 1, bd%jsd - 1:bd%jed + 1) :: k2e_coef_u, k2e_coef_v

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

    npz = flagstruct%npz

    k2e_nord = dg%k2e_nord

    ull = -99999.
    vll = -99999.
    ullp1 = -99999.
    vllp1 = -99999.

    c2l = -99999.
    l2c = -99999.

    u_N = -99899.
    v_N = -99899.
    u_S = -99899.
    v_S = -99899.
    u_E = -99899.
    v_E = -99899.
    u_W = -99899.
    v_W = -99899.

!necessary
    u_N1 = -99899.
    v_N1 = -99899.
    u_S1 = -99899.
    v_S1 = -99899.
    u_E1 = -99899.
    v_E1 = -99899.
    u_W1 = -99899.
    v_W1 = -99899.

    ! Since we are transforming C/D variables to A using c2l, we only use the A-grid
    ! remapping coeff. C/D extensions coeff are available, using them however with the current
    ! logic produced more noised compared to the first method.
    ! To CHECK:  do C/D on the fly without passing by A.

    call timing_on('extvector')

    do j = jsd - 1, jed + 1
    do i = isd - 1, ied + 1
      k2e_loc_u(i, j) = dg%k2e_loc(i, j)
      k2e_coef_u(:, i, j) = dg%k2e_coef(:, i, j)
      k2e_loc_v(i, j) = dg%k2e_loc(i, j)
      k2e_coef_v(:, i, j) = dg%k2e_coef(:, i, j)
    end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform cubed to latlon winds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0) then  !AGRID
      do j = jsd, jed
      do i = isd, ied
        c2l(:, :, i, j) = dg%a_c2l(:, :, i, j)
        l2c(:, :, i, j) = dg%a_l2c(:, :, i, j)
      end do
      end do
      call mpp_update_domains(u_in, v_in, domain, gridtype=AGRID) ! update interior pe first
      call c2l_agrid(u_in, v_in, ull, vll, npz, bd, c2l) !zonal-merdional vel
    end if

    if (ieu_stag == 0 .and. jeu_stag == 1 .and. iev_stag == 1 .and. jev_stag == 0) then  !DGRID
      call mpp_update_domains(u_in, v_in, domain, gridtype=DGRID_NE) ! update interior pe first
      call c2l_ord2(u_in, v_in, ull, vll, gridstruct, npz, 0, bd, .true.) !zonal-merdional vel
    end if

    if (ieu_stag == 1 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 1) then  !CGRID
      call mpp_update_domains(u_in, v_in, domain, gridtype=CGRID_NE) !update interior pe first
      call c2l_ord2_cgrid(u_in, v_in, ull, vll, gridstruct, npz, 0, bd, .true.) !zonal-merdional vel
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Update latlonwinds haloes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0)) then  ! AGRID
      call mpp_update_domains(ull, domain, complete=.false.)
      call mpp_update_domains(vll, domain, complete=.true.)
      call create_buff(ull, dg, bd, u_N, u_E, u_S, u_W, 0, 0)
      call create_buff(vll, dg, bd, v_N, v_E, v_S, v_W, 0, 0)
    else !CDGRID
      !$OMP parallel do default(none) shared(npz,is,ie,js,je,ull,vll,ullp1,vllp1)
      do k = 1, npz
        do j = js, je + 1
          do i = is, ie + 1
            ullp1(i, j, k) = ull(i, j, k)
            vllp1(i, j, k) = vll(i, j, k)
          end do
        end do
      end do
      call mpp_update_domains(ullp1, dg%domain_for_duo, complete=.false.)
      call mpp_update_domains(vllp1, dg%domain_for_duo, complete=.true.)
      call create_buff(ullp1, dg, dg%bd, u_N1, u_E1, u_S1, u_W1, 0, 0)
      call create_buff(vllp1, dg, dg%bd, v_N1, v_E1, v_S1, v_W1, 0, 0)
    end if

    !!!!!!!!!!!!!!!!!!
    !--- remap to ext
    !!!!!!!!!!!!!!!!!!
    if ((ieu_stag == 0 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 0)) then   !AGRID

      call timing_on('cube_buffrmp&cornerbuff')

      !$OMP parallel do default(none) shared(npz,dg,bd,domain, u_N,u_E,u_S,u_W,k2e_loc_u,k2e_coef_u, v_N,v_E,v_S,v_W,k2e_loc_v,k2e_coef_v, isd, ied, jsd, jed, k2e_nord)
      do k = 1, npz
        call cube_rmp_buff(u_N(:, :, k), u_E(:, :, k), u_S(:, :, k), u_W(:, :, k), dg, bd, domain, 0, 0, k2e_loc_u (isd:ied, jsd:jed), k2e_coef_u(1:k2e_nord, isd:ied, jsd:jed) )
        call cube_rmp_buff(v_N(:, :, k), v_E(:, :, k), v_S(:, :, k), v_W(:, :, k), dg, bd, domain, 0, 0, k2e_loc_v (isd:ied, jsd:jed), k2e_coef_v(1:k2e_nord, isd:ied, jsd:jed) )
        call fill_corner_buffer(u_N(:, :, k), u_E(:, :, k), u_S(:, :, k), u_W(:, :, k), bd, dg, 0, 0)
        call fill_corner_buffer(v_N(:, :, k), v_E(:, :, k), v_S(:, :, k), v_W(:, :, k), bd, dg, 0, 0)
      end do

      call timing_off('cube_buffrmp&cornerbuff')

      call timing_on('fillcornerdomaindecomp')
      call fill_buff_corners_domain_decomp(u_N, u_E, u_S, u_W, dg, bd, domain, 0, 0)
      call fill_buff_corners_domain_decomp(v_N, v_E, v_S, v_W, dg, bd, domain, 0, 0)
      call timing_off('fillcornerdomaindecomp')

      call timing_on('l2c_halo_buff')
      call l2c_halo_buff_only(npz, u_N, v_N, bd, bd%isd, bd%ied, bd%je, bd%jed, l2c)
      call l2c_halo_buff_only(npz, u_W, v_W, bd, bd%isd, bd%is, bd%jsd, bd%jed, l2c)
      call l2c_halo_buff_only(npz, u_S, v_S, bd, bd%isd, bd%ied, bd%jsd, bd%js, l2c)
      call l2c_halo_buff_only(npz, u_E, v_E, bd, bd%ie, bd%ied, bd%jsd, bd%jed, l2c)
      call timing_off('l2c_halo_buff')

      call apply_buff(u_in, dg, bd, u_N, u_E, u_S, u_W, 0, 0)
      call apply_buff(v_in, dg, bd, v_N, v_E, v_S, v_W, 0, 0)

    else !CDGRID

      call timing_on('cube_buffrmp&cornerbuff')
      !$OMP parallel do default(none) shared(npz,u_N1,u_E1,u_S1,u_W1,dg,k2e_loc_u,k2e_coef_u, v_N1,v_E1,v_S1,v_W1,k2e_loc_v,k2e_coef_v)
      do k = 1, npz
        call cube_rmp_buff(u_N1(:, :, k), u_E1(:, :, k), u_S1(:, :, k), u_W1(:, :, k), \
        dg, dg%bd, dg%domain_for_duo, 0, 0, k2e_loc_u, k2e_coef_u)
        call cube_rmp_buff(v_N1(:, :, k), v_E1(:, :, k), v_S1(:, :, k), v_W1(:, :, k), \
        dg, dg%bd, dg%domain_for_duo, 0, 0, k2e_loc_v, k2e_coef_v)
        call fill_corner_buffer(u_N1(:, :, k), u_E1(:, :, k), u_S1(:, :, k), u_W1(:, :, k), dg%bd, dg, 0, 0)
        call fill_corner_buffer(v_N1(:, :, k), v_E1(:, :, k), v_S1(:, :, k), v_W1(:, :, k), dg%bd, dg, 0, 0)
      end do

      call timing_off('cube_buffrmp&cornerbuff')

      call timing_on('fillcornerdomaindecomp')
      call fill_buff_corners_domain_decomp(u_N1, u_E1, u_S1, u_W1, dg, dg%bd, dg%domain_for_duo, 0, 0)
      call fill_buff_corners_domain_decomp(v_N1, v_E1, v_S1, v_W1, dg, dg%bd, dg%domain_for_duo, 0, 0)
      call timing_off('fillcornerdomaindecomp')

      if (ieu_stag == 1 .and. jeu_stag == 0 .and. iev_stag == 0 .and. jev_stag == 1) then  !CGRID
      call timing_on('cubed_a2c_halo_buff')
        call cubed_a2c_halo_buff_only(npz, u_N1, v_N1, bd, dg, dg%bd%isd, dg%bd%ied, dg%bd%je, dg%bd%jed)
        call cubed_a2c_halo_buff_only(npz, u_W1, v_W1, bd, dg, dg%bd%isd, dg%bd%is, dg%bd%jsd, dg%bd%jed)
        call cubed_a2c_halo_buff_only(npz, u_S1, v_S1, bd, dg, dg%bd%isd, dg%bd%ied, dg%bd%jsd, dg%bd%js)
        call cubed_a2c_halo_buff_only(npz, u_E1, v_E1, bd, dg, dg%bd%ie, dg%bd%ied, dg%bd%jsd, dg%bd%jed)
      call timing_off('cubed_a2c_halo_buff')
      else !DGRID
      call timing_on('cubed_a2d_halo_buff')
        call cubed_a2d_halo_buff_only(npz, u_N1, v_N1, bd, dg, dg%bd%isd, dg%bd%ied, dg%bd%je, dg%bd%jed)
        call cubed_a2d_halo_buff_only(npz, u_W1, v_W1, bd, dg, dg%bd%isd, dg%bd%is, dg%bd%jsd, dg%bd%jed)
        call cubed_a2d_halo_buff_only(npz, u_S1, v_S1, bd, dg, dg%bd%isd, dg%bd%ied, dg%bd%jsd, dg%bd%js)
        call cubed_a2d_halo_buff_only(npz, u_E1, v_E1, bd, dg, dg%bd%ie, dg%bd%ied, dg%bd%jsd, dg%bd%jed)
      call timing_off('cubed_a2d_halo_buff')
      end if
      call apply_buff_stag(u_in, dg, bd, u_N1, u_E1, u_S1, u_W1, ieu_stag, jeu_stag)
      call apply_buff_stag(v_in, dg, bd, v_N1, v_E1, v_S1, v_W1, iev_stag, jev_Stag)
    end if

    call timing_off('extvector')
  end subroutine ext_vector
  !===========================================================================
  !===========================================================================

  subroutine cube_rmp_buff(N_buff, E_buff, S_buff, W_buff, dg, bd, domain, istag, jstag, k2e_loc, k2e_coef)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(domain2d), intent(INOUT) :: domain
    integer, intent(in) :: istag, jstag
    real, intent(inout), dimension(bd%isd:, bd%je + jstag:) :: N_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: S_buff
    real, intent(inout), dimension(bd%ie + istag:, bd%jsd:) :: E_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: W_buff
    ! local
!    real, dimension( bd%isd:,         bd%je+jstag:,       1:) :: N_buff_local
!    real, dimension( bd%isd:,         bd%jsd:,            1:) :: S_buff_local
!    real, dimension( bd%ie+istag:,    bd%jsd:,             1:) :: E_buff_local
!    real, dimension( bd%isd:,         bd%jsd:,            1:) :: W_buff_local

    real, dimension(bd%isd:bd%ied + istag, bd%je + jstag:bd%jed + jstag) :: N_buff_local
    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%js) :: S_buff_local
    real, dimension(bd%ie + istag:bd%ied + istag, bd%jsd:bd%jed + jstag) :: E_buff_local
    real, dimension(bd%isd:bd%is, bd%jsd:bd%jed + jstag) :: W_buff_local

    logical :: rmp_w, rmp_e, rmp_s, rmp_n, rmp_sw, rmp_se, rmp_ne, rmp_nw
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: loc, lo, offset, jmin, jmax
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

    call timing_on('cube_buffrmp_routine')
    do ii = 1, ng
      do i = isd, ied + istag
        j = js - ii
        S_buff_local(i, j) = S_buff(i, j)
        j = je + jstag + ii
        N_buff_local(i, j) = N_buff(i, j)
      end do
      do j = jsd, jed + jstag
        i = is - ii
        W_buff_local(i, j) = W_buff(i, j)
        i = ie + istag + ii
        E_buff_local(i, j) = E_buff(i, j)
      end do
    end do

    !--- south
    if (rmp_s) then
!      S_buff_local=S_buff
      do ii = 1, ng

        j = js - ii

        do i = is, ie + istag

          loc = k2e_loc(i, j)
          lo = loc - offset

          !var(i, j) = 0.
          S_buff(i, j) = 0
          do n = 1, dg%k2e_nord
            S_buff(i, j) = S_buff(i, j) + S_buff_local(lo + n, j)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !--- north
    if (rmp_n) then
!      N_buff_local=N_buff
      do ii = 1, ng

        j = je + jstag + ii

        do i = is, ie + istag

          loc = k2e_loc(i, j)
          lo = loc - offset

          !var(i, j) = 0.
          N_buff(i, j) = 0
          do n = 1, dg%k2e_nord
            N_buff(i, j) = N_buff(i, j) + N_buff_local(lo + n, j)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    !--- west
    if (rmp_w) then
!      W_buff_local=W_buff
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

          !var(i, j) = 0.
          W_buff(i, j) = 0
          do n = 1, dg%k2e_nord
            W_buff(i, j) = W_buff(i, j) + W_buff_local(i, lo + n)*k2e_coef(n, i, j)
!write(34000+mpp_pe()) 'ii i j n var', ii,i,j,n, w_buff(i,j)
          end do

        end do

      end do
    end if

    !--- east
    if (rmp_e) then
!      E_buff_local=E_buff
      do ii = 1, ng

        i = ie + istag + ii

        do j = js, je + jstag

          loc = k2e_loc(i, j)
          lo = loc - offset

          !var(i, j) = 0.
          E_buff(i, j) = 0
          do n = 1, dg%k2e_nord
            E_buff(i, j) = E_buff(i, j) + E_buff_local(i, lo + n)*k2e_coef(n, i, j)
          end do

        end do

      end do
    end if

    call timing_off('cube_buffrmp_routine')
    !move this outside of cube_rmp to do the comms 3D
    !call fill_corners_domain_decomp(var, dg, bd, domain, istag, jstag)

  end subroutine cube_rmp_buff

  subroutine create_buff_2d(var, dg, bd, N_buff, E_buff, S_buff, W_buff, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:) :: var
    real, intent(out), dimension(bd%isd:, bd%je + jstag:) :: N_buff
    real, intent(out), dimension(bd%isd:, bd%jsd:) :: S_buff
    real, intent(out), dimension(bd%ie + istag:, bd%jsd:) :: E_buff
    real, intent(out), dimension(bd%isd:, bd%jsd:) :: W_buff

    ! local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: k, dims(3), npz

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

    if (dg%rmp_n) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = je + jstag + ii
          N_buff(i, j) = var(i, j)
        end do
      end do
    end if

    if (dg%rmp_e) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = ie + istag + ii
          E_buff(i, j) = var(i, j)
        end do
      end do
    end if

    if (dg%rmp_s) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = js - ii
          S_buff(i, j) = var(i, j)
        end do
      end do
    end if

    if (dg%rmp_w) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = is - ii
          W_buff(i, j) = var(i, j)
        end do
      end do
    end if

  end subroutine create_buff_2d

  subroutine create_buff(var, dg, bd, N_buff, E_buff, S_buff, W_buff, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, 1:) :: var
    real, intent(out), dimension(bd%isd:, bd%je + jstag:, 1:) :: N_buff
    real, intent(out), dimension(bd%isd:, bd%jsd:, 1:) :: S_buff
    real, intent(out), dimension(bd%ie + istag:, bd%jsd:, 1:) :: E_buff
    real, intent(out), dimension(bd%isd:, bd%jsd:, 1:) :: W_buff

    ! local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: k, dims(3), npz

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
    dims = shape(var)
    npz = dims(3)

!$OMP parallel do default(none) shared(npz,dg,ng,is,ie,js,je,isd,ied,jsd,jed,istag,jstag,var,N_buff,S_buff,W_buff,E_buff)
    do k = 1, npz
    if (dg%rmp_n) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = je + jstag + ii
          N_buff(i, j, k) = var(i, j, k)
!if (k==1) write(mpp_pe()+321100,*) 'ijnbuff', i,j,n_buff(i,j,k), istag, jstag
        end do
      end do
    end if

    if (dg%rmp_e) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = ie + istag + ii
          E_buff(i, j, k) = var(i, j, k)
        end do
      end do
    end if

    if (dg%rmp_s) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = js - ii
          S_buff(i, j, k) = var(i, j, k)
        end do
      end do
    end if

    if (dg%rmp_w) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = is - ii
          W_buff(i, j, k) = var(i, j, k)
        end do
      end do
    end if
    end do

  end subroutine create_buff

  subroutine apply_buff_2d(var, dg, bd, N_buff, E_buff, S_buff, W_buff, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:) :: var
    real, intent(in), dimension(bd%isd:, bd%je + jstag:) :: N_buff
    real, intent(in), dimension(bd%isd:, bd%jsd:) :: S_buff
    real, intent(in), dimension(bd%ie + istag:, bd%jsd:) :: E_buff
    real, intent(in), dimension(bd%isd:, bd%jsd:) :: W_buff
    ! local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: k, dims(3), npz

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

    if (dg%rmp_n) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = je + jstag + ii
          var(i, j) = N_buff(i, j)
!if (k==1) write(mpp_pe()+323200,*) 'ijnbuff', i,j,n_buff(i,j,k), istag, jstag
        end do
      end do
    end if

    if (dg%rmp_e) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = ie + istag + ii
          var(i, j) = E_buff(i, j)
!if (k==1) write(mpp_pe()+323210,*) 'ijebuff', i,j,e_buff(i,j,k)
        end do
      end do
    end if

    if (dg%rmp_s) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = js - ii
          var(i, j) = S_buff(i, j)
        end do
      end do
    end if

    if (dg%rmp_w) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = is - ii
          var(i, j) = W_buff(i, j)
        end do
      end do
    end if

  end subroutine apply_buff_2d

!! APPLY 4th haloes buffs to normal three haloes u,v,uc,vc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine apply_buff_stag(var, dg, bd, N_buff, E_buff, S_buff, W_buff, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, 1:) :: var
    real, intent(in), dimension(dg%bd%isd:, dg%bd%je:, 1:) :: N_buff
    real, intent(in), dimension(dg%bd%isd:, dg%bd%jsd:, 1:) :: S_buff
    real, intent(in), dimension(dg%bd%ie:, dg%bd%jsd:, 1:) :: E_buff
    real, intent(in), dimension(dg%bd%isd:, dg%bd%jsd:, 1:) :: W_buff
    ! local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: k, dims(3), npz

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
    dims = shape(var)
    npz = dims(3)

    if (dg%rmp_s) then
      ! south
!$OMP parallel do default(none) shared(npz,jsd,js,isd,ied,istag,var,S_buff)
      do k = 1, npz
        do j = jsd, js - 1
          do i = isd, ied + istag
            var(i, j, k) = S_buff(i, j, k)
          end do
        end do
      end do
    end if

    if (dg%rmp_n) then
      ! north
!$OMP parallel do default(none) shared(npz,je,jed,isd,ied,istag,jstag,var,N_buff)
      do k = 1, npz
        do j = je + jstag + 1, jed + jstag
          do i = isd, ied + istag
            var(i, j, k) = N_buff(i, j, k)
          end do
        end do
      end do
    end if

    if (dg%rmp_w) then
      ! west
!$OMP parallel do default(none) shared(npz,jsd,jed,isd,is,jstag,var,W_buff)
      do k = 1, npz
        do j = jsd, jed + jstag
          do i = isd, is - 1
            var(i, j, k) = W_buff(i, j, k)
          end do
        end do
      end do
    end if

    if (dg%rmp_e) then
      ! east
!$OMP parallel do default(none) shared(npz,jsd,jed,ie,ied,jstag,istag,var,E_buff)
      do k = 1, npz
        do j = jsd, jed + jstag
          do i = ie + istag + 1, ied + istag
            var(i, j, k) = E_buff(i, j, k)
          end do
        end do
      end do
    end if

  end subroutine apply_buff_stag

  subroutine apply_buff(var, dg, bd, N_buff, E_buff, S_buff, W_buff, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag
    real, dimension(bd%isd:, bd%jsd:, 1:) :: var
    real, intent(in), dimension(bd%isd:, bd%je + jstag:, 1:) :: N_buff
    real, intent(in), dimension(bd%isd:, bd%jsd:, 1:) :: S_buff
    real, intent(in), dimension(bd%ie + istag:, bd%jsd:, 1:) :: E_buff
    real, intent(in), dimension(bd%isd:, bd%jsd:, 1:) :: W_buff
    ! local
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, ii, n
    integer :: k, dims(3), npz

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
    dims = shape(var)
    npz = dims(3)

!$OMP parallel do default(none) shared(npz,dg,ng,is,ie,js,je,isd,ied,jsd,jed,istag,jstag,var,N_buff,S_buff,W_buff,E_buff)
    do k = 1, npz
    if (dg%rmp_n) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = je + jstag + ii
          var(i, j, k) = N_buff(i, j, k)
!if (k==1) write(mpp_pe()+323200,*) 'ijnbuff', i,j,n_buff(i,j,k), istag, jstag
        end do
      end do
    end if

    if (dg%rmp_e) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = ie + istag + ii
          var(i, j, k) = E_buff(i, j, k)
!if (k==1) write(mpp_pe()+323210,*) 'ijebuff', i,j,e_buff(i,j,k)
        end do
      end do
    end if

    if (dg%rmp_s) then
      do ii = 1, ng
        do i = isd, ied + istag
          j = js - ii
          var(i, j, k) = S_buff(i, j, k)
        end do
      end do
    end if

    if (dg%rmp_w) then
      do ii = 1, ng
        !do j = js, je + jstag
        do j = jsd, jed + jstag
          i = is - ii
          var(i, j, k) = W_buff(i, j, k)
        end do
      end do
    end if
    end do

  end subroutine apply_buff

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

  subroutine fill_buff_corners_domain_decomp_2d(N_buff, E_buff, S_buff, W_buff, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(domain2d), intent(IN)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag

    real, intent(inout), dimension(bd%isd:, bd%je + jstag:) :: N_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: S_buff
    real, intent(inout), dimension(bd%ie + istag:, bd%jsd:) :: E_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: W_buff
!    real, dimension(bd%isd:, bd%jsd:) :: var
    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, layoutx, layouty, gid, kk
    integer :: tile(1), npes_per_tile
    real, dimension(1:bd%ng, 1:bd%ng) :: corner
!    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: var
    integer :: upperpelist_y(dg%layout(2)-1)
    integer :: lowerpelist_y(dg%layout(2)-1)
    integer :: upperpelist_x(dg%layout(1)-1)
    integer :: lowerpelist_x(dg%layout(1)-1)

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

        upperpelist_y = (/ (kk * layoutx + (tile(1) - 1) * npes_per_tile, kk = 1, layouty - 1) /)
        lowerpelist_y = (/ ((kk - 1) * layoutx + (tile(1) - 1) * npes_per_tile, kk = 1, layouty - 1) /)

          if (any(gid == lowerpelist_y)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = w_buff(isd + i - 1, je - ng + j)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=(gid+layoutx))
          end if

          if (any(gid == upperpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=(gid - layoutx))
            do j = 1, ng
              do i = 1, ng
                w_buff(isd + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          if (any(gid == upperpelist_y)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = w_buff(isd + i - 1, js + jstag + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-layoutx)
          end if

          if (any(gid == lowerpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=(gid+layoutx))
            do j = 1, ng
              do i = 1, ng
                w_buff(isd + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

      end if !endif rmp_w

!!!!!!!!!!!!!!!!!!
!!!!!! EAST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_e .and. layouty > 1) then

        upperpelist_y = (/ (kk * layoutx + (tile - 1) * npes_per_tile +layoutx-1, kk = 1, layouty - 1) /)
        lowerpelist_y = (/ ((kk - 1) * layoutx + (tile - 1) * npes_per_tile+layoutx-1, kk = 1, layouty - 1) /)

          if (any(gid == lowerpelist_y)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = e_buff(ie + istag + 1 + i - 1, je - ng + j)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=(gid+layoutx))
          end if

          if (any(gid == upperpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=(gid - layoutx))
            do j = 1, ng
              do i = 1, ng
                e_buff(ie + istag + 1 + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          if (any(gid == upperpelist_y)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = e_buff(ie + istag + 1 + i - 1, js + jstag + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=(gid-layoutx))
          end if

          if (any(gid == lowerpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=(gid + layoutx))
            do j = 1, ng
              do i = 1, ng
                e_buff(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!SOUTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_s .and. layoutx > 1) then

        upperpelist_x = (/ (kk + (tile(1)-1)*npes_per_tile, kk = 1, layoutx - 1) /)
        lowerpelist_x = (/ ((kk - 1) + (tile(1) - 1) * npes_per_tile, kk = 1, layoutx - 1) /)

          if (any(gid == lowerpelist_x)) then
          do j = 1, ng
            do i = 1, ng
              corner(i, j) = s_buff(ie - ng + i, jsd + j - 1)
            end do
          end do
          call mpp_send(corner, size(corner), to_pe=gid+1)
          end if

          if (any(gid == upperpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid-1)
            do j = 1, ng
              do i = 1, ng
                s_buff(isd + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

          if (any(gid == upperpelist_x)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = s_buff(is + istag + i - 1, jsd + j - 1)
              end do
          end do
          call mpp_send(corner, size(corner), to_pe=gid-1)
          end if

          if (any(gid == lowerpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid+1)
            do j = 1, ng
              do i = 1, ng
                s_buff(ie + istag + 1 + i - 1, jsd + j - 1) = corner(i, j)
              end do
            end do
          end if

      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!NORTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_n .and. layoutx > 1) then

        upperpelist_x = (/ (kk + (layoutx*(layouty - 1)) + (tile(1)-1)*npes_per_tile, kk = 1, layoutx - 1) /)
        lowerpelist_x = (/ ((kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1) * npes_per_tile, kk = 1, layoutx - 1) /)

          if (any(gid == lowerpelist_x)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = n_buff(ie - ng + i, je + jstag + 1 + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid+1)
          end if

          if (any(gid==upperpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid-1)
            do j = 1, ng
              do i = 1, ng
                n_buff(isd + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

          if (any(gid==upperpelist_x)) then
            do j = 1, ng
              do i = 1, ng
                corner(i, j) = n_buff(is + istag + i - 1, je + jstag + 1 + j - 1)
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-1)
          end if

          if (any(gid==lowerpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid+1)
            do j = 1, ng
              do i = 1, ng
                n_buff(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1) = corner(i, j)
              end do
            end do
          end if

      end if ! rmpn

    end if ! layoux+layouty>2

  end subroutine fill_buff_corners_domain_decomp_2d

  subroutine fill_buff_corners_domain_decomp_3d(N_buff, E_buff, S_buff, W_buff, dg, bd, domain, istag, jstag)
    type(duogrid_type), intent(in) :: dg
    type(domain2d), intent(IN)         :: domain
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: istag, jstag

    !real, intent(inout), dimension( bd%isd:bd%ied + istag, 1:bd%ng, 1:) :: N_buff
    !real, intent(inout), dimension( bd%isd:bd%ied + istag, 1:bd%ng, 1:) :: S_buff
    !real, intent(inout), dimension( 1:bd%ng, bd%jsd:bd%jed + jstag, 1:) :: E_buff
    !real, intent(inout), dimension( 1:bd%ng, bd%jsd:bd%jed + jstag, 1:) :: W_buff
    real, intent(inout), dimension(bd%isd:, bd%je + jstag:, 1:) :: N_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:, 1:) :: S_buff
    real, intent(inout), dimension(bd%ie + istag:, bd%jsd:, 1:) :: E_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:, 1:) :: W_buff
    !real, dimension(bd%isd:, bd%jsd:, 1:) :: var
    real, allocatable, dimension(:, :, :) :: corner
    ! real, dimension(1:bd%ng, 1:bd%ng) :: corner

    integer :: is, ie, js, je, isd, ied, jsd, jed, ng
    integer :: i, j, layoutx, layouty, gid, k, kk
    integer :: tile(1), npes_per_tile, npz, dims(3)
    integer :: upperpelist_y(dg%layout(2)-1)
    integer :: lowerpelist_y(dg%layout(2)-1)
    integer :: upperpelist_x(dg%layout(1)-1)
    integer :: lowerpelist_x(dg%layout(1)-1)

    dims = shape(N_buff)
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

        upperpelist_y = (/ (kk * layoutx + (tile(1) - 1) * npes_per_tile, kk = 1, layouty - 1) /)
        lowerpelist_y = (/ ((kk - 1) * layoutx + (tile(1) - 1) * npes_per_tile, kk = 1, layouty - 1) /)

! Lower to upper
!!!!!!!!!!!!!!!!
          if (any(gid == lowerpelist_y)) then
!$OMP parallel do default(none) shared(npz,ng,corner,w_buff,isd,je)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = w_buff(isd + i - 1, je - ng + j, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=(gid+layoutx))
          end if

          if (any(gid == upperpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=(gid - layoutx))
!$OMP parallel do default(none) shared(npz,ng,corner,w_buff,isd,jsd)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  w_buff(isd + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

! Upper to Lower
!!!!!!!!!!!!!!!!
          if (any(gid == upperpelist_y)) then
!$OMP parallel do default(none) shared(npz,ng,corner,w_buff,isd,js,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = w_buff(isd + i - 1, js + jstag + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-layoutx)
          end if

          if (any(gid == lowerpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=gid+layoutx)
!$OMP parallel do default(none) shared(npz,ng,corner,w_buff,isd,je,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  w_buff(isd + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

      end if !endif rmp_w

!!!!!!!!!!!!!!!!!!
!!!!!! EAST !!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_e .and. layouty > 1) then

        upperpelist_y = (/ (kk * layoutx + (tile - 1) * npes_per_tile +layoutx-1, kk = 1, layouty - 1) /)
        lowerpelist_y = (/ ((kk - 1) * layoutx + (tile - 1) * npes_per_tile+layoutx-1, kk = 1, layouty - 1) /)

          if (any(gid == lowerpelist_y)) then
!$OMP parallel do default(none) shared(npz,ng,corner,e_buff,ie,je,istag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = e_buff(ie + istag + 1 + i - 1, je - ng + j, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid+layoutx)
          end if

          if (any(gid == upperpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=gid-layoutx)
!$OMP parallel do default(none) shared(npz,ng,corner,e_buff,ie,jsd,istag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  e_buff(ie + istag + 1 + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          if (any(gid == upperpelist_y)) then
!$OMP parallel do default(none) shared(npz,ng,corner,e_buff,ie,js,istag,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = e_buff(ie + istag + 1 + i - 1, js + jstag + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-layoutx)
          end if

          if (any(gid == lowerpelist_y)) then
            call mpp_recv(corner, size(corner), from_pe=gid+layoutx)
!$OMP parallel do default(none) shared(npz,ng,corner,e_buff,ie,je,istag,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  e_buff(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!SOUTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_s .and. layoutx > 1) then

        upperpelist_x = (/ (kk + (tile(1)-1)*npes_per_tile, kk = 1, layoutx - 1) /)
        lowerpelist_x = (/ ((kk - 1) + (tile(1) - 1) * npes_per_tile, kk = 1, layoutx - 1) /)

          if (any(gid == lowerpelist_x)) then
!$OMP parallel do default(none) shared(npz,ng,corner,s_buff,ie,jsd)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = s_buff(ie - ng + i, jsd + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid+1)
          end if

          if (any(gid == upperpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid-1)
!$OMP parallel do default(none) shared(npz,ng,corner,s_buff,isd,jsd)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  s_buff(isd + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          if (any(gid == upperpelist_x)) then
!$OMP parallel do default(none) shared(npz,ng,corner,s_buff,is,jsd,istag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = s_buff(is + istag + i - 1, jsd + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-1)
          end if

          if (any(gid == lowerpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid+1)
!$OMP parallel do default(none) shared(npz,ng,corner,s_buff,ie,istag,jsd)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  s_buff(ie + istag + 1 + i - 1, jsd + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

      end if

!!!!!!!!!!!!!!!!!!
!!!!!!!NORTH!!!!!!
!!!!!!!!!!!!!!!!!!

      if (dg%rmp_n .and. layoutx > 1) then

        upperpelist_x = (/ (kk + (layoutx*(layouty - 1)) + (tile(1)-1)*npes_per_tile, kk = 1, layoutx - 1) /)
        lowerpelist_x = (/ ((kk - 1) + (layoutx*(layouty - 1)) + (tile(1) - 1) * npes_per_tile, kk = 1, layoutx - 1) /)

          if (any(gid == lowerpelist_x)) then
!$OMP parallel do default(none) shared(npz,ng,corner,n_buff,ie,je,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = n_buff(ie - ng + i, je + jstag + 1 + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid+1)
          end if

          if (any (gid==upperpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid-1)
!$OMP parallel do default(none) shared(npz,ng,corner,n_buff,isd,je,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  n_buff(isd + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

          if (any(gid==upperpelist_x)) then
!$OMP parallel do default(none) shared(npz,ng,corner,n_buff,is,je,istag,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  corner(i, j, k) = n_buff(is + istag + i - 1, je + jstag + 1 + j - 1, k)
                end do
              end do
            end do
            call mpp_send(corner, size(corner), to_pe=gid-1)
          end if

          if (any(gid==lowerpelist_x)) then
            call mpp_recv(corner, size(corner), from_pe=gid+1)
!$OMP parallel do default(none) shared(npz,ng,corner,n_buff,ie,je,istag,jstag)
            do k = 1, npz
              do j = 1, ng
                do i = 1, ng
                  n_buff(ie + istag + 1 + i - 1, je + jstag + 1 + j - 1, k) = corner(i, j, k)
                end do
              end do
            end do
          end if

      end if ! rmpn

    end if ! layoux+layouty>2

    deallocate (corner)

  end subroutine fill_buff_corners_domain_decomp_3d

  subroutine fill_corner_buffer(N_buff, E_buff, S_buff, W_buff, bd, dg, istag, jstag)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), target, intent(in) :: dg
    integer, intent(IN) :: istag, jstag
    real, intent(inout), dimension(bd%isd:, bd%je + jstag:) :: N_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: S_buff
    real, intent(inout), dimension(bd%ie + istag:, bd%jsd:) :: E_buff
    real, intent(inout), dimension(bd%isd:, bd%jsd:) :: W_buff

    integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, ng
    !integer :: interporder = 3 !order is interporder+1

    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: veltemp
    real, dimension(bd%isd:bd%ied + istag, bd%jsd:bd%jed + jstag) :: veltempp
    real :: temp1, temp2

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xm
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: yp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ym

    call timing_on('lagbuff')

    if (istag == 0 .and. jstag == 0) then
      xp => dg%xp
      xm => dg%xm
      yp => dg%yp
      ym => dg%ym
    end if
    if (istag == 1 .and. jstag == 1) then
      xp => dg%xp3
      xm => dg%xm3
      yp => dg%yp3
      ym => dg%ym3
    end if

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    ng= bd%ng

    !lastpoint in compute domain
    ie = ie + istag
    je = je + jstag

      call timing_on('fillcornerbuffer_routine')
    ! NE corner
!     |   1--2--3
!     |   |     |
!     |   4  5  6
!     |   |     |
!     |   7--8--9
!     -------------

! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! TAKING THE AVERAGE OF BOTH SIDES INTERPOLATIONS
! IS GIVING LESS ERRORS COMPARED TO ONE INTERPOLATION FORM EACH BUFFER
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    if (dg%rmp_ne) then

!     do i=ie+1,ie+bd%ng
!       do j=je+1,je+bd%ng
!         if (j>i) then ! 1 2 4
!           call lagrange_poly_buff(N_buff, i, j, bd, dg, 'X+', istag, jstag, interporder, isd, je)
!           E_buff(i,j)=N_buff(i,j)
!         elseif (i>j) then !8 9 6
!           call lagrange_poly_buff(E_buff, i, j, bd, dg, 'Y+', istag, jstag, interporder, ie, jsd)
!           N_buff(i,j)=E_buff(i,j)
!         endif
!       enddo
!     enddo

!diag
      do i = ie + 1, ie + ng
        do j = je + 1, je + ng
!         if (j==i) then ! 7 5 3
          call lagrange_poly_buff_fast(N_buff(:,j), i, xp(:,i,j,2*istag+jstag+1), istag, jstag, interporder, isd, bd%ie)
!          call lagrange_poly_buff(N_buff, i, j, bd, dg, 'X+', istag, jstag, interporder, isd, je)
          temp1 = N_buff(i, j)
          call lagrange_poly_buff_fast(E_buff(i,:), j, yp(:,i,j,2*istag+jstag+1), istag, jstag, interporder, jsd, bd%je)

         ! call lagrange_poly_buff(E_buff, i, j, bd, dg, 'Y+', istag, jstag, interporder, ie, jsd)
          temp2 = E_buff(i, j)
          N_buff(i, j) = (temp1 + temp2)*0.5
!           if (i==ie+1) then
!           N_buff(i,j)=0.5*(N_buff(i+1, j) + N_buff(i, j+1))
!           elseif (i==ie+bd%ng) then
!            N_buff(i,j)=0.5*(N_buff(i-1, j) + N_buff(i, j-1))
!           else
!            N_buff(i,j)=0.25*(N_buff(i-1, j) + N_buff(i, j-1) + N_buff(i+1,j) + N_buff(i,j+1))
!           endif
!         endif
          E_buff(i, j) = N_buff(i, j)
        end do
      end do

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

!     do i=is-1,is-bd%ng,-1
!       do j=je+1,je+bd%ng
!         if (j-(je+1)>abs(i)) then ! 1 2 4
!           call lagrange_poly_buff(N_buff, i, j, bd, dg, 'X-', istag, jstag, interporder, isd, je)
!           W_buff(i,j)=N_buff(i,j)
!         elseif (abs(i)>j-(je+1)) then !8 9 6
!           call lagrange_poly_buff(W_buff, i, j, bd, dg, 'Y+', istag, jstag, interporder, isd, jsd)
!           N_buff(i,j)=W_buff(i,j)
!         endif
!       enddo
!     enddo

!diag
      do i = is - 1, is - ng, -1
        do j = je + 1, je + ng
          !if (j-(je+1)==abs(i)) then ! 7 5 3
          !call lagrange_poly_buff(N_buff, i, j, bd, dg, 'X-', istag, jstag, interporder, isd, je)
          call lagrange_poly_buff_fast_minus(N_buff(:,j), i, xm(:,i,j,2*istag+jstag+1), istag, jstag, interporder, isd, bd%is)
          temp1 = N_buff(i, j)
          !call lagrange_poly_buff(W_buff, i, j, bd, dg, 'Y+', istag, jstag, interporder, isd, jsd)
          call lagrange_poly_buff_fast(W_buff(i,:), j, yp(:,i,j,2*istag+jstag+1), istag, jstag, interporder, jsd, bd%je)
          temp2 = W_buff(i, j)
          N_buff(i, j) = (temp1 + temp2)*0.5
          !   if (i==is-1) then
          !    N_buff(i,j)=0.5*(N_buff(i-1, j) + N_buff(i, j+1))
          !   elseif (i==is-bd%ng) then
          !    N_buff(i,j)=0.5*(N_buff(i+1, j) + N_buff(i, j-1))
          !   else
          !    N_buff(i,j)=0.25*(N_buff(i-1, j) + N_buff(i, j-1) + N_buff(i+1,j) + N_buff(i,j+1))
          !   endif
          ! endif
          W_buff(i, j) = N_buff(i, j)
        end do
      end do

    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! SE corner
!!     --------------
!!     |   1--2--3
!!     |   |     |
!!     |   4  5  6
!!     |   |     |
!!     |   7--8--9
!!
    if (dg%rmp_se) then

!     do i=ie+1,ie+bd%ng
!       do j=js-1,js-bd%ng,-1
!         if (i-(ie+1)>abs(j)) then !2 3 6
!           call lagrange_poly_buff(E_buff, i, j, bd, dg, 'Y-', istag, jstag, interporder, ie, jsd)
!           S_buff(i,j)=E_buff(i,j)
!         elseif (abs(j)>i-(ie+1)) then !4 7 8
!           call lagrange_poly_buff(S_buff, i, j, bd, dg, 'X+', istag, jstag, interporder, isd, jsd)
!           E_buff(i,j)=S_buff(i,j)
!         endif
!       enddo
!     enddo

!diag
      do i = ie + 1, ie + ng
        do j = js - 1, js - ng, -1
          !if (i-(ie+1)==abs(j)) then ! 7 5 3
          !call lagrange_poly_buff(E_buff, i, j, bd, dg, 'Y-', istag, jstag, interporder, ie, jsd)
          call lagrange_poly_buff_fast_minus(E_buff(i,:), j, ym(:,i,j,2*istag+jstag+1), istag, jstag, interporder, jsd, bd%js)
          temp1 = E_buff(i, j)
          !call lagrange_poly_buff(S_buff, i, j, bd, dg, 'X+', istag, jstag, interporder, isd, jsd)
          call lagrange_poly_buff_fast(S_buff(:,j), i, xp(:,i,j,2*istag+jstag+1), istag, jstag, interporder, isd, bd%ie)
          temp2 = S_buff(i, j)
          S_buff(i, j) = (temp1 + temp2)*0.5
          !   if (j==js-1) then
          !    S_buff(i,j)=0.5*(S_buff(i+1, j) + S_buff(i, j-1))
          !   elseif (j==js-bd%ng) then
          !    S_buff(i,j)=0.5*(S_buff(i-1, j) + S_buff(i, j+1))
          !   else
          !    S_buff(i,j)=0.25*(S_buff(i-1, j) + S_buff(i, j-1) + S_buff(i+1,j) + S_buff(i,j+1))
          !   endif
          ! endif
          E_buff(i, j) = S_buff(i, j)
        end do
      end do

    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ! SW corner
!!     -------------
!!        1--2--3  |
!!        |     |  |
!!        4  5  6  |
!!        |     |  |
!!        7--8--9  |
!!
!
    if (dg%rmp_sw) then

!     do i=is-1,isd,-1
!       do j=js-1,jsd,-1
!         if (abs(j)>abs(i)) then ! 1 2 4
!           call lagrange_poly_buff(S_buff, i, j, bd, dg, 'X-', istag, jstag, interporder, isd, jsd)
!           W_buff(i,j)=S_buff(i,j)
!         elseif (abs(i)>abs(j)) then !8 9 6
!           call lagrange_poly_buff(W_buff, i, j, bd, dg, 'Y-', istag, jstag, interporder, isd, jsd)
!           S_buff(i,j)=W_buff(i,j)
!         endif
!       enddo
!     enddo

!diag
      do i = is - 1, isd, -1
        do j = js - 1, jsd, -1
          ! if (j==i) then ! 7 5 3
          !call lagrange_poly_buff(S_buff, i, j, bd, dg, 'X-', istag, jstag, interporder, isd, jsd)
          call lagrange_poly_buff_fast_minus(S_buff(:,j), i, xm(:,i,j,2*istag+jstag+1), istag, jstag, interporder, isd, bd%is)
          temp1 = S_buff(i, j)
          !call lagrange_poly_buff(W_buff, i, j, bd, dg, 'Y-', istag, jstag, interporder, isd, jsd)
          call lagrange_poly_buff_fast_minus(W_buff(i,:), j, ym(:,i,j,2*istag+jstag+1), istag, jstag, interporder, jsd, bd%js)
          temp2 = W_buff(i, j)
          S_buff(i, j) = (temp1 + temp2)*0.5
          !if (i==isd) then
          ! S_buff(i,j)=0.5*(S_buff(i+1, j) + S_buff(i, j+1))
          !elseif (i==is-1) then
          ! S_buff(i,j)=0.5*(S_buff(i-1, j) + S_buff(i, j-1))
          !else
          ! S_buff(i,j)=0.25*(S_buff(i-1, j) + S_buff(i, j-1) + S_buff(i+1,j) + S_buff(i,j+1))
          !endif
          ! endif
          W_buff(i, j) = S_buff(i, j)
        end do
      end do

    end if

      call timing_off('fillcornerbuffer_routine')
  end subroutine fill_corner_buffer

! subroutine fill_corner_region will fill the corner regions of pes lying on
! the corners of the cubes sphere nw, ne, se, sw.
! the data is extrapolated from the haloe regions present on the pe following
! a logic similar to ZA22's using a lagrangian interpolation function.

! Ideally, values should be interoplated (no extrapolation needed) using the exact
! ZA22 formulation; however, this is not done here yet and requires more complicated work.

! In addition, using the lagrangian interpolation function for extrapolation
! purposes is not highly recommended in general; however, it is tolerable in this case given that
! the extrapolation is performed for couple points just outside of the domain on which the data reside.

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
  subroutine compute_lagrange_coeff(xp, xm, yp, ym, bd, dg, gg)

    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(inout) :: dg
    type(global_grid_type), intent(in) :: gg
    real(kind=R_GRID), intent(out) :: xp(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: xm(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: yp(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: ym(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4)

    !local
    integer ii, iii, it, jstag, istag
    integer is, js, ie, je, i, j
    integer isd, jsd, ied, jed
    integer istagup, jstagup
    real(kind=R_GRID) :: ss, dist1, dist2, S

    if (bd%isd == dg%bd%isd) then
      istagup = 0
      jstagup = 0
      if (mpp_pe() == mpp_root_pe()) print *, "Generating interp coeff for ng4"
    else
      istagup = 1
      jstagup = 1
      if (mpp_pe() == mpp_root_pe()) print *, "Generating interp coeff for ng3"
    end if
    is = bd%is
    js = bd%js
    !xp(:,:,:,:)=-3232.
    !do istag = 0, 1
    !do jstag = 0, 1
    do istag = 0, istagup
    do jstag = 0, jstagup
    do i = bd%ie - interporder, bd%ied + istag
    do j = bd%jsd, bd%jed + jstag

      do ii = bd%ie - interporder, bd%ie
        ss = 1.0
        do iii = bd%ie - interporder, bd%ie
          if (ii < iii) it = -1
          if (ii > iii) it = 1
          if (ii .eq. iii) cycle
          if (j > bd%je - 1) then
            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_kik(:, iii, j))
            dist2 = great_circle_dist(dg%a_kik(:, ii, j), dg%a_kik(:, iii, j))
            if (great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_kik(:, bd%ie - interporder, j))  <  great_circle_dist(dg%a_kik(:, iii, j), dg%a_kik(:, bd%ie - interporder, j)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, ii, j), dg%a_kik(:, bd%ie - interporder, j))        <  great_circle_dist(dg%a_kik(:, iii, j), dg%a_kik(:, bd%ie - interporder, j)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, iii, j - jstag))
            dist2 = great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))
            if (great_circle_dist(dg%a_pt(:, i - istag, j-jstag), dg%a_kik(:, bd%ie - interporder, j-jstag))  <  great_circle_dist(dg%a_kik(:, iii, j-jstag), dg%a_kik(:, bd%ie - interporder, j-jstag)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, ii, j-jstag), dg%a_kik(:, bd%ie - interporder, j-jstag))        <  great_circle_dist(dg%a_kik(:, iii, j-jstag), dg%a_kik(:, bd%ie - interporder, j-jstag)) ) dist2=dist2*(-1)
!if (mpp_pe()==0) print*, i,j, ii,iii,dist1,dist2
            ss = ss*dist1/dist2

            ! ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
            !             great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, iii, j - jstag))/ &
            !            (great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))+0.000000001))
          end if
        end do
!write(1200+mpp_pe(),*) istag, jstag, i, j, ii, ss
        !xp(bd%ie - ii + 1, i, j, 2*istag + jstag + 1) = ss
        xp(ii - bd%ie + interporder + 1, i, j, 2*istag + jstag + 1) = ss
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j), dg%a_kik(:, iii, j))
            dist2 = great_circle_dist(dg%a_kik(:, ii, j), dg%a_kik(:, iii, j))
            if (great_circle_dist(dg%a_pt(:, i , j), dg%a_kik(:, bd%is + interporder, j))  <  great_circle_dist(dg%a_kik(:, iii, j), dg%a_kik(:, bd%is + interporder, j)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, ii, j), dg%a_kik(:, bd%is + interporder, j))        <  great_circle_dist(dg%a_kik(:, iii, j), dg%a_kik(:, bd%is + interporder, j)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !  ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, iii, j))/ &
            !              great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j)))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, iii, j - jstag))
            dist2 = great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))
            if (great_circle_dist(dg%a_pt(:, i , j-jstag), dg%a_kik(:, bd%is + interporder, j-jstag))  <  great_circle_dist(dg%a_kik(:, iii, j-jstag), dg%a_kik(:, bd%is + interporder, j-jstag)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, ii, j-jstag), dg%a_kik(:, bd%is + interporder, j-jstag))        <  great_circle_dist(dg%a_kik(:, iii, j-jstag), dg%a_kik(:, bd%is + interporder, j-jstag)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !  ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
            !              great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, iii, j - jstag))/ &
            !            (great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))+0.000001))
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, i, iii))
            dist2 = great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))
            if (great_circle_dist(dg%a_pt(:, i , j-jstag), dg%a_kik(:, i, bd%je - interporder))  <  great_circle_dist(dg%a_kik(:, i, iii), dg%a_kik(:, i, bd%je - interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, bd%je - interporder))        <  great_circle_dist(dg%a_kik(:, i, iii), dg%a_kik(:, i, bd%je - interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            ! ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, i, iii))/ &
            !              great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, i, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))+0.000001))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, i - istag, iii))
            dist2 = great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))
            if (great_circle_dist(dg%a_pt(:, i -istag, j-jstag), dg%a_kik(:, i-istag, bd%je - interporder))  <  great_circle_dist(dg%a_kik(:, i-istag, iii), dg%a_kik(:, i-istag, bd%je - interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, i-istag, ii), dg%a_kik(:, i-istag, bd%je - interporder))        <  great_circle_dist(dg%a_kik(:, i-istag, iii), dg%a_kik(:, i-istag, bd%je - interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, i - istag, iii))/ &
            !            great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, i - istag, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))+0.000001))

          end if

        end do
        !yp(bd%je - ii + 1, i, j, 2*istag + jstag + 1) = ss
        yp(ii - bd%je + interporder + 1, i, j, 2*istag + jstag + 1) = ss
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j), dg%a_kik(:, i, iii))
            dist2 = great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))
            if (great_circle_dist(dg%a_pt(:, i , j), dg%a_kik(:, i, bd%js + interporder))  <  great_circle_dist(dg%a_kik(:, i, iii), dg%a_kik(:, i, bd%js + interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, bd%js + interporder))        <  great_circle_dist(dg%a_kik(:, i, iii), dg%a_kik(:, i, bd%js + interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

!            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, i, iii))/ &
!                        great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_kik(:, i, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))+0.000001))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_kik(:, i - istag, iii))
            dist2 = great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))
            if (great_circle_dist(dg%a_pt(:, i -istag, j), dg%a_kik(:, i-istag, bd%js + interporder))  <  great_circle_dist(dg%a_kik(:, i-istag, iii), dg%a_kik(:, i-istag, bd%js + interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_kik(:, i-istag, ii), dg%a_kik(:, i-istag, bd%js + interporder))        <  great_circle_dist(dg%a_kik(:, i-istag, iii), dg%a_kik(:, i-istag, bd%js + interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

!            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, i - istag, iii))/ &
!                        great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_kik(:, i - istag, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))+0.000001))
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

  subroutine compute_lagrange_coeff_extra(xp, xm, yp, ym, bd, dg, gg)

    type(fv_grid_bounds_type), intent(IN) :: bd
    type(duogrid_type), intent(inout) :: dg
    type(global_grid_type), intent(in) :: gg
    real(kind=R_GRID), intent(out) :: xp(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: xm(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: yp(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4)
    real(kind=R_GRID), intent(out) :: ym(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4)

    !local
    integer ii, iii, it, jstag, istag
    integer is, js, ie, je, i, j
    integer isd, jsd, ied, jed
    integer istagup, jstagup
    real(kind=R_GRID) :: ss, dist1, dist2, S

    if (bd%isd == dg%bd%isd) then
      istagup = 0
      jstagup = 0
      if (mpp_pe() == mpp_root_pe()) print *, "Generating interp coeff for ng4"
    else
      istagup = 1
      jstagup = 1
      if (mpp_pe() == mpp_root_pe()) print *, "Generating interp coeff for ng3"
    end if
    is = bd%is
    js = bd%js
    !xp(:,:,:,:)=-3232.
    !do istag = 0, 1
    !do jstag = 0, 1
    do istag = 0, istagup
    do jstag = 0, jstagup
    do i = bd%ie - interporder, bd%ied + istag
    do j = bd%jsd, bd%jed + jstag

      do ii = bd%ie - interporder, bd%ie
        ss = 1.0
        do iii = bd%ie - interporder, bd%ie
          if (ii < iii) it = -1
          if (ii > iii) it = 1
          if (ii .eq. iii) cycle
          if (j > bd%je - 1) then
            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, iii, j))
            dist2 = great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j))
            if (great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, bd%ie - interporder, j))  <  great_circle_dist(dg%a_pt(:, iii, j), dg%a_pt(:, bd%ie - interporder, j)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, bd%ie - interporder, j))        <  great_circle_dist(dg%a_pt(:, iii, j), dg%a_pt(:, bd%ie - interporder, j)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, iii, j - jstag))
            dist2 = great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag))
            if (great_circle_dist(dg%a_pt(:, i - istag, j-jstag), dg%a_pt(:, bd%ie - interporder, j-jstag))  <  great_circle_dist(dg%a_pt(:, iii, j-jstag), dg%a_pt(:, bd%ie - interporder, j-jstag)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, ii, j-jstag), dg%a_pt(:, bd%ie - interporder, j-jstag))        <  great_circle_dist(dg%a_pt(:, iii, j-jstag), dg%a_pt(:, bd%ie - interporder, j-jstag)) ) dist2=dist2*(-1)
!if (mpp_pe()==0) print*, i,j, ii,iii,dist1,dist2
            ss = ss*dist1/dist2

            ! ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
            !             great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, iii, j - jstag))/ &
            !            (great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))+0.000000001))
          end if
        end do
!write(1200+mpp_pe(),*) istag, jstag, i, j, ii, ss
        !xp(bd%ie - ii + 1, i, j, 2*istag + jstag + 1) = ss
        xp(ii - bd%ie + interporder + 1, i, j, 2*istag + jstag + 1) = ss
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, iii, j))
            dist2 = great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j))
            if (great_circle_dist(dg%a_pt(:, i , j), dg%a_pt(:, bd%is + interporder, j))  <  great_circle_dist(dg%a_pt(:, iii, j), dg%a_pt(:, bd%is + interporder, j)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, bd%is + interporder, j))        <  great_circle_dist(dg%a_pt(:, iii, j), dg%a_pt(:, bd%is + interporder, j)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !  ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, iii, j))/ &
            !              great_circle_dist(dg%a_pt(:, ii, j), dg%a_pt(:, iii, j)))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, iii, j - jstag))
            dist2 = great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag))
            if (great_circle_dist(dg%a_pt(:, i , j-jstag), dg%a_pt(:, bd%is + interporder, j-jstag))  <  great_circle_dist(dg%a_pt(:, iii, j-jstag), dg%a_pt(:, bd%is + interporder, j-jstag)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, ii, j-jstag), dg%a_pt(:, bd%is + interporder, j-jstag))        <  great_circle_dist(dg%a_pt(:, iii, j-jstag), dg%a_pt(:, bd%is + interporder, j-jstag)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !  ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, iii, j - jstag))/ &
            !              great_circle_dist(dg%a_pt(:, ii, j - jstag), dg%a_pt(:, iii, j - jstag)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, iii, j - jstag))/ &
            !            (great_circle_dist(dg%a_kik(:, ii, j - jstag), dg%a_kik(:, iii, j - jstag))+0.000001))
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, i, iii))
            dist2 = great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii))
            if (great_circle_dist(dg%a_pt(:, i , j-jstag), dg%a_pt(:, i, bd%je - interporder))  <  great_circle_dist(dg%a_pt(:, i, iii), dg%a_pt(:, i, bd%je - interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, bd%je - interporder))        <  great_circle_dist(dg%a_pt(:, i, iii), dg%a_pt(:, i, bd%je - interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            ! ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_pt(:, i, iii))/ &
            !              great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j - jstag), dg%a_kik(:, i, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))+0.000001))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, i - istag, iii))
            dist2 = great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii))
            if (great_circle_dist(dg%a_pt(:, i -istag, j-jstag), dg%a_pt(:, i-istag, bd%je - interporder))  <  great_circle_dist(dg%a_pt(:, i-istag, iii), dg%a_pt(:, i-istag, bd%je - interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, i-istag, ii), dg%a_pt(:, i-istag, bd%je - interporder))        <  great_circle_dist(dg%a_pt(:, i-istag, iii), dg%a_pt(:, i-istag, bd%je - interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_pt(:, i - istag, iii))/ &
            !            great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j - jstag), dg%a_kik(:, i - istag, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))+0.000001))

          end if

        end do
        !yp(bd%je - ii + 1, i, j, 2*istag + jstag + 1) = ss
        yp(ii - bd%je + interporder + 1, i, j, 2*istag + jstag + 1) = ss
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

            dist1 = great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, i, iii))
            dist2 = great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii))
            if (great_circle_dist(dg%a_pt(:, i , j), dg%a_pt(:, i, bd%js + interporder))  <  great_circle_dist(dg%a_pt(:, i, iii), dg%a_pt(:, i, bd%js + interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, bd%js + interporder))        <  great_circle_dist(dg%a_pt(:, i, iii), dg%a_pt(:, i, bd%js + interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

!            ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_pt(:, i, iii))/ &
!                        great_circle_dist(dg%a_pt(:, i, ii), dg%a_pt(:, i, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i, j), dg%a_kik(:, i, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i, ii), dg%a_kik(:, i, iii))+0.000001))
          else

            dist1 = great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, i - istag, iii))
            dist2 = great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii))
            if (great_circle_dist(dg%a_pt(:, i -istag, j), dg%a_pt(:, i-istag, bd%js + interporder))  <  great_circle_dist(dg%a_pt(:, i-istag, iii), dg%a_pt(:, i-istag, bd%js + interporder)) ) dist1=dist1*(-1)
            if (great_circle_dist(dg%a_pt(:, i-istag, ii), dg%a_pt(:, i-istag, bd%js + interporder))        <  great_circle_dist(dg%a_pt(:, i-istag, iii), dg%a_pt(:, i-istag, bd%js + interporder)) ) dist2=dist2*(-1)
            ss = ss*dist1/dist2

!            ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_pt(:, i - istag, iii))/ &
!                        great_circle_dist(dg%a_pt(:, i - istag, ii), dg%a_pt(:, i - istag, iii)))
            !ss = ss*it*(great_circle_dist(dg%a_pt(:, i - istag, j), dg%a_kik(:, i - istag, iii))/ &
            !            (great_circle_dist(dg%a_kik(:, i - istag, ii), dg%a_kik(:, i - istag, iii))+0.000001))
          end if

        end do
        ym(js + interporder - ii + 1, i, j, 2*istag + jstag + 1) = ss
      end do
    end do
    end do
    end do
    end do

    if (mpp_pe() == mpp_root_pe()) print *, "DUO lag coeff: DONE!"
  end subroutine compute_lagrange_coeff_extra



  subroutine lagrange_poly_buff_fast(field, loc, arr, istag, jstag, order, ilbound, upperbound)

    integer, intent(in) :: loc, istag, jstag, order, ilbound, upperbound
    real, dimension(ilbound:) :: field
    real(kind=R_GRID), dimension(1:) :: arr

    !local
    integer ii
    real(kind=R_GRID) :: S

    call timing_on('lagbuff')

    S = 0. ! here we sum all the polynomial*solution

!OK with field 2D
!!!!!!!!!!!!!!!!      do ii = ie-order,ie
!!!!!!!!!!!!!!!!        !S = S + xp(ii - ie + order + 1, i, j, 2*istag + jstag + 1)*field(ii + istag, j)
!!!!!!!!!!!!!!!!        S = S + arr(ii - ie + order + 1)*field(ii + istag, j)
!!!!!!!!!!!!!!!!      end do
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!      field(i, j) = S

     do ii = upperbound-order,upperbound
       S = S + arr(ii - upperbound + order + 1)*field(ii + istag)
     end do

     field(loc) = S

    call timing_off('lagbuff')
  end subroutine lagrange_poly_buff_fast



  subroutine lagrange_poly_buff_fast_minus(field, loc, arr, istag, jstag, order, ilbound, upperbound)

    integer, intent(in) :: loc, istag, jstag, order, ilbound, upperbound
    real, dimension(ilbound:) :: field
    real(kind=R_GRID), dimension(1:) :: arr

    !local
    integer ii
    real(kind=R_GRID) :: S

    call timing_on('lagbuff')

    S = 0. ! here we sum all the polynomial*solution

!OK with field 2D
!!!!!!!!!!!!!!!!      do ii = ie-order,ie
!!!!!!!!!!!!!!!!        !S = S + xp(ii - ie + order + 1, i, j, 2*istag + jstag + 1)*field(ii + istag, j)
!!!!!!!!!!!!!!!!        S = S + arr(ii - ie + order + 1)*field(ii + istag, j)
!!!!!!!!!!!!!!!!      end do
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!      field(i, j) = S


     do ii = upperbound+order,upperbound,-1
       S = S + arr(upperbound + order -ii + 1)*field(ii)
     end do

     field(loc) = S

    call timing_off('lagbuff')
  end subroutine lagrange_poly_buff_fast_minus



  subroutine lagrange_poly_buff(field, i, j, bd, dg, direction, istag, jstag, order, ilbound, jlbound, corner)

    integer, intent(in) :: i, j, istag, jstag, order, ilbound, jlbound
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, dimension(ilbound:, jlbound:) :: field
    type(duogrid_type), target, intent(in) :: dg
    character(len=*), intent(in) :: direction
    logical, optional, intent(in) :: corner

    !local
    integer ii, iii, it
    integer is, js, ie, je, isd, jsd
    real(kind=R_GRID) :: ss, S
    !real(kind=R_GRID) :: ss_store(5, 1:interporder+1, bd%isd:bd%ied+1, bd%jsd:bd%jed+1, 0:1,0:1)

    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: xm
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: yp
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ym

    call timing_on('lagbuff')

    if (istag == 0 .and. jstag == 0) then
      xp => dg%xp
      xm => dg%xm
      yp => dg%yp
      ym => dg%ym
    end if
    if (istag == 1 .and. jstag == 1) then
      xp => dg%xp3
      xm => dg%xm3
      yp => dg%yp3
      ym => dg%ym3
    end if

    S = 0. ! here we sum all the polynomial*solution

    is = bd%is
    js = bd%js
    ie = bd%ie
    je = bd%je
    it = 1

    !call compute_lagrange_poly_interp_2d(ss_store, i, j, bd, dg, direction, istag, jstag, order, corner)

!TO check: compute weight coefficients a priori?!

    if (direction == 'X+') then

      do ii = ie-order,ie
        S = S + xp(ii - ie + order + 1, i, j, 2*istag + jstag + 1)*field(ii + istag, j)
      end do

      field(i, j) = S

    end if

    if (direction == 'X-') then

      do ii = is+order,is,-1
        S = S + xm(is + order - ii + 1, i, j, 2*istag + jstag + 1)*field(ii, j)
      end do

      field(i, j) = S

    end if

    if (direction == 'Y+') then

      do ii = je-order, je
        S = S + yp(ii - je + order + 1, i, j, 2*istag + jstag + 1)*field(i, ii + jstag)
      end do

      field(i, j) = S

    end if

    if (direction == 'Y-') then

      do ii = js+order,js,-1
        S = S + ym(js + order - ii + 1, i, j, 2*istag + jstag + 1)*field(i, ii)
      end do

      field(i, j) = S

    end if
    nullify (xp, xm, yp, ym)

    call timing_off('lagbuff')
  end subroutine lagrange_poly_buff

  subroutine cubed_a2d_halo_buff_only(npz, ubuff, vbuff, bd, dg, ilbound, iubound, jlbound, jubound)

    integer, intent(in) ::  npz, ilbound, jlbound, iubound, jubound
    real, intent(inout), dimension(ilbound:iubound, jlbound:jubound, npz) :: ubuff, vbuff
    type(fv_grid_bounds_type), intent(IN) :: bd

    type(duogrid_type), target, intent(in) :: dg

    !local
    integer i, j, k

    real v3(3, ilbound:iubound, jlbound:jubound)
    real ue(3, ilbound:iubound, jlbound:jubound)
    real ve(3, ilbound:iubound, jlbound:jubound)
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ew, es
    real(kind=R_GRID), pointer, dimension(:, :, :)   :: vlon, vlat
      call timing_on('a2d_halo_buff_routine')
    vlon => dg%vlon_ext
    vlat => dg%vlat_ext
    ew => dg%ew_ext
    es => dg%es_ext

!$OMP parallel do default(none) shared(npz,jlbound,jubound,ilbound,iubound,ubuff,vbuff,vlon,vlat,es,ew) &
!$OMP                          private(v3,ue,ve)
    do k = 1, npz

! Compute 3D wind on A grid
      do j = jlbound, jubound
        do i = ilbound, iubound
          v3(1, i, j) = ubuff(i, j, k)*vlon(i, j, 1) + vbuff(i, j, k)*vlat(i, j, 1)
          v3(2, i, j) = ubuff(i, j, k)*vlon(i, j, 2) + vbuff(i, j, k)*vlat(i, j, 2)
          v3(3, i, j) = ubuff(i, j, k)*vlon(i, j, 3) + vbuff(i, j, k)*vlat(i, j, 3)
        end do
      end do

! A --> D
! Interpolate to cell edges
      do j = jlbound + 1, jubound
        do i = ilbound, iubound
          ue(1, i, j) = 0.5*(v3(1, i, j - 1) + v3(1, i, j))
          ue(2, i, j) = 0.5*(v3(2, i, j - 1) + v3(2, i, j))
          ue(3, i, j) = 0.5*(v3(3, i, j - 1) + v3(3, i, j))
        end do
      end do

      do j = jlbound, jubound
        do i = ilbound + 1, iubound
          ve(1, i, j) = 0.5*(v3(1, i - 1, j) + v3(1, i, j))
          ve(2, i, j) = 0.5*(v3(2, i - 1, j) + v3(2, i, j))
          ve(3, i, j) = 0.5*(v3(3, i - 1, j) + v3(3, i, j))
        end do
      end do

      do j = jlbound + 1, jubound
        do i = ilbound, iubound
          ubuff(i, j, k) = ue(1, i, j)*es(1, i, j, 1) + &
                           ue(2, i, j)*es(2, i, j, 1) + &
                           ue(3, i, j)*es(3, i, j, 1)
        end do
      end do
      do j = jlbound, jubound
        do i = ilbound + 1, iubound
          vbuff(i, j, k) = ve(1, i, j)*ew(1, i, j, 2) + &
                           ve(2, i, j)*ew(2, i, j, 2) + &
                           ve(3, i, j)*ew(3, i, j, 2)
        end do
      end do

    end do

    nullify (vlon,vlat,ew,es)
      call timing_off('a2d_halo_buff_routine')
  end subroutine cubed_a2d_halo_buff_only

  subroutine cubed_a2c_halo_buff_only(npz, ubuff, vbuff, bd, dg, ilbound, iubound, jlbound, jubound)

    integer, intent(in) :: npz, ilbound, jlbound, iubound, jubound
    real, intent(inout), dimension(ilbound:iubound, jlbound:jubound, npz) :: ubuff, vbuff
    type(fv_grid_bounds_type), intent(IN) :: bd

    type(duogrid_type), target, intent(in) :: dg

    !local
    integer i, j, k

    real v3(3, ilbound:iubound, jlbound:jubound)
    real ue(3, ilbound:iubound, jlbound:jubound)
    real ve(3, ilbound:iubound, jlbound:jubound)
    real(kind=R_GRID), pointer, dimension(:, :, :, :) :: ew, es
    real(kind=R_GRID), pointer, dimension(:, :, :)   :: vlon, vlat

      call timing_on('a2c_halo_buff_routine')
    vlon => dg%vlon_ext
    vlat => dg%vlat_ext
    ew => dg%ew_ext
    es => dg%es_ext

!$OMP parallel do default(none) shared(npz,jlbound,jubound,ilbound,iubound,ubuff,vbuff,vlon,vlat,es,ew) &
!$OMP                          private(v3,ue,ve)
    do k = 1, npz

! Compute 3D wind on A grid
      do j = jlbound, jubound
        do i = ilbound, iubound
          v3(1, i, j) = ubuff(i, j, k)*vlon(i, j, 1) + vbuff(i, j, k)*vlat(i, j, 1)
          v3(2, i, j) = ubuff(i, j, k)*vlon(i, j, 2) + vbuff(i, j, k)*vlat(i, j, 2)
          v3(3, i, j) = ubuff(i, j, k)*vlon(i, j, 3) + vbuff(i, j, k)*vlat(i, j, 3)
        end do
!      end do

! A --> C
! Interpolate to cell edges
!      do j = jsd, jed
        do i = ilbound + 1, iubound
          ue(1, i, j) = 0.5*(v3(1, i - 1, j) + v3(1, i, j))
          ue(2, i, j) = 0.5*(v3(2, i - 1, j) + v3(2, i, j))
          ue(3, i, j) = 0.5*(v3(3, i - 1, j) + v3(3, i, j))
        end do
      end do

      do j = jlbound + 1, jubound
        do i = ilbound, iubound
          ve(1, i, j) = 0.5*(v3(1, i, j - 1) + v3(1, i, j))
          ve(2, i, j) = 0.5*(v3(2, i, j - 1) + v3(2, i, j))
          ve(3, i, j) = 0.5*(v3(3, i, j - 1) + v3(3, i, j))
        end do
      end do

      do j = jlbound + 1, jubound
        do i = ilbound, iubound
          ubuff(i, j, k) = ue(1, i, j)*ew(1, i, j, 1) + &
                           ue(2, i, j)*ew(2, i, j, 1) + &
                           ue(3, i, j)*ew(3, i, j, 1)
        end do
      end do
      do j = jlbound + 1, jubound
        do i = ilbound, iubound
          vbuff(i, j, k) = ve(1, i, j)*es(1, i, j, 2) + &
                           ve(2, i, j)*es(2, i, j, 2) + &
                           ve(3, i, j)*es(3, i, j, 2)
        end do
      end do

    end do

    nullify (vlon,vlat,ew,es)
      call timing_on('a2c_halo_buff_routine')
  end subroutine cubed_a2c_halo_buff_only

  subroutine l2c_halo_buff_only(npz, ubuff, vbuff, bd, ilbound, iubound, jlbound, jubound, l2c)
    integer, intent(in) ::  npz, ilbound, jlbound, iubound, jubound
    type(fv_grid_bounds_type), intent(IN) :: bd
    real(kind=R_GRID), intent(in), dimension(2, 2, bd%isd:bd%ied, bd%jsd:bd%jed) :: l2c
    real, intent(inout), dimension(ilbound:iubound, jlbound:jubound, npz) :: ubuff, vbuff
    !local
    integer i, j, k
    real tempu

    ! Compute 3D wind on A grid
!$OMP parallel do default(none) shared(npz,jlbound,jubound,ilbound,iubound,ubuff,vbuff,l2c) &
!$OMP                          private(tempu)
    do k = 1, npz
      do j = jlbound, jubound
        do i = ilbound, iubound
          tempu = l2c(1, 1, i, j)*ubuff(i, j, k) + l2c(1, 2, i, j)*vbuff(i, j, k)
          vbuff(i, j, k) = l2c(2, 1, i, j)*ubuff(i, j, k) + l2c(2, 2, i, j)*vbuff(i, j, k)
          ubuff(i, j, k) = tempu
        end do
      end do
    end do
  end subroutine l2c_halo_buff_only

  subroutine c2l_agrid(u_in, v_in, ull, vll, npz, bd, c2l)
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: npz
    real(kind=R_GRID), intent(in) ::  c2l(2, 2, bd%isd:bd%ied, bd%jsd:bd%jed)
    real, intent(in), dimension(bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: u_in, v_in
    real, intent(out), dimension(bd%isd:bd%ied, bd%jsd:bd%jed, npz) :: ull, vll
    ! Local
    integer i, j, k
    integer :: is, ie, js, je

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    !$OMP parallel do default(none) shared(npz,is,ie,js,je,ull,vll,c2l,u_in,v_in)
    do k = 1, npz
      do j = js, je
        do i = is, ie
          ull(i, j, k) = c2l(1, 1, i, j)*u_in(i, j, k) + c2l(1, 2, i, j)*v_in(i, j, k)
          vll(i, j, k) = c2l(2, 1, i, j)*u_in(i, j, k) + c2l(2, 2, i, j)*v_in(i, j, k)
        end do
      end do
    end do

  end subroutine c2l_agrid

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

!############# from fv_grid_utils to avoid circular modules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine normalize_vect(e)
!                              Make e an unit vector
    real(kind=R_GRID), intent(inout):: e(3)
    real(f_p):: pdot
    integer k

    pdot = e(1)**2 + e(2)**2 + e(3)**2
    pdot = sqrt(pdot)

    do k = 1, 3
      e(k) = e(k)/pdot
    end do

  end subroutine normalize_vect

  subroutine mid_pt_sphere(p1, p2, pm)
    real(kind=R_GRID), intent(IN)  :: p1(2), p2(2)
    real(kind=R_GRID), intent(OUT) :: pm(2)
!------------------------------------------
    real(kind=R_GRID) e1(3), e2(3), e3(3)

    call latlon2xyz(p1, e1)
    call latlon2xyz(p2, e2)
    call mid_pt3_cart(e1, e2, e3)
    call cart_to_latlon(1, e3, pm(1), pm(2))

  end subroutine mid_pt_sphere

  subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
    integer, intent(in):: np
    real(kind=R_GRID), intent(inout):: q(3, np)
    real(kind=R_GRID), intent(inout):: xs(np), ys(np)
! local
    real(kind=R_GRID), parameter:: esl = 1.d-10
    real(f_p):: p(3)
    real(f_p):: dist, lat, lon
    integer i, k

    do i = 1, np
      do k = 1, 3
        p(k) = q(k, i)
      end do
      dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
      do k = 1, 3
        p(k) = p(k)/dist
      end do

      if ((abs(p(1)) + abs(p(2))) < esl) then
        lon = real(0., kind=f_p)
      else
        lon = atan2(p(2), p(1))   ! range [-pi,pi]
      end if

      if (lon < 0.) lon = real(2., kind=f_p)*pi + lon
! RIGHT_HAND system:
      lat = asin(p(3))

      xs(i) = lon
      ys(i) = lat
! q Normalized:
      do k = 1, 3
        q(k, i) = p(k)
      end do
    end do

  end subroutine cart_to_latlon

  subroutine mid_pt3_cart(p1, p2, e)
    real(kind=R_GRID), intent(IN)  :: p1(3), p2(3)
    real(kind=R_GRID), intent(OUT) :: e(3)
!
    real(f_p):: q1(3), q2(3)
    real(f_p):: dd, e1, e2, e3
    integer k

    do k = 1, 3
      q1(k) = p1(k)
      q2(k) = p2(k)
    end do

    e1 = q1(1) + q2(1)
    e2 = q1(2) + q2(2)
    e3 = q1(3) + q2(3)

    dd = sqrt(e1**2 + e2**2 + e3**2)
    e1 = e1/dd
    e2 = e2/dd
    e3 = e3/dd

    e(1) = e1
    e(2) = e2
    e(3) = e3

  end subroutine mid_pt3_cart

  subroutine mid_pt_cart(p1, p2, e3)
    real(kind=R_GRID), intent(IN)  :: p1(2), p2(2)
    real(kind=R_GRID), intent(OUT) :: e3(3)
!-------------------------------------
    real(kind=R_GRID) e1(3), e2(3)

    call latlon2xyz(p1, e1)
    call latlon2xyz(p2, e2)
    call mid_pt3_cart(e1, e2, e3)

  end subroutine mid_pt_cart

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

  subroutine latlon2xyz2(lon, lat, p3)
    real(kind=R_GRID), intent(in):: lon, lat
    real(kind=R_GRID), intent(out):: p3(3)
    real(kind=R_GRID) e(2)

    e(1) = lon; e(2) = lat
    call latlon2xyz(e, p3)

  end subroutine latlon2xyz2

  subroutine latlon2xyz(p, e, id)
!
! Routine to map (lon, lat) to (x,y,z)
!
    real(kind=R_GRID), intent(in) :: p(2)
    real(kind=R_GRID), intent(out):: e(3)
    integer, optional, intent(in):: id   ! id=0 do nothing; id=1, right_hand

    integer n
    real(f_p):: q(2)
    real(f_p):: e1, e2, e3

    do n = 1, 2
      q(n) = p(n)
    end do

    e1 = cos(q(2))*cos(q(1))
    e2 = cos(q(2))*sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

  end subroutine latlon2xyz

  real function great_circle_dist(q1, q2, radius)
    real(kind=R_GRID), intent(IN)           :: q1(2), q2(2)
    real(kind=R_GRID), intent(IN), optional :: radius

    real(f_p):: p1(2), p2(2)
    real(f_p):: beta
    integer n

    do n = 1, 2
      p1(n) = q1(n)
      p2(n) = q2(n)
    end do

    beta = asin(sqrt(sin((p1(2) - p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))* &
                     sin((p1(1) - p2(1))/2.)**2))*2.

    if (present(radius)) then
      great_circle_dist = radius*beta
    else
      great_circle_dist = beta   ! Returns the angle
    end if

  end function great_circle_dist

  subroutine c2l_ord2(u, v, ua, va, gridstruct, km, grid_type, bd, do_halo)
    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(in) :: km, grid_type
    real, intent(in) ::  u(bd%isd:bd%ied, bd%jsd:bd%jed + 1, km)
    real, intent(in) ::  v(bd%isd:bd%ied + 1, bd%jsd:bd%jed, km)
    type(fv_grid_type), intent(IN), target :: gridstruct
    logical, intent(in) :: do_halo
!
    real, intent(out):: ua(bd%isd:bd%ied, bd%jsd:bd%jed, km)
    real, intent(out):: va(bd%isd:bd%ied, bd%jsd:bd%jed, km)
!--------------------------------------------------------------
! Local
    real wu(bd%is - 1:bd%ie + 1, bd%js - 1:bd%je + 2)
    real wv(bd%is - 1:bd%ie + 2, bd%js - 1:bd%je + 1)
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
        do j = js, je + 1
          do i = is, ie
            wu(i, j) = u(i, j, k)*dx(i, j)
          end do
        end do
        do j = js, je
          do i = is, ie + 1
            wv(i, j) = v(i, j, k)*dy(i, j)
          end do
        end do

        do j = js, je
          do i = is, ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
            u1(i) = 2.*(wu(i, j) + wu(i, j + 1))/(dx(i, j) + dx(i, j + 1))
            v1(i) = 2.*(wv(i, j) + wv(i + 1, j))/(dy(i, j) + dy(i + 1, j))
!!!          u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
!!!          v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
            ua(i, j, k) = a11(i, j)*u1(i) + a12(i, j)*v1(i)
            va(i, j, k) = a21(i, j)*u1(i) + a22(i, j)*v1(i)
          end do
        end do
      else
! 2nd order:
        do j = js, je
          do i = is, ie
            ua(i, j, k) = 0.5*(u(i, j, k) + u(i, j + 1, k))
            va(i, j, k) = 0.5*(v(i, j, k) + v(i + 1, j, k))
          end do
        end do
      end if
    end do

  end subroutine c2l_ord2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    allocate (dg%a_kik(2, isd:ied, jsd:jed))
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
    !allocate (dg%xp(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4))
    !allocate (dg%xm(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4))
    !allocate (dg%yp(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4))
    !allocate (dg%ym(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4))

    allocate (dg%xp3(1:interporder + 1, bd%ie - interporder:bd%ied + 1, bd%jsd:bd%jed + 1, 4))
    allocate (dg%xm3(1:interporder + 1, bd%isd:bd%is + interporder, bd%jsd:bd%jed + 1, 4))
    allocate (dg%yp3(1:interporder + 1, bd%isd:bd%ied + 1, bd%je - interporder:bd%jed + 1, 4))
    allocate (dg%ym3(1:interporder + 1, bd%isd:bd%ied + 1, bd%jsd:bd%js + interporder, 4))

    allocate (dg%xp(1:interporder + 1, dg%bd%ie - interporder:dg%bd%ied + 1, dg%bd%jsd:dg%bd%jed + 1, 4))
    allocate (dg%xm(1:interporder + 1, dg%bd%isd:dg%bd%is + interporder, dg%bd%jsd:dg%bd%jed + 1, 4))
    allocate (dg%yp(1:interporder + 1, dg%bd%isd:dg%bd%ied + 1, dg%bd%je - interporder:dg%bd%jed + 1, 4))
    allocate (dg%ym(1:interporder + 1, dg%bd%isd:dg%bd%ied + 1, dg%bd%jsd:dg%bd%js + interporder, 4))
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
    deallocate (dg%a_kik)
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

end module duogrid_mod

