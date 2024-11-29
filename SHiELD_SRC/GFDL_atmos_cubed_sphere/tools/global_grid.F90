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

!-------------------------------------------------------------------------------
!> @author Xi.Chen <xic@princeton.edu>
!> @date 01/10/2020
!
!  REVISION HISTORY:
!  12/01/2022 - Joseph.Mouallem @FV3 GFDL/Princeton
!  Combine all routines, 4halo updates, stag k2e coef
!-------------------------------------------------------------------------------

module global_grid_mod

  use platform_mod, only: r8_kind
  use constants_mod, only: pi => pi_8, GRAV, OMEGA
  use mpp_mod, only: mpp_pe, mpp_root_pe
  !------ AC modules
  implicit none
  private

  public :: global_grid_type
  public :: global_grid_init, global_grid_end

#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter :: R_GRID = r8_kind
  real(kind=R_GRID), parameter :: RADIUS = 6.3712d+6 !< Radius of the Earth [m]

  real, parameter :: M12vm = 7./12.
  real, parameter :: M22vm = -1./12.

  real, parameter :: M13vm = 37./60.
  real, parameter :: M23vm = -2./15.
  real, parameter :: M33vm = 1./60.

  type global_grid_type

    ! grid data

    !--- dim parameters
    integer :: grid_type = 0
    integer :: res = -999, ng = -999

    real, dimension(2, 6) :: west_pole = -1., north_pole = -1., tile_center = -1.

    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: pt_ext, pt_kik

    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_x, kik_x
    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_y, kik_y

    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_dx, kik_dx
    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_dy, kik_dy
    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_da, kik_da

    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: ext_e1co, kik_e1co
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: ext_e2co, kik_e2co
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: ext_elon, kik_elon
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: ext_elat, kik_elat

    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_sina, kik_sina
    real(kind=R_GRID), dimension(:, :, :), allocatable :: ext_cosa, kik_cosa

    real(kind=R_GRID), dimension(:, :, :, :, :), allocatable :: &
      ext_g_co, ext_g_ctr, &
      kik_g_co, kik_g_ctr, &
      ext_ct2ort_x, ext_ort2ct_x, &
      ext_ct2ort_y, ext_ort2ct_y, &
      ext_ct2rll, ext_rll2ct, &
      ext_c2l, ext_l2c, & ! co
      kik_ct2rll, kik_rll2ct, &
      kik_c2l, kik_l2c

    !--- k2e rmp parameter
    integer :: k2e_nord = 4
    integer, dimension(:, :, :), allocatable :: k2e_loc
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef
    integer, dimension(:, :, :), allocatable :: k2e_loc_b
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef_b
    integer, dimension(:, :, :), allocatable :: k2e_loc_c_x
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef_c_x
    integer, dimension(:, :, :), allocatable :: k2e_loc_c_y
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef_c_y
    integer, dimension(:, :, :), allocatable :: k2e_loc_d_x
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef_d_x
    integer, dimension(:, :, :), allocatable :: k2e_loc_d_y
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: k2e_coef_d_y

    !--- type status
    logical :: is_initialized = .false.

  end type global_grid_type

contains

!!!BASE!!!
!##########
  subroutine global_grid_init(gg, res, ng, grid_type, do_schmidt, target_lon, target_lat)
    type(global_grid_type), intent(inout) :: gg
    integer, intent(in) :: res, ng, grid_type
    real(kind=R_GRID), intent(inout) :: target_lon, target_lat
    logical, intent(in) :: do_schmidt

!        target_lon = target_lon * pi/180.
!        target_lat = target_lat * pi/180.

    !--- check sgatus
    if (gg%is_initialized) return

    !--- set dim parameters and grid_type
    gg%res = res
    gg%ng = ng
    gg%grid_type = grid_type

    !--- alloc
    call global_grid_alloc(gg)

    !--- init lonlat values
    call global_grid_gen_lonlat_equal_edge(gg, do_schmidt, target_lon, target_lat)

    !--- init remap coords
    call global_grid_gen_coords(gg)

    !--- xic: do stretching here, need to update tropical belt rmp coords

    !--- init dx, dy
    call global_grid_gen_ds(gg)

    !--- init da
    call global_grid_gen_da(gg)

    !--- init vec
    call global_grid_gen_eco(gg)
    call global_grid_gen_elonlat(gg)
    call global_grid_gen_mat(gg)

    !--- init k2e rmp
    call global_grid_gen_k2e(gg)

    !--- set status
    gg%is_initialized = .true.

  end subroutine global_grid_init
  !===========================================================================
  !> @brief end method for global_grid_type
  subroutine global_grid_end(gg)
    type(global_grid_type), intent(inout) :: gg

    !--- check status
    if (.not. gg%is_initialized) return

    !--- dealloc
    call global_grid_dealloc(gg)

    !--- set status
    gg%is_initialized = .false.

  end subroutine global_grid_end

!ALLOC!
!########

  subroutine global_grid_alloc(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: isd, ied, jsd, jed

    !--- assign dims
    isd = 1 - 2*gg%ng
    ied = 2*gg%res + 1 + 2*gg%ng
    jsd = isd
    jed = ied

    !--- allocate arrays
    allocate (gg%pt_ext(2, isd:ied, jsd:jed, 6))
    allocate (gg%pt_kik(2, isd:ied, jsd:jed, 6))

    allocate (gg%ext_x(isd:ied, jsd:jed, 6))
    allocate (gg%kik_x(isd:ied, jsd:jed, 6))

    allocate (gg%ext_y(isd:ied, jsd:jed, 6))
    allocate (gg%kik_y(isd:ied, jsd:jed, 6))

    allocate (gg%ext_dx(isd:ied, jsd:jed, 6))
    allocate (gg%kik_dx(isd:ied, jsd:jed, 6))

    allocate (gg%ext_dy(isd:ied, jsd:jed, 6))
    allocate (gg%kik_dy(isd:ied, jsd:jed, 6))

    allocate (gg%ext_da(isd:ied, jsd:jed, 6))
    allocate (gg%kik_da(isd:ied, jsd:jed, 6))

    allocate (gg%ext_e1co(3, isd:ied, jsd:jed, 6))
    allocate (gg%kik_e1co(3, isd:ied, jsd:jed, 6))

    allocate (gg%ext_e2co(3, isd:ied, jsd:jed, 6))
    allocate (gg%kik_e2co(3, isd:ied, jsd:jed, 6))

    allocate (gg%ext_elon(3, isd:ied, jsd:jed, 6))
    allocate (gg%kik_elon(3, isd:ied, jsd:jed, 6))

    allocate (gg%ext_elat(3, isd:ied, jsd:jed, 6))
    allocate (gg%kik_elat(3, isd:ied, jsd:jed, 6))

    allocate (gg%ext_sina(isd:ied, jsd:jed, 6))
    allocate (gg%kik_sina(isd:ied, jsd:jed, 6))

    allocate (gg%ext_cosa(isd:ied, jsd:jed, 6))
    allocate (gg%kik_cosa(isd:ied, jsd:jed, 6))

    allocate (gg%ext_g_co(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_g_ctr(2, 2, isd:ied, jsd:jed, 6))

    allocate (gg%kik_g_co(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%kik_g_ctr(2, 2, isd:ied, jsd:jed, 6))

    allocate (gg%ext_ct2ort_x(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_ort2ct_x(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_ct2ort_y(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_ort2ct_y(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_ct2rll(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_rll2ct(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%kik_ct2rll(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%kik_rll2ct(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_c2l(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%ext_l2c(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%kik_c2l(2, 2, isd:ied, jsd:jed, 6))
    allocate (gg%kik_l2c(2, 2, isd:ied, jsd:jed, 6))

    !--- reassign a-pt dims for k2e
    isd = 1 - gg%ng
    ied = gg%res + gg%ng
    jsd = isd
    jed = ied

    allocate (gg%k2e_loc(isd:ied, jsd:jed, 6))
    allocate (gg%k2e_coef(gg%k2e_nord, isd:ied, jsd:jed, 6))
    allocate (gg%k2e_loc_b(isd:ied + 1, jsd:jed + 1, 6))
    allocate (gg%k2e_coef_b(gg%k2e_nord, isd:ied + 1, jsd:jed + 1, 6))
    allocate (gg%k2e_loc_c_x(isd:ied + 1, jsd:jed, 6))
    allocate (gg%k2e_coef_c_x(gg%k2e_nord, isd:ied + 1, jsd:jed, 6))
    allocate (gg%k2e_loc_c_y(isd:ied, jsd:jed + 1, 6))
    allocate (gg%k2e_coef_c_y(gg%k2e_nord, isd:ied, jsd:jed + 1, 6))
    allocate (gg%k2e_loc_d_x(isd:ied + 1, jsd:jed, 6))
    allocate (gg%k2e_coef_d_x(gg%k2e_nord, isd:ied + 1, jsd:jed, 6))
    allocate (gg%k2e_loc_d_y(isd:ied, jsd:jed + 1, 6))
    allocate (gg%k2e_coef_d_y(gg%k2e_nord, isd:ied, jsd:jed + 1, 6))

  end subroutine global_grid_alloc
  !===========================================================================
  subroutine global_grid_dealloc(gg)
    type(global_grid_type), intent(inout) :: gg

    !--- deallocate arrays
    deallocate (gg%k2e_loc)
    deallocate (gg%k2e_coef)
    deallocate (gg%k2e_loc_b)
    deallocate (gg%k2e_coef_b)
    deallocate (gg%k2e_loc_c_x)
    deallocate (gg%k2e_coef_c_x)
    deallocate (gg%k2e_loc_c_y)
    deallocate (gg%k2e_coef_c_y)
    deallocate (gg%k2e_loc_d_x)
    deallocate (gg%k2e_coef_d_x)
    deallocate (gg%k2e_loc_d_y)
    deallocate (gg%k2e_coef_d_y)

    deallocate (gg%ext_ct2ort_x)
    deallocate (gg%ext_ort2ct_x)
    deallocate (gg%ext_ct2ort_y)
    deallocate (gg%ext_ort2ct_y)
    deallocate (gg%ext_ct2rll)
    deallocate (gg%ext_rll2ct)
    deallocate (gg%kik_ct2rll)
    deallocate (gg%kik_rll2ct)
    deallocate (gg%ext_c2l)
    deallocate (gg%ext_l2c)
    deallocate (gg%kik_c2l)
    deallocate (gg%kik_l2c)

    deallocate (gg%ext_g_co)
    deallocate (gg%ext_g_ctr)
    deallocate (gg%kik_g_co)
    deallocate (gg%kik_g_ctr)

    deallocate (gg%ext_sina)
    deallocate (gg%kik_sina)

    deallocate (gg%ext_cosa)
    deallocate (gg%kik_cosa)

    deallocate (gg%ext_elon)
    deallocate (gg%kik_elon)

    deallocate (gg%ext_elat)
    deallocate (gg%kik_elat)

    deallocate (gg%ext_e1co)
    deallocate (gg%kik_e1co)

    deallocate (gg%ext_e2co)
    deallocate (gg%kik_e2co)

    deallocate (gg%ext_dx)
    deallocate (gg%kik_dx)

    deallocate (gg%ext_dy)
    deallocate (gg%kik_dy)

    deallocate (gg%ext_da)
    deallocate (gg%kik_da)

    deallocate (gg%ext_x)
    deallocate (gg%kik_x)

    deallocate (gg%ext_y)
    deallocate (gg%kik_y)

    deallocate (gg%pt_ext)
    deallocate (gg%pt_kik)

  end subroutine global_grid_dealloc

!GEN_CELL
  subroutine global_grid_gen_ds(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: ii, jj
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(R_GRID), dimension(:, :, :), allocatable :: pt
    real(R_GRID), dimension(:, :), allocatable :: dy

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- allocate arrays
    allocate (pt(2, isd:ied, jsd:jed))
    allocate (dy(isd:nc, jsd:nc - 1)) !--- intentional dims to reveal errors

    !--- assign tile 1 ext
    pt(:, :, :) = gg%pt_ext(:, :, :, 1)

    !--- gen ds ext
    do j = jsd, nc - 1
      do i = isd, nc
        dy(i, j) = (pt(2, i, j + 1) - pt(2, i, j))*RADIUS
      end do
    end do

    !--- populate to gg (both ext and kik, replace kik halo ds later)
    do n = 1, 6
      do j = jsd, nc - 1
        do i = isd, nc
          !--- dy
          ii = ie - i + 1
          jj = je - j
          gg%ext_dy(i, j, n) = dy(i, j) ! sw
          gg%ext_dy(ii, j, n) = dy(i, j) ! se
          gg%ext_dy(i, jj, n) = dy(i, j) ! nw
          gg%ext_dy(ii, jj, n) = dy(i, j) ! ne
          gg%kik_dy(i, j, n) = dy(i, j) ! sw
          gg%kik_dy(ii, j, n) = dy(i, j) ! se
          gg%kik_dy(i, jj, n) = dy(i, j) ! nw
          gg%kik_dy(ii, jj, n) = dy(i, j) ! ne
          !--- dx
          gg%ext_dx(j, i, n) = dy(i, j) ! sw
          gg%ext_dx(j, ii, n) = dy(i, j) ! nw
          gg%ext_dx(jj, i, n) = dy(i, j) ! se
          gg%ext_dx(jj, ii, n) = dy(i, j) ! ne
          gg%kik_dx(j, i, n) = dy(i, j) ! sw
          gg%kik_dx(j, ii, n) = dy(i, j) ! nw
          gg%kik_dx(jj, i, n) = dy(i, j) ! se
          gg%kik_dx(jj, ii, n) = dy(i, j) ! ne
        end do
      end do
    end do

    !--- assign tile1 kik for halo
    pt(:, :, :) = gg%pt_kik(:, :, :, 1)

    !--- gen ds kik
    do j = js, nc - 1
      do i = isd, is - 1
        dy(i, j) = (pt(2, i, j + 1) - pt(2, i, j))*RADIUS
      end do
    end do

    !--- populate to gg kik halo ds
    do n = 1, 6
      do j = js, nc - 1
        do i = isd, is - 1
          !--- dy
          ii = ie - i + 1
          jj = je - j
          gg%kik_dy(i, j, n) = dy(i, j) ! sw
          gg%kik_dy(ii, j, n) = dy(i, j) ! se
          gg%kik_dy(i, jj, n) = dy(i, j) ! nw
          gg%kik_dy(ii, jj, n) = dy(i, j) ! ne
          !--- dx
          gg%kik_dx(j, i, n) = dy(i, j) ! sw
          gg%kik_dx(j, ii, n) = dy(i, j) ! nw
          gg%kik_dx(jj, i, n) = dy(i, j) ! se
          gg%kik_dx(jj, ii, n) = dy(i, j) ! ne
        end do
      end do
    end do

    !--- deallocate arrays
    deallocate (pt)
    deallocate (dy)

  end subroutine global_grid_gen_ds
  !===========================================================================
  !> @brief generate da
  subroutine global_grid_gen_da(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: ii, jj
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(R_GRID), dimension(:, :, :), allocatable :: pt
    real(R_GRID), dimension(:, :), allocatable :: da

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- allocate arrays
    allocate (pt(2, isd:ied, jsd:jed))
    allocate (da(isd:nc - 1, jsd:nc - 1)) !--- intentional dims to reveal errors

    !--- assign tile 1 ext
    pt(:, :, :) = gg%pt_ext(:, :, :, 1)

    !--- gen da ext
    do j = jsd, nc - 1
      do i = isd, nc - 1
        da(i, j) = lib_4pt_area( &
                   pt(:, i, j), pt(:, i + 1, j), pt(:, i + 1, j + 1), pt(:, i, j + 1), &
                   RADIUS)
      end do
    end do

    !--- populate to gg (both ext and kik, replace kik halo ds later)
    do n = 1, 6
      do j = jsd, nc - 1
        do i = isd, nc - 1
          ii = ie - i
          jj = je - j
          gg%ext_da(i, j, n) = da(i, j) ! sw
          gg%ext_da(ii, j, n) = da(i, j) ! se
          gg%ext_da(i, jj, n) = da(i, j) ! nw
          gg%ext_da(ii, jj, n) = da(i, j) ! ne
          gg%kik_da(i, j, n) = da(i, j) ! sw
          gg%kik_da(ii, j, n) = da(i, j) ! se
          gg%kik_da(i, jj, n) = da(i, j) ! nw
          gg%kik_da(ii, jj, n) = da(i, j) ! ne
        end do
      end do
    end do

    !--- assign tile1 kik for halo
    pt(:, :, :) = gg%pt_kik(:, :, :, 1)

    !--- gen ds ext
    do j = js, nc - 1
      do i = isd, is - 1
        da(i, j) = lib_4pt_area( &
                   pt(:, i, j), pt(:, i + 1, j), pt(:, i + 1, j + 1), pt(:, i, j + 1), &
                   RADIUS)
      end do
    end do

    !--- populate to gg (both ext and kik, replace kik halo ds later)
    do n = 1, 6
      do j = js, nc - 1
        do i = isd, is - 1
          !--- y-dir
          ii = ie - i
          jj = je - j
          gg%kik_da(i, j, n) = da(i, j) ! sw
          gg%kik_da(ii, j, n) = da(i, j) ! se
          gg%kik_da(i, jj, n) = da(i, j) ! nw
          gg%kik_da(ii, jj, n) = da(i, j) ! ne
          !--- x-dir
          gg%kik_da(j, i, n) = da(i, j) ! sw
          gg%kik_da(j, ii, n) = da(i, j) ! nw
          gg%kik_da(jj, i, n) = da(i, j) ! se
          gg%kik_da(jj, ii, n) = da(i, j) ! ne
        end do
      end do
    end do

    !--- deallocate arrays
    deallocate (pt)
    deallocate (da)

  end subroutine global_grid_gen_da

!GENCOORD
!#######################
  subroutine global_grid_gen_coords(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(kind=R_GRID), dimension(:), allocatable :: line
    real(kind=R_GRID), dimension(:, :), allocatable :: pt

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied
!!!!!!
! should we create global extended kik_x/y for staggered fields or this ok?
!!!!!!
    nc = 1 + gg%res

    !--- allocate arrays
    allocate (line(isd:ied))
    allocate (pt(2, isd:ied))

    line(:) = -999.

    !do k = 1, ng
    do k = 0, ng

      !--- assign one layer of rmp kik pts
      pt(:, :) = gg%pt_kik(:, is - k, :, 1)

      !--- get coords
      line(nc) = 0.
      !do j = js, nc-1
      do j = js - 1, nc
!                line(j) = - lib_great_circ_dist( &
!                    pt(:,j), pt(:,nc))
        line(j) = pt(2, j) - pt(2, nc) !--- take advantage of meridians
        line(je - j + 1) = -line(j)
      end do

      !--- populate the coords to cubed sphere
      do n = 1, 6
        !do j = js,je
        do j = js - 1, je + 1
          !--- west and south
          i = is - k
          gg%kik_y(i, j, n) = line(j)
          gg%kik_x(j, i, n) = line(j)
          !--- east and north
          i = ie + k
          gg%kik_y(i, j, n) = line(j)
          gg%kik_x(j, i, n) = line(j)
        end do
      end do

      !--- assign one layer of rmp ext pts
      pt(:, :) = gg%pt_ext(:, is - k, :, 1)

      !--- get coords
      line(nc) = 0.
      !do j = jsd, nc-1
      do j = jsd, nc
!                line(j) = - lib_great_circ_dist( &
!                    pt(:,j), pt(:,nc))
        line(j) = pt(2, j) - pt(2, nc) !--- take advantage of meridians
        line(je - j + 1) = -line(j)
      end do

      !--- populate the coords to cubed sphere
      do n = 1, 6
        do j = jsd, jed
          !--- west and south
          i = is - k
          gg%ext_y(i, j, n) = line(j)
          gg%ext_x(j, i, n) = line(j)
          !--- east and north
          i = ie + k
          gg%ext_y(i, j, n) = line(j)
          gg%ext_x(j, i, n) = line(j)
        end do
      end do

    end do

    !--- deallocate arrays
    deallocate (line)
    deallocate (pt)

  end subroutine global_grid_gen_coords

!genlonlat

  subroutine global_grid_gen_lonlat_equal_edge(gg, do_schmidt, target_lon, target_lat)
    type(global_grid_type), intent(inout) :: gg
    real(kind=R_GRID), intent(in) :: target_lon, target_lat
    logical, intent(in) :: do_schmidt

    !--- local
    integer :: i, j, ii, jj, nc, n, ng
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(kind=R_GRID) :: rsq3, rsq2, alpha, dela
    real(kind=R_GRID) :: x, y, angl_x, angl_y, angl_g, angl, sina, cosa
    real(kind=R_GRID), dimension(:, :, :, :), allocatable :: cart
    real(kind=R_GRID), dimension(:), allocatable :: line
    real(kind=R_GRID), parameter :: shift_fac = -18.
    real(kind=R_GRID), dimension(3, 3) :: rot_z
    real(kind=R_GRID) :: c = 1.            ! Stretching factor

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- allocate arrays
    allocate (cart(3, isd:ied, jsd:jed, 6))
    allocate (line(isd:ied))

    !------ xic: equal_angular part
    if (gg%grid_type == 2) then
      !--- get the angular step
      alpha = pi/4.0
      dela = 2.d0*alpha/real(ie - is, kind=R_GRID)

      !--- form half line: computational domain
      line(:) = -999.
      line(js) = -1.
      line(nc) = 0.
      do j = js + 1, nc - 1
        angl_y = (j - 1)*dela - alpha
        y = tan(angl_y)
        line(j) = y
      end do
    end if
    !------ xic: end equal_edge part

    !------ xic: equal_edge part
    if (gg%grid_type == 0) then
      !--- get the angular step
      rsq3 = 1.d0/sqrt(3.d0)
      rsq2 = 1.d0/sqrt(2.d0)
      alpha = asin(rsq3)
      dela = 2.d0*alpha/real(ie - is, kind=R_GRID)

      !--- form half line: computational domain
      line(:) = -999.
      line(js) = -1.
      line(nc) = 0.
      do j = js + 1, nc - 1
        angl_y = (j - 1)*dela - alpha
        y = tan(angl_y)*sqrt(2.d0)
        line(j) = y
      end do
    end if
    !------ xic: end equal_edge part

!recheck the ng logic here for low res:
!for a C8, isd=1-2*5=-9
!j -9<-->0
!jj 2<-->11 which is > 8
!so line goes overboard

    !--- form half line: ghost cells
    do j = jsd, js - 1
      jj = 2*js - j
      angl_g = -0.5*pi - atan(line(jj))
      line(j) = tan(angl_g)
    end do

    !--- form full line
    do j = jsd, nc - 1
      jj = je - j + 1
      line(jj) = -line(j)
    end do

    !--- populate to 6 tiles cart
    do j = jsd, jed
      do i = isd, ied
        x = line(i)
        y = line(j)
        cart(:, i, j, 1) = [1.d0, x, y]
        cart(:, i, j, 2) = [-x, 1.d0, y]
        cart(:, i, j, 3) = [-x, -y, 1.d0]
        cart(:, i, j, 4) = [-1.d0, -y, -x]
        cart(:, i, j, 5) = [y, -1.d0, -x]
        cart(:, i, j, 6) = [y, x, -1.d0]
      end do
    end do

    !--- Shift the corner away from Japan
!        angl = pi/shift_fac
!        sina = sin(angl)
!        cosa = cos(angl)
!        rot_z(:,:) = 0.
!        rot_z(1,1) = cosa
!        rot_z(2,2) = cosa
!        rot_z(3,3) = 1.
!        rot_z(1,2) = -sina
!        rot_z(2,1) = sina
!
!        do n = 1,6
!            do j = jsd, jed
!                do i = isd, ied
!
!                    cart(:,i,j,n) = matmul(rot_z, cart(:,i,j,n))
!
!                enddo
!            enddo
!        enddo

    !--- populate to 6 tiles lonlat
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied
          gg%pt_ext(:, i, j, n) = lib_cart2lonlat(cart(:, i, j, n))
        end do
      end do
    end do

!        call direct_transform(1, isd, ied, jsd, jed, &
!            Atm%flagstruct%target_lon, Atm%flagstruct%target_lat, &
!            n, gg%pt_ext(1,:,:,n),gg%pt_ext(2,:,:,n))
    if (do_schmidt) then
      do n = 1, 6
        call direct_transform(c, isd, ied, jsd, jed, &
                              target_lon, target_lat, &
                              n, gg%pt_ext(1, :, :, n), gg%pt_ext(2, :, :, n))
        ! call direct_transform(1., isd, ied, jsd, jed, &
        !     3.14159265358979  ,     -1.57079632679490, &
        !     n, gg%pt_ext(1,:,:,n),gg%pt_ext(2,:,:,n))
      end do
    end if

    call global_grid_update_kik(gg%pt_ext, gg%pt_kik, isd, ied, ng, 2)

    !--- deallocate arrays
    deallocate (cart)
    deallocate (line)

  end subroutine global_grid_gen_lonlat_equal_edge
  !===========================================================================

  subroutine direct_transform(c, i1, i2, j1, j2, lon_p, lat_p, n, lon, lat)
!
! This is a direct transformation of the standard (symmetrical) cubic grid
! to a locally enhanced high-res grid on the sphere; it is an application
! of the Schmidt transformation at the south pole followed by a
! pole_shift_to_target (rotation) operation
!
    real(kind=R_GRID), intent(in):: c              ! Stretching factor
    real(kind=R_GRID), intent(in):: lon_p, lat_p   ! center location of the target face, radian
    integer, intent(in):: n              ! grid face number
    integer, intent(in):: i1, i2, j1, j2
!  0 <= lon <= 2*pi ;    -pi/2 <= lat <= pi/2
    real(kind=R_GRID), intent(inout), dimension(i1:i2, j1:j2):: lon, lat
!
    real(f_p):: lat_t, sin_p, cos_p, sin_lat, cos_lat, sin_o, p2, two_pi
    real(f_p):: c2p1, c2m1
    integer:: i, j

    p2 = 0.5d0*pi
    two_pi = 2.d0*pi

    c2p1 = 1.d0 + c*c
    c2m1 = 1.d0 - c*c

    sin_p = sin(lat_p)
    cos_p = cos(lat_p)

    do j = j1, j2
      do i = i1, i2
        if (abs(c2m1) > 1.d-7) then
          sin_lat = sin(lat(i, j))
          lat_t = asin((c2m1 + c2p1*sin_lat)/(c2p1 + c2m1*sin_lat))
        else         ! no stretching
          lat_t = lat(i, j)
        end if
        sin_lat = sin(lat_t)
        cos_lat = cos(lat_t)
        sin_o = -(sin_p*sin_lat + cos_p*cos_lat*cos(lon(i, j)))
        if ((1.-abs(sin_o)) < 1.d-7) then    ! poles
          lon(i, j) = 0.d0
          lat(i, j) = sign(p2, sin_o)
        else
          lat(i, j) = asin(sin_o)
          lon(i, j) = lon_p + atan2(-cos_lat*sin(lon(i, j)), &
                                    -sin_lat*cos_p + cos_lat*sin_p*cos(lon(i, j)))
          if (lon(i, j) < 0.d0) then
            lon(i, j) = lon(i, j) + two_pi
          elseif (lon(i, j) >= two_pi) then
            lon(i, j) = lon(i, j) - two_pi
          end if
        end if
      end do
    end do

  end subroutine direct_transform

!GENVEC

  subroutine global_grid_gen_eco(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: ii, jj
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(R_GRID), dimension(2) :: p1, p2

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- get ext eco
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          !--- e1co
          if (i == ied) then
            p1 = gg%pt_ext(:, i, j, n)
            p2 = gg%pt_ext(:, i - 1, j, n)
            gg%ext_e1co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_ext(:, i, j, n)
            p2 = gg%pt_ext(:, i + 1, j, n)
            gg%ext_e1co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if

          !--- e2co
          if (j == jed) then
            p1 = gg%pt_ext(:, i, j, n)
            p2 = gg%pt_ext(:, i, j - 1, n)
            gg%ext_e2co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_ext(:, i, j, n)
            p2 = gg%pt_ext(:, i, j + 1, n)
            gg%ext_e2co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if

        end do
      end do
    end do

    !--- get kik eco (a little more work)
    !------ default null
    gg%kik_e1co(:, :, :, :) = -999.
    gg%kik_e2co(:, :, :, :) = -999.

    !------ interior
    do n = 1, 6
      do j = js, je
        do i = is, ie
          gg%kik_e1co(:, i, j, n) = gg%ext_e1co(:, i, j, n)
          gg%kik_e2co(:, i, j, n) = gg%ext_e2co(:, i, j, n)
        end do
      end do
    end do

    !------ e1co
    do n = 1, 6

      do j = js, je
        ! west
        do i = isd, is
          if (i == is) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i - 1, j, n)
            gg%kik_e1co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i + 1, j, n)
            gg%kik_e1co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
        ! east --- repeat the code to enhance readability
        do i = ie, ied
          if (i == ied) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i - 1, j, n)
            gg%kik_e1co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i + 1, j, n)
            gg%kik_e1co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do

      end do

      !--- south --- tile edge equals to ext
      do j = jsd, js - 1
        do i = is, ie
          if (i == ie) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i - 1, j, n)
            gg%kik_e1co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i + 1, j, n)
            gg%kik_e1co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
      end do
      !--- north
      do j = je + 1, jed
        do i = is, ie
          if (i == ie) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i - 1, j, n)
            gg%kik_e1co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i + 1, j, n)
            gg%kik_e1co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
      end do

    end do ! n=1,6

    !------ e2co
    do n = 1, 6

      do j = js, je
        ! west --- tile edge equals to ext
        do i = isd, is - 1
          if (j == je) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j - 1, n)
            gg%kik_e2co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j + 1, n)
            gg%kik_e2co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
        ! east --- tile edge equals to ext
        do i = ie + 1, ied
          if (j == je) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j - 1, n)
            gg%kik_e2co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j + 1, n)
            gg%kik_e2co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
      end do

      ! south
      do j = jsd, js
        do i = is, ie
          if (j == js) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j - 1, n)
            gg%kik_e2co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j + 1, n)
            gg%kik_e2co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
      end do

      ! north
      do j = je, jed
        do i = is, ie
          if (j == jed) then
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j - 1, n)
            gg%kik_e2co(:, i, j, n) = -lib_2pt_unit_vec(p1, p2)
          else
            p1 = gg%pt_kik(:, i, j, n)
            p2 = gg%pt_kik(:, i, j + 1, n)
            gg%kik_e2co(:, i, j, n) = lib_2pt_unit_vec(p1, p2)
          end if
        end do
      end do

    end do ! n=1,6

    !--- nullify the tile corners (shared by different halos)
    do n = 1, 6
      gg%kik_e1co(:, is, js, n) = -999.
      gg%kik_e1co(:, ie, js, n) = -999.
      gg%kik_e1co(:, ie, je, n) = -999.
      gg%kik_e1co(:, is, je, n) = -999.
      gg%kik_e2co(:, is, js, n) = -999.
      gg%kik_e2co(:, ie, js, n) = -999.
      gg%kik_e2co(:, ie, je, n) = -999.
      gg%kik_e2co(:, is, je, n) = -999.
    end do

  end subroutine global_grid_gen_eco
  !===========================================================================
  !> @brief generate regular elon and elat on the cubed-sphere
  subroutine global_grid_gen_elonlat(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: ii, jj
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(R_GRID), dimension(2) :: p1, p2
    real(R_GRID) :: lon, lat

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- get ext elonlat
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          lon = gg%pt_ext(1, i, j, n)
          lat = gg%pt_ext(2, i, j, n)
          gg%ext_elon(1, i, j, n) = -sin(lon)
          gg%ext_elon(2, i, j, n) = cos(lon); 
          gg%ext_elon(3, i, j, n) = 0.
          gg%ext_elat(1, i, j, n) = -sin(lat)*cos(lon); 
          gg%ext_elat(2, i, j, n) = -sin(lat)*sin(lon); 
          gg%ext_elat(3, i, j, n) = cos(lat)

        end do
      end do
    end do

    !--- get kik elonlat
    !------ default null
    gg%kik_elon(:, :, :, :) = -999.
    gg%kik_elat(:, :, :, :) = -999.

    call global_grid_update_kik(gg%ext_elon, gg%kik_elon, isd, ied, ng, 3)
    call global_grid_update_kik(gg%ext_elat, gg%kik_elat, isd, ied, ng, 3)

  end subroutine global_grid_gen_elonlat
  !===========================================================================
  !> @brief generate sina cosa on the cubed-sphere
  subroutine global_grid_gen_mat(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, nc, n, ng
    integer :: ii, jj
    integer :: isd, ied, jsd, jed, is, ie, js, je
    real(R_GRID), dimension(2) :: p1, p2
    real(R_GRID) :: lon, lat, sina, cosa, coef
    logical, dimension(:, :), allocatable :: assigned
    real(R_GRID) :: gct11, gco12, gco22
    real(R_GRID) :: gct22, gco11, gco21
    real(R_GRID), dimension(3) :: e1, e2, elon, elat
    real(R_GRID), dimension(2, 2) :: mat

    !--- assign dims
    is = 1
    ie = 2*gg%res + 1
    js = is
    je = ie

    ng = 2*gg%ng

    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    nc = 1 + gg%res

    !--- allocate working array, and make corners true
    allocate (assigned(isd:ied, jsd:jed))
    assigned(:, :) = .true.
    assigned(is:ie, jsd:jed) = .false.
    assigned(isd:ied, js:je) = .false.
    assigned(is, js) = .true.
    assigned(is, je) = .true.
    assigned(ie, je) = .true.
    assigned(ie, js) = .true.

    !--- get ext sina and cosa
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          gg%ext_cosa(i, j, n) = dot_product( &
                                 gg%ext_e1co(:, i, j, n), gg%ext_e2co(:, i, j, n))
          gg%ext_sina(i, j, n) = sqrt(1.- &
                                      gg%ext_cosa(i, j, n)*gg%ext_cosa(i, j, n))

        end do
      end do
    end do

    !--- get kik elonlat
    !------ default cosa = 0, sina = 1
    gg%kik_cosa(:, :, :) = 0.
    gg%kik_sina(:, :, :) = 1.

    !------ interior
    do n = 1, 6
      do j = js + 1, je - 1
        do i = is + 1, ie - 1

          gg%kik_cosa(i, j, n) = gg%ext_cosa(i, j, n)
          gg%kik_sina(i, j, n) = gg%ext_sina(i, j, n)
          assigned(i, j) = .true.

        end do
      end do
    end do

    !------ halo
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          if (.not. assigned(i, j)) then

            gg%kik_cosa(i, j, n) = dot_product( &
                                   gg%kik_e1co(:, i, j, n), gg%kik_e2co(:, i, j, n))
            gg%kik_sina(i, j, n) = sqrt(1.- &
                                        gg%kik_cosa(i, j, n)*gg%kik_cosa(i, j, n))

            if (n == 6) then
              assigned(i, j) = .true.
            end if

          end if

        end do
      end do
    end do

    !--- 2x2 matrix init
    assigned(:, :) = .false.

    do j = js, je
      do i = isd, ied
        assigned(i, j) = .true.
      end do
    end do

    do j = jsd, jed
      do i = is, ie
        assigned(i, j) = .true.
      end do
    end do

    do j = js + 1, je - 1
      do i = is + 1, ie - 1
        assigned(i, j) = .false.
      end do
    end do

    assigned(is, js) = .false.
    assigned(is, je) = .false.
    assigned(ie, je) = .false.
    assigned(ie, js) = .false.

    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          sina = gg%ext_sina(i, j, n)
          cosa = gg%ext_cosa(i, j, n)

          ! ext g_co g_ctr
          gg%ext_g_co(1, 1, i, j, n) = 1.
          gg%ext_g_co(1, 2, i, j, n) = cosa
          gg%ext_g_co(2, 1, i, j, n) = cosa
          gg%ext_g_co(2, 2, i, j, n) = 1.

          coef = 1./sina/sina
          gg%ext_g_ctr(1, 1, i, j, n) = coef
          gg%ext_g_ctr(1, 2, i, j, n) = -coef*cosa
          gg%ext_g_ctr(2, 1, i, j, n) = -coef*cosa
          gg%ext_g_ctr(2, 2, i, j, n) = coef

          ! ext ct2ort_x
          gct11 = gg%ext_g_ctr(1, 1, i, j, n)
          gco12 = gg%ext_g_co(1, 2, i, j, n)
          gco22 = gg%ext_g_co(2, 2, i, j, n)
          gct22 = gg%ext_g_ctr(2, 2, i, j, n)
          gco11 = gg%ext_g_co(1, 1, i, j, n)
          gco21 = gg%ext_g_co(2, 1, i, j, n)

          gg%ext_ct2ort_x(1, 1, i, j, n) = 1./sqrt(gct11)
          gg%ext_ct2ort_x(1, 2, i, j, n) = 0.
          gg%ext_ct2ort_x(2, 1, i, j, n) = gco12/sqrt(gco22)
          gg%ext_ct2ort_x(2, 2, i, j, n) = sqrt(gco22)

          ! ext ort2ct_x
          gg%ext_ort2ct_x(1, 1, i, j, n) = sqrt(gct11)
          gg%ext_ort2ct_x(1, 2, i, j, n) = 0.
          gg%ext_ort2ct_x(2, 1, i, j, n) = -gco12/gco22*sqrt(gct11)
          gg%ext_ort2ct_x(2, 2, i, j, n) = 1./sqrt(gco22)

          ! ext ct2ort_y
          gg%ext_ct2ort_y(1, 1, i, j, n) = 0.
          gg%ext_ct2ort_y(1, 2, i, j, n) = 1./sqrt(gct22)
          gg%ext_ct2ort_y(2, 1, i, j, n) = sqrt(gco11)
          gg%ext_ct2ort_y(2, 2, i, j, n) = gco21/sqrt(gco11)

          ! ext ort2ct_y
          gg%ext_ort2ct_y(1, 1, i, j, n) = -gco21/gco11*sqrt(gct22)
          gg%ext_ort2ct_y(1, 2, i, j, n) = 1./sqrt(gco11)
          gg%ext_ort2ct_y(2, 1, i, j, n) = sqrt(gct22)
          gg%ext_ort2ct_y(2, 2, i, j, n) = 0.

          ! ext c2l
          e1 = gg%ext_e1co(:, i, j, n)
          e2 = gg%ext_e2co(:, i, j, n)
          elon = gg%ext_elon(:, i, j, n)
          elat = gg%ext_elat(:, i, j, n)
          mat(1, 1) = dot_product(e1, elon)
          mat(1, 2) = dot_product(e2, elon)
          mat(2, 1) = dot_product(e1, elat)
          mat(2, 2) = dot_product(e2, elat)

          gg%ext_ct2rll(:, :, i, j, n) = mat(:, :)

          coef = 1./(mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1))

          gg%ext_rll2ct(1, 1, i, j, n) = coef*mat(2, 2)
          gg%ext_rll2ct(1, 2, i, j, n) = -coef*mat(1, 2)
          gg%ext_rll2ct(2, 1, i, j, n) = -coef*mat(2, 1)
          gg%ext_rll2ct(2, 2, i, j, n) = coef*mat(1, 1)

          gg%ext_c2l(:, :, i, j, n) = matmul( &
                                      gg%ext_ct2rll(:, :, i, j, n), gg%ext_g_ctr(:, :, i, j, n))
          gg%ext_l2c(:, :, i, j, n) = matmul( &
                                      gg%ext_g_co(:, :, i, j, n), gg%ext_rll2ct(:, :, i, j, n))

        end do
      end do
    end do

    ! kik ct2rll
    assigned(:, :) = .false.

    do j = js, je
      do i = isd, ied
        assigned(i, j) = .true.
      end do
    end do

    do j = jsd, jed
      do i = is, ie
        assigned(i, j) = .true.
      end do
    end do

    do j = js, je
      do i = is, ie
        assigned(i, j) = .false.
      end do
    end do

    gg%kik_ct2rll(:, :, :, :, :) = -999.
    gg%kik_rll2ct(:, :, :, :, :) = -999.
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          if (assigned(i, j)) then
            e1 = gg%kik_e1co(:, i, j, n)
            e2 = gg%kik_e2co(:, i, j, n)
            elon = gg%kik_elon(:, i, j, n)
            elat = gg%kik_elat(:, i, j, n)
            mat(1, 1) = dot_product(e1, elon)
            mat(1, 2) = dot_product(e2, elon)
            mat(2, 1) = dot_product(e1, elat)
            mat(2, 2) = dot_product(e2, elat)

            gg%kik_ct2rll(:, :, i, j, n) = mat(:, :)

            coef = 1./(mat(1, 1)*mat(2, 2) - mat(1, 2)*mat(2, 1))

            gg%kik_rll2ct(1, 1, i, j, n) = coef*mat(2, 2)
            gg%kik_rll2ct(1, 2, i, j, n) = -coef*mat(1, 2)
            gg%kik_rll2ct(2, 1, i, j, n) = -coef*mat(2, 1)
            gg%kik_rll2ct(2, 2, i, j, n) = coef*mat(1, 1)

          end if

        end do
      end do
    end do

    assigned(:, :) = .false.

    do j = js, je
      do i = isd, ied
        assigned(i, j) = .true.
      end do
    end do

    do j = jsd, jed
      do i = is, ie
        assigned(i, j) = .true.
      end do
    end do

    do j = js + 1, je - 1
      do i = is + 1, ie - 1
        assigned(i, j) = .false.
      end do
    end do

    assigned(is, js) = .false.
    assigned(is, je) = .false.
    assigned(ie, je) = .false.
    assigned(ie, js) = .false.

    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied

          if (assigned(i, j)) then

            sina = gg%kik_sina(i, j, n)
            cosa = gg%kik_cosa(i, j, n)

            ! ext g_co g_ctr
            gg%kik_g_co(1, 1, i, j, n) = 1.
            gg%kik_g_co(1, 2, i, j, n) = cosa
            gg%kik_g_co(2, 1, i, j, n) = cosa
            gg%kik_g_co(2, 2, i, j, n) = 1.

            coef = 1./sina/sina
            gg%kik_g_ctr(1, 1, i, j, n) = coef
            gg%kik_g_ctr(1, 2, i, j, n) = -coef*cosa
            gg%kik_g_ctr(2, 1, i, j, n) = -coef*cosa
            gg%kik_g_ctr(2, 2, i, j, n) = coef

            gg%kik_c2l(:, :, i, j, n) = matmul( &
                                        gg%kik_ct2rll(:, :, i, j, n), gg%kik_g_ctr(:, :, i, j, n))
            gg%kik_l2c(:, :, i, j, n) = matmul( &
                                        gg%kik_g_co(:, :, i, j, n), gg%kik_rll2ct(:, :, i, j, n))

          end if

        end do
      end do
    end do

    !--- deallocate working array
    deallocate (assigned)

  end subroutine global_grid_gen_mat
  !=====

!UTIL
!##########

  subroutine global_grid_update_kik(pt_ext, pt_kik, isd, ied, ng, np)
    integer, intent(in) :: isd, ied, ng, np
    real(kind=R_GRID), dimension(np, isd:ied, isd:ied, 6), intent(in) :: pt_ext
    real(kind=R_GRID), dimension(np, isd:ied, isd:ied, 6), intent(out):: pt_kik
    !--- local
    integer :: is, ie, js, je, jsd, jed
    integer :: i, j, k, n, n_src, ii, jj
    integer :: nw, ne, ns, nn, n_neighbor(4)
    integer :: isk, iek, jsk, jek

    !--- assign dims
    jsd = isd
    jed = ied

    is = isd + ng
    js = jsd + ng
    ie = ied - ng
    je = jed - ng

    !--- copy inner values
    pt_kik(:, :, :, :) = -999.
    do n = 1, 6
      do j = js - 2, je + 2
        do i = is - 2, ie + 2

          pt_kik(:, i, j, n) = pt_ext(:, i, j, n)

        end do
      end do
    end do

    do n = 1, 6

      call get_neighbor_tile_num(n, nw, ne, ns, nn)
      n_neighbor = [nw, ne, ns, nn]
      do k = 1, 4
        n_src = n_neighbor(k)

        call get_neighbor_bounds(isk, iek, jsk, jek, n, n_src, ie, je, ng)
        do j = jsk, jek
          do i = isk, iek
            call get_neighbor_index(ii, jj, i, j, n, n_src, ie, je)
            pt_kik(:, i, j, n) = pt_kik(:, ii, jj, n_src)
          end do
        end do

      end do

    end do

  end subroutine global_grid_update_kik
  !=============================================================================
  subroutine get_neighbor_index(ii, jj, i, j, n, n_src, npx, npy)
    integer, intent(out) :: ii, jj
    integer, intent(in) :: i, j, n, n_src, npx, npy
    ! local
    integer :: nw, ne, ns, nn, isc = 0, jsc = 0, iec = 0, jec = 0
    logical :: is_even_tile_num

    is_even_tile_num = mod(n, 2) .eq. 0

    isc = 1; jsc = 1; iec = npx; jec = npy
    call get_neighbor_tile_num(n, nw, ne, ns, nn)

    if (.not. is_even_tile_num) then
      if (n_src .eq. nw) then
        ii = iec - (j - jsc)
        jj = jec + (i - isc)
      end if
      if (n_src .eq. ne) then
        ii = isc + (i - iec)
        jj = jsc + (j - jsc)
      end if
      if (n_src .eq. ns) then
        ii = isc + (i - isc)
        jj = jec + (j - jsc)
      end if
      if (n_src .eq. nn) then
        ii = isc + (j - jec)
        jj = jec - (i - isc)
      end if
    else  ! is_even_tile_num
      if (n_src .eq. nw) then
        ii = iec + (i - isc)
        jj = jsc + (j - jsc)
      end if
      if (n_src .eq. ne) then
        ii = iec - (j - jsc)
        jj = jsc + (i - iec)
      end if
      if (n_src .eq. ns) then
        ii = iec + (j - jsc)
        jj = jec - (i - isc)
      end if
      if (n_src .eq. nn) then
        ii = isc + (i - isc)
        jj = jsc + (j - jec)
      end if
    end if ! is_even_tile_num

  end subroutine get_neighbor_index
  !===========================================================================
  subroutine get_neighbor_bounds(is, ie, js, je, n, n_src, npx, npy, ng)
    ! n refer to the tile to get kik edge assigned
    ! n_src refer to the neighbor tile of the source pt info
    integer, intent(out) :: is, ie, js, je
    integer, intent(in) :: n, n_src, npx, npy, ng
    ! local
    integer :: nw, ne, ns, nn, isc, jsc, iec, jec

    isc = 1; jsc = 1; iec = npx; jec = npy
    call get_neighbor_tile_num(n, nw, ne, ns, nn)

    if (n_src .eq. nw) then
      is = isc - ng; ie = isc - 1; js = jsc; je = jec
    end if
    if (n_src .eq. ne) then
      is = iec + 1; ie = iec + ng; js = jsc; je = jec
    end if
    if (n_src .eq. ns) then
      is = isc; ie = iec; js = jsc - ng; je = jsc - 1
    end if
    if (n_src .eq. nn) then
      is = isc; ie = iec; js = jec + 1; je = jec + ng
    end if

  end subroutine get_neighbor_bounds
  !===========================================================================
  subroutine get_neighbor_tile_num(n, nw, ne, ns, nn)
    ! call get_neighbor_tile_num((in)n,(out)nw,(out)ne,(out)ns,(out)nn)
    ! get the neighbor tile number of the cubed sphere tile of number n
    ! nw => west tile; ne => east tile;
    ! ns => south tile; nn => north tile
    integer, intent(in) :: n
    integer, intent(out) :: nw, ne, ns, nn

    if (mod(n, 2) .eq. 0) then
      nn = mod(n + 0, 6) + 1
      ne = mod(n + 1, 6) + 1
      ns = mod(n + 3, 6) + 1
      nw = mod(n + 4, 6) + 1
    else
      ne = mod(n + 0, 6) + 1
      nn = mod(n + 1, 6) + 1
      nw = mod(n + 3, 6) + 1
      ns = mod(n + 4, 6) + 1
    end if

  end subroutine get_neighbor_tile_num
  !=======

  pure function lib_cart2lonlat(pt_cart) result(pt_lonlat)
    real(kind=R_GRID), dimension(2) :: pt_lonlat
    real(kind=R_GRID), dimension(3), intent(in) :: pt_cart
    ! local
    real(kind=R_GRID) :: p(3)
    real(kind=R_GRID) :: lat, lon
    real(kind=R_GRID), parameter :: esl = 1.e-5

    p(:) = pt_cart(:)/sqrt(pt_cart(1)*pt_cart(1) + &
                           pt_cart(2)*pt_cart(2) + pt_cart(3)*pt_cart(3))

    if ((abs(p(1)) + abs(p(2))) < esl) then
      lon = 0.
    else
      lon = atan2(p(2), p(1))   ! range [-pi,pi]
    end if

    if (lon < 0.) lon = 2.*pi + lon
    lat = asin(p(3))

    pt_lonlat = [lon, lat]

  end function lib_cart2lonlat
  !===========================================================================
  pure function lib_lonlat2cart(pt_lonlat) result(e_cart)
    ! pt_lonlat is dim(2) with lon and lat
    ! e_cart is dim(3) with normalized cartisian unit basis (x, y, z)
    real(kind=R_GRID), dimension(3) :: e_cart
    real(kind=R_GRID), dimension(2), intent(in) :: pt_lonlat
    ! local
    real(kind=R_GRID) :: lat, lon

    lon = pt_lonlat(1)
    lat = pt_lonlat(2)

    e_cart(1) = cos(lat)*cos(lon)
    e_cart(2) = cos(lat)*sin(lon)
    e_cart(3) = sin(lat)

  end function lib_lonlat2cart
  !===========================================================================
  pure function lib_great_circ_dist(lonlat1, lonlat2, radiuss) result(dist)
    real(kind=R_GRID), dimension(2), intent(in) :: lonlat1, lonlat2
    real(kind=R_GRID), intent(in), optional :: radiuss
    real(kind=R_GRID) :: dist
    ! local
    real(kind=R_GRID) :: lat1, lat2, lon1, lon2, dlat, dlon
    real(kind=R_GRID) :: beta, tmp1, tmp2

    lon1 = lonlat1(1); lon2 = lonlat2(1)
    lat1 = lonlat1(2); lat2 = lonlat2(2)
    dlat = lat2 - lat1; dlon = lon2 - lon1

    tmp1 = sqrt((cos(lat2)*sin(dlon))**2 + &
                (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))**2)
    tmp2 = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)

    beta = atan2(tmp1, tmp2)

    if (present(radiuss)) then
      dist = radiuss*beta
    else
      dist = beta
    end if

  end function lib_great_circ_dist
  !=============================================================================
  real function lib_4pt_area(p1_ll, p2_ll, p3_ll, p4_ll, r) result(area)
    real(kind=R_GRID), dimension(2), intent(in):: p1_ll, p2_ll, p3_ll, p4_ll
    real(kind=R_GRID), intent(in) :: r
    ! local
    real(kind=R_GRID), dimension(3) :: p1, p2, p3, p4
    real(kind=R_GRID), dimension(3) :: p12, p23, p34, p41
    real(kind=R_GRID) :: c123, c234, c341, c412
    real(kind=R_GRID) :: a123, a234, a341, a412

    p1 = lib_lonlat2cart(p1_ll)
    p2 = lib_lonlat2cart(p2_ll)
    p3 = lib_lonlat2cart(p3_ll)
    p4 = lib_lonlat2cart(p4_ll)

    p12 = cross_prod(p1, p2)
    p23 = cross_prod(p2, p3)
    p34 = cross_prod(p3, p4)
    p41 = cross_prod(p4, p1)

    p12 = norm_vec_3d(p12)
    p23 = norm_vec_3d(p23)
    p34 = norm_vec_3d(p34)
    p41 = norm_vec_3d(p41)

    c123 = dot_product(p12, p23)
    c234 = dot_product(p23, p34)
    c341 = dot_product(p34, p41)
    c412 = dot_product(p41, p12)

    a123 = acos(-c123)
    a234 = acos(-c234)
    a341 = acos(-c341)
    a412 = acos(-c412)

    area = (a123 + a234 + a341 + a412 - 2*pi)*r**2

  end function lib_4pt_area
  !===========================================================================
  function lib_2pt_unit_vec(p1_ll, p2_ll) result(e)
    real(kind=R_GRID), dimension(3) :: e
    real(kind=R_GRID), dimension(2), intent(in):: p1_ll, p2_ll
    ! local
    real(kind=R_GRID), dimension(3) :: p1, p2
    real(kind=R_GRID), dimension(3) :: p12

    p1 = lib_lonlat2cart(p1_ll)
    p2 = lib_lonlat2cart(p2_ll)

    p12 = cross_prod(p1, p2)
    e = cross_prod(p12, p1)

    e = norm_vec_3d(e)

  end function lib_2pt_unit_vec
  !===========================================================================
  pure function norm_vec_3d(r) result(e)
    ! e = norm_vec_3d(r)
    ! Return a vector that is the same vector as r except the length is 1
    real(kind=R_GRID), dimension(3) :: e
    real(kind=R_GRID), intent(in) :: r(3)
    ! local
    real(kind=R_GRID) dist

    dist = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
    e = r/dist

  end function norm_vec_3d
  !===========================================================================
  !pure function cross_prod_r(a,b) result(rlt)
  pure function cross_prod(a, b) result(rlt)
    real(kind=R_GRID), dimension(3) :: rlt
    real(kind=R_GRID), dimension(3), intent(in) :: a, b

    rlt(1) = a(2)*b(3) - a(3)*b(2)
    rlt(2) = a(3)*b(1) - a(1)*b(3)
    rlt(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_prod
  !end function cross_prod_r
  !---------------------------------------------------------------------------
!    pure function cross_prod_d(a,b) result(rlt)
!        integer, dimension(3) :: rlt
!        integer, dimension(3), intent(in) :: a,b
!
!        rlt(1) = a(2)*b(3) - a(3)*b(2)
!        rlt(2) = a(3)*b(1) - a(1)*b(3)
!        rlt(3) = a(1)*b(2) - a(2)*b(1)
!
!    end function cross_prod_d
  !=============================================================================
  ! Lagrangian polynomial interpolation only make sense when remapping on fixed
  ! locations. Otherwise, use the lib_interp function.
  function lib_interp_lag_get_coef(x, vx) result(coef)
    ! coef = lib_interp_lag_get_coef(x,vx)
    real(kind=R_GRID), dimension(:), intent(in) :: vx
    real(kind=R_GRID), intent(in) :: x
    real(kind=R_GRID), dimension(size(vx)) :: coef
    ! local
    integer :: i, j, n

    n = size(vx)
    do j = 1, n
      coef(j) = 1.
      do i = 1, n
        if (i .ne. j) then
          coef(j) = coef(j)*(x - vx(i))/(vx(j) - vx(i))
        end if
      end do
    end do

  end function lib_interp_lag_get_coef
  !===========================================================================
  pure function lib_interp_lag_eval(coef, vy) result(y)
    ! y = lib_interp_lag_eval(coef,vy)
    real :: y
    real(kind=R_GRID), dimension(:), intent(in) :: coef
    real, dimension(size(coef)), intent(in) :: vy
    ! local
    integer :: i

    y = 0.
    do i = 1, size(vy)
      y = y + vy(i)*coef(i)
    end do

  end function lib_interp_lag_eval

! GEN k2e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine global_grid_gen_k2e(gg)
    type(global_grid_type), intent(inout) :: gg
    !--- local
    integer :: i, j, k, n, ng, k2e_nord
    integer :: ii, jj, np
    integer :: isd, ied, jsd, jed, is, ie, js, je, khi, klo
    real(R_GRID), dimension(:, :, :), allocatable :: &
      ext_x, ext_y, kik_x, kik_y ! A-pt copies
    real(R_GRID), dimension(:), allocatable :: x, y ! x: kik(org), y: ext(tar)
    real(R_GRID) :: yy
    real(R_GRID), dimension(:, :, :), allocatable :: &
      ext_x_b, ext_y_b, kik_x_b, kik_y_b, & ! B-pt copies
      ext_x_c, ext_y_c, kik_x_c, kik_y_c, & ! c-pt copies
      ext_x_d, ext_y_d, kik_x_d, kik_y_d ! d-pt copies
    real(R_GRID), dimension(:), allocatable :: x_b, y_b ! x: kik(org), y: ext(tar)
    real(R_GRID), dimension(:), allocatable :: x_c, y_c, x_c_one ! x: kik(org), y: ext(tar)
    real(R_GRID), dimension(:), allocatable :: x_d, y_d, x_d_one ! x: kik(org), y: ext(tar)

    !--- assign dims
    is = 1
    ie = gg%res
    js = is
    je = ie

    ng = gg%ng - 2

    print *, 'gen_k2e ng', ng
    isd = is - ng
    ied = ie + ng
    jsd = isd
    jed = ied

    k2e_nord = gg%k2e_nord
    np = k2e_nord/2 - 1

    !--- allocate a-pt copies
    allocate (ext_x(isd:ied, jsd:jed, 6))
    allocate (ext_y(isd:ied, jsd:jed, 6))
    allocate (kik_x(isd:ied, jsd:jed, 6))
    allocate (kik_y(isd:ied, jsd:jed, 6))
    allocate (x(is:ie))
    allocate (y(isd:ied))

    !--- assign a-pt copies
    do n = 1, 6
      do j = jsd, jed
        do i = isd, ied
          ! point
          ii = i*2; jj = j*2

          ext_x(i, j, n) = gg%ext_x(ii, jj, n)
          ext_y(i, j, n) = gg%ext_y(ii, jj, n)
          kik_x(i, j, n) = gg%kik_x(ii, jj, n)
          kik_y(i, j, n) = gg%kik_y(ii, jj, n)

        end do
      end do
    end do

    !--- get k2e_loc
    gg%k2e_loc(:, :, :) = -999
    gg%k2e_coef(:, :, :, :) = -999.
    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x
        j = je + ii
        call get_loc_x

        !--- rmp w e
        i = is - ii
        call get_loc_y
        i = ie + ii
        call get_loc_y

      end do
    end do

    !--- deallocate a-pt copies
    deallocate (x)
    deallocate (y)

    deallocate (ext_x)
    deallocate (ext_y)
    deallocate (kik_x)
    deallocate (kik_y)

!############################################################
    !--- allocate b-pt copies
    allocate (ext_x_b(isd:ied + 1, jsd:jed + 1, 6))
    allocate (ext_y_b(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_x_b(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_y_b(isd:ied + 1, jsd:jed + 1, 6))
    allocate (x_b(is:ie + 1))
    allocate (y_b(isd:ied + 1))

    print *, 'genk2e c', is, ie, ied, jed, ng
    !--- assign b-pt copies
    do n = 1, 6
      do j = jsd, jed + 1
        do i = isd, ied + 1
          ! point
          ii = i*2 - 1; jj = j*2 - 1

          ext_x_b(i, j, n) = gg%ext_x(ii, jj, n)
          ext_y_b(i, j, n) = gg%ext_y(ii, jj, n)
          kik_x_b(i, j, n) = gg%kik_x(ii, jj, n)
          kik_y_b(i, j, n) = gg%kik_y(ii, jj, n)
        end do
      end do
    end do

    !--- get k2e_loc_b
    gg%k2e_loc_b(:, :, :) = -999
    gg%k2e_coef_b(:, :, :, :) = -999.
    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x_b
        j = je + 1 + ii
        call get_loc_x_b

        !--- rmp w e
        i = is - ii
        call get_loc_y_b
        i = ie + 1 + ii
        call get_loc_y_b

      end do
    end do

    !--- deallocate b-pt copies
    deallocate (x_b)
    deallocate (y_b)

    deallocate (ext_x_b)
    deallocate (ext_y_b)
    deallocate (kik_x_b)
    deallocate (kik_y_b)

!############################################################
    !--- allocate c-pt copies
    ! +1 in both x and y to create the extfields for u and v
    allocate (ext_x_c(isd:ied + 1, jsd:jed + 1, 6))
    allocate (ext_y_c(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_x_c(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_y_c(isd:ied + 1, jsd:jed + 1, 6))
    allocate (x_c(is:ie))
    allocate (x_c_one(is:ie + 1))
    allocate (y_c(isd:ied))

    print *, 'genk2e c', is, ie, ied, jed, ng
    !--- assign c-pt copies
    do n = 1, 6
      do j = jsd, jed + 1
        do i = isd, ied + 1
          ! point
          ii = i*2 - 1; jj = j*2

          ext_x_c(i, j, n) = gg%ext_x(ii, jj, n)
          ext_y_c(i, j, n) = gg%ext_y(ii, jj, n)
          kik_x_c(i, j, n) = gg%kik_x(ii, jj, n)
          kik_y_c(i, j, n) = gg%kik_y(ii, jj, n)

        end do
      end do
    end do

    !--- get k2e_loc_c
    gg%k2e_loc_c_x(:, :, :) = -999
    gg%k2e_coef_c_x(:, :, :, :) = -999.
    gg%k2e_loc_c_y(:, :, :) = -999
    gg%k2e_coef_c_y(:, :, :, :) = -999.

    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x_c_x
        j = je + ii
        call get_loc_x_c_x

        !--- rmp w e
        i = is - ii
        call get_loc_y_c_x
        i = ie + 1 + ii
        call get_loc_y_c_x

      end do
    end do

    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x_c_y
        j = je + 1 + ii
        call get_loc_x_c_y

        !--- rmp w e
        i = is - ii
        call get_loc_y_c_y
        i = ie + ii
        call get_loc_y_c_y

      end do
    end do

    !--- deallocate c-pt copies
    deallocate (x_c)
    deallocate (x_c_one)
    deallocate (y_c)

    deallocate (ext_x_c)
    deallocate (ext_y_c)
    deallocate (kik_x_c)
    deallocate (kik_y_c)

!############################################################
    !--- allocate D-pt copies
    ! +1 in both x and y to create the extfields for u and v
    allocate (ext_x_d(isd:ied + 1, jsd:jed + 1, 6))
    allocate (ext_y_d(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_x_d(isd:ied + 1, jsd:jed + 1, 6))
    allocate (kik_y_d(isd:ied + 1, jsd:jed + 1, 6))
    allocate (x_d(is:ie))
    allocate (x_d_one(is:ie + 1))
    allocate (y_d(isd:ied))

    !--- assign D-pt copies
    do n = 1, 6
      do j = jsd, jed + 1
        do i = isd, ied + 1
          ! point
          ii = i*2; jj = j*2 - 1

          ext_x_d(i, j, n) = gg%ext_x(ii, jj, n)
          ext_y_d(i, j, n) = gg%ext_y(ii, jj, n)
          kik_x_d(i, j, n) = gg%kik_x(ii, jj, n)
          kik_y_d(i, j, n) = gg%kik_y(ii, jj, n)

        end do
      end do
    end do

    !--- get k2e_loc_d
    gg%k2e_loc_d_x(:, :, :) = -999
    gg%k2e_coef_d_x(:, :, :, :) = -999.
    gg%k2e_loc_d_y(:, :, :) = -999
    gg%k2e_coef_d_y(:, :, :, :) = -999.

    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x_d_x
        j = je + ii
        call get_loc_x_d_x

        !--- rmp w e
        i = is - ii
        call get_loc_y_d_x
        i = ie + 1 + ii
        call get_loc_y_d_x

      end do
    end do

    do n = 1, 6
      do ii = 1, ng

        !--- rmp s n
        j = js - ii
        call get_loc_x_d_y
        j = je + 1 + ii
        call get_loc_x_d_y

        !--- rmp w e
        i = is - ii
        call get_loc_y_d_y
        i = ie + ii
        call get_loc_y_d_y

      end do
    end do

    !--- deallocate d-pt copies
    deallocate (x_d)
    deallocate (x_d_one)
    deallocate (y_d)

    deallocate (ext_x_d)
    deallocate (ext_y_d)
    deallocate (kik_x_d)
    deallocate (kik_y_d)

  contains
    !-----------------------------------------------------------------------
    subroutine get_loc_x
      x(:) = kik_x(is:ie, j, n)
      y(:) = ext_x(isd:ied, j, n)

      do i = is - ii + 1, ie + ii - 1

        yy = y(i)
        klo = is
        khi = ie

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie - np - 1)
        khi = klo + 1

        gg%k2e_loc(i, j, n) = klo

        gg%k2e_coef(:, i, j, n) = lib_interp_lag_get_coef(yy, x(klo - np:khi + np))

      end do

    end subroutine get_loc_x
    !-----------------------------------------------------------------------
    subroutine get_loc_y
      x(:) = kik_y(i, js:je, n)
      y(:) = ext_y(i, jsd:jed, n)

      do j = js - ii + 1, je + ii - 1

        yy = y(j)
        klo = js
        khi = je

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je - np - 1)
        khi = klo + 1

        gg%k2e_loc(i, j, n) = klo

        gg%k2e_coef(:, i, j, n) = lib_interp_lag_get_coef(yy, x(klo - np:khi + np))

      end do

    end subroutine get_loc_y
    !-----------------------------------------------------------------------

    subroutine get_loc_x_b
      x_b(:) = kik_x_b(is:ie + 1, j, n)
      y_b(:) = ext_x_b(isd:ied + 1, j, n)

      do i = is - ii + 1, ie + 1 + ii - 1

        yy = y_b(i)
        klo = is
        khi = ie + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_b(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_b(i, j, n) = klo

        gg%k2e_coef_b(:, i, j, n) = lib_interp_lag_get_coef(yy, x_b(klo - np:khi + np))

      end do

    end subroutine get_loc_x_b
    !-----------------------------------------------------------------------
    subroutine get_loc_y_b
      x_b(:) = kik_y_b(i, js:je + 1, n)
      y_b(:) = ext_y_b(i, jsd:jed + 1, n)

      do j = js - ii + 1, je + 1 + ii - 1

        yy = y_b(j)
        klo = js
        khi = je + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_b(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_b(i, j, n) = klo

        gg%k2e_coef_b(:, i, j, n) = lib_interp_lag_get_coef(yy, x_b(klo - np:khi + np))

      end do

    end subroutine get_loc_y_b

!C grid
    subroutine get_loc_x_c_x
      x_c_one(:) = kik_x_c(is:ie + 1, j, n)
      y_c(:) = ext_x_c(isd:ied, j, n)

      do i = is - ii + 1, ie + 1 + ii - 1

        yy = y_c(i)
        klo = is
        khi = ie + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_c_one(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_c_x(i, j, n) = klo

        gg%k2e_coef_c_x(:, i, j, n) = lib_interp_lag_get_coef(yy, x_c_one(klo - np:khi + np))

      end do

    end subroutine get_loc_x_c_x
    !-----------------------------------------------------------------------
    subroutine get_loc_y_c_x
      x_c(:) = kik_y_c(i, js:je, n)
      y_c(:) = ext_y_c(i, jsd:jed, n)

      do j = js - ii + 1, je + ii - 1

        yy = y_c(j)
        klo = js
        khi = je

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_c(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je - np - 1)
        khi = klo + 1

        gg%k2e_loc_c_x(i, j, n) = klo

        gg%k2e_coef_c_x(:, i, j, n) = lib_interp_lag_get_coef(yy, x_c(klo - np:khi + np))

      end do

    end subroutine get_loc_y_c_x

    subroutine get_loc_x_c_y
      x_c(:) = kik_x_c(is:ie, j, n)
      y_c(:) = ext_x_c(isd:ied, j, n)

      do i = is - ii + 1, ie + ii - 1

        yy = y_c(i)
        klo = is
        khi = ie

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_c(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie - np - 1)
        khi = klo + 1

        gg%k2e_loc_c_y(i, j, n) = klo

        gg%k2e_coef_c_y(:, i, j, n) = lib_interp_lag_get_coef(yy, x_c(klo - np:khi + np))

      end do

    end subroutine get_loc_x_c_y
    !-----------------------------------------------------------------------
    subroutine get_loc_y_c_y
      x_c_one(:) = kik_y_c(i, js:je + 1, n)
      y_c(:) = ext_y_c(i, jsd:jed, n)

      do j = js - ii + 1, je + 1 + ii - 1

        yy = y_c(j)
        klo = js
        khi = je + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_c_one(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_c_y(i, j, n) = klo

        gg%k2e_coef_c_y(:, i, j, n) = lib_interp_lag_get_coef(yy, x_c_one(klo - np:khi + np))

      end do

    end subroutine get_loc_y_c_y

!same as before, change kik_d_c and ext_d_c
    subroutine get_loc_x_d_x
      x_d_one(:) = kik_x_d(is:ie + 1, j, n)
      y_d(:) = ext_x_d(isd:ied, j, n)

      do i = is - ii + 1, ie + 1 + ii - 1

        yy = y_d(i)
        klo = is
        khi = ie + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_d_one(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_d_x(i, j, n) = klo

        gg%k2e_coef_d_x(:, i, j, n) = lib_interp_lag_get_coef(yy, x_d_one(klo - np:khi + np))

      end do

    end subroutine get_loc_x_d_x
    !-----------------------------------------------------------------------
    subroutine get_loc_y_d_x
      x_d(:) = kik_y_d(i, js:je, n)
      y_d(:) = ext_y_d(i, jsd:jed, n)

      do j = js - ii + 1, je + ii - 1

        yy = y_d(j)
        klo = js
        khi = je

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_d(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je - np - 1)
        khi = klo + 1

        gg%k2e_loc_d_x(i, j, n) = klo

        gg%k2e_coef_d_x(:, i, j, n) = lib_interp_lag_get_coef(yy, x_d(klo - np:khi + np))

      end do

    end subroutine get_loc_y_d_x

    subroutine get_loc_x_d_y
      x_d(:) = kik_x_d(is:ie, j, n)
      y_d(:) = ext_x_d(isd:ied, j, n)

      do i = is - ii + 1, ie + ii - 1

        yy = y_d(i)
        klo = is
        khi = ie

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_d(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, is + np)
        klo = min(klo, ie - np - 1)
        khi = klo + 1

        gg%k2e_loc_d_y(i, j, n) = klo

        gg%k2e_coef_d_y(:, i, j, n) = lib_interp_lag_get_coef(yy, x_d(klo - np:khi + np))

      end do

    end subroutine get_loc_x_d_y
    !-----------------------------------------------------------------------
    subroutine get_loc_y_d_y
      x_d_one(:) = kik_y_d(i, js:je + 1, n)
      y_d(:) = ext_y_d(i, jsd:jed, n)

      do j = js - ii + 1, je + 1 + ii - 1

        yy = y_d(j)
        klo = js
        khi = je + 1

        do while (khi - klo > 1)
          k = (khi + klo)/2
          if (x_d_one(k) > yy) then
            khi = k
          else
            klo = k
          end if
        end do
        klo = max(klo, js + np)
        klo = min(klo, je + 1 - np - 1)
        khi = klo + 1

        gg%k2e_loc_d_y(i, j, n) = klo

        gg%k2e_coef_d_y(:, i, j, n) = lib_interp_lag_get_coef(yy, x_d_one(klo - np:khi + np))

      end do

    end subroutine get_loc_y_d_y

  end subroutine global_grid_gen_k2e

end module global_grid_mod

