!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @defgroup mosaic2_mod mosaic2_mod
!> @ingroup mosaic2
!> @brief Implements some utility routines to read mosaic information.
!!
!> @author Zhi Liang
!!
!! Implements some utility routines to read mosaic information.
!! The information includes number of tiles and contacts in the mosaic,
!! mosaic grid resolution of each tile, mosaic contact information, mosaic exchange
!! grid information. Each routine will call a C-version routine to get these information.

!> @file
!> @brief File for @ref mosaic2_mod

!> @addtogroup mosaic2_mod
!> @{
module mosaic2_mod

!use fms_mod,    only : write_version_number

use mpp_mod,    only : mpp_error, FATAL, mpp_pe, mpp_root_pe
use mpp_domains_mod, only : domain2D, mpp_get_current_ntile, mpp_get_tile_id
use constants_mod, only : PI, RADIUS
use fms2_io_mod,   only : FmsNetcdfFile_t, open_file, close_file, get_dimension_size
use fms2_io_mod,   only : read_data, variable_exists
use platform_mod

implicit none
private

character(len=*), parameter :: &
     grid_dir  = 'INPUT/'      !> root directory for all grid files

integer, parameter :: &
     MAX_NAME = 256,  & !> max length of the variable names
     MAX_FILE = 1024, & !> max length of the file names
     X_REFINE = 2,    & !> supergrid size/model grid size in x-direction
     Y_REFINE = 2       !> supergrid size/model grid size in y-direction

! --- public interface


public :: get_mosaic_ntiles
public :: get_mosaic_ncontacts
public :: get_mosaic_grid_sizes
public :: get_mosaic_contact
public :: get_mosaic_xgrid_size
public :: get_mosaic_xgrid
public :: get_mosaic_tile_grid
public :: calc_mosaic_grid_area
public :: calc_mosaic_grid_great_circle_area
public :: is_inside_polygon

logical :: module_is_initialized = .true.
! Include variable "version" to be written to log file.
#include<file_version.h>

contains

!#######################################################################

!> @brief Initialize the mosaic_mod.
!> Initialization routine for the mosaic module. It writes the
!! version information to the log file.
!! <br>Example usage:
!!     call mosaic_init ( )
subroutine mosaic_init()

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

!--------- write version number and namelist ------------------
!  call write_version_number("MOSAIC_MOD", version)

end subroutine mosaic_init

!#######################################################################
!> @brief return exchange grid size of mosaic xgrid file.
!!
!> <br>Example usage:
!!              nxgrid = get_mosaic_xgrid_size(xgrid_file)
  function get_mosaic_xgrid_size(fileobj)
    type(FmsNetcdfFile_t), intent(in) :: fileobj !> File that contains exchange grid information.
    integer                           :: get_mosaic_xgrid_size

    call get_dimension_size(fileobj, "ncells", get_mosaic_xgrid_size)

    return

  end function get_mosaic_xgrid_size
!#######################################################################
!> @brief Get exchange grid information from mosaic xgrid file.
!> Example usage:
!!
!!           call get_mosaic_xgrid(fileobj, nxgrid, i1, j1, i2, j2, area)
!!
  subroutine get_mosaic_xgrid(fileobj, i1, j1, i2, j2, area, ibegin, iend)
    type(FmsNetcdfFile_t), intent(in) :: fileobj !> The file that contains exchange grid information.
    integer,       intent(inout) :: i1(:), j1(:), i2(:), j2(:) !> i and j indices for grids 1 and 2
    class(*),          intent(inout) :: area(:) !> area of the exchange grid. The area is scaled to
                                            !! represent unit earth area
    integer, optional, intent(in) :: ibegin, iend

    integer                            :: start(4), nread(4), istart
    real,    dimension(2, size(i1(:))) :: tile1_cell, tile2_cell
    integer                            :: nxgrid, n
    real                               :: garea
    real                               :: get_global_area

    garea = get_global_area()

    ! When start and nread present, make sure nread(1) is the same as the size of the data
    if(present(ibegin) .and. present(iend)) then
       istart = ibegin
       nxgrid = iend - ibegin + 1
       if(nxgrid .NE. size(i1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i1(:))")
       if(nxgrid .NE. size(j1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j1(:))")
       if(nxgrid .NE. size(i2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i2(:))")
       if(nxgrid .NE. size(j2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j2(:))")
       if(nxgrid .NE. size(area(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(area(:))")
    else
       istart = 1
       nxgrid = size(i1(:))
    endif

    start  = 1; nread = 1
    start(1) = istart; nread(1) = nxgrid

    select type(area)
    type is (real(r4_kind))
      call read_data(fileobj, 'xgrid_area', area, corner=start, edge_lengths=nread)
    type is (real(r8_kind))
      call read_data(fileobj, 'xgrid_area', area, corner=start, edge_lengths=nread)
    class default
      call mpp_error(FATAL,"get_mosaic_xgrid: invalid data type for area, must be real(r4_kind) or real(r8_kind)")
    end select

    start = 1; nread = 1
    nread(1) = 2
    start(2) = istart; nread(2) = nxgrid

    select type(area)
    type is (real(r4_kind))
      call read_data(fileobj, 'tile1_cell', tile1_cell, corner=start, edge_lengths=nread)
      call read_data(fileobj, 'tile2_cell', tile2_cell, corner=start, edge_lengths=nread)
    type is (real(r8_kind))
      call read_data(fileobj, 'tile1_cell', tile1_cell, corner=start, edge_lengths=nread)
      call read_data(fileobj, 'tile2_cell', tile2_cell, corner=start, edge_lengths=nread)
    end select

     do n = 1, nxgrid
       i1(n) = int(tile1_cell(1,n))
       j1(n) = int(tile1_cell(2,n))
       i2(n) = int(tile2_cell(1,n))
       j2(n) = int(tile2_cell(2,n))
       select type(area)
       type is (real(r4_kind))
         area(n) = real(area(n)/garea, r4_kind)
       type is (real(r8_kind))
         area(n) = real(area(n)/garea, r8_kind)
       end select
    end do

    return

  end subroutine get_mosaic_xgrid

  !###############################################################################
  !> Get number of tiles in the mosaic_file.
  !> @param fileobj mosaic file object
  !> @return Number of tiles in given file
  !> <br>Example usage:
  !!     ntiles = get_mosaic_ntiles( mosaic_file)
  !!
  function get_mosaic_ntiles(fileobj)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    integer                           :: get_mosaic_ntiles

    call get_dimension_size(fileobj, "ntiles", get_mosaic_ntiles)

    return

  end function get_mosaic_ntiles

  !###############################################################################
  !> Get number of contacts in the mosaic_file.
  !> @param fileobj mosaic file object
  !> @return number of contacts in a given file
  !> <br>Example usage:
  !!         ntiles = get_mosaic_ncontacts( mosaic_file)
  function get_mosaic_ncontacts(fileobj)
    type(FmsNetcdfFile_t), intent(in) :: fileobj
    integer                           :: get_mosaic_ncontacts

    if(variable_exists(fileobj, "contacts") ) then
      call get_dimension_size(fileobj, "ncontact", get_mosaic_ncontacts)
    else
      get_mosaic_ncontacts = 0
    endif

    return

  end function get_mosaic_ncontacts


  !###############################################################################
  !> Get grid size of each tile from mosaic_file
  !> @param fileobj mosaic file object
  !> @param[inout] nx List of grid size in x-direction of each tile
  !> @param[inout] ny List of grid size in y-direction of each tile
  !> <br>Example usage:
  !!            call get_mosaic_grid_sizes(mosaic_file, nx, ny)
  subroutine get_mosaic_grid_sizes( fileobj, nx, ny)
    type(FmsNetcdfFile_t), intent(in)    :: fileobj
    integer, dimension(:), intent(inout) :: nx, ny

    character(len=MAX_FILE) :: gridfile
    integer                 :: ntiles, n
    type(FmsNetcdfFile_t)   :: gridobj

    ntiles = get_mosaic_ntiles(fileobj)
    if(ntiles .NE. size(nx(:)) .OR. ntiles .NE. size(ny(:)) ) then
      call mpp_error(FATAL, "get_mosaic_grid_sizes: size of nx/ny does not equal to ntiles")
    endif

    do n = 1, ntiles
      call read_data(fileobj, 'gridfiles', gridfile, corner=n)
      gridfile = grid_dir//trim(gridfile)

      if(.not. open_file(gridobj, gridfile, 'read')) then
         call mpp_error(FATAL, 'mosaic_mod(get_mosaic_grid_sizes):Error in opening file '//trim(gridfile))
      endif
      call get_dimension_size(gridobj, "nx", nx(n))
      call get_dimension_size(gridobj, "ny", ny(n))
      call close_file(gridobj)
      if(mod(nx(n),x_refine) .NE. 0) call mpp_error(FATAL, "get_mosaic_grid_sizes: nx is not divided by x_refine");
      if(mod(ny(n),y_refine) .NE. 0) call mpp_error(FATAL, "get_mosaic_grid_sizes: ny is not divided by y_refine");
      nx(n) = nx(n)/x_refine;
      ny(n) = ny(n)/y_refine;
    enddo

    return

  end subroutine get_mosaic_grid_sizes

  !###############################################################################
  !> Get contact information from mosaic_file
  !> <br>Example usage:
  !!         call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1,
  !!                                 istart2, iend2, jstart2, jend2)
  !> @param fileobj mosaic file object
  !> @param[inout] tile1 list of tile numbers in tile 1 of each contact
  !> @param[inout] tile2 list of tile numbers in tile 2 of each contact
  !> @param[inout] istart1 list of starting i-index in tile 1 of each contact
  !> @param[inout] iend1 list of ending i-index in tile 1 of each contact
  !> @param[inout] jstart1 list of starting j-index in tile 1 of each contact
  !> @param[inout] jend1 list of ending j-index in tile 1 of each contact
  !> @param[inout] istart2 list of starting i-index in tile 2 of each contact
  !> @param[inout] iend2 list of ending i-index in tile 2 of each contact
  !> @param[inout] jstart2 list of starting j-index in tile 2 of each contact
  !> @param[inout] jend2 list of ending j-index in tile 2 of each contact
  subroutine get_mosaic_contact( fileobj, tile1, tile2, istart1, iend1, jstart1, jend1, &
                                   istart2, iend2, jstart2, jend2)
    type(FmsNetcdfFile_t),    intent(in) :: fileobj
    integer, dimension(:), intent(inout) :: tile1, tile2
    integer, dimension(:), intent(inout) :: istart1, iend1, jstart1, jend1
    integer, dimension(:), intent(inout) :: istart2, iend2, jstart2, jend2
    character(len=MAX_NAME), allocatable :: gridtiles(:)
    character(len=MAX_NAME), allocatable :: contacts(:)
    character(len=MAX_NAME), allocatable :: contacts_index(:)
    character(len=MAX_NAME)              :: strlist(8)
    integer :: ntiles, n, m, ncontacts, nstr, ios
    integer :: i1_type, j1_type, i2_type, j2_type
    logical :: found

    ntiles = get_mosaic_ntiles(fileobj)
    allocate(gridtiles(ntiles))
    call read_data(fileobj, 'gridtiles', gridtiles)

    ncontacts = get_mosaic_ncontacts(fileobj)
    if(ncontacts>0) then
       allocate(contacts(ncontacts), contacts_index(ncontacts))
       call read_data(fileobj, "contacts", contacts)
       call read_data(fileobj, "contact_index", contacts_index)
    endif
    do n = 1, ncontacts
      nstr = parse_string(contacts(n), ":", strlist)
      if(nstr .NE. 4) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): number of elements in contact seperated by :/:: should be 4")
      found = .false.
      do m = 1, ntiles
        if(trim(gridtiles(m)) == trim(strlist(2)) ) then !found the tile name
          found = .true.
          tile1(n) = m
          exit
        endif
      enddo

      if(.not.found) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact):the first tile name specified in contact is not found in tile list")

      found = .false.
      do m = 1, ntiles
        if(trim(gridtiles(m)) == trim(strlist(4)) ) then !found the tile name
          found = .true.
          tile2(n) = m
          exit
        endif
      enddo

      if(.not.found) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact):the second tile name specified in contact is not found in tile list")

      nstr = parse_string(contacts_index(n), ":,", strlist)
      if(nstr .NE. 8) then
        if(mpp_pe()==mpp_root_pe()) then
          print*, "nstr is ", nstr
          print*, "contacts is ", contacts_index(n)
          do m = 1, nstr
            print*, "strlist is ", trim(strlist(m))
          enddo
        endif
        call mpp_error(FATAL, &
               "mosaic_mod(get_mosaic_contact): number of elements in contact_index seperated by :/, should be 8")
      endif
      read(strlist(1), *, iostat=ios) istart1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading istart1")
      read(strlist(2), *, iostat=ios) iend1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading iend1")
      read(strlist(3), *, iostat=ios) jstart1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jstart1")
      read(strlist(4), *, iostat=ios) jend1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jend1")
      read(strlist(5), *, iostat=ios) istart2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading istart2")
      read(strlist(6), *, iostat=ios) iend2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading iend2")
      read(strlist(7), *, iostat=ios) jstart2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jstart2")
      read(strlist(8), *, iostat=ios) jend2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jend2")

      i1_type = transfer_to_model_index(istart1(n), iend1(n), x_refine)
      j1_type = transfer_to_model_index(jstart1(n), jend1(n), y_refine)
      i2_type = transfer_to_model_index(istart2(n), iend2(n), x_refine)
      j2_type = transfer_to_model_index(jstart2(n), jend2(n), y_refine)

      if( i1_type == 0 .AND. j1_type == 0 ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): istart1==iend1 and jstart1==jend1")
      if( i2_type == 0 .AND. j2_type == 0 ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): istart2==iend2 and jstart2==jend2")
      if( i1_type + j1_type .NE. i2_type + j2_type ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): It is not a line or overlap contact")

   enddo

   deallocate(gridtiles)
   if(ncontacts>0) deallocate(contacts, contacts_index)

  end subroutine get_mosaic_contact


function transfer_to_model_index(istart, iend, refine_ratio)
   integer, intent(inout) :: istart, iend
   integer                :: refine_ratio
   integer                :: transfer_to_model_index
   integer                :: istart_in, iend_in

   istart_in = istart
   iend_in = iend

   if( istart_in == iend_in ) then
      transfer_to_model_index = 0
      istart = (istart_in + 1)/refine_ratio
      iend   = istart
   else
      transfer_to_model_index = 1
      if( iend_in > istart_in ) then
        istart = istart_in + 1
        iend   = iend_in
      else
        istart = istart_in
        iend   = iend_in + 1
      endif
      if( mod(istart, refine_ratio) .NE. 0 .OR. mod(iend,refine_ratio) .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(transfer_to_model_index): mismatch between refine_ratio and istart/iend")
      istart = istart/refine_ratio
      iend = iend/refine_ratio

   endif

   return

end function transfer_to_model_index

  !###############################################################################
  !> @brief Calculate grid cell area.
  !> Calculate the grid cell area. The purpose of this routine is to make
  !! sure the consistency between model grid area and exchange grid area.
  !> @param lon geographical longitude of grid cell vertices.
  !> @param lat geographical latitude of grid cell vertices.
  !> @param[inout] area grid cell area.
  !> <br>Example usage:
  !!            call calc_mosaic_grid_area(lon, lat, area)
  subroutine calc_mosaic_grid_area(lon, lat, area)
     class(*), dimension(:,:), intent(in)    :: lon
     class(*), dimension(:,:), intent(in)    :: lat
     class(*), dimension(:,:), intent(inout) :: area
     integer                             :: nlon, nlat
     logical                                 :: valid_types = .false.

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")
     select type (lon)
     type is (real(r4_kind))
       select type (lat)
       type is (real(r4_kind))
         select type(area)
         type is (real(r4_kind))
           call get_grid_area( nlon, nlat, real(lon, r8_kind), real(lat, r8_kind), real(area, r8_kind))
           valid_types = .true.
         end select
       end select
     type is (real(r8_kind))
       select type (lat)
       type is (real(r8_kind))
         select type(area)
         type is (real(r8_kind))
           call get_grid_area( nlon, nlat, lon, lat, area)
           valid_types = .true.
         end select
       end select
     end select

     if(.not. valid_types) call mpp_error(FATAL, "calc_mosaic_grid_area: invalid types given." &
                                                 //" Arguments must be all r4_kind or r8_kind")

  end subroutine calc_mosaic_grid_area

  !###############################################################################
  !> @brief Calculate grid cell area using great cirlce algorithm
  !> Calculate the grid cell area. The purpose of this routine is to make
  !! sure the consistency between model grid area and exchange grid area.
  !> @param lon geographical longitude of grid cell vertices.
  !> @param lat geographical latitude of grid cell vertices.
  !> @param[inout] area grid cell area.
  !> <br>Example usage:
  !!            call calc_mosaic_grid_great_circle_area(lon, lat, area)
  subroutine calc_mosaic_grid_great_circle_area(lon, lat, area)
     class(*), dimension(:,:), intent(in)    :: lon
     class(*), dimension(:,:), intent(in)    :: lat
     class(*), dimension(:,:), intent(inout) :: area
     integer                             :: nlon, nlat
     logical                             :: valid_types = .false.

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     select type (lon)
     type is (real(r4_kind))
       select type (lat)
       type is (real(r4_kind))
         select type(area)
         type is (real(r4_kind))
           call get_grid_great_circle_area( nlon, nlat, real(lon, r8_kind), real(lat, r8_kind), real(area, r8_kind))
           valid_types = .true.
         end select
       end select
     type is (real(r8_kind))
       select type (lat)
       type is (real(r8_kind))
         select type(area)
         type is (real(r8_kind))
           call get_grid_great_circle_area( nlon, nlat, lon, lat, area)
           valid_types = .true.
         end select
       end select
     end select

     if(.not. valid_types) call mpp_error(FATAL, "calc_mosaic_grid_area: invalid types given." &
                                                 //" Arguments must be all r4_kind or r8_kind")

  end subroutine calc_mosaic_grid_great_circle_area

  !#####################################################################
  !> This function check if a point (lon1,lat1) is inside a polygon (lon2(:), lat2(:))
  !! lon1, lat1, lon2, lat2 are in radians.
  function is_inside_polygon(lon1, lat1, lon2, lat2 )
     real, intent(in) :: lon1, lat1
     real, intent(in) :: lon2(:), lat2(:)
     logical          :: is_inside_polygon
     integer                        :: npts, isinside
     integer                        :: inside_a_polygon

     npts = size(lon2(:))

     isinside = inside_a_polygon(lon1, lat1, npts, lon2, lat2)
     if(isinside == 1) then
        is_inside_polygon = .TRUE.
     else
        is_inside_polygon = .FALSE.
     endif

     return

  end function is_inside_polygon

  function parse_string(string, set, value)
  character(len=*),  intent(in) :: string
  character(len=*),  intent(in) :: set
  character(len=*), intent(out) :: value(:)
  integer                       :: parse_string
  integer :: nelem, length, first, last

     nelem = size(value(:))
     length = len_trim(string)

     first = 1; last = 0
     parse_string = 0

     do while(first .LE. length)
       parse_string = parse_string + 1
       if(parse_string>nelem) then
         call mpp_error(FATAL, "mosaic_mod(parse_string) : number of element is greater than size(value(:))")
       endif
       last = first - 1 + scan(string(first:length), set)
       if(last == first-1 ) then  ! not found, end of string
         value(parse_string) = string(first:length)
         exit
       else
         if(last <= first) then
           call mpp_error(FATAL, "mosaic_mod(parse_string) : last <= first")
         endif
         value(parse_string) = string(first:(last-1))
         first = last + 1
         ! scan to make sure the next is not the character in the set
         do while (first == last+1)
           last = first - 1 + scan(string(first:length), set)
           if(last == first) then
             first = first+1
           else
             exit
           endif
         end do
       endif
     enddo

     return

  end function parse_string

  !#############################################################################
  !> Gets the name of a mosaic tile grid file
  !> @param[out] grid_file name of grid file
  !> @param fileobj mosaic file object
  !> @param domain current domain
  !> @param tile_count optional count of tiles
  subroutine get_mosaic_tile_grid(grid_file, fileobj, domain, tile_count)
    character(len=*), intent(out)          :: grid_file
    type(FmsNetcdfFile_t), intent(in)      :: fileobj
    type(domain2D),   intent(in)           :: domain
    integer,          intent(in), optional :: tile_count
    integer                                :: tile, ntileMe
    integer, dimension(:), allocatable     :: tile_id
    character(len=256), allocatable        :: filelist(:)
    integer :: ntiles

    ntiles = get_mosaic_ntiles(fileobj)
    allocate(filelist(ntiles))
    tile = 1
    if(present(tile_count)) tile = tile_count
    ntileMe = mpp_get_current_ntile(domain)
    allocate(tile_id(ntileMe))
    tile_id = mpp_get_tile_id(domain)
    call read_data(fileobj, "gridfiles", filelist)
    grid_file = 'INPUT/'//trim(filelist(tile_id(tile)))
    deallocate(tile_id, filelist)

  end subroutine get_mosaic_tile_grid


end module mosaic2_mod
!> @}
! close documentation grouping
