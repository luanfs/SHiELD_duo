!-------------------------------------------------------------------------------
!> @brief geometry lib for grid functions
!> @author Xi.Chen <xic@princeton.edu>
!> @date 01/17/2020
!
!  REVISION HISTORY:
!  01/17/2020 - Initial Version
!-------------------------------------------------------------------------------

module lib_grid_mod

    !------ fms modules
    use platform_mod,       only: r8_kind
    use constants_mod,      only: pi=>pi_8, GRAV, OMEGA
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    !------ AC modules
    implicit none
    private

    integer, parameter :: R_GRID = r8_kind
    real(kind=R_GRID), parameter :: RADIUS = 6.3712d+6 !< Radius of the Earth [m]

    real, parameter :: M12vm = 7./12.
    real, parameter :: M22vm = -1./12.

    real, parameter :: M13vm = 37./60.
    real, parameter :: M23vm = -2./15.
    real, parameter :: M33vm = 1./60.

    public :: R_GRID
    public :: pi
    public :: RADIUS
    public :: GRAV
    public :: OMEGA

    public :: M12vm, M22vm, M13vm, M23vm, M33vm

    public :: lib_cart2lonlat
    public :: lib_lonlat2cart
    public :: lib_great_circ_dist
    public :: lib_4pt_area
    public :: lib_2pt_unit_vec

    public :: norm_vec_3d
    public :: cross_prod

    public :: lib_interp_lag_get_coef
    public :: lib_interp_lag_eval

    !===========================================================================
    interface cross_prod
      module procedure cross_prod_r, cross_prod_d
    endinterface cross_prod
    !===========================================================================
contains
    !===========================================================================
    pure function lib_cart2lonlat(pt_cart) result(pt_lonlat)
        real(kind=R_GRID), dimension(2) :: pt_lonlat
        real(kind=R_GRID), dimension(3), intent(in) :: pt_cart
        ! local
        real(kind=R_GRID) :: p(3)
        real(kind=R_GRID) :: lat, lon
        real(kind=R_GRID), parameter :: esl = 1.e-5

        p(:) = pt_cart(:)/sqrt(pt_cart(1)*pt_cart(1) + &
            pt_cart(2)*pt_cart(2) + pt_cart(3)*pt_cart(3))

        if ( (abs(p(1))+abs(p(2)))  < esl ) then 
            lon = 0.
        else
            lon = atan2( p(2), p(1) )   ! range [-pi,pi]
        endif

        if ( lon < 0.) lon = 2.*pi + lon
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
        real(kind=R_GRID) :: lat,lon

        lon = pt_lonlat(1)
        lat = pt_lonlat(2)

        e_cart(1) = cos(lat)*cos(lon)
        e_cart(2) = cos(lat)*sin(lon)
        e_cart(3) = sin(lat)
    
    end function lib_lonlat2cart
    !===========================================================================
    pure function lib_great_circ_dist(lonlat1, lonlat2, radius) result(dist)
        real(kind=R_GRID), dimension(2), intent(in) :: lonlat1, lonlat2
        real(kind=R_GRID), intent(in), optional :: radius
        real(kind=R_GRID) :: dist
        ! local
        real(kind=R_GRID) :: lat1,lat2,lon1,lon2,dlat,dlon
        real(kind=R_GRID) :: beta, tmp1, tmp2
        
        lon1 = lonlat1(1); lon2 = lonlat2(1)
        lat1 = lonlat1(2); lat2 = lonlat2(2)
        dlat = lat2-lat1; dlon = lon2-lon1

        tmp1 = sqrt((cos(lat2)*sin(dlon))**2 + &
               (cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon))**2)
        tmp2 = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)

        beta = atan2(tmp1,tmp2)

        if (present(radius)) then
            dist = radius * beta
        else
            dist = beta
        endif

    end function lib_great_circ_dist
    !=============================================================================
    real function lib_4pt_area(p1_ll,p2_ll,p3_ll,p4_ll,r) result(area)
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

        p12 = cross_prod(p1,p2)
        p23 = cross_prod(p2,p3)
        p34 = cross_prod(p3,p4)
        p41 = cross_prod(p4,p1)

        p12 = norm_vec_3d(p12)
        p23 = norm_vec_3d(p23)
        p34 = norm_vec_3d(p34)
        p41 = norm_vec_3d(p41)

        c123 = dot_product(p12,p23)
        c234 = dot_product(p23,p34)
        c341 = dot_product(p34,p41)
        c412 = dot_product(p41,p12)

        a123 = acos(-c123)
        a234 = acos(-c234)
        a341 = acos(-c341)
        a412 = acos(-c412)

        area = ( a123+a234+a341+a412-2*pi ) * r**2

    end function lib_4pt_area
    !===========================================================================
    function lib_2pt_unit_vec(p1_ll,p2_ll) result(e)
        real(kind=R_GRID), dimension(3) :: e
        real(kind=R_GRID), dimension(2), intent(in):: p1_ll, p2_ll
        ! local
        real(kind=R_GRID), dimension(3) :: p1, p2
        real(kind=R_GRID), dimension(3) :: p12

        p1 = lib_lonlat2cart(p1_ll)
        p2 = lib_lonlat2cart(p2_ll)

        p12 = cross_prod(p1,p2)
        e = cross_prod(p12,p1)

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

        dist = sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
        e = r/dist

    end function norm_vec_3d
    !===========================================================================
    pure function cross_prod_r(a,b) result(rlt)
        real(kind=R_GRID), dimension(3) :: rlt
        real(kind=R_GRID), dimension(3), intent(in) :: a,b

        rlt(1) = a(2)*b(3) - a(3)*b(2)
        rlt(2) = a(3)*b(1) - a(1)*b(3)
        rlt(3) = a(1)*b(2) - a(2)*b(1)

    end function cross_prod_r
    !---------------------------------------------------------------------------
    pure function cross_prod_d(a,b) result(rlt)
        integer, dimension(3) :: rlt
        integer, dimension(3), intent(in) :: a,b

        rlt(1) = a(2)*b(3) - a(3)*b(2)
        rlt(2) = a(3)*b(1) - a(1)*b(3)
        rlt(3) = a(1)*b(2) - a(2)*b(1)

    end function cross_prod_d
    !=============================================================================
    ! Lagrangian polynomial interpolation only make sense when remapping on fixed
    ! locations. Otherwise, use the lib_interp function.
    function lib_interp_lag_get_coef(x,vx) result(coef)
        ! coef = lib_interp_lag_get_coef(x,vx)
        real(kind=R_GRID), dimension(:), intent(in) :: vx
        real(kind=R_GRID), intent(in) :: x
        real(kind=R_GRID), dimension(size(vx)) :: coef
        ! local
        integer :: i,j,n

        n = size(vx)
        do j = 1,n
            coef(j) = 1.
            do i = 1,n
                if (i.ne.j) then
                    coef(j) = coef(j)*(x-vx(i))/(vx(j)-vx(i))
                endif
            enddo
        enddo

    end function lib_interp_lag_get_coef
    !===========================================================================
    pure function lib_interp_lag_eval(coef,vy) result(y)
        ! y = lib_interp_lag_eval(coef,vy)
        real :: y
        real(kind=R_GRID), dimension(:), intent(in) :: coef
        real, dimension(size(coef)), intent(in) :: vy
        ! local
        integer :: i

        y = 0.
        do i = 1,size(vy)
            y = y+vy(i)*coef(i)
        enddo

    end function lib_interp_lag_eval
    !=============================================================================
end module lib_grid_mod

