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
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module fv_tracer2d_mod
   use tp_core_mod,       only: fv_tp_2d, copy_corners
   use fv_mp_mod,         only: mp_reduce_max
   use fv_mp_mod,         only: mp_gather, is_master
   use fv_mp_mod,         only: group_halo_update_type
   use fv_mp_mod,         only: start_group_halo_update, complete_group_halo_update
   use mpp_domains_mod,   only: mpp_update_domains, CGRID_NE, domain2d, mpp_get_boundary
   use fv_timing_mod,     only: timing_on, timing_off
   use boundary_mod,      only: nested_grid_BC_apply_intT
   use fv_regional_mod,   only: regional_boundary_update
   use fv_regional_mod,   only: current_time_in_seconds
   use fv_arrays_mod,     only: fv_grid_type, fv_nest_type, fv_atmos_type, fv_grid_bounds_type
   use mpp_mod,           only: mpp_error, FATAL, mpp_broadcast, mpp_send, mpp_recv, mpp_sum, mpp_max, mpp_pe
   use duogrid_mod,       only: ext_scalar

implicit none
private

public :: tracer_2d, tracer_2d_nested, tracer_2d_1L

real, allocatable, dimension(:,:,:) :: nest_fx_west_accum, nest_fx_east_accum, nest_fx_south_accum, nest_fx_north_accum

contains

  subroutine flux_adj(fx, fy, bd, domain, npx, npy, npz)
    type(fv_grid_bounds_type), intent(IN) :: bd
    real, intent(inout) :: fx(bd%is:bd%ie + 1, bd%js:bd%je)
    real, intent(inout) :: fy(bd%is:bd%ie, bd%js:bd%je + 1)
    type(domain2d), intent(INOUT) :: domain
    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz

    real :: fxx(bd%isd:bd%ied + 1, bd%jsd:bd%jed) ! working array to average fluxes
    real :: fyy(bd%isd:bd%ied, bd%jsd:bd%jed + 1) ! working array to average fluxes
    real wbuffer(npy + 2)
    real ebuffer(npy + 2)
    real nbuffer(npx + 2)
    real sbuffer(npx + 2)
    integer :: i, j, is, ie, js, je
    integer :: isd, ied, jsd, jed

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    do j = js, je
      do i = is, ie + 1
        fxx(i, j) = fx(i, j)
      end do
    end do
    do j = js, je + 1
      do i = is, ie
        fyy(i, j) = fy(i, j)
      end do
    end do

    call mpp_get_boundary(fxx, fyy, domain, &
                          wbufferx=wbuffer, ebufferx=ebuffer, &
                          sbuffery=sbuffer, nbuffery=nbuffer, gridtype=CGRID_NE)

    do i = is, ie
      fyy(i, js) = 0.5*(fyy(i, js) + sbuffer(i - is + 1))
      fyy(i, je + 1) = 0.5*(fyy(i, je + 1) + nbuffer(i - is + 1))
    end do
    do j = js, je
      fxx(is, j) = 0.5*(fxx(is, j) + wbuffer(j - js + 1))
      fxx(ie + 1, j) = 0.5*(fxx(ie + 1, j) + ebuffer(j - js + 1))
    end do

    do j = js, je
      do i = is, ie + 1
        fx(i, j) = fxx(i, j)
      end do
    end do
    do j = js, je + 1
      do i = is, ie
        fy(i, j) = fyy(i, j)
      end do
    end do

  end subroutine flux_adj

!-----------------------------------------------------------------------
! !ROUTINE: Perform 2D horizontal-to-lagrangian transport
!-----------------------------------------------------------------------

  subroutine tracer_2d_1L(q, dp1, delp, mfx, mfy, cx, cy, cx_dp2, cy_dp2, &
                          gridstruct, bd, domain, npx, npy, npz, &
                          nq, hord, q_split, dt, id_divg, q_pack, dp1_pack, nord_tr, trdm, lim_fac)

    type(fv_grid_bounds_type), intent(IN) :: bd
    integer, intent(IN) :: npx
    integer, intent(IN) :: npy
    integer, intent(IN) :: npz
    integer, intent(IN) :: nq    ! number of tracers to be advected
    integer, intent(IN) :: hord, nord_tr
    integer, intent(IN) :: q_split
    integer, intent(IN) :: id_divg
    real, intent(IN) :: dt, trdm
    real, intent(IN) :: lim_fac
    type(group_halo_update_type), intent(inout) :: q_pack, dp1_pack
    real, intent(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, npz, nq)   ! Tracers
    real, intent(INOUT) :: dp1(bd%isd:bd%ied, bd%jsd:bd%jed, npz)        ! DELP before dyn_core
    real, intent(INOUT) :: delp(bd%isd:bd%ied, bd%jsd:bd%jed, npz)       ! DELP after  dyn_core
    real, intent(INOUT) :: mfx(bd%is:bd%ie + 1, bd%js:bd%je, npz)    ! Mass Flux X-Dir
    real, intent(INOUT) :: mfy(bd%is:bd%ie, bd%js:bd%je + 1, npz)    ! Mass Flux Y-Dir
    real, intent(INOUT) ::  cx(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)  ! Courant Number X-Dir
    real, intent(INOUT) ::  cy(bd%isd:bd%ied, bd%js:bd%je + 1, npz)  ! Courant Number Y-Dir
    real, intent(INOUT) ::  cx_dp2(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)  ! Courant Number X-Dir
    real, intent(INOUT) ::  cy_dp2(bd%isd:bd%ied, bd%js:bd%je + 1, npz)  ! Courant Number Y-Dir
    type(fv_grid_type), intent(INOUT), target :: gridstruct
    type(domain2d), intent(INOUT) :: domain

! Local Arrays
    real :: qn2(bd%isd:bd%ied, bd%jsd:bd%jed, nq)   ! 3D tracers
    real :: dp2(bd%is:bd%ie, bd%js:bd%je)
    real :: fx(bd%is:bd%ie + 1, bd%js:bd%je)
    real :: fy(bd%is:bd%ie, bd%js:bd%je + 1)

    real :: fxx(bd%isd:bd%ied + 1, bd%jsd:bd%jed) ! working array to average fluxes
    real :: fyy(bd%isd:bd%ied, bd%jsd:bd%jed + 1) ! working array to average fluxes
    ! real wbuffer(npy+2,npz)
    ! real ebuffer(npy+2,npz)
    ! real nbuffer(npx+2,npz)
    ! real sbuffer(npx+2,npz)
    real wbuffer(npy + 2)
    real ebuffer(npy + 2)
    real nbuffer(npx + 2)
    real sbuffer(npx + 2)

    real :: ra_x(bd%is:bd%ie, bd%jsd:bd%jed)
    real :: ra_y(bd%isd:bd%ied, bd%js:bd%je)
    real :: xfx(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)
    real :: yfx(bd%isd:bd%ied, bd%js:bd%je + 1, npz)
    real :: xfx_dp2(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)
    real :: yfx_dp2(bd%isd:bd%ied, bd%js:bd%je + 1, npz)
    real :: cmax(npz)
    real :: frac
    integer :: nsplt
    integer :: adv_scheme
    integer :: i, j, k, it, iq

    real, pointer, dimension(:, :) :: area, rarea
    real, pointer, dimension(:, :, :) :: sin_sg
    real, pointer, dimension(:, :) :: dxa, dya, dx, dy

    real, pointer, dimension(:)   :: dxa_cs, dya_cs
    real, pointer, dimension(:)   :: dxc_cs, dyc_cs

    integer :: is, ie, js, je
    integer :: isd, ied, jsd, jed

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed
    adv_scheme = gridstruct%adv_scheme

    area => gridstruct%area
    rarea => gridstruct%rarea

    sin_sg => gridstruct%sin_sg
    dxa => gridstruct%dxa
    dya => gridstruct%dya
    dx => gridstruct%dx
    dy => gridstruct%dy

    dxa_cs => gridstruct%dxa_cs
    dya_cs => gridstruct%dya_cs
    dxc_cs => gridstruct%dx_cs
    dyc_cs => gridstruct%dy_cs
 
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  adv_scheme,sin_sg,cy,yfx,dya,dx,cmax, &
!$OMP                                  xfx_dp2,yfx_dp2,cx_dp2,cy_dp2, &
!$OMP                                  dxa_cs,dya_cs,dxc_cs,dyc_cs)
    do k = 1, npz
      do j = jsd, jed
        do i = is, ie + 1
          if (cx(i, j, k) > 0.) then
            xfx(i, j, k) = cx(i, j, k)*dxa(i - 1, j)*dy(i, j)*sin_sg(i - 1, j, 3)
          else
            xfx(i, j, k) = cx(i, j, k)*dxa(i, j)*dy(i, j)*sin_sg(i, j, 1)
          end if
        end do
      end do
      do j = js, je + 1
        do i = isd, ied
          if (cy(i, j, k) > 0.) then
            yfx(i, j, k) = cy(i, j, k)*dya(i, j - 1)*dx(i, j)*sin_sg(i, j - 1, 4)
          else
            yfx(i, j, k) = cy(i, j, k)*dya(i, j)*dx(i, j)*sin_sg(i, j, 2)
          end if
        end do
      end do

      if(adv_scheme==2)then
        do j = jsd, jed
          do i = is, ie + 1
            xfx_dp2(i, j, k) = cx_dp2(i, j, k)*dxc_cs(i)*dya_cs(j)
          end do
        end do
        do j = js, je + 1
          do i = isd, ied
            yfx_dp2(i, j, k) = cy_dp2(i, j, k)*dxa_cs(i)*dyc_cs(j)
          end do
        end do
      endif

      cmax(k) = 0.
      if (k < npz/6) then
        do j = js, je
          do i = is, ie
            cmax(k) = max(cmax(k), abs(cx(i, j, k)), abs(cy(i, j, k)))
          end do
        end do
      else
        do j = js, je
          do i = is, ie
            cmax(k) = max(cmax(k), max(abs(cx(i, j, k)), abs(cy(i, j, k))) + 1.-sin_sg(i, j, 5))
          end do
        end do
      end if
    end do  ! k-loop

    if (trdm > 1.e-4) then
      call timing_on('COMM_TOTAL')
      call timing_on('COMM_TRACER')
      if (.not. gridstruct%dg%is_initialized) then
        call complete_group_halo_update(dp1_pack, domain)
      else
        call ext_scalar(dp1, gridstruct%dg, bd, domain, 0, 0)
      end if
      call timing_off('COMM_TRACER')
      call timing_off('COMM_TOTAL')

    end if
    call mp_reduce_max(cmax, npz)

    if(adv_scheme==2)then
      call timing_on('COMM_TOTAL')
      call timing_on('COMM_TRACER')
      call ext_scalar(dp1, gridstruct%dg, bd, domain, 0, 0)
      call timing_off('COMM_TRACER')
      call timing_off('COMM_TOTAL')
    endif

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx, &
!$OMP                                  cy,yfx,mfx,mfy,cmax)   &
!$OMP                          private(nsplt, frac)
    do k = 1, npz

      nsplt = int(1.+cmax(k))
      if (nsplt > 1) then
        frac = 1./real(nsplt)
        do j = jsd, jed
          do i = is, ie + 1
            cx(i, j, k) = cx(i, j, k)*frac
            xfx(i, j, k) = xfx(i, j, k)*frac
          end do
        end do
        do j = js, je
          do i = is, ie + 1
            mfx(i, j, k) = mfx(i, j, k)*frac
          end do
        end do
        do j = js, je + 1
          do i = isd, ied
            cy(i, j, k) = cy(i, j, k)*frac
            yfx(i, j, k) = yfx(i, j, k)*frac
          end do
        end do
        do j = js, je + 1
          do i = is, ie
            mfy(i, j, k) = mfy(i, j, k)*frac
          end do
        end do
      end if

    end do
    call timing_on('COMM_TOTAL')
    call timing_on('COMM_TRACER')
    if (.not. gridstruct%dg%is_initialized) then
      call complete_group_halo_update(q_pack, domain)
    else
      call ext_scalar(q, gridstruct%dg, bd, domain, 0, 0)
    end if
    call timing_off('COMM_TRACER')
    call timing_off('COMM_TOTAL')

! Begin k-independent tracer transport; can not be OpenMPed because the mpp_update call.
    do k = 1, npz

!$OMP parallel do default(none) shared(k,is,ie,js,je,isd,ied,jsd,jed,xfx,area,yfx,ra_x,ra_y)
      do j = jsd, jed
        do i = is, ie
          ra_x(i, j) = area(i, j) + xfx(i, j, k) - xfx(i + 1, j, k)
        end do
        if (j >= js .and. j <= je) then
          do i = isd, ied
            ra_y(i, j) = area(i, j) + yfx(i, j, k) - yfx(i, j + 1, k)
          end do
        end if
      end do

      nsplt = int(1.+cmax(k))

      if (nsplt == 1) then
        if(adv_scheme==1) then
!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx,mfy,dp1,dp2)
          do j = js, je
            do i = is, ie
              dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k) - mfx(i + 1, j, k) + mfy(i, j, k) - mfy(i, j + 1, k))*rarea(i, j)
            end do
          end do
          do iq = 1, nq
            call fv_tp_2d(q(isd, jsd, k, iq), cx(is, jsd, k), cy(isd, js, k), &
                          npx, npy, hord, fx, fy, xfx(is, jsd, k), yfx(isd, js, k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is, js, k), mfy=mfy(is, js, k))
!          if (gridstruct%dg%is_initialized) then
!            call flux_adj(fx, fy, bd, domain, npx, npy, npz)
!          end if !if duo

            do j = js, je
              do i = is, ie
            q(i, j, k, iq) = (q(i, j, k, iq)*dp1(i, j, k) + (fx(i, j) - fx(i + 1, j) + fy(i, j) - fy(i, j + 1))*rarea(i, j))/dp2(i, j)
              end do
            end do
          end do! iq
        else if(adv_scheme==2) then
          do iq = 1, nq
            !$OMP parallel do default(none) shared(k,isd,ied,jsd,jed,iq,dp1,q)
            do j = jsd, jed
              do i = isd, ied
                q(i, j, k, iq) = dp1(i, j, k)*q(i, j, k, iq)
              end do
            end do
            call fv_tp_2d(q(isd, jsd, k, iq), cx_dp2, cy_dp2, npx, npy, hord, fx, fy, &
                          xfx_dp2, yfx_dp2, gridstruct, bd, ra_x, ra_y, lim_fac, &
                          advscheme=adv_scheme)
!          if (gridstruct%dg%is_initialized) then
!            call flux_adj(fx, fy, bd, domain, npx, npy, npz)
!          end if !if duo

            do j = js, je
              do i = is, ie
                q(i, j, k, iq) = (q(i, j, k, iq) + (fx(i, j) - fx(i + 1, j) + fy(i, j) - fy(i, j + 1))*rarea(i, j))/delp(i, j, k)
              end do
            end do
          end do! iq

        end if
      end if

      if (nsplt /= 1) then

        do it = 1, nsplt

!$OMP parallel do default(none) shared(k,is,ie,js,je,rarea,mfx,mfy,dp1,dp2)
          do j = js, je
            do i = is, ie
              dp2(i, j) = dp1(i, j, k) + (mfx(i, j, k) - mfx(i + 1, j, k) + mfy(i, j, k) - mfy(i, j + 1, k))*rarea(i, j)
            end do
          end do

!$OMP parallel do default(none) shared(k,nsplt,it,is,ie,js,je,isd,ied,jsd,jed,npx,npy,cx,xfx,hord,trdm, &
!$OMP                                  nord_tr,nq,gridstruct,bd,cy,yfx,mfx,mfy,qn2,q,ra_x,ra_y,dp1,dp2,rarea,lim_fac) &
!$OMP                          private(fx,fy)
          do iq = 1, nq
            if (it == 1) then
              do j = jsd, jed
                do i = isd, ied
                  qn2(i, j, iq) = q(i, j, k, iq)
                end do
              end do
            end if
            call fv_tp_2d(qn2(isd, jsd, iq), cx(is, jsd, k), cy(isd, js, k), &
                          npx, npy, hord, fx, fy, xfx(is, jsd, k), yfx(isd, js, k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is, js, k), mfy=mfy(is, js, k))
!            if (gridstruct%dg%is_initialized) then
!              call flux_adj(fx, fy, bd, domain, npx, npy, npz)
!            end if !if duo

            if (it < nsplt) then   ! not last call
              do j = js, je
              do i = is, ie
            qn2(i, j, iq) = (qn2(i, j, iq)*dp1(i, j, k) + (fx(i, j) - fx(i + 1, j) + fy(i, j) - fy(i, j + 1))*rarea(i, j))/dp2(i, j)
              end do
              end do
            else
              do j = js, je
              do i = is, ie
           q(i, j, k, iq) = (qn2(i, j, iq)*dp1(i, j, k) + (fx(i, j) - fx(i + 1, j) + fy(i, j) - fy(i, j + 1))*rarea(i, j))/dp2(i, j)
              end do
              end do
            end if

          end do   !  tracer-loop

          if (it < nsplt) then   ! not last call
            do j = js, je
              do i = is, ie
                dp1(i, j, k) = dp2(i, j)
              end do
            end do
            call timing_on('COMM_TOTAL')
            call timing_on('COMM_TRACER')
            if (.not. gridstruct%dg%is_initialized) then
              call mpp_update_domains(qn2, domain)
            else
              call ext_scalar(qn2, gridstruct%dg, bd, domain, 0, 0)
            end if
            call timing_off('COMM_TRACER')
            call timing_off('COMM_TOTAL')
          end if
        end do  ! time-split loop
      end if
    end do    ! k-loop

  end subroutine tracer_2d_1L




subroutine tracer_2d(q, dp1, mfx, mfy, cx, cy, cx_dp2, cy_dp2, &
                     gridstruct,bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, dp1_pack, nord_tr, trdm, lim_fac)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack, dp1_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      real   , intent(INOUT) ::  cx_dp2(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy_dp2(bd%isd:bd%ied, bd%js:bd%je + 1, npz)  ! Courant Number Y-Dir
      type(fv_grid_type), intent(INOUT), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: xfx_dp2(bd%is:bd%ie + 1, bd%jsd:bd%jed, npz)
      real :: yfx_dp2(bd%isd:bd%ied, bd%js:bd%je + 1, npz)
      real :: cmax(npz)
      real :: c_global
      real :: frac, rdt
      integer :: ksplt(npz)
      integer :: nsplt
      integer :: adv_scheme
      integer :: i,j,k,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      real, pointer, dimension(:)   :: dxa_cs, dya_cs
      real, pointer, dimension(:)   :: dxc_cs, dyc_cs

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed
      adv_scheme = gridstruct%adv_scheme

      if(mpp_pe()==0) print*, 'hitracer2d'
       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa
      dya    => gridstruct%dya
      dx     => gridstruct%dx
      dy     => gridstruct%dy

      dxa_cs => gridstruct%dxa_cs
      dya_cs => gridstruct%dya_cs
      dxc_cs => gridstruct%dx_cs
      dyc_cs => gridstruct%dy_cs

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  adv_scheme,sin_sg,cy,yfx,dya,dx,cmax,q_split,ksplt, &
!$OMP                                  xfx_dp2,yfx_dp2,cx_dp2,cy_dp2, &
!$OMP                                  dxa_cs,dya_cs,dxc_cs,dyc_cs)
    do k=1,npz
       do j=jsd,jed
          do i=is,ie+1
             if (cx(i,j,k) > 0.) then
                 xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
             else
                 xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
             endif
          enddo
       enddo
       do j=js,je+1
          do i=isd,ied
              if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
              else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
              endif
          enddo
       enddo

       if(adv_scheme==2)then
         do j = jsd, jed
           do i = is, ie + 1
             xfx_dp2(i, j, k) = cx_dp2(i, j, k)*dxc_cs(i)*dya_cs(j)
           end do
         end do
         do j = js, je + 1
           do i = isd, ied
             yfx_dp2(i, j, k) = cy_dp2(i, j, k)*dxa_cs(i)*dyc_cs(j)
           end do
         end do
       endif

       if ( q_split == 0 ) then
         cmax(k) = 0.
         if ( k < npz/6 ) then
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), abs(cx(i,j,k)), abs(cy(i,j,k)) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax(k) = max( cmax(k), max(abs(cx(i,j,k)),abs(cy(i,j,k)))+1.-sin_sg(i,j,5) )
               enddo
            enddo
         endif
       endif
       ksplt(k) = 1

    enddo

!--------------------------------------------------------------------------------

! Determine global nsplt:
  if ( q_split == 0 ) then
      call mp_reduce_max(cmax,npz)
! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)
      if ( is_master() .and. nsplt > 4 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global
   else
      nsplt = q_split
   endif

!--------------------------------------------------------------------------------

    if( nsplt /= 1 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,mfx,cy,yfx,mfy,cmax,nsplt,ksplt) &
!$OMP                          private( frac )
        do k=1,npz

#ifdef GLOBAL_CFL
           ksplt(k) = nsplt
#else
           ksplt(k) = int(1. + cmax(k))
#endif
           frac  = 1. / real(ksplt(k))

           do j=jsd,jed
              do i=is,ie+1
                 cx(i,j,k) =   cx(i,j,k) * frac
                 xfx(i,j,k) = xfx(i,j,k) * frac
              enddo
           enddo
           do j=js,je
              do i=is,ie+1
                 mfx(i,j,k) = mfx(i,j,k) * frac
              enddo
           enddo

           do j=js,je+1
              do i=isd,ied
                 cy(i,j,k) =  cy(i,j,k) * frac
                yfx(i,j,k) = yfx(i,j,k) * frac
              enddo
           enddo
           do j=js,je+1
              do i=is,ie
                mfy(i,j,k) = mfy(i,j,k) * frac
              enddo
           enddo

        enddo
    endif

    if (trdm>1.e-4) then
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      if (.not. gridstruct%dg%is_initialized) then
        call complete_group_halo_update(dp1_pack, domain)
      else
        call ext_scalar(dp1, gridstruct%dg, bd, domain, 0, 0)
      end if
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

    endif
    do it=1,nsplt
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      if (.not. gridstruct%dg%is_initialized) then
        call complete_group_halo_update(q_pack, domain)
      else
        call ext_scalar(q, gridstruct%dg, bd, domain, 0, 0)
      end if
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq,ksplt,&
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm,lim_fac) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
     do k=1,npz

       if ( it .le. ksplt(k) ) then

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j,k) + (mfx(i,j,k)-mfx(i+1,j,k)+mfy(i,j,k)-mfy(i,j+1,k))*rarea(i,j)
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
            enddo
         enddo

         do iq=1,nq
         if ( it==1 .and. trdm>1.e-4 ) then
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
         else
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
         endif
!          if (gridstruct%dg%is_initialized) then
!            call flux_adj(fx, fy, bd, domain, npx, npy, npz)
!          end if !if duo
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j) )/dp2(i,j)
               enddo
               enddo
            enddo

         if ( it /= nsplt ) then
              do j=js,je
                 do i=is,ie
                    dp1(i,j,k) = dp2(i,j)
                 enddo
              enddo
         endif

       endif   ! ksplt

     enddo ! npz

      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
      endif

   enddo  ! nsplt


end subroutine tracer_2d


subroutine tracer_2d_nested(q, dp1, mfx, mfy, cx, cy, gridstruct, bd, domain, npx, npy, npz,   &
                     nq,  hord, q_split, dt, id_divg, q_pack, dp1_pack, nord_tr, trdm, &
                     k_split, neststruct, parent_grid, n_map, lim_fac)

      type(fv_grid_bounds_type), intent(IN) :: bd
      integer, intent(IN) :: npx
      integer, intent(IN) :: npy
      integer, intent(IN) :: npz
      integer, intent(IN) :: nq    ! number of tracers to be advected
      integer, intent(IN) :: hord, nord_tr
      integer, intent(IN) :: q_split, k_split, n_map
      integer, intent(IN) :: id_divg
      real   , intent(IN) :: dt, trdm
      real   , intent(IN) :: lim_fac
      type(group_halo_update_type), intent(inout) :: q_pack, dp1_pack
      real   , intent(INOUT) :: q(bd%isd:bd%ied,bd%jsd:bd%jed,npz,nq)   ! Tracers
      real   , intent(INOUT) :: dp1(bd%isd:bd%ied,bd%jsd:bd%jed,npz)        ! DELP before dyn_core
      real   , intent(INOUT) :: mfx(bd%is:bd%ie+1,bd%js:bd%je,  npz)    ! Mass Flux X-Dir
      real   , intent(INOUT) :: mfy(bd%is:bd%ie  ,bd%js:bd%je+1,npz)    ! Mass Flux Y-Dir
      real   , intent(INOUT) ::  cx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)  ! Courant Number X-Dir
      real   , intent(INOUT) ::  cy(bd%isd:bd%ied,bd%js :bd%je +1,npz)  ! Courant Number Y-Dir
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_nest_type), intent(INOUT) :: neststruct
      type(fv_atmos_type), pointer, intent(IN) :: parent_grid
      type(domain2d), intent(INOUT) :: domain

! Local Arrays
      real :: dp2(bd%is:bd%ie,bd%js:bd%je)
      real :: fx(bd%is:bd%ie+1,bd%js:bd%je )
      real :: fy(bd%is:bd%ie , bd%js:bd%je+1)
      real :: ra_x(bd%is:bd%ie,bd%jsd:bd%jed)
      real :: ra_y(bd%isd:bd%ied,bd%js:bd%je)
      real :: xfx(bd%is:bd%ie+1,bd%jsd:bd%jed  ,npz)
      real :: yfx(bd%isd:bd%ied,bd%js: bd%je+1, npz)
      real :: cmax(npz)
      real :: cmax_t
      real :: c_global
      real :: frac, rdt
      real :: reg_bc_update_time
      integer :: nsplt, nsplt_parent, msg_split_steps = 1
      integer :: i,j,k,it,iq

      real, pointer, dimension(:,:) :: area, rarea
      real, pointer, dimension(:,:,:) :: sin_sg
      real, pointer, dimension(:,:) :: dxa, dya, dx, dy

      integer :: is,  ie,  js,  je
      integer :: isd, ied, jsd, jed

      is  = bd%is
      ie  = bd%ie
      js  = bd%js
      je  = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      if(mpp_pe()==0) print*, '2d-nested'
       area => gridstruct%area
      rarea => gridstruct%rarea

      sin_sg => gridstruct%sin_sg
      dxa    => gridstruct%dxa
      dya    => gridstruct%dya
      dx     => gridstruct%dx
      dy     => gridstruct%dy

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,xfx,dxa,dy, &
!$OMP                                  sin_sg,cy,yfx,dya,dx)
      do k=1,npz
         do j=jsd,jed
            do i=is,ie+1
               if (cx(i,j,k) > 0.) then
                  xfx(i,j,k) = cx(i,j,k)*dxa(i-1,j)*dy(i,j)*sin_sg(i-1,j,3)
               else
                  xfx(i,j,k) = cx(i,j,k)*dxa(i,j)*dy(i,j)*sin_sg(i,j,1)
               endif
            enddo
         enddo
         do j=js,je+1
            do i=isd,ied
               if (cy(i,j,k) > 0.) then
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j-1)*dx(i,j)*sin_sg(i,j-1,4)
               else
                  yfx(i,j,k) = cy(i,j,k)*dya(i,j)*dx(i,j)*sin_sg(i,j,2)
               endif
            enddo
         enddo
      enddo

!--------------------------------------------------------------------------------
  if ( q_split == 0 ) then
! Determine nsplt

!$OMP parallel do default(none) shared(is,ie,js,je,npz,cmax,cx,cy,sin_sg) &
!$OMP                          private(cmax_t )
      do k=1,npz
         cmax(k) = 0.
         if ( k < 4 ) then
! Top layers: C < max( abs(c_x), abs(c_y) )
            do j=js,je
               do i=is,ie
                  cmax_t  = max( abs(cx(i,j,k)), abs(cy(i,j,k)) )
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         else
            do j=js,je
               do i=is,ie
                  cmax_t  = max(abs(cx(i,j,k)), abs(cy(i,j,k))) + 1.-sin_sg(i,j,5)
                  cmax(k) = max( cmax_t, cmax(k) )
               enddo
            enddo
         endif
      enddo
      call mp_reduce_max(cmax,npz)

! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
      c_global = cmax(1)
      if ( npz /= 1 ) then                ! if NOT shallow water test case
         do k=2,npz
            c_global = max(cmax(k), c_global)
         enddo
      endif
      nsplt = int(1. + c_global)
      if ( is_master() .and. nsplt > 3 )  write(*,*) 'Tracer_2d_split=', nsplt, c_global
   else
      nsplt = q_split
      if (gridstruct%nested .and. neststruct%nestbctype > 1) msg_split_steps = max(q_split/parent_grid%flagstruct%q_split,1)
   endif

!--------------------------------------------------------------------------------

   frac  = 1. / real(nsplt)

      if( nsplt /= 1 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,cx,frac,xfx,mfx,cy,yfx,mfy)
          do k=1,npz
             do j=jsd,jed
                do i=is,ie+1
                   cx(i,j,k) =  cx(i,j,k) * frac
                   xfx(i,j,k) = xfx(i,j,k) * frac
                enddo
             enddo
             do j=js,je
                do i=is,ie+1
                   mfx(i,j,k) = mfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=isd,ied
                   cy(i,j,k) =  cy(i,j,k) * frac
                  yfx(i,j,k) = yfx(i,j,k) * frac
                enddo
             enddo

             do j=js,je+1
                do i=is,ie
                  mfy(i,j,k) = mfy(i,j,k) * frac
                enddo
             enddo
          enddo
      endif


    do it=1,nsplt
       if ( gridstruct%nested ) then
          neststruct%tracer_nest_timestep = neststruct%tracer_nest_timestep + 1
       end if
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
      call complete_group_halo_update(q_pack, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

      if (gridstruct%nested) then
            do iq=1,nq
                 call nested_grid_BC_apply_intT(q(isd:ied,jsd:jed,:,iq), &
                      0, 0, npx, npy, npz, bd, &
                      real(neststruct%tracer_nest_timestep)+real(nsplt*k_split), real(nsplt*k_split), &
                 neststruct%q_BC(iq), bctype=neststruct%nestbctype  )
           enddo
      endif

      if (gridstruct%regional) then
            !This is more accurate than the nested BC calculation
            ! since it takes into account varying nsplit
            reg_bc_update_time=current_time_in_seconds+(real(n_map-1) + real(it-1)*frac)*dt
            do iq=1,nq
                 call regional_boundary_update(q(:,:,:,iq), 'q', &
                                               isd, ied, jsd, jed, npz, &
                                               is,  ie,  js,  je,       &
                                               isd, ied, jsd, jed,      &
                                               reg_bc_update_time,      &
                                               it, iq )
            enddo
      endif

      if (trdm>1.e-4) then
                        call timing_on('COMM_TOTAL')
                            call timing_on('COMM_TRACER')
         call complete_group_halo_update(dp1_pack, domain)
                           call timing_off('COMM_TRACER')
                       call timing_off('COMM_TOTAL')

      endif


!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,npz,dp1,mfx,mfy,rarea,nq, &
!$OMP                                  area,xfx,yfx,q,cx,cy,npx,npy,hord,gridstruct,bd,it,nsplt,nord_tr,trdm,lim_fac) &
!$OMP                          private(dp2, ra_x, ra_y, fx, fy)
      do k=1,npz

         do j=js,je
            do i=is,ie
               dp2(i,j) = dp1(i,j,k) + (mfx(i,j,k)-mfx(i+1,j,k)+mfy(i,j,k)-mfy(i,j+1,k))*rarea(i,j)
            enddo
         enddo

         do j=jsd,jed
            do i=is,ie
               ra_x(i,j) = area(i,j) + xfx(i,j,k) - xfx(i+1,j,k)
            enddo
         enddo
         do j=js,je
            do i=isd,ied
               ra_y(i,j) = area(i,j) + yfx(i,j,k) - yfx(i,j+1,k)
            enddo
         enddo

         do iq=1,nq
         if ( it==1 .and. trdm>1.e-4 ) then
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k),   &
                          mass=dp1(isd,jsd,k), nord=nord_tr, damp_c=trdm)
         else
            call fv_tp_2d(q(isd,jsd,k,iq), cx(is,jsd,k), cy(isd,js,k), &
                          npx, npy, hord, fx, fy, xfx(is,jsd,k), yfx(isd,js,k), &
                          gridstruct, bd, ra_x, ra_y, lim_fac, mfx=mfx(is,js,k), mfy=mfy(is,js,k))
         endif
            do j=js,je
               do i=is,ie
                  q(i,j,k,iq) = ( q(i,j,k,iq)*dp1(i,j,k) + &
                                (fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))*rarea(i,j) )/dp2(i,j)
               enddo
               enddo
          enddo
      enddo ! npz

      if ( it /= nsplt ) then
                      call timing_on('COMM_TOTAL')
                          call timing_on('COMM_TRACER')
           call start_group_halo_update(q_pack, q, domain)
                          call timing_off('COMM_TRACER')
                      call timing_off('COMM_TOTAL')
      endif

   enddo  ! nsplt

   if ( id_divg > 0 ) then
        rdt = 1./(frac*dt)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,dp1,xfx,yfx,rarea,rdt)
        do k=1,npz
        do j=js,je
           do i=is,ie
              dp1(i,j,k) = (xfx(i+1,j,k)-xfx(i,j,k) + yfx(i,j+1,k)-yfx(i,j,k))*rarea(i,j)*rdt
           enddo
        enddo
        enddo
   endif

 end subroutine tracer_2d_nested

end module fv_tracer2d_mod
