! -*-f90-*-
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
    subroutine MPP_DO_UPDATE_3D_V_(f_addrsx,f_addrsy, domain, update_x, update_y, &
                                   d_type, ke, gridtype, flags)
!updates data domain of 3D field whose computational domains have been computed
      integer(i8_kind),  intent(in)        :: f_addrsx(:,:), f_addrsy(:,:)
      type(domain2d),      intent(in)        :: domain
      type(overlapSpec),   intent(in)        :: update_x, update_y
      integer,             intent(in)        :: ke
      MPP_TYPE_, intent(in)                  :: d_type  ! creates unique interface
      integer, intent(in)                    :: gridtype
      integer, intent(in),          optional :: flags

      MPP_TYPE_ :: fieldx(update_x%xbegin:update_x%xend, update_x%ybegin:update_x%yend,ke)
      MPP_TYPE_ :: fieldy(update_y%xbegin:update_y%xend, update_y%ybegin:update_y%yend,ke)
      pointer(ptr_fieldx, fieldx)
      pointer(ptr_fieldy, fieldy)

      integer :: update_flags
      integer :: l_size, l, i, j, k, is, ie, js, je, n, m
      integer :: pos, nlist, msgsize, isd, ied, jsd, jed
      integer :: to_pe, from_pe, midpoint
      integer :: tMe, dir

      integer :: send_start_pos, nsend
      integer :: send_msgsize(2*MAXLIST)
      integer :: send_pe(2*MAXLIST)
      integer,    allocatable :: msg1(:), msg2(:), msg3(:)
      logical :: send(8), recv(8), update_edge_only
      MPP_TYPE_ :: buffer(size(mpp_domains_stack(:)))
      pointer(ptr,buffer )
      integer :: buffer_pos
      character(len=8) :: text
      integer :: buffer_recv_size, shift
      integer :: rank_x, rank_y, ind_x, ind_y, cur_rank
      integer :: nsend_x, nsend_y, nrecv_x, nrecv_y, outunit

      outunit = stdout()
      update_flags = XUPDATE+YUPDATE   !default
      if( PRESENT(flags) ) then
          update_flags = flags
          ! The following test is so that SCALAR_PAIR can be used alone with the
          ! same default update pattern as without.
          if (BTEST(update_flags,SCALAR_BIT)) then
            if (.NOT.(BTEST(update_flags,WEST) .OR. BTEST(update_flags,EAST) &
                 .OR. BTEST(update_flags,NORTH) .OR. BTEST(update_flags,SOUTH))) &
              update_flags = update_flags + XUPDATE+YUPDATE   !default with SCALAR_PAIR
          end if
      end if

      if( BTEST(update_flags,NORTH) .AND. BTEST(domain%fold,NORTH) .AND. BTEST(gridtype,SOUTH) ) &
           call mpp_error( FATAL, 'MPP_DO_UPDATE_V: Incompatible grid offset and fold.' )

      update_edge_only = BTEST(update_flags, EDGEONLY)
      recv = .false.
      recv(1) = BTEST(update_flags,EAST)
      recv(3) = BTEST(update_flags,SOUTH)
      recv(5) = BTEST(update_flags,WEST)
      recv(7) = BTEST(update_flags,NORTH)
      if( update_edge_only ) then
         if( .NOT. (recv(1) .OR. recv(3) .OR. recv(5) .OR. recv(7)) ) then
            recv(1) = .true.
            recv(3) = .true.
            recv(5) = .true.
            recv(7) = .true.
         endif
      else
         recv(2) = recv(1) .AND. recv(3)
         recv(4) = recv(3) .AND. recv(5)
         recv(6) = recv(5) .AND. recv(7)
         recv(8) = recv(7) .AND. recv(1)
      endif

      send    = recv

      l_size = size(f_addrsx,1)
      nlist = size(domain%list(:))
      ptr = LOC(mpp_domains_stack)

!recv
      nsend_x = update_x%nsend
      nsend_y = update_y%nsend
      nrecv_x = update_x%nrecv
      nrecv_y = update_y%nrecv

      if(debug_message_passing) then
         allocate(msg1(0:nlist-1), msg2(0:nlist-1), msg3(0:nlist-1) )
         msg1 = 0
         msg2 = 0
         msg3 = 0
         cur_rank = get_rank_recv(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y)

         do while (ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y)
            msgsize = 0
            if(cur_rank == rank_x) then
               from_pe = update_x%recv(ind_x)%pe
               do n = 1, update_x%recv(ind_x)%count
                  dir = update_x%recv(ind_x)%dir(n)
                  if(recv(dir)) then
                     is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                     js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_x = ind_x+1
               if(ind_x .LE. nrecv_x) then
                  rank_x = update_x%recv(ind_x)%pe - domain%pe
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = -1
               endif
            endif
            if(cur_rank == rank_y) then
               from_pe = update_y%recv(ind_y)%pe
               do n = 1, update_y%recv(ind_y)%count
                  dir = update_y%recv(ind_y)%dir(n)
                  if(recv(dir)) then
                     is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                     js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_y = ind_y+1
               if(ind_y .LE. nrecv_y) then
                  rank_y = update_y%recv(ind_y)%pe - domain%pe
                  if(rank_y .LE.0) rank_y = rank_y + nlist
               else
                  rank_y = -1
               endif
            endif
            cur_rank = max(rank_x, rank_y)
            m = from_pe-mpp_root_pe()
            msg2(m) = msgsize
         end do

         cur_rank = get_rank_send(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y)
         do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
            msgsize = 0
            if(cur_rank == rank_x) then
               to_pe = update_x%send(ind_x)%pe
               do n = 1, update_x%send(ind_x)%count
                  dir = update_x%send(ind_x)%dir(n)
                  if( send(dir) ) then
                     is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                     js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_x = ind_x+1
               if(ind_x .LE. nsend_x) then
                  rank_x = update_x%send(ind_x)%pe - domain%pe
                  if(rank_x .LT.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
            endif
            if(cur_rank == rank_y) then
               to_pe = update_y%send(ind_y)%pe
               do n = 1, update_y%send(ind_y)%count
                  dir = update_y%send(ind_y)%dir(n)
                  if( send(dir) ) then
                     is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                     js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_y = ind_y+1
               if(ind_y .LE. nsend_y) then
                  rank_y = update_y%send(ind_y)%pe - domain%pe
                  if(rank_y .LT.0) rank_y = rank_y + nlist
               else
                  rank_y = nlist+1
               endif
            endif
            m = to_pe-mpp_root_pe()
            msg3(m) = msgsize
            cur_rank = min(rank_x, rank_y)
         enddo
         call mpp_alltoall(msg3, 1, msg1, 1)
!         call mpp_sync_self(check=EVENT_RECV)
         do m = 0, nlist-1
            if(msg1(m) .NE. msg2(m)) then
               print*, "My pe = ", mpp_pe(), ",domain name =", trim(domain%name), ",from pe=", &
                    domain%list(m)%pe, ":send size = ", msg1(m), ", recv size = ", msg2(m)
               call mpp_error(FATAL, "mpp_do_updateV: mismatch on send and recv size")
            endif
         enddo

!         call mpp_sync_self()
         write(outunit,*)"NOTE from mpp_do_updateV: message sizes are matched between send and recv for domain " &
              //trim(domain%name)
         deallocate(msg1, msg2, msg3)
      endif

      !--- recv
      buffer_pos = 0
      cur_rank = get_rank_recv(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y)
      call mpp_clock_begin(recv_clock)
      do while (ind_x .LE. nrecv_x .OR. ind_y .LE. nrecv_y)
         msgsize = 0
         select case(gridtype)
         case(BGRID_NE, BGRID_SW, AGRID)
            if(cur_rank == rank_x) then
               from_pe = update_x%recv(ind_x)%pe
               do n = 1, update_x%recv(ind_x)%count
                  dir = update_x%recv(ind_x)%dir(n)
                  if(recv(dir)) then
                     tMe = update_x%recv(ind_x)%tileMe(n)
                     is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                     js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               msgsize = msgsize*2
               ind_x = ind_x+1
               ind_y = ind_x
               if(ind_x .LE. nrecv_x) then
                  rank_x = update_x%recv(ind_x)%pe - domain%pe
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = -1
               endif
               rank_y = rank_x
            endif
         case(CGRID_NE, CGRID_SW)
            if(cur_rank == rank_x) then
               from_pe = update_x%recv(ind_x)%pe
               do n = 1, update_x%recv(ind_x)%count
                  dir = update_x%recv(ind_x)%dir(n)
                  if(recv(dir)) then
                     is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                     js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_x = ind_x+1
               if(ind_x .LE. nrecv_x) then
                  rank_x = update_x%recv(ind_x)%pe - domain%pe
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = -1
               endif
            endif
            if(cur_rank == rank_y) then
               from_pe = update_y%recv(ind_y)%pe
               do n = 1, update_y%recv(ind_y)%count
                  dir = update_y%recv(ind_y)%dir(n)
                  if(recv(dir)) then
                     is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                     js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                     msgsize = msgsize + (ie-is+1)*(je-js+1)
                  end if
               end do
               ind_y = ind_y+1
               if(ind_y .LE. nrecv_y) then
                  rank_y = update_y%recv(ind_y)%pe - domain%pe
                  if(rank_y .LE.0) rank_y = rank_y + nlist
               else
                  rank_y = -1
               endif
            endif
         end select
         cur_rank = max(rank_x, rank_y)
         msgsize = msgsize*ke*l_size

         if( msgsize.GT.0 )then
             mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, buffer_pos+msgsize )
             if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
                write( text,'(i8)' )mpp_domains_stack_hwm
                call mpp_error( FATAL, 'MPP_DO_UPDATE_V: mpp_domains_stack overflow, '// &
                     'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
             end if
             call mpp_recv( buffer(buffer_pos+1), glen=msgsize, from_pe=from_pe, block=.false., tag=COMM_TAG_2 )
             buffer_pos = buffer_pos + msgsize
         end if
      end do
      call mpp_clock_end(recv_clock)
      buffer_recv_size = buffer_pos
      send_start_pos = buffer_pos

      !--- send
      cur_rank = get_rank_send(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y)
      nsend = 0
      call mpp_clock_begin(pack_clock)
      do while (ind_x .LE. nsend_x .OR. ind_y .LE. nsend_y)
         pos = buffer_pos
         !--- make sure the domain stack size is big enough
         msgsize = 0
         if(cur_rank == rank_x) then
            do n = 1, update_x%send(ind_x)%count
               dir = update_x%send(ind_x)%dir(n)
               if( send(dir) ) msgsize = msgsize +  update_x%send(ind_x)%msgsize(n)
            enddo
         endif
         if(cur_rank == rank_y) then
            do n = 1, update_y%send(ind_y)%count
               dir = update_y%send(ind_y)%dir(n)
               if( send(dir) ) msgsize = msgsize +  update_y%send(ind_y)%msgsize(n)
            enddo
         endif

         if( msgsize.GT.0 )then
            msgsize = msgsize*ke*l_size
            mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, pos+msgsize )
            if( mpp_domains_stack_hwm.GT.mpp_domains_stack_size )then
               write( text,'(i8)' )mpp_domains_stack_hwm
               call mpp_error( FATAL, 'MPP_DO_UPDATE_V: mpp_domains_stack overflow, ' // &
                    'call mpp_domains_set_stack_size('//trim(text)//') from all PEs.')
            end if
         end if
         select case( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            if(cur_rank == rank_x) then
               to_pe = update_x%send(ind_x)%pe
               do n = 1, update_x%send(ind_x)%count
                  dir = update_x%send(ind_x)%dir(n)
                  if( send(dir) ) then
                     tMe = update_x%send(ind_x)%tileMe(n)

                     is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                     js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                     select case( update_x%send(ind_x)%rotation(n) )
                     case(ZERO)
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l,tMe)
                           ptr_fieldy = f_addrsy(l,tMe)
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 2
                                    buffer(pos-1) = fieldx(i,j,k)
                                    buffer(pos)   = fieldy(i,j,k)
                                 end do
                              end do
                           end do
                        end do
                     case( MINUS_NINETY )
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 2
                                       buffer(pos-1) = fieldy(i,j,k)
                                       buffer(pos)   = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 2
                                       buffer(pos-1) = -fieldy(i,j,k)
                                       buffer(pos)   =  fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     case( NINETY )
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 2
                                       buffer(pos-1) = fieldy(i,j,k)
                                       buffer(pos)   = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 2
                                       buffer(pos-1) = fieldy(i,j,k)
                                       buffer(pos)   = -fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     case( ONE_HUNDRED_EIGHTY )
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 2
                                       buffer(pos-1) =  fieldx(i,j,k)
                                       buffer(pos)   =  fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l,tMe)
                              ptr_fieldy = f_addrsy(l,tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 2
                                       buffer(pos-1) =  -fieldx(i,j,k)
                                       buffer(pos)   =  -fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     end select ! select case( rotation(n) )
                  end if ! if( send(dir) )
               end do ! do n = 1, update_x%send(ind_x)%count
               ind_x = ind_x+1
               ind_y = ind_x
               if(ind_x .LE. nsend_x) then
                  rank_x = update_x%send(ind_x)%pe - domain%pe
                  if(rank_x .LT.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
               rank_y = rank_x
            endif
         case(CGRID_NE, CGRID_SW)
            if(cur_rank == rank_x) then
               to_pe = update_x%send(ind_x)%pe
               do n = 1, update_x%send(ind_x)%count
                  dir = update_x%send(ind_x)%dir(n)
                  if( send(dir) ) then
                     tMe = update_x%send(ind_x)%tileMe(n)
                     is = update_x%send(ind_x)%is(n); ie = update_x%send(ind_x)%ie(n)
                     js = update_x%send(ind_x)%js(n); je = update_x%send(ind_x)%je(n)
                     select case( update_x%send(ind_x)%rotation(n) )
                     case(ZERO)
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                    buffer(pos) = fieldx(i,j,k)
                                 end do
                              end do
                           end do
                        end do
                     case(MINUS_NINETY)
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 1
                                       buffer(pos) = fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do i = is, ie
                                    do j = je, js, -1
                                       pos = pos + 1
                                       buffer(pos) = -fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     case(NINETY)
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1, ke
                              do i = ie, is, -1
                                 do j = js, je
                                    pos = pos + 1
                                    buffer(pos) = fieldy(i,j,k)
                                 end do
                              end do
                           end do
                        end do
                     case(ONE_HUNDRED_EIGHTY)
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 1
                                       buffer(pos) = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 1
                                       buffer(pos) = -fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     end select
                  end if
               end do
               ind_x = ind_x+1
               if(ind_x .LE. nsend_x) then
                  rank_x = update_x%send(ind_x)%pe - domain%pe
                  if(rank_x .LT.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
            endif
            if(cur_rank == rank_y) then
               to_pe = update_y%send(ind_y)%pe
               do n = 1, update_y%send(ind_y)%count
                  dir = update_y%send(ind_y)%dir(n)
                  if( send(dir) ) then
                     tMe = update_y%send(ind_y)%tileMe(n)
                     is = update_y%send(ind_y)%is(n); ie = update_y%send(ind_y)%ie(n)
                     js = update_y%send(ind_y)%js(n); je = update_y%send(ind_y)%je(n)
                     select case( update_y%send(ind_y)%rotation(n) )
                     case(ZERO)
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
                              do j = js, je
                                 do i = is, ie
                                    pos = pos + 1
                                    buffer(pos) = fieldy(i,j,k)
                                 end do
                              end do
                           end do
                        end do
                     case(MINUS_NINETY)
                        do l=1,l_size  ! loop over number of fields
                           ptr_fieldx = f_addrsx(l, tMe)
                           ptr_fieldy = f_addrsy(l, tMe)
                           do k = 1,ke
                              do i = is, ie
                                 do j = je, js, -1
                                    pos = pos + 1
                                    buffer(pos) = fieldx(i,j,k)
                                 end do
                              end do
                           end do
                        end do
                     case(NINETY)
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 1
                                       buffer(pos) = fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do i = ie, is, -1
                                    do j = js, je
                                       pos = pos + 1
                                       buffer(pos) = -fieldx(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     case(ONE_HUNDRED_EIGHTY)
                        if( BTEST(update_flags,SCALAR_BIT) ) then
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 1
                                       buffer(pos) = fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        else
                           do l=1,l_size  ! loop over number of fields
                              ptr_fieldx = f_addrsx(l, tMe)
                              ptr_fieldy = f_addrsy(l, tMe)
                              do k = 1,ke
                                 do j = je, js, -1
                                    do i = ie, is, -1
                                       pos = pos + 1
                                       buffer(pos) = -fieldy(i,j,k)
                                    end do
                                 end do
                              end do
                           end do
                        end if
                     end select
                  endif
               enddo
               ind_y = ind_y+1
               if(ind_y .LE. nsend_y) then
                  rank_y = update_y%send(ind_y)%pe - domain%pe
                  if(rank_y .LT.0) rank_y = rank_y + nlist
               else
                  rank_y = nlist+1
               endif
            endif
         end select
         cur_rank = min(rank_x, rank_y)
         nsend = nsend + 1
         send_pe(nsend) = to_pe
         send_msgsize(nsend) = pos - buffer_pos
         buffer_pos = pos
      end do

      buffer_pos = send_start_pos
      call mpp_clock_end(pack_clock)
      call mpp_clock_begin(send_clock)
      do m = 1, nsend
         msgsize = send_msgsize(m)
         if( msgsize.GT.0 )then
            call mpp_send( buffer(buffer_pos+1), plen=msgsize, to_pe=send_pe(m), tag=COMM_TAG_2 )
            buffer_pos = buffer_pos + msgsize
         end if
      end do
      call mpp_clock_end(send_clock)

!unpack recv
!unpack halos in reverse order
      call mpp_clock_begin(wait_clock)
      call mpp_sync_self(check=EVENT_RECV)
      call mpp_clock_end(wait_clock)
      buffer_pos = buffer_recv_size
      cur_rank = get_rank_unpack(domain, update_x, update_y, rank_x, rank_y, ind_x, ind_y)

      call mpp_clock_begin(unpk_clock)
      do while (ind_x > 0 .OR. ind_y > 0)
         pos = buffer_pos
         select case ( gridtype )
         case(BGRID_NE, BGRID_SW, AGRID)
            if(cur_rank == rank_x) then
               do n = update_x%recv(ind_x)%count, 1, -1
                  dir = update_x%recv(ind_x)%dir(n)
                  if( recv(dir) ) then
                     tMe = update_x%recv(ind_x)%tileMe(n)
                     is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                     js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                     msgsize = (ie-is+1)*(je-js+1)*ke*2*l_size
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do l=1, l_size  ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 2
                                 fieldx(i,j,k) = buffer(pos-1)
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end do
                  end if ! end if( recv(dir) )
               end do  ! do dir=8,1,-1
               ind_x = ind_x-1
               ind_y = ind_x
               if(ind_x .GT. 0) then
                  rank_x = update_x%recv(ind_x)%pe - domain%pe
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
               rank_y = rank_x
            endif
         case(CGRID_NE, CGRID_SW)
            if(cur_rank == rank_y) then
               do n = update_y%recv(ind_y)%count, 1, -1
                  dir = update_y%recv(ind_y)%dir(n)
                  if( recv(dir) ) then
                     tMe = update_y%recv(ind_y)%tileMe(n)
                     is = update_y%recv(ind_y)%is(n); ie = update_y%recv(ind_y)%ie(n)
                     js = update_y%recv(ind_y)%js(n); je = update_y%recv(ind_y)%je(n)
                     msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do l=1,l_size  ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldy(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end do
                  end if
               end do
               ind_y = ind_y-1
               if(ind_y .GT. 0) then
                  rank_y = update_y%recv(ind_y)%pe - domain%pe
                  if(rank_y .LE.0) rank_y = rank_y + nlist
               else
                  rank_y = nlist+1
               endif
            endif
            if(cur_rank == rank_x) then
               do n = update_x%recv(ind_x)%count, 1, -1
                  dir = update_x%recv(ind_x)%dir(n)
                  if( recv(dir) ) then
                     tMe = update_x%recv(ind_x)%tileMe(n)
                     is = update_x%recv(ind_x)%is(n); ie = update_x%recv(ind_x)%ie(n)
                     js = update_x%recv(ind_x)%js(n); je = update_x%recv(ind_x)%je(n)
                     msgsize = (ie-is+1)*(je-js+1)*ke*l_size
                     pos = buffer_pos - msgsize
                     buffer_pos = pos
                     do l=1,l_size  ! loop over number of fields
                        ptr_fieldx = f_addrsx(l, tMe)
                        ptr_fieldy = f_addrsy(l, tMe)
                        do k = 1,ke
                           do j = js, je
                              do i = is, ie
                                 pos = pos + 1
                                 fieldx(i,j,k) = buffer(pos)
                              end do
                           end do
                        end do
                     end do
                  end if
               end do
               ind_x = ind_x-1
               if(ind_x .GT. 0) then
                  rank_x = update_x%recv(ind_x)%pe - domain%pe
                  if(rank_x .LE.0) rank_x = rank_x + nlist
               else
                  rank_x = nlist+1
               endif
            endif
         end select
         cur_rank = min(rank_x, rank_y)
      end do
      call mpp_clock_end(unpk_clock)

     ! ---northern boundary fold
      shift = 0
      if(domain%symmetry) shift = 1
      if( BTEST(domain%fold,NORTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then
         isd = domain%x(1)%compute%begin - update_x%whalo;
         ied = domain%x(1)%compute%end + update_x%ehalo;
         jsd = domain%y(1)%compute%begin - update_y%shalo;
         jed = domain%y(1)%compute%end + update_y%nhalo;

         j = domain%y(1)%global%end+shift
         if( jsd .LE. j .AND. j.LE.jed+shift )then !fold is within domain
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
               j  = domain%y(1)%global%end+shift
               is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
               if( .NOT. domain%symmetry ) is = is - 1
               do i = is ,ie, midpoint
                  if( isd.LE.i .AND. i.LE. ied+shift )then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from
            ! other half may have the wrong sign
            !off west edge, when update north or west direction
            j = domain%y(1)%global%end+shift
            if ( recv(7) .OR. recv(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  if(domain%symmetry) then
                     is = domain%x(1)%global%begin
                  else
                     is = domain%x(1)%global%begin - 1
                  end if
                  if( is.GT.isd )then
                     if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift  ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do i = isd,is-1
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.isd )then
                     if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-north CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do i = isd,is-1
                              fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domain%x(1)%global%end
            if(domain%x(1)%cyclic .AND. is.LT.ied )then
               ie = ied
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  is = is + shift
                  ie = ie + shift
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do i = is, ie
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,SOUTH) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---southern
                                                                                               !! boundary fold
         ! NOTE: symmetry is assumed for fold-south boundary
         j = domain%y(1)%global%begin
         if( domain%y(1)%data%begin.LE.j .AND. j.LE.domain%y(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%x(1)%global%begin+domain%x(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               j  = domain%y(1)%global%begin
               is = domain%x(1)%global%begin; ie = domain%x(1)%global%end+shift
               do i = is ,ie, midpoint
                  if( domain%x(1)%data%begin.LE.i .AND. i.LE. domain%x(1)%data%end+shift )then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from
            ! other half may have the wrong sign
            !off west edge, when update north or west direction
            j = domain%y(1)%global%begin
            if ( recv(3) .OR. recv(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.domain%x(1)%data%begin )then

                     if( 2*is-domain%x(1)%data%begin.GT.domain%x(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-south BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do i = domain%x(1)%data%begin,is-1
                              fieldx(i,j,k) = fieldx(2*is-i,j,k)
                              fieldy(i,j,k) = fieldy(2*is-i,j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  is = domain%x(1)%global%begin
                  if( is.GT.domain%x(1)%data%begin )then
                     if( 2*is-domain%x(1)%data%begin-1.GT.domain%x(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-south CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do i = domain%x(1)%data%begin,is-1
                              fieldy(i,j,k) = fieldy(2*is-i-1,j,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off east edge
            is = domain%x(1)%global%end
            if(domain%x(1)%cyclic .AND. is.LT.domain%x(1)%data%end )then
               ie = domain%x(1)%data%end
               is = is + 1
               select case(gridtype)
               case(BGRID_NE)
                  is = is + shift
                  ie = ie + shift
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do i = is,ie
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do i = is, ie
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,WEST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern
                                                                                              !! boundary fold
         ! NOTE: symmetry is assumed for fold-west boundary
         i = domain%x(1)%global%begin
         if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               i  = domain%x(1)%global%begin
               js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
               do j = js ,je, midpoint
                  if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from
            ! other half may have the wrong sign
            !off south edge, when update south or west direction
            i = domain%x(1)%global%begin
            if ( recv(3) .OR. recv(5) ) then
               select case(gridtype)
               case(BGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then

                     if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-west BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              fieldx(i,j,k) = fieldx(i,2*js-j,k)
                              fieldy(i,j,k) = fieldy(i,2*js-j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then
                     if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-west CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off north edge
            js = domain%y(1)%global%end
            if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
               je = domain%y(1)%data%end
               js = js + 1
               select case(gridtype)
               case(BGRID_NE)
                  js = js + shift
                  je = je + shift
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do j = js,je
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     do k = 1,ke
                        do j = js, je
                           fieldx(i,j,k) = -fieldx(i,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      else if( BTEST(domain%fold,EAST) .AND. (.NOT.BTEST(update_flags,SCALAR_BIT)) )then      ! ---eastern
                                                                                              !! boundary fold
         ! NOTE: symmetry is assumed for fold-west boundary
         i = domain%x(1)%global%end+shift
         if( domain%x(1)%data%begin.LE.i .AND. i.LE.domain%x(1)%data%end+shift )then !fold is within domain
            midpoint = (domain%y(1)%global%begin+domain%y(1)%global%end-1+shift)/2
            !poles set to 0: BGRID only
            if( gridtype.EQ.BGRID_NE )then
               i  = domain%x(1)%global%end+shift
               js = domain%y(1)%global%begin; je = domain%y(1)%global%end+shift
               do j = js ,je, midpoint
                  if( domain%y(1)%data%begin.LE.j .AND. j.LE. domain%y(1)%data%end+shift )then
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           fieldx(i,j,k) = 0.
                           fieldy(i,j,k) = 0.
                        end do
                     end do
                  end if
               end do
            endif

            ! the following code code block correct an error where the data in your halo coming from
            ! other half may have the wrong sign
            !off south edge, when update south or west direction
            i = domain%x(1)%global%end+shift
            if ( recv(3) .OR. recv(1) ) then
               select case(gridtype)
               case(BGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then

                     if( 2*js-domain%y(1)%data%begin.GT.domain%y(1)%data%end+shift ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-east BGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        ptr_fieldy = f_addrsy(l, 1)
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              fieldx(i,j,k) = fieldx(i,2*js-j,k)
                              fieldy(i,j,k) = fieldy(i,2*js-j,k)
                           end do
                        end do
                     end do
                  end if
               case(CGRID_NE)
                  js = domain%y(1)%global%begin
                  if( js.GT.domain%y(1)%data%begin )then
                     if( 2*js-domain%y(1)%data%begin-1.GT.domain%y(1)%data%end ) &
                          call mpp_error( FATAL, 'MPP_DO_UPDATE_V: folded-east CGRID_NE west edge ubound error.' )
                     do l=1,l_size
                        ptr_fieldx = f_addrsx(l, 1)
                        do k = 1,ke
                           do j = domain%y(1)%data%begin,js-1
                              fieldx(i,j,k) = fieldx(i, 2*js-j-1,k)
                           end do
                        end do
                     end do
                  end if
               end select
            end if

            !off north edge
            js = domain%y(1)%global%end
            if(domain%y(1)%cyclic .AND. js.LT.domain%y(1)%data%end )then
               je = domain%y(1)%data%end
               js = js + 1
               select case(gridtype)
               case(BGRID_NE)
                  js = js + shift
                  je = je + shift
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     ptr_fieldy = f_addrsy(l, 1)
                     do k = 1,ke
                        do j = js,je
                           fieldx(i,j,k) = -fieldx(i,j,k)
                           fieldy(i,j,k) = -fieldy(i,j,k)
                        end do
                     end do
                  end do
               case(CGRID_NE)
                  do l=1,l_size
                     ptr_fieldx = f_addrsx(l, 1)
                     do k = 1,ke
                        do j = js, je
                           fieldx(i,j,k) = -fieldx(i,j,k)
                        end do
                     end do
                  end do
               end select
            end if
         end if
      end if

      call mpp_clock_begin(wait_clock)
      call mpp_sync_self( )
      call mpp_clock_end(wait_clock)

      return

    end subroutine MPP_DO_UPDATE_3D_V_
