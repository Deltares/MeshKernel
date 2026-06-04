      subroutine msho29( kelem, nelem, npoint, iwork, ibuurp,
     +                   nipnt, coor, icancel )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.4  date 04-02-2010 Check for correct use of both
!                                      cases ( 3 and 4 neighbours)
!        version  2.3  date 24-03-2009 Use only case 3 neighbours
!        version  2.2  date 11-01-2005 Adjust kelem (2dim array)
!        version  2.1  date 03-10-2000 Declaration of coor
!        version  2.0  date 22-02-1994 New norms
!        version  1.0  date 09-06-1989
!
!   copyright (c) 1989-2010  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard copy or soft copy granted only by written license
!   obtained from "Ingenieursbureau SEPRA".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system ( e.g., in memory, disk, or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
! **********************************************************************
!
!                       DESCRIPTION
!
!     Determine and remove internal nodal points with either three or
!     four neighbour nodes.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     point_removal
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kelem(3,*), nelem, npoint, iwork(*), ibuurp(*),
     +        nipnt, icancel
      double precision coor(*)

!     coor     i/o   coordinate array
!     ibuurp    i    array with nodal point numbers of neighbour points
!                    for each point
!     icancel  i/o   counter
!     iwork     i    helparray with start addresses of neighbour points
!     kelem    i/o   element array
!     nelem    i/o   number of elements in mesh
!     nipnt     i    number of boundary points
!     npoint    i    total number of points in the mesh
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision opp
      integer i, j, np, nel, ia, ib, ic, id, iex, i1, i2, i3

!     i          loop variable
!     i1         node triangle
!     i2         node triangle
!     i3         node triangle
!     ia         neighbour node
!     ib         neighbour node
!     ic         neighbour node
!     id         neighbour node
!     iex        indicator
!     j          loop variable
!     nel        element number
!     np         number of neighbours of point considered
!     opp        area of triangle
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     MSHO19   compute area of triangle
!     MSHO27   determine best diagonal in quadrilateral
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!     run through all points and cancel the point if it has only three
!     neighbours
!     if point has four neighbours then cancel the four triangles
!     and make two new triangles
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Run through all non boundary points and cancel points with only
!         three or four neighbours

      icancel = 0

      do i = nipnt + 1 , npoint

!     --- Determine number of neighbours

          np = iwork(i) - iwork(i-1)

          if ( np==3 ) then
!         --- Neighbours

               ia = ibuurp(iwork(i-1)+1)
               ib = ibuurp(iwork(i-1)+2)
               ic = ibuurp(iwork(i-1)+3)
!              --- Cancel this point

               icancel = icancel + 1

               nel = 0

               do j = 1 , nelem

                   i1 = kelem(1,j)
                   i2 = kelem(2,j)
                   i3 = kelem(3,j)

                   if ( i1.ne.i .and. i2.ne.i .and. i3.ne.i ) then
!                  --- Triangle has to be kept

                        nel = nel + 1

                        kelem(1,nel) = i1
                        kelem(2,nel) = i2
                        kelem(3,nel) = i3

                   end if

               end do

               nel = nel + 1

               call msho19( coor, ia, ib, ic, opp )

               if ( opp.gt.0 ) then
                    kelem(1,nel) = ia
                    kelem(2,nel) = ib
                    kelem(3,nel) = ic
               else
                    kelem(1,nel) = ia
                    kelem(2,nel) = ic
                    kelem(3,nel) = ib
               end if

               nelem = nel

          else if ( np==4 ) then
!         --- Start considering the four neighbours
               ia = ibuurp(iwork(i-1)+1)
               ib = ibuurp(iwork(i-1)+2)
               ic = ibuurp(iwork(i-1)+3)
               id = ibuurp(iwork(i-1)+4)
!              --- Check existence of neighbouring points
               if ( ia>nipnt ) then
                    iex = iwork(ia)-iwork(ia-1)
                    if ( iex==3 .or. ia<i .and. iex==4 ) goto 200
               end if

               if ( ib>nipnt ) then
                    iex = iwork(ib)-iwork(ib-1)
                    if ( iex==3 .or. ib<i .and. iex==4 ) goto 200
               end if

               if ( ic>nipnt ) then
                    iex = iwork(ic)-iwork(ic-1)
                    if ( iex==3 .or. ic<i .and. iex==4 ) goto 200
               end if

               if ( id>nipnt ) then
                    iex = iwork(id)-iwork(id-1)
                    if ( iex==3 .or. id<i .and. iex==4 ) goto 200
               end if
!              --- All neighbours will remain during this cycle, so cancel
!                  point i
               ia = 0
               ib = 0
               ic = 0
               id = 0

               icancel = icancel + 1000

!              --- Cancel the four triangles of this point
               nel = 0
               iex = 0

               do j = 1 , nelem

                  i1 = kelem(1,j)
                  i2 = kelem(2,j)
                  i3 = kelem(3,j)

                  if ( i1/=i .and. i2/=i .and. i3/=i ) then
!                 --- Triangle has no relation, keep it:
                     nel = nel + 1

                     kelem(1,nel) = i1
                     kelem(2,nel) = i2
                     kelem(3,nel) = i3
                  else
                     if ( iex==0 ) then
                        if ( i1==i ) then
                           ia  = i2
                           ib  = i3
                        else if ( i2==i ) then
                           ia = i3
                           ib = i1
                        else if ( i3==i ) then
                           ia = i1
                           ib = i2
                        end if
                        iex = 1
                     else if ( iex==1 ) then
                        if ( i1==i ) then
                           if ( i2/=ib .and. i3/=ia ) then
                              ic  = i2
                              id  = i3
                              iex = iex + 1
                            end if
                        else if ( i2==i ) then
                           if ( i1/=ia .and. i3/=ib ) then
                              ic  = i3
                              id  = i1
                              iex = iex + 1
                           end if
                        else if ( i3==i ) then
                           if ( i1/=ib .and. i2/=ia ) then
                              ic  = i1
                              id  = i2
                              iex = iex + 1
                           end if
                        end if
                     end if
                  end if
               end do

!              --- Create two new triangles

               if ( ia==0 .or. ib==0 .or. ic==0 .or. id==0 ) then
                  write(irefwr,*) 'Unexpected situation in msho29'
                  write(irefwr,*) 'execution stopped'
                  write(irefwr,*) 'If your input is correct'
                  write(irefwr,*) 'Please inform SEPRA'
!AvD                  call instop
               end if

               call msho27(coor,ia,ib,id,ic)

               nel = nel + 1

               kelem(1,nel) = ia
               kelem(2,nel) = ib
               kelem(3,nel) = id

               nel = nel + 1

               kelem(1,nel) = ib
               kelem(2,nel) = ic
               kelem(3,nel) = id

               nelem = nel

          end if

200   end do

      end
