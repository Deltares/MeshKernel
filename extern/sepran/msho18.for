      subroutine msho18 ( kelem, nelem, coor, npoint, nipnt, iwork,
     +                    ibuurp, coars )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.4  date 07-12-2010 New method Laplace
!        version  2.3  date 27-08-2009 Enlarge stopvalue of repositioning
!        version  2.2  date 24-03-2009 Check eventually correctness Laplacian
!                                      repositioning (Plaxis)
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
!     Perform a Laplacian repositioning of the nodes
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     repositioning
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
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cmcdpr'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*), coars
      integer nelem, npoint, nipnt, iwork(*), ibuurp(*), kelem(3,nelem)

!     coars   i     coarseness, reference length
!     coor   i/o    coordinate array
!     ibuurp  i     array with nodal point numbers of neighbour points
!                   for each point
!     iwork   i     helparray with startaddresses of IBUURP
!     kelem   i     element array
!     nelem   i     number of elements in mesh
!     nipnt   i     number of boundary points
!     npoint  i     total number of points in the mesh
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision afstnd, xz, yz, xh, yh, xn, yn, af,
     +                 surf, surf1, surfre
      integer jr, ikn, jstart, jeind, nbuur, i, jkn,
     +        iallow, i1, i2, i3
      logical debug, finished

!     af            extreme value for distances
!     afstnd        distance between old and new position of point
!     debug         if true debug statements are carried out otherwise
!                   they are not
!     finished      logical to indicate whether loop is finished or not
!     i             loop variable
!     iallow        allowance indicator
!     ikn           loop variable
!     jeind         endaddress neighbours
!     jkn           local node number
!     jr            loop variable
!     jstart        startaddress neighbours
!     nbuur         number of neighbours
!     surf          surface of all elements containing specified node
!     surf1         surface one special element
!     surfre        reflection value for surface
!     xh            old x-coordinate
!     xz            new x-coordinate
!     yh            old y-coordinate
!     yz            new y-coordinate
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!     Compute stop value for coarseness
!     while not finished do
!         compute max repositioning distance
!         check whether distance small enough
!         if so finished else reposition again
!     end while
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Compute stop value for coarseness

      debug = .false.
      coars = 0.02 * coars
      finished = .false.
      jr = 1
      do while ( .not. finished )
         jr = jr + 1
         afstnd = 0d0
         do ikn = nipnt + 1 , npoint

!        --- Node ikn, compute surface

            surf  = 0d0
            surf1 = 0d0

!           --- Start values for ikn:

            xh = coor(1,ikn)
            yh = coor(2,ikn)

!           --- Add surf of all elements ikn is node of:

            do i = 1, nelem
               i1 = kelem( 1,i )
               i2 = kelem( 2,i )
               i3 = kelem( 3,i )

!              --- Check triangle and add surface:

               if ( i1==ikn.or.i2==ikn.or.i3==ikn ) then
                    call msho19(coor,i1,i2,i3,surf1)
                    surf = surf + surf1
               end if
            end do

!           --- Set reference value for surface

            surfre = surf

!           --- Run through all neighbours and compute (xz,yz):

            jstart = 0
            if ( ikn/=1 ) jstart = iwork( ikn - 1 )
            jstart = jstart + 1
            jeind  = iwork( ikn )
            if ( jstart<=jeind ) then

!           --- jstart<=jeind, there are neighbours, initialize:

               xz  =  0d0
               yz  =  0d0

!              --- Set number of neighbours:

               nbuur = jeind - jstart + 1

!              --- Run through all neighbours and compute barycenter:

               do i = jstart , jeind
                  jkn = ibuurp(i)
                  if ( jkn>0 ) then
                     xz = xz + coor(1,jkn)
                     yz = yz + coor(2,jkn)
                  else
                     nbuur = nbuur - 1
                  end if
               end do

!              --- New(?) coordinates

               xz = xz / nbuur
               yz = yz / nbuur

!              --- Two possibilities: No relaxation or Overrelaxation:

               xn = xh + 1.62 * ( xz - xh )
               yn = yh + 1.62 * ( yz - yh )

!              --- Try overrelaxation: Set coordinates:

               coor(1,ikn) = xn
               coor(2,ikn) = yn

!              --- Check surface:

               surf = 0d0

!              --- Run through all elements:

               do i = 1, nelem
                  i1 = kelem( 1,i )
                  i2 = kelem( 2,i )
                  i3 = kelem( 3,i )

!                 --- Check triangle and add surface:

                  if ( i1==ikn.or.i2==ikn.or.i3==ikn ) then
                       call msho19(coor,i1,i2,i3,surf1)
                       surf = surf + abs(surf1)
                  end if
               end do
               iallow = 0
               if ( abs(surf - surfre)<100 * epsmac ) iallow = 2

!              --- Allowed? Then ready:

               if ( iallow/=2 ) then

!              --- Not allowed: Check for barycenter:

                  coor(1,ikn) = xz
                  coor(2,ikn) = yz

!                 --- Check surface:

                  surf = 0d0

!                 --- Run through all elements:

                  do i = 1, nelem
                     i1 = kelem( 1,i )
                     i2 = kelem( 2,i )
                     i3 = kelem( 3,i )

!                    --- Check triangle and add surface:

                     if ( i1==ikn.or.i2==ikn.or.i3==ikn ) then
                        call msho19(coor,i1,i2,i3,surf1)
                        surf = surf + surf1
                     end if
                  end do
                  if ( abs(surf - surfre)<100 * epsmac ) iallow = 1

!                 --- Allowed? Then ready:

                  if ( iallow==0 ) then

!                 --- Nothing allowed use old values:

                     coor(1,ikn) = xh
                     coor(2,ikn) = yh
                  end if
               end if

!              --- Check whether "new" point is different from "old" point:

               xz = coor(1,ikn)
               yz = coor(2,ikn)

!              --- Compute distance:

               af = (xh-xz)*(xh-xz)+(yh-yz)*(yh-yz)  ! squared distance
               if ( af>afstnd ) afstnd=af
            end if  ! ( jstart<=jeind )
         end do

!        --- Check whether wanted accuracy is found:

         if ( debug )
     +      write(irefwr,*) 'jr =',jr,'afstand=',afstnd,'coars=',coars
         if ( jr>2  .and. sqrt(afstnd)<coars ) finished = .true.
         if ( jr==20 ) finished = .true.
      end do ! while ( .not. finished )
      end
