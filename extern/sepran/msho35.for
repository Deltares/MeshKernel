      subroutine msho35 ( npoint, coor, xstart, ystart, dismax,
     +                    coar, ncoar, icube, nx, kelem, nelem,
     +                    kstapl, kstap, itri, isurnr, userpoints,
     +                    kbndpt, nbndpt )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  1.2  date 04-07-2008 Extra internal lines added
!        version  1.1  date 18-02-2005 Accuracy added
!        version  1.0  date 11-01-2005
!
!   copyright (c) 2005-2008  "Ingenieursbureau SEPRA"
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
!     Place special fixed nodes and elements around these special nodes
!     and place extra internal lines nodes in array kbndpt.
!
! **********************************************************************
!
!                       KEYWORDS
!
!     2d
!     mesh
!     mesh_generation
!     node
!     plaxis
!     user_point
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use msherror
      use mshdummymethods
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cconst'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*), coar(3,*), xstart, ystart, dismax
      integer npoint, ncoar, kelem(3,*), icube(*), isurnr,
     +        nelem, nx, kstap, kstapl(2,*), nbndpt,
     +        kbndpt(*), itri(*), userpoints(*)

!     coar           i    array with coordinates and coarseness of fixed points
!     coor           i    coordinate array
!     dismax         i    help distance for blockdetermination
!     icube          i    array with for each node the number of the
!                         quadrilateral that node belongs to
!     isurnr         i    Surface sequence number
!     itri          i/o   help array containing information indicating whether
!                         a point is in the (actual) boundary or not
!     kelem         i/o   element array (created triangles)
!     kstap         i/o   length of boundary elements array
!     kstapl        i/o   array with nodes of boundary elements
!     ncoar          i    number of fixed points
!     nelem         i/o   number of created elements in interior in total
!     npoint         i    number of points in line elements of contour
!     nx             i    number of quadrilaterals in x-direction considered
!     userpoints     i    integer array of length ncoar containing the
!                         user point numbers of the extra points
!     xstart         i    x-coordinate left under enveloping quadrilateral
!     ystart         i    y-coordinate left under enveloping quadrilateral
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision dist, xp, yp, xpoin, ypoin, eps
      integer i, j, ja, n1, n2, nc, npoin, npoine
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     dist           help variable to store local meshsize
!     eps            accuracy parameter
!     i              loop variable
!     j              loop variable
!     ja             indicator whether point is inside area or not
!     n1             help variable x-number of quadrilateral
!     n2             help variable y-number of quadrilateral
!     nc             quadrilateral number
!     npoin          number of points at start
!     npoine         number of points npoin plus number of coars points
!     xp             x-coordinate
!     xpoin          x-coordinate new point to be placed
!     yp             y-coordinate
!     ypoin          y-coordinate new point to be placed
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERREAL         Put real in error message
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSHO07         Check whether point belongs to given area
!     PRININ1        print 2d integer array
!     PRINRL1        Print 2d real vector
! **********************************************************************
!
!                       I/O
!
!    None
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    1344   User point not inside surface
!    1346   User point has incorrect coarseness inside surface
! **********************************************************************
!
!                       PSEUDO CODE
!
!    Run through all (ncoar) special nodes
!
!    Check whether triangles with given coarseness are permitted.
!    If not: stop, error
!    Else: place nodes and elements and adjust arrys.
!
! **********************************************************************
!
!                       DATA STATEMENTS
!
!    None
! ======================================================================
!
      call eropen ( 'msho35' )
      debug = .false.

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msho35'
         write(irefwr,100) 'kstap', kstap
         write(irefwr,110) 'xstart, ystart', xstart, ystart
100      format ( a, 1x, (10i6) )
110      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Set help values

      npoin  = npoint
      npoint = npoint + ncoar

!     --- Run through all coarse nodes:

      do i = 1, ncoar

!     --- Set helpvariable

         npoine = npoin + ncoar + (i-1) * 6

!        --- Check position of nodes:

         dist = coar(3,i)

         xp = coar(1,i)
         yp = coar(2,i)

         eps = 1d-10 * dist

         call msho07 ( xp+eps, yp+eps, xstart, coor, kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted point that is not inside region

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call erreal ( xp, 1 )
            call erreal ( yp, 2 )
            call errsub ( 1344, 2, 2, 0 )
            go to 1000

         end if

         nbndpt = nbndpt + 1

         kbndpt( nbndpt ) = npoin + i

         coor( 1, npoin+i ) = xp
         coor( 2, npoin+i ) = yp

         n1 = int ( (xp-xstart)/dismax )
         n2 = int ( (yp-ystart)/dismax )

!        --- This point belongs to quadrilateral with number

         nc = 1 + n1 + n2*nx

         icube( npoin+i ) = nc

!        --- Six other points:

         dist = coar(3,i)
         if ( debug ) write(irefwr,110) 'dist', dist

!        --- Point 1:

         xpoin = xp - 0.40 * dist
         ypoin = yp + 0.67 * dist

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 1, 3 )
            call erreal ( xp, 1 )
            call erreal ( yp, 2 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Point 2:

         xpoin = xp + 0.40 * dist
         ypoin = yp + 0.67 * dist

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 2, 3 )
            call erreal ( dist, 1 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Point 3:

         xpoin = xp - 0.8 * dist
         ypoin = yp

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 3, 3 )
            call erreal ( dist, 1 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Point 4:

         xpoin = xp + 0.8 * dist
         ypoin = yp

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 4, 3 )
            call erreal ( dist, 1 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Point 5:

         xpoin = xp - 0.40 * dist
         ypoin = yp - 0.67 * dist

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 5, 3 )
            call erreal ( dist, 1 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Point 6:

         xpoin = xp + 0.40 * dist
         ypoin = yp - 0.67 * dist

         call msho07 ( xpoin+eps, ypoin+eps, xstart, coor,
     +                 kstapl, kstap, ja )

         if ( ja==0 ) then

!        --- Stop: User submitted wrong coarseness:

            call errint ( userpoints(i), 1 )
            call errint ( isurnr, 2 )
            call errint ( 6, 3 )
            call erreal ( dist, 1 )
            call erreal ( xpoin, 2 )
            call erreal ( ypoin, 3 )
            call errsub ( 1346, 3, 3, 0 )
            go to 1000

         end if

         npoint = npoint + 1

         coor(1,npoint) = xpoin
         coor(2,npoint) = ypoin

         icube( npoint ) = nc

!        --- Add elements:

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 1
         kelem(3,nelem) = npoine + 3

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 3
         kelem(3,nelem) = npoine + 5

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 5
         kelem(3,nelem) = npoine + 6

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 6
         kelem(3,nelem) = npoine + 4

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 4
         kelem(3,nelem) = npoine + 2

         nelem          = nelem  + 1
         kelem(1,nelem) = npoin  + i
         kelem(2,nelem) = npoine + 2
         kelem(3,nelem) = npoine + 1

!        --- Adjust boundary array

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 3
         kstapl(2,kstap) = npoine + 1

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 5
         kstapl(2,kstap) = npoine + 3

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 6
         kstapl(2,kstap) = npoine + 5

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 4
         kstapl(2,kstap) = npoine + 6

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 2
         kstapl(2,kstap) = npoine + 4

         kstap           = kstap  + 1
         kstapl(1,kstap) = npoine + 1
         kstapl(2,kstap) = npoine + 2

         do j = 1, 6

            itri(npoine + j ) = 2

         end do

      end do

1000  call erclos ( 'msho35' )

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End msho35'

      end if

      end
