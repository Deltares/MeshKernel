      subroutine msho24 ( kstapl, kstap, coor, i1, i2, icheck )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.1  date 23-02-2010 Adjustment for plaxis pnts, use
!                                      MSHO75 in stead of MSHG03
!        version  3.0  date 14-02-2005 Update
!        version  2.2  date 25-06-1999 For savety leave routine if second point
!                                      is plaxis
!        version  2.1  date 26-05-1997 Check coincidence of double
!                                      Plaxis points
!        version  2.0  date 21-02-1994  New norms
!        version  1.1  date 14-10-1990  use routine MSHG03
!        version  1.0  date 17-04-1989
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
!     Check whether line i1 - i2 has a common point with one of the
!     boundary lines
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     crossing_point
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
!      include 'SPcommon/cmcdpr'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*)
      integer kstap, kstapl(2,kstap), i1, i2, icheck

!     coor      i   coordinate array
!     i1        i   first node of line to be considered
!     i2        i   second node of line to be considered
!     icheck    o   if a line crosses the new element then
!                   icheck # 0  and equal to number of line
!                   if no point is "inside" the new element then
!                   icheck = -1 else icheck is nodenumber of point
!                   inside
!     kstap     i   number of lines in kstapl
!     kstapl    i   actual array with boundary lines
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision x1, y1, x2, y2, xmin, xmax, ymin, ymax, xa, ya,
     +                 xb, yb, xlm, ylm, eps, dis
      integer ia, ib, ih, il

!     dis        distance
!     eps        accuracy
!     ia         node number
!     ib         node number
!     ih         indicator
!     il         loop variable
!     x1         x-coordinate
!     x2         x-coordinate
!     xa         x-coordinate
!     xb         x-coordinate
!     xlm        extreme x-value
!     xmax       extreme x-value
!     xmin       extreme x-value
!     y1         y-coordinate
!     y2         y-coordinate
!     ya         y-coordinate
!     yb         y-coordinate
!     ylm        extreme y-value
!     ymax       extreme y-value
!     ymin       extreme y-value
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS    Resets old name of previous subroutine of higher level
!     EROPEN    Produces concatenated name of local subroutine
!     MSHO75    check whether two line segments have a common point
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
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen( 'msho24' )

      eps    =  10 * epsmac

      icheck = -1

      x1 = coor(1,i1)
      y1 = coor(2,i1)

      x2 = coor(1,i2)
      y2 = coor(2,i2)

      xmin = min(x1,x2)
      xmax = max(x1,x2)
      ymin = min(y1,y2)
      ymax = max(y1,y2)

!     --- Check whether line i1 - i2 has point in common
!         with old boundary lines.

      do il = 1 , kstap

         ia = kstapl(1,il)
         ib = kstapl(2,il)

!        --- First check whether node numbers of boundary element
!            are either equal i1 or i2

         if ( ia==i1 .or. ia==i2 ) goto 200
         if ( ib==i1 .or. ib==i2 ) goto 200

         xa = coor(1,ia)
         ya = coor(2,ia)

         xb = coor(1,ib)
         yb = coor(2,ib)

!        --- Check whether coordinates of nodes are exactly the same, i.e.
!            whether it concerns Plaxis points
!            Point ia and i1

         dis = (x1-xa)*(x1-xa) + (y1-ya)*(y1-ya)

         if ( dis<eps ) goto 200

!        --- Point ia and i2

         dis = (x2-xa)*(x2-xa) + (y2-ya)*(y2-ya)

         if ( dis<eps ) then

!        --- plaxis point ?
!            for savety jump out of routine

            icheck = il
            go to 300

         end if

!        --- Point ib and i1

         dis = (x1-xb)*(x1-xb) + (y1-yb)*(y1-yb)

         if ( dis<eps ) goto 200

!        --- Point ib and i2

         dis = (x2-xb)*(x2-xb) + (y2-yb)*(y2-yb)

         if ( dis<eps ) then

!        --- plaxis point ?
!            for savety jump out of routine

            icheck = il
            go to 300

         end if

         xlm = min(xa,xb)
         if ( xlm>xmax ) goto 200
         xlm = max(xa,xb)
         if ( xlm<xmin ) goto 200
         ylm = min(ya,yb)
         if ( ylm>ymax ) goto 200
         ylm = max(ya,yb)
         if ( ylm<ymin ) goto 200

!        --- Check

         call msho75( xa, ya, xb, yb, x1, y1, x2, y2, ih )
         if ( ih==0 ) then

!        --- Common point found, ready

              icheck = il
              goto 300
         end if

200   end do

      icheck = 0

300   call erclos( 'msho24' )

      end
