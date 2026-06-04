      subroutine msho28 ( coor, i1, i2, i3, xn, yn, npoint, itri,
     +                    iperm )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.1  date 23-03-2009 Accuracy adjusted using norm
!        version  3.0  date 15-02-2005 Update
!        version  2.2  date 29-09-1998 Correction for double line
!        version  2.1  date 21-02-1997 Adjust for double line and hence
!                                      for "double" points
!        version  2.0  date 22-02-1994 New norms
!        version  1.0  date 07-06-1989
!
!   copyright (c) 1989-2009  "Ingenieursbureau SEPRA"
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
!     Check whether there are boundary points inside the
!     triangle with nodes i1 - i2 - i3.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     check
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpr'
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*), xn, yn
      integer i1, i2, i3, npoint, itri(*), iperm

!     coor        i   coordinate array
!     i1          i   first node
!     i2          i   second node
!     i3          i   last node triangle, if 0 then new coordinates
!                     are given by xn,yn
!     iperm       o   result , IPERM = 1 there are internal points
!                              IPERM = 0 no internal points
!     itri        i   indicator whether point i is in boundary or not
!     npoint      i   number of points in COOR
!     xn          i   coordinate new point if i3 = 0
!     yn          i   coordinate new point if i3 = 0
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision x1, y1, x2, y2, x3, y3, xm, ym, opp, det, eps,
     +                 dist, norm
      integer i

!     det        area of triangle
!     dist       distance
!     eps        accuracy
!     i          loop variable
!     norm       normed value for length
!     opp        area of triangle
!     x1         x-coordinate
!     x2         x-coordinate
!     x3         x-coordinate
!     xm         x-coordinate
!     y1         y-coordinate
!     y2         y-coordinate
!     y3         y-coordinate
!     ym         y-coordinate
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
!     Check for all boundary points whether they are inside
!     If point inside then check whether point is a "double"
!     point
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      iperm = 1
      eps   = 1d-2 * epsmac

!     --- Check whether boundary points are inside triangle

      x1 = coor(1,i1)
      y1 = coor(2,i1)

      x2 = coor(1,i2)
      y2 = coor(2,i2)

      if ( i3==0 ) then
           x3 = xn
           y3 = yn
      else
           x3 = coor(1,i3)
           y3 = coor(2,i3)
      end if

!     --- Compute relative value for length of "vector":
      norm = abs(x1)+abs(x2)+abs(x3)+abs(y1)+abs(y2)+abs(y3)

!     --- Compute area of triangle
      det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)

      do i = 1, npoint
!     --- If not point of triangle and point i in boundary then check
         if (i/=i1.and.i/=i2.and.i/=i3.and.itri(i)/=0) then

            xm = coor(1,i)
            ym = coor(2,i)

            opp = (x1*(y2-ym) + x2*(ym-y1) + xm*(y1-y2))/det
            if ( opp<-eps .or. opp>1d0+eps ) goto 100

            opp = (x2*(y3-ym) + x3*(ym-y2) + xm*(y2-y3))/det
            if ( opp<-eps .or. opp>1d0+eps ) goto 100

            opp = (x3*(y1-ym) + x1*(ym-y3) + xm*(y3-y1))/det
            if ( opp<-eps .or. opp>1d0+eps ) goto 100

!           --- Point i is inside, check whether point coincides
!               with one of the cornerpoints (and hence is a "double"
!               point)

            dist = abs(x1-xm) + abs(y1-ym)

            if ( dist<eps*norm ) goto 100

            dist = abs(x2-xm) + abs(y2-ym)

            if ( dist<eps*norm ) goto 100
!              --- Point i is inside and not a "double" point, return
               goto 500
         end if

100   end do

!     --- None of the points inside

      iperm = 0

500   end
