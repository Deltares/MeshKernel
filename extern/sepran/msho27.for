      subroutine msho27( coor, i1, i2, i3, i4 )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 22-02-1994 New norms
!        version  2.0  date 07-09-1989 New method utilizing a circle-concept
!        version  1.0  date 17-04-1989
!
!   copyright (c) 1989-1999  "Ingenieursbureau SEPRA"
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
!     Determine best diagonal in quadrangle. At input diagonal is
!     line 2-3. Check whether 1-4 is better.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     diagonal
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

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(*)
      integer i1, i2, i3, i4

!     coor       i   coordinate array
!     i1         i   node
!     i2         i   node
!     i3         i   node
!     i4         i   node
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision af1, af2, xa, xb, xc, xd, ya, yb, yc, yd,
     +                 xm, ym, det
      integer j1, j2, j3, j4

!     af1       distance
!     af2       distance
!     det       area of triangle
!     j1        helpnode
!     j2        helpnode
!     j3        helpnode
!     j4        helpnode
!     xa        x-coordinate
!     xb        x-coordinate
!     xc        x-coordinate
!     xd        x-coordinate
!     xm        x-coordinate
!     ya        y-coordinate
!     yb        y-coordinate
!     yc        y-coordinate
!     yd        y-coordinate
!     ym        y-coordinate
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
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- consider quadrangle i1 - i2 - i3 - i4

      j1 = i1
      j2 = i2
      j3 = i3
      j4 = i4

!     --- Determine best diagonal

      xa = coor(2*i1-1)
      xb = coor(2*i2-1)
      xc = coor(2*i3-1)
      xd = coor(2*i4-1)
      ya = coor(2*i1  )
      yb = coor(2*i2  )
      yc = coor(2*i3  )
      yd = coor(2*i4  )
      det = (xa-xb)*(ya-yc) - (xa-xc)*(ya-yb)
      if ( abs(det).lt.10 * epsmac ) then

!     --- Triangle 1-2-3 is not good, take diagonal 1-4

           i1 = j2
           i2 = j4
           i3 = j1
           i4 = j3
      else
         xm = ( ( xa*xa + ya*ya ) * ( yb-yc ) +
     +          ( xb*xb + yb*yb ) * ( yc-ya ) +
     +          ( xc*xc + yc*yc ) * ( ya-yb ) ) / (2*det)
         ym = ( ( xa*xa + ya*ya ) * ( xb-xc ) +
     +          ( xb*xb + yb*yb ) * ( xc-xa ) +
     +          ( xc*xc + yc*yc ) * ( xa-xb ) ) / (-2*det)
         af1 = ( xa-xm ) * ( xa-xm ) + ( ya-ym ) * ( ya-ym )
         af2 = ( xd-xm ) * ( xd-xm ) + ( yd-ym ) * ( yd-ym )
         if ( af2.lt.af1 ) then

!        --- Change diagonal : (4-1 instead of 2-3)

              i1 = j2
              i2 = j4
              i3 = j1
              i4 = j3
         end if
      end if
      end
