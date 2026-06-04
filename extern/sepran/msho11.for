      subroutine msho11 ( i1, i2, i3, coor, angle )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 11-02-2005 Alternative accuracy treatment
!        version  2.1  date 28-02-1997 Improvements
!        version  2.0  date 14-02-1994 New norms
!        version  1.0  date 13-04-1989
!
!   copyright (c) 1989-2005  "Ingenieursbureau SEPRA"
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
!     Compute the angle between line i1 - i2  and line i2 - i3 using the
!     dot product
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     angle
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
      integer i1, i2, i3
      double precision coor(2,*), angle

!     angle      o   cosine of angle between lines 1-2 and 2-3
!     coor       i   coordinate array
!     i1         i   first node
!     i2         i   second node
!     i3         i   third node
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision x1, y1, x2, y2, x3, y3, e1, e2, e3, e4, eleng,
     +                 h1, h2, hnorm

!     e1        first component direction line 3-2
!     e2        second component direction line 3-2
!     e3        first component direction line 1-2
!     e4        second component direction line 1-2
!     eleng     scaling parameter
!     h1        helpvariable for coincidences
!     h2        helpvariable for coincidences
!     hnorm     norm of helpvariables
!     x1        x-coordinate node 1
!     x2        x-coordinate node 2
!     x3        x-coordinate node 3
!     y1        y-coordinate node 1
!     y2        y-coordinate node 2
!     y3        y-coordinate node 3
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS   deconcatenate name subroutine from calling string
!     EROPEN   concatenate name subroutine to calling string
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
      call eropen( 'msho11' )
      x1 = coor(1,i1)
      y1 = coor(2,i1)
      x2 = coor(1,i2)
      y2 = coor(2,i2)
      x3 = coor(1,i3)
      y3 = coor(2,i3)
      e1  = x3 - x2
      e2  = y3 - y2
      eleng = sqrt ( e1*e1 + e2*e2 )
      h1  = ( x2 + x3 ) / 2d0
      h2  = ( y2 + y3 ) / 2d0
      hnorm = sqrt( h1*h1 + h2*h2 )
      if ( eleng<1.d2 * epsmac * hnorm ) then

!     --- Points coincide, skip this possibility

           angle = -1.5d0
           goto 1000
      end if
      e1  = e1 / eleng
      e2  = e2 / eleng
      e3  = x1 - x2
      e4  = y1 - y2
      eleng = sqrt ( e3*e3 + e4*e4 )
      h1  =  ( x1 + x2 ) / 2d0
      h2  =  ( y1 + y2 ) / 2d0
      hnorm = sqrt( h1*h1 + h2*h2 )
      if ( eleng<1.d2 * epsmac * hnorm ) then

!     --- Points coincide, skip this possibility

           angle = -1.5d0
           goto 1000
      end if
      e3  = e3 / eleng
      e4  = e4 / eleng
      angle = e1*e3 + e2*e4
1000  call erclos( 'msho11' )
      end
