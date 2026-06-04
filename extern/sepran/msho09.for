      subroutine msho09( coor, i, j, e1, e2, xm, ym )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.0  date 14-02-1994 New norms
!        version  1.0  date 13-04-1989
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
!     Determine the vector ( e1, e2 ) which is perpendicular to the line
!     of the points i-j. Also point (xm,ym) in the middle of line i-j is
!     computed
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     perpendicular
! **********************************************************************
!
!                       MODULES USED
!
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer i, j
      double precision coor(2,*), e1, e2, xm, ym

!     coor        i    coordinate array
!     e1          o    first coordinate vector perpendicular to line
!     e2          o    second coordinate vector perpendicular to line
!     i           i    first node of line
!     j           i    second node of line
!     xm          o    first coordinate midpoint of line
!     ym          o    second coordinate midpoint of line
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision xi, yi, xj, yj, eleng

!     eleng        normalisation parameter
!     xi           x-coordinate point i
!     xj           x-coordinate point j
!     yi           y-coordinate point i
!     yj           y-coordinate point j
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
!     Determine coordinates two line endpoints
!     Compute normalizing length and the unit vector
!     Compute mid point
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Determine coordinates

      xi = coor( 1, i )
      yi = coor( 2, i )
      xj = coor( 1, j )
      yj = coor( 2, j )
      e1 = - ( yj - yi )
      e2 =   ( xj - xi )

!     --- Normalize

      eleng = sqrt( e1*e1 + e2*e2 )
      e1 = e1 / eleng
      e2 = e2 / eleng

!     --- Mid point

      xm  =  0.5000001234 * xi + 0.4999998766 * xj
      ym  =  0.5000001234 * yi + 0.4999998766 * yj
      end
