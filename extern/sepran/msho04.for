      subroutine msho04 ( coor, npoint, xmin, xmax, ymin, ymax )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 11-02-2005 New for min and max determination
!        version  2.0  date 03-02-1994 New norms
!        version  1.0  date 12-04-1989
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
!     Subroutine to determine extreme values of region considered
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     extreme
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
      double precision coor(2,*), xmin, xmax, ymin, ymax
      integer npoint

!     coor        i   coordinate array  ( Standard SEPRAN )
!     npoint      i   number of points in coor
!     xmax        o   max x-coordinate
!     xmin        o   min x-coordinate
!     ymax        o   max y-coordinate
!     ymin        o   min y-coordinate
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i
      double precision x, y

!     i         loop variable
!     x         x-coordinate
!     y         y-coordinate
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
!     Run through all points and compare with temporary extreme
!     values; if necessary replace values by better new extremes
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Set starting values

      xmax = -rinfin
      ymax = -rinfin

      xmin = -xmax
      ymin = -ymax

!     --- Run through all npoint points and compute

      do i = 1, npoint

          x = coor( 1,i )
          y = coor( 2,i )

          xmin = min( xmin, x)
          xmax = max( xmax, x)
          ymin = min( ymin, y)
          ymax = max( ymax, y)

      end do

      end
