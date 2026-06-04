      subroutine msho07 ( xcub, ycub, xmini, coor, kbound, nbound, ja )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.2  date 09-01-2008 Replace mshg03 by msho75
!        version  2.1  date 26-01-2005 Layout
!        version  2.0  date 03-04-1994 new norms
!        version  1.1  date 10-10-1990 msho23 replaced by mshg03
!        version  1.0  date 13-04-1989
!
!   copyright (c) 1989-2008  "Ingenieursbureau SEPRA"
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
!     Subroutine to determine whether point (xcub,ycub) belongs to
!     the region given by kbound.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
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
      double precision xcub, ycub, xmini, coor(2,*)
      integer kbound(*), nbound, ja

!     coor       i  coordinate array
!     ja         o  output parameter
!                     if point belongs to area ( ja = 1 )
!                     else ( ja = 0 )
!     kbound     i  boundary element array
!     nbound     i  number of boundary elements
!     xcub       i  x-coordinate of point considered
!     xmini      i  smallest x-coordinate of all nodes
!     ycub       i  y-coordinate of point considered
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision xleft, yleft, xmin, ymin, ymax, x1, y1, x2, y2
      integer i, isnij, ih, i1, i2

!     i         loop variable
!     i1        first node boundary segment
!     i2        second node boundary segment
!     ih        help indicator for crossing point
!     isnij     number of crossings with boundary segments
!     x1        x-coordinate of node i1
!     x2        x-coordinate of node i2
!     xleft     reference x-coordinate (far away)
!     xmin      min value of x-values of boundary line
!     y1        y-coordinate of node i1
!     y2        y-coordinate of node i2
!     yleft     reference y-coordinate
!     ymax      max value of y-values of boundary line
!     ymin      min value of y-values of boundary line
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     MSHO75  Check whether two line segments have a point in common
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
!     Suppose that a line is drawn from point (xcub,ycub) to reference
!     point (xleft,yleft). Count all crossings of this line with boundary
!     lines (i1,i2) from array KBOUND. If the number of crossings is odd
!     than point (xcub,ycub) belongs to area given by KBOUND.
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      xleft = xmini - 10d0
      yleft = ycub

      isnij = 0
      ja    = 0

      do i = 1, nbound

         i1 = kbound(2*i-1)
         i2 = kbound(2*i  )

         x1 = coor( 1,i1 )
         y1 = coor( 2,i1 )

         x2 = coor( 1,i2 )
         y2 = coor( 2,i2 )

         xmin = min(x1,x2)

         ymin = min(y1,y2)
         ymax = max(y1,y2)

         if ( xcub<xmin ) goto 100
         if ( ycub<ymin ) goto 100
         if ( ycub>ymax ) goto 100

!        --- Potential point

         call msho75( x1, y1, x2, y2, xleft, yleft, xcub, ycub, ih )

         if ( ih==1 ) goto 100

!        --- Real intersection point

         isnij = isnij + 1

100   end do

      if ( (isnij - (isnij/2) * 2)/=0 ) ja = 1

      end
