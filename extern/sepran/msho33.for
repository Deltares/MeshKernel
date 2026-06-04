      subroutine msho33 ( coor, ratio, i1, i2, i3 )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  1.0  date 17-01-2007
!
!   copyright (c) 2007-2007  "Ingenieursbureau SEPRA"
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
!     Determine the ratio (2*Rout/Rin) for all triangles of the mesh
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     radius
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
      integer i1, i2, i3
      double precision ratio, coor(2,*)

!     coor       i    array with coordinates of nodes
!     i1         i    first node of triangle
!     i2         i    second node of triangle
!     i3         i    third node of triangle
!     ratio      o    2 * radius divided by radius outer circle
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision opp, ri, ro, s, s1, s2, s3,
     +                 x1, x2, x3, y1, y2, y3

!     opp      surface
!     ri       radius inner circle
!     ro       radius outer circle
!     s        helpvariable to store sum of length of sides
!     s1       euclidian distance side 1
!     s2       euclidian distance side 2
!     s3       euclidian distance side 3
!     x1       x-coordinate node 1
!     x2       x-coordinate node 2
!     x3       x-coordinate node 3
!     y1       y-coordinate node 1
!     y2       y-coordinate node 2
!     y3       y-coordinate node 3
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
!     trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Determine coordinates

      x1 = coor(1,i1)
      y1 = coor(2,i1)

      x2 = coor(1,i2)
      y2 = coor(2,i2)

      x3 = coor(1,i3)
      y3 = coor(2,i3)

!     --- Determine length of the three sides:

      s1 = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
      s2 = sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) )
      s3 = sqrt( (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) )

      s = 0.5 * ( s1 + s2 + s3 )

!     --- Calculate the surfacearea of the triangle:

      opp = sqrt( s * ( s - s1 ) * ( s - s2 ) * ( s - s3 ) )

!     --- Calculate the values of the radii of outer and inner circle:

      ro = ( s1 * s2 * s3 ) / ( 4 * opp )

      ri = opp / s

!     --- Finally determine the ratio value:

      ratio = 2 * ri / ro

      end
