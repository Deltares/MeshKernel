      subroutine msho19 ( coor, i1, i2, i3, surf )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 14-02-2005 Update
!        version  2.0  date 21-02-1994 New norms
!        version  1.0  date 17-04-1989
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
!     Compute area of triangle i1 - i2 - i3
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     area
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
      double precision coor(2,*), surf
      integer i1, i2, i3

!     coor      i   coordinate array
!     i1        i   node one of triangle
!     i2        i   node two of triangle
!     i3        i   node three of triangle
!     surf      o   area of triangle
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision x1, y1, x2, y2, x3, y3

!     x1       x-coordinate node one
!     x2       x-coordinate node two
!     x3       x-coordinate node three
!     y1       y-coordinate node one
!     y2       y-coordinate node two
!     y3       y-coordinate node three
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
      x1 = coor(1,i1)
      y1 = coor(2,i1)
      x2 = coor(1,i2)
      y2 = coor(2,i2)
      x3 = coor(1,i3)
      y3 = coor(2,i3)
      surf  = ( x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2) ) / 2d0
      end
