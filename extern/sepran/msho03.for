      subroutine msho03( i, j, coor, afst)

! =======================================================================
!
!
!
!         programmer   niek praagman
!
!         version 2.0  date  03-02-1994 New norms
!         version 1.0  date  12-04-1989
!
!
!
!   copyright (c) 1989-1994  "ingenieursbureau sepra"
!   permission to copy or distribute this software or documentation
!   in hard or soft copy granted only by written license obtained
!   from "ingenieursbureau sepra".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system (e.g. , in memory, disk or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
!
! ***********************************************************************
!
!     DESCRIPTION
!
!     Subroutine to compute the Euclidian distance of two nodes
!
! ***********************************************************************
!
!     KEYWORDS
!
!     mesh_generation
!     2d
!     distance
!
! ***********************************************************************
!
!
!     INPUT / OUTPUT  PARAMETERS
!
      implicit none
      double precision afst, coor(2,*)
      integer i, j

!     afst        o   Euclidian distance
!     coor        i   coordinate array  ( Standard SEPRAN )
!     i           i   first node to be considered
!     j           i   second node to be considered
!
! ***********************************************************************
!
!     COMMON BLOCKS
!
! ***********************************************************************
!
!     LOCAL PARAMETERS
!
      double precision dx, dy

!     dx     difference in x-coordinates of i and j
!     dy     difference in y-coordinates of i and j
!
! ***********************************************************************
!
!     SUBROUTINES CALLED
!
! ***********************************************************************
!
!     ERROR MESSAGES
!
! ***********************************************************************
!
!     METHOD
!
!     trivial
!
! =======================================================================

!     Compute distance

      dx = coor( 1,i ) - coor( 1,j )
      dy = coor( 2,i ) - coor( 2,j )

      afst = sqrt( dx*dx + dy*dy )

      end
