      subroutine msh406( kelem, j, i1, i2, i3, i4, i5, i6 )

! =======================================================================
!
!
!         programmer   niek praagman
!
!         version 2.0  date  03-03-1993 adjustment layout
!         version 1.1  date  22-08-1991 Zdenek: no dimension, cdc, ce
!         version 1    date  19-5-1986
!
!
!
!   copyright (c) 1986-1993  "ingenieursbureau sepra"
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
!     Subroutine to place nodes in element array in case of a quadratic
!     triangle
!
! ***********************************************************************
!
!     KEYWORDS
!
!     mesh_generation
!     quadratic
!     triangle
!
! ***********************************************************************
!
!     INPUT / OUTPUT PARAMETERS
!
      implicit none

      integer i1, i2, i3, i4, i5, i6, j, kelem( 6,j )

!     i1      first node
!     i2      second node
!     i3      third node
!     i4      fourth node
!     i5      fifth node
!     i6      sixth node
!     j       triangle number to be considered
!     kelem   element array
!
! ***********************************************************************
!
!     COMMON BLOCKS
!
! ***********************************************************************
!
!     LOCAL PARAMETERS
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
!     PSEUDO CODE
!
! =======================================================================

      kelem( 1,j ) = i1
      kelem( 2,j ) = i2
      kelem( 3,j ) = i3
      kelem( 4,j ) = i4
      kelem( 5,j ) = i5
      kelem( 6,j ) = i6

      end

