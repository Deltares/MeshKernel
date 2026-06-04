      subroutine msho21 ( cube, ncube, refvol )
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
!     Calculate reference surf values for all cubes
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     reference_value
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
      double precision  cube(*), refvol(*)
      integer ncube

!     cube     i    array containing the coarsenesses in the
!                   blocks
!     ncube    i    number of blocks
!     refvol   o    reference volume value in block
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i

!     i       loop variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     None
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     None
! **********************************************************************
!
!                       PSEUDO CODE
!
!     In array cube the coarseness is found for each cube
!     Next a referencevalue for the surface (of a triangle) is computed
!     and stored
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      do i = 1, ncube
         refvol(i) = 0.5 * cube(i) * cube(i)
      end do

      end
