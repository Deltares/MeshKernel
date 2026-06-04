      subroutine msho42( kbndpt, lenbnd, i, j, ja )
! ======================================================================
!
!        programmer    niek praagman
!        version  1.0  date 02-02-2010
!
!   copyright (c) 2010-2010  "Ingenieursbureau SEPRA"
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
!     Check whether line i - j is a piece in one of the internal curves
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
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
      integer          :: i, j, ja, lenbnd
      integer          :: kbndpt(lenbnd)

!     i            i   first node of linepiece considered
!     j            i   second node of linepiece considered
!     ja           o   indicator that i-j is internal (=1) or not (=0)
!     kbndpt       i   linepieces of internal boundaries
!     lenbnd       i   length of array kbndpt
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer       :: ia, ib, k

!     ia          node number
!     ib          node number
!     k           loop variable
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
!     Trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!

      ja = 0

!     --- Run through all pieces of kbndpt and check:

      do k = 1, lenbnd, 2

         ia = kbndpt(  k     )
         ib = kbndpt(  k + 1 )
!        --- Check for both directions:
         if ( ia==i .and. ib==j .or.
     +        ia==j .and. ib==i ) ja = 1

      end do

      end
