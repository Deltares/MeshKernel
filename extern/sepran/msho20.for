      subroutine msho20 ( kstapl, kstap, itri )
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
!     Reorder boundary elements (first to last)
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     reordering
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
      integer kstap, kstapl(2,kstap), itri(*)

!     itri      i,o   array with for each point actual number of
!                     occurrences in boundary triangles
!     kstap      i    number of triangles in kstapl
!     kstapl    i,o   array with boundary triangles
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i1, i2, ikl

!     i1      node number
!     i2      node number
!     ikl     loop variable
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
!     --- Change array kstapl : place first element on last positions

      i1 = kstapl(1,1)
      i2 = kstapl(2,1)
      do ikl = 2, kstap
           kstapl (1,ikl-1) = kstapl (1,ikl)
           kstapl (2,ikl-1) = kstapl (2,ikl)
      end do
      kstapl (1,kstap) = i1
      kstapl (2,kstap) = i2

!     --- Refill itri accordingly

      itri (i1) = itri (i1) + 1
      itri (i2) = itri (i2) + 1
      end
