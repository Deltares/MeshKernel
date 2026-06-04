      subroutine msho32 ( kelem, nelem, inpelm, npunt )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.0  date 15-02-2005 Update
!        version  2.0  date 21-02-1994 New norms
!        version  1.0  date 12-04-1991
!
!   copyright (c) 1991-2005  "Ingenieursbureau SEPRA"
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
!     Recreate a triangular mesh from quadratic triangles.
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     recreation
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
      integer kelem(*), nelem, inpelm, npunt

!     inpelm     i    number of nodes in one element
!     kelem     i,o   array containing the triangular quadratic elements
!                     at input and the linear ones at output
!     nelem      i    number of elements
!     npunt      o    highest point number in triangular mesh
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j1, j2

!     i        loop variable
!     j1       local pointer
!     j2       local pointer
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
!     --- Run through array kelem to make linear triangles

      do i = 1 , nelem
         j1 = 3 * ( i - 1 )
         j2 = inpelm * ( i - 1 )
         kelem( j1 + 1 ) = kelem( j2 + 1 )
         kelem( j1 + 2 ) = kelem( j2 + 3 )
         kelem( j1 + 3 ) = kelem( j2 + 5 )
      end do

!     --- Find highest node number

      npunt = 0
      do i = 1 , 3 * nelem
         if ( kelem(i)>npunt ) npunt = kelem(i)
      end do
      end
