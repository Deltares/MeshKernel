      subroutine msho15 ( coor, npoint, kelem, nelem, kstapl, kstap,
     +                    itri, i1, i2, xnx, yny )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 14-02-2005 Update
!        version  2.0  date 16-02-1994 New norms
!        version  1.0  date 14-04-1989
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
!     Fill new element and point in respectively arrays
!     coor and kelem and adjust arrays kstapl and itri.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     element
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
      integer npoint, kelem(3,*), nelem, kstapl(*), kstap, itri(*),
     +        i1, i2
      double precision coor(2,*), xnx, yny

!     coor        i,o   coordinate array
!     i1           i    first node of basissegment
!     i2           i    second node basissegment
!     itri        i,o   array to indicate whether point is in boundary
!                       or not
!     kelem       i,o   element array
!     kstap       i,o   number of elements in kstapl
!     kstapl      i,o   array with current boundary segments
!     nelem       i,o   number of elements in kelem
!     npoint      i,o   number of nodes in mesh
!     xnx          i    x-coordinate new point
!     yny          i    y-coordinate of new point
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer ib

!     ib      :  loopvariable
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
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Fill array coor

      npoint = npoint + 1
      coor( 1, npoint ) = xnx
      coor( 2, npoint ) = yny

!     --- Fill array kelem

      nelem = nelem + 1
      kelem(1,nelem) = i1
      kelem(2,nelem) = i2
      kelem(3,nelem) = npoint

!     --- Adjust kstapl and itri

      do ib = 1, 2 * kstap - 2
          kstapl(ib) = kstapl(ib+2)
      end do
      ib    = 2 * kstap - 2
      kstap = kstap + 1
      kstapl(ib+1) = i1
      kstapl(ib+2) = npoint
      kstapl(ib+3) = npoint
      kstapl(ib+4) = i2
      itri(npoint) = 2
      itri(i1)     = itri(i1) + 1
      itri(i2)     = itri(i2) + 1
      end
