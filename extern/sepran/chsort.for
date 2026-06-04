      subroutine chsort ( keysrt, kgrad, npoint )
! ======================================================================
!
!       programmer     Jos van Kan/Onno Hoitinga
!       version 3.1    date 26-12-1998 Remove continues
!       version 3.0    date 25-01-1995 Complete new subroutine, programmed
!                                      by Onno Hoitinga
!       version 2.3    date 20-06-1994 Adaptation comments
!       version 2.2    date 03-02-1994 New comments
!       version 2.1    date 22-11-1990 Change declarations to prevent errors
!       version 2.0    date 12-05-1989 Complete revision, with extra parameters
!       version 1.0    date 21-06-1987
!
!
!   copyright (c) 1987-1998  "Ingenieursbureau SEPRA"
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
!   This routine gets as input array keysrt of length npoint
!   on output array kgrad contains the indices to sort array keysrt
!   The Heap sort algorithm is used
! **********************************************************************
!
!                       KEYWORDS
!
!     sort
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer npoint, keysrt(npoint), kgrad(npoint)

!     keysrt   i    sortkey (integer) on array to be sorted.
!     kgrad    o    contains permutation index that sorts keysrt
!                   in ascending order on output.
!     npoint   i    number of entries in keysrt & kgrad
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j, l, ir, kgradt, q

!     i       Counting variable
!     ir      ?
!     j       Counting variable
!     kgradt  ?
!     l       General loop variable
!     q       ?
! **********************************************************************
!
!                       I/O
!
!     none
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!  This routine originates from Numerical Recipes (fortran edition)
!     Cambridge University Press 1989
! **********************************************************************
!
!                       MODULES USED
!
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      do j = 1, npoint
         kgrad(j) = j
      end do
      if ( npoint.le.1 ) go to 1000
      l = npoint/2+1
      ir = npoint

200     if ( l.gt.1 ) then
          l = l-1
          kgradt = kgrad(l)
          q = keysrt(kgradt)
        else
          kgradt = kgrad(ir)
          q = keysrt(kgradt)
          kgrad(ir) = kgrad(1)
          ir = ir-1
          if ( ir.eq.1 ) then
            kgrad(1) = kgradt
            go to 1000
          end if
        end if
        i = l
        j = l+l
300     if ( j.le.ir ) then
          if ( j.lt.ir ) then
            if ( keysrt(kgrad(j)).lt.keysrt(kgrad(j+1)) ) j = j+1
          end if
          if ( q.lt.keysrt(kgrad(j)) ) then
            kgrad(i) = kgrad(j)
            i = j
            j = j+j
          else
            j = ir+1
          end if
        go to 300
        end if
        kgrad(i) = kgradt
      go to 200

1000  end
