      subroutine mshchkstapl ( kstapl, lenkstapl, text )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 10-02-2003
!
!   copyright (c) 2003-2003  "Ingenieursbureau SEPRA"
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
!     Check the contents of array kstapl
!     It is checked if a node that is on the left-hand side is also the
!     right-hand side
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
!     check
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!AvD      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kstapl(2,*), lenkstapl
      character *(*) text

!     kstapl         i    Array containing the actual boundary build by
!                         pairs of nodes
!     lenkstapl      i    Number of pairs in kstapl
!     text           i    String variable containing some text to be
!                         printed as part of the debugging
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j, iseq1(lenkstapl), iseq2(lenkstapl),
     +        iwork(lenkstapl), ihelp(lenkstapl)
      logical debug, found

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     found          If true an error has been found
!     i              Counting variable
!     ihelp          integer work array to be used for the sorting
!     iseq1          Array containing the nodes of the left-hand side
!                    of the pairs after sorting
!     iseq2          Array containing the node after sortings of the right-hand
!                    side of the pairs
!     iwork          integer work array to store the nodes of either the
!                    left-hand side or the right-hand side of the pairs
!     j              Counting variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     CHSORT         Sort integer array for increasing sequence
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     INSTOP         Stop the program
!     PRININ1        print 2d integer array
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
      call eropen ( 'mshchkstapl' )
      debug = .false.

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshchkstapl'

      end if
      if ( lenkstapl.eq.0 ) go to 1000

!     --- Check array by sorting each row of pairs

      do i = 1, lenkstapl
         iwork(i) = kstapl(1,i)
      end do
      call chsort ( iwork, ihelp, lenkstapl )
      do i = 1, lenkstapl
         iseq1(i) = iwork(ihelp(i))
      end do

      do i = 1, lenkstapl
         iwork(i) = kstapl(2,i)
      end do
      call chsort ( iwork, ihelp, lenkstapl )
      do i = 1, lenkstapl
         iseq2(i) = iwork(ihelp(i))
      end do

      found = .false.
      do i = 1, lenkstapl
         if ( iseq1(i).ne.iseq2(i) ) then

!        --- first point found in second column

            found = .true.
            j = i
            go to 100

         end if

      end do

100   if ( found ) then

!     --- Both rows do not contain the same node numbers

         write(irefwr,*) text
         write(irefwr,*) 'Both columns of kstapl do not contain the ',
     +                   'same node numbers'
         write(irefwr,*) 'The ', j, '-th nodes are different'
         do i = 1, lenkstapl
            write(irefwr,110) i, iseq1(i), iseq2(i)
110         format ( i5, ':  ', 2i6 )
         end do
         call prinin1 ( kstapl, lenkstapl, 2, 'kstapl' )
         !AvD: call instop

      end if

1000  call erclos ( 'mshchkstapl' )
      if ( debug ) write(irefwr,*) 'End mshchkstapl'

      end
