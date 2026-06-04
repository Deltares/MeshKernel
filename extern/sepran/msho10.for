      subroutine msho10 ( kelem, nelem, i1, i2, jpn, kstapl,
     +                    kstap, itri )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.1  date 21-02-2003 Layout
!        version  2.0  date 14-02-1994 New norms
!        version  1.0  date 13-04-1989
!
!   copyright (c) 1989-2003  "Ingenieursbureau SEPRA"
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
!     Fill new element in array KELEM and adjust arrays KSTAPL and ITRI
!
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
      use mshconstants
      use mshdummymethods
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kelem(3,*), nelem, i1, i2, jpn, kstapl(2,*), kstap,
     +        itri(*)

!     i1          i    first node of new element
!     i2          i    second node of new element
!     itri       i,o   array indicating whether a node n is on the boundary
!                      itri(n) = 2 or not : itri(n) = 0
!     jpn         i    third node of element just created
!     kelem      i,o   element array
!     kstap      i,o   number of elements in kstapl
!     kstapl     i,o   array containing the actual boundary
!     nelem      i,o   actual number of elements
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer idrie, ieen, itwe, ks, ia, ib
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     ia             node of boundary element
!     ib             node of boundary element
!     idrie          countvariable for new number of elements in kstapl
!     ieen           check for line i2  - jpn
!     itwe           check for line jpn - i1
!     ks             loop variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
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
!     Fill new created element in in array kelem
!     Adjust kstapl and itri
!     Run through all elements of kstapl
!       if element is one of the sides of the new element : cancel
!     Check which elements of the three sides have already been cancelled
!     If sides are left place them in kstapl
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'msho10' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msho10'
         write(irefwr,*) 'nelem, i1, i2, jpn ', nelem, i1, i2, jpn

      end if

!     --- Fill new created element in in array kelem

      nelem             = nelem + 1
      kelem( 1, nelem ) = i1
      kelem( 2, nelem ) = i2
      kelem( 3, nelem ) = jpn

!     --- Adjust kstapl and itri

      idrie  = 0
      ieen   = 1
      itwe   = 1
      do ks = 2, kstap
         ia = kstapl(1,ks)
         ib = kstapl(2,ks)
         if ( ia.eq.i2 .and. ib.eq.jpn ) then
            ieen = 0
         else if ( ia.eq.jpn .and. ib.eq.i1 ) then
            itwe = 0
         else
            idrie = idrie + 1
            kstapl(1,idrie) = ia
            kstapl(2,idrie) = ib
         end if
      end do
      if ( debug )
     +    write(irefwr,*) 'ieen, itwe, idrie ',ieen, itwe, idrie

!     --- Check which elements have already been used

      if ( ieen.eq.1 ) then
         idrie = idrie + 1
         kstapl(1,idrie) = jpn
         kstapl(2,idrie) = i2
         itri(i2 ) = itri(i2 ) + 1
         itri(jpn) = itri(jpn) + 1
      else
         itri(i2 ) = itri(i2 ) - 1
         itri(jpn) = itri(jpn) - 1
      end if

      if ( itwe.eq.1 ) then
         idrie = idrie + 1
         kstapl(1,idrie) = i1
         kstapl(2,idrie) = jpn
         itri(jpn) = itri(jpn) + 1
         itri(i1 ) = itri(i1 ) + 1
      else
         itri(jpn) = itri(jpn) - 1
         itri(i1 ) = itri(i1 ) - 1
      end if
      kstap = idrie

      call erclos ( 'msho10' )

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End msho10'

      end if

      end
