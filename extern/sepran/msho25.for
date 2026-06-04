      subroutine msho25 ( kstapl, kstap, coor )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.0  date 14-02-2005 Update: Smaller (0.98) testvalue
!        version  2.2  date 06-03-2003 Make start less arbitrary
!        version  2.1  date 14-02-2003 New norms
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
!     Place smallest element at top of kstapl
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     element_size
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpr'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*)
      integer kstap, kstapl(2,kstap)

!     coor      i   coordinate array
!     kstap     i   number of triangles in kstapl
!     kstapl   i/o  actual array with boundary triangles
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision  afref, xa, ya, xb, yb, dx, dy, afst
      integer ielem, il, ia, ib, i1, i2

!     afref         distance
!     afst          distance
!     dx            x-distance
!     dy            y-distance
!     i1            node number
!     i2            node number
!     ia            node number
!     ib            node number
!     ielem         element number
!     il            loop variable
!     xa            x-coordinate
!     xb            x-coordinate
!     ya            y-coordinate
!     yb            y-coordinate
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
!     run through all line segments and compare the length of each
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      ielem = 0

      afref = rinfin

!     --- Find smallest element in kstapl

      do il = 1 , kstap

          ia = kstapl(1,il)
          ib = kstapl(2,il)

          xa = coor(1,ia)
          ya = coor(2,ia)

          xb = coor(1,ib)
          yb = coor(2,ib)

          dx = xb - xa
          dy = yb - ya

          afst = dx*dx + dy*dy

          if ( afst<0.98d0 * afref ) then

!         --- Only if the element is really smaller than the other elements
!             it is marked. This is done to avoid too much coincidence

             afref = afst
             ielem = il

          end if

      end do

!     --- Interchange two elements

      if ( ielem>1 ) then

           i1 = kstapl(1,ielem)
           i2 = kstapl(2,ielem)

           kstapl(1,ielem) = kstapl(1,1)
           kstapl(2,ielem) = kstapl(2,1)

           kstapl(1,1) = i1
           kstapl(2,1) = i2

      end if

      end
