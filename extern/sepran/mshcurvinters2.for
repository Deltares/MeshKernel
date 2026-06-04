      subroutine mshcurvinters2 ( kbound, nbound, coor, isurnr )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 30-11-2005
!
!   copyright (c) 2005-2005  "Ingenieursbureau SEPRA"
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
!     Check if edges in boundary of surface do not intersect
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
!     surface
!     boundary
!     check
!     intersection
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use msherror
      use mshdummymethods
      
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kbound(2,*), nbound, isurnr
      double precision coor(2,*)

!     coor           i    contains the coordinates of the nodes on the boundary
!     isurnr         i    Surface sequence number
!     kbound         i    Contains the node numbers of the boundary edges
!                         kbound(1,i) first node of edge i
!                         kbound(2,i) last node of edge i
!     nbound         i    Number of edges on the boundary
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j
      double precision eps, fact1, fact2, x1(2), x2(2), x3(2), x4(2)
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     eps            Local accuracy
!     fact1          multiplication factor
!     fact2          multiplication factor
!     i              Counting variable
!     j              Counting variable
!     x1             coordinates of first point of edge in icurnr
!     x2             coordinates of last point of edge in icurnr
!     x3             coordinates of first point of edge in jcurnr
!     x4             coordinates of last point of edge in jcurnr
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERREAL         Put real in error message
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSHCROSSLINE1  Find the intersection of the lines x1-x2 and x3-x4
!     PRININ1        print 2d integer array
!     PRINRL1        Print 2d real vector
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2788   boundary intersects itself
! **********************************************************************
!
!                       PSEUDO CODE
!
!     for all edges on boundary do
!        for all edges with higher number on boundary do
!           if edges intersect give error message
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      character(len=260) localName
      localName = 'mshcurvinters2'
      call eropen( localName )
      
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshcurvinters2'
         write(irefwr,100) 'isurnr', isurnr
100      format ( a, 1x, (10i6) )
110      format ( a, 1x, (5d15.7) )

      end if
      if ( ierror/=0 ) go to 1000
      eps = 1d-3

!     --- for all edges on boundary do

      do i = 1, nbound-1

         x1(1) = coor(1,kbound(1,i)) ! first node of edge
         x1(2) = coor(2,kbound(1,i))
         x2(1) = coor(1,kbound(2,i)) ! last node of edge
         x2(2) = coor(2,kbound(2,i))

!        --- for all edges with higher number on boundary do

         do j = i+1, nbound

            x3(1) = coor(1,kbound(1,j)) ! first node of edge
            x3(2) = coor(2,kbound(1,j))
            x4(1) = coor(1,kbound(2,j)) ! last node of edge
            x4(2) = coor(2,kbound(2,j))

!           --- if edges intersect give error message

            call mshcrossline1 ( x1, x2, x3, x4, fact1, fact2, eps )
            if ( debug ) then
               write(irefwr,100) 'i, j', i, j
               write(irefwr,100) 'kbound(i)', kbound(1,i), kbound(2,i)
               write(irefwr,100) 'kbound(j)', kbound(1,j), kbound(2,j)
               write(irefwr,110) 'x1, x2', x1, x2
               write(irefwr,110) 'x3, x4', x3, x4
               write(irefwr,110) 'fact1, fact2', fact1, fact2
            end if  ! ( debug )

            if ( fact1>eps .and. fact1<1d0-eps .and.
     +           fact2>eps .and. fact2<1d0-eps ) then

!           --- curves intersect, give error message and leave subroutine

               call errint ( isurnr, 1 )
               call errint ( i, 2 )
               call errint ( j, 3 )
               call erreal ( x1(1) + fact1 * (x2(1)-x1(1)), 1 )
               call erreal ( x1(2) + fact1 * (x2(2)-x1(2)), 2 )
               call errsub ( 2788, 3, 2, 0 )

            end if  ! ( fact1>sqreps .and. fact1<1d0-sqreps ... )

         end do  ! j = i+1, nbound

      end do  !  i = 1, nbound-1

1000  call erclos ( 'mshcurvinters2' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End mshcurvinters2'

      end if

      end

