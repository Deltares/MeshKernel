      subroutine mshcurvinters1 ( icurnr, jcurnr, icurvs, curves )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 28-11-2005
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
!     Give an error if 2 curves intersect
!     2d only
! **********************************************************************
!
!                       KEYWORDS
!
!     intersection
!     curve
!     2d
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use msherror
      
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer icurnr, jcurnr, icurvs(*)
      double precision curves(2,*)

!     curves         i    array containing the coordinates of the curves
!     icurnr         i    Curve sequence number of first curve
!     icurvs         i    array containing the number of points in the curves
!                         accumulated
!     jcurnr         i    Curve sequence number of second curve
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer k, j, istart, jstart, inodes, jnodes
      double precision x1(2), x2(2), x3(2), x4(2), eps, fact1, fact2
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     eps            Local accuracy
!     fact1          multiplication factor
!     fact2          multiplication factor
!     inodes         number of nodes on icurnr
!     istart         Last position used at icurnr
!     j              Counting variable
!     jnodes         number of nodes on jcurnr
!     jstart         Last position used at jcurnr
!     k              Counting variable
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
!     ERRWAR         Warnings
!     MSHCROSSLINE1  Find the intersection of the lines x1-x2 and x3-x4
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2787   curves intersect
! **********************************************************************
!
!                       PSEUDO CODE
!
!     for all edges on curve 1 do
!        for all edges on curve 2 do
!           if edges intersect give error message
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'mshcurvinters1' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshcurvinters1'
         write(irefwr,100) 'icurnr, jcurnr', icurnr, jcurnr
100      format ( a, 1x, (10i6) )
110      format ( a, 1x, (5d12.4) )

      end if
      if ( ierror/=0 ) go to 1000
      eps = sqreps

!     --- for all edges on curve 1 do

      if ( icurnr==1 ) then
         istart = 0
      else
         istart = icurvs(icurnr-1)
      end if  ! ( icurnr==1 )
      inodes = icurvs(icurnr)-istart
      x1(1) = curves(1,istart+1)  ! first node of edge
      x1(2) = curves(2,istart+1)

      do j = 2, inodes

         x2(1) = curves(1,istart+j) ! last node of edge
         x2(2) = curves(2,istart+j)

!        --- for all edges on curve 2 do

         if ( jcurnr==1 ) then
            jstart = 0
         else
            jstart = icurvs(jcurnr-1)
         end if  ! ( icurnr==1 )
         jnodes = icurvs(jcurnr)-jstart
         x3(1) = curves(1,jstart+1)  ! first node of edge
         x3(2) = curves(2,jstart+1)

         do k = 2, jnodes

            x4(1) = curves(1,jstart+k) ! last node of edge
            x4(2) = curves(2,jstart+k)

!           --- if edges intersect give error message

            call mshcrossline1 ( x1, x2, x3, x4, fact1, fact2, eps )
            if ( debug ) then
               write(irefwr,100) 'j, k', j, k
               write(irefwr,110) 'x1, x2', x1, x2
               write(irefwr,110) 'x3, x4', x3, x4
               write(irefwr,110) 'fact1, fact2', fact1, fact2
            end if  ! ( debug )

            if ( fact1>sqreps .and. fact1<1d0-sqreps .and.
     +           fact2>sqreps .and. fact2<1d0-sqreps ) then

!           --- curves intersect, give error message and leave subroutine

               call errint ( icurnr, 1 )
               call errint ( jcurnr, 2 )
               call erreal ( x1(1) + fact1 * (x2(1)-x1(1)), 1 )
               call erreal ( x1(2) + fact1 * (x2(2)-x1(2)), 2 )
               call errwar ( 2787, 2, 2, 0 )
               go to 1000

            end if  ! ( fact1>sqreps .and. fact1<1d0-sqreps ... )
            x3(1) = x4(1)
            x3(2) = x4(2)

         end do  ! k = 2, jnodes
         x1(1) = x2(1)
         x1(2) = x2(2)

      end do  ! j = 1, nnodes-1

1000  call erclos ( 'mshcurvinters1' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End mshcurvinters1'

      end if

      end

