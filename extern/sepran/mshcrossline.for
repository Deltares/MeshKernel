      subroutine mshcrossline ( x1, x2, x3, x4, fact1, fact2, eps )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.2  date 16-07-2008 adaptation for zero determinant
!        version  1.1  date 29-06-2006 Debug statements
!        version  1.0  date 30-11-2000
!
!   copyright (c) 2000-2008  "Ingenieursbureau SEPRA"
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
!    Find the intersection of the lines x1-x2 and x3-x4
!    The result is the parameter fact, where the intersection point is defined
!    as x1 + fact*(x2-x1)
!    If no intersection point is found fact gets the value -infinity
!    Two-dimensions only
! **********************************************************************
!
!                       KEYWORDS
!
!     2d
!     intersection
!     line
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
      double precision x1(2), x2(2), x3(2), x4(2), fact1, fact2, eps

!     eps            i    local accuracy
!     fact1          o    Defines the factor for the intersection point
!                         The intersection point on the face is defined by
!                         x1 + fact1 * (x2-x1)
!     fact2          o    Defines the factor for the intersection point
!                         The intersection point on the face is defined by
!                         x3 + fact2 * (x4-x3)
!     x1             i    Contains the coordinates of the point x1
!     x2             i    Contains the coordinates of the point x2
!     x3             i    Contains the coordinates of the point x3
!     x4             i    Contains the coordinates of the point x4
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision det, amax
      logical debug

!     amax           maximum value of coordinates
!     debug          If true debug statements are carried out otherwise
!                    they are not
!     det            determinant of system of equations to be solved
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRINRL         Print 1d real vector
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
!     The line through x1-x2 is defined as
!     x = x1 + fact1 * (x2-x1)
!     The line through x3-x4 is defined as
!     x = x3 + fact2 * (x4-x3)
!     If an intersection is present then
!        x1 + fact * (x2-x1) = x3 + alpha1 * (x4-x3)
!        Hence:
!        [ (x2-x1) (x4-x3) ] [fact -alpha1]^T = x3-x1
!        This system can be solved if det( [ (x2-x1) (x4-x3) ] ) # 0
!        The solution is:
!
!        [fact1 -fact2]^T =  [ (x2-x1) (x4-x3) ]^(-1) (x3-x1)
!
!        or:
!
!        det = (x2-x1)*(y4-y3)-(x4-x3)*(y2-y1)
!        fact1 = ((y4-y3)*(x3-x1)+(x3-x4)*(y3-y1))/det
!        fact2 = -((y1-y2)*(x3-x1)+(x2-x1)*(y3-y1))/det
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      character(len=260) localName
      localName = 'mshcrossline'
      call eropen( localName )
      
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshcrossline'
  2      format ( a, 1x, (5d12.4) )
!         call prinrl ( x1, 2, 'x1' )
!         call prinrl ( x2, 2, 'x2' )
!         call prinrl ( x3, 2, 'x3' )
!         call prinrl ( x4, 2, 'x4' )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

      det = (x2(1)-x1(1))*(x4(2)-x3(2))-(x4(1)-x3(1))*(x2(2)-x1(2))
      amax = max(abs(x1(1)), abs(x1(2)), abs(x2(1)), abs(x2(2)),
     +           abs(x3(1)), abs(x3(2)), abs(x4(1)), abs(x4(2)) )
      if ( debug ) write(irefwr,2) 'det', det
      if ( abs(det)<=epsmac*amax*100d0 ) then
         fact1 = -rinfin
         fact2 = -rinfin
      else
         fact1 = ((x4(2)-x3(2))*(x3(1)-x1(1))+
     +           (x3(1)-x4(1))*(x3(2)-x1(2)))/det
         if ( fact1>=-eps .and. fact1<=1d0+eps ) then
            fact2 = -((x1(2)-x2(2))*(x3(1)-x1(1))+
     +                (x2(1)-x1(1))*(x3(2)-x1(2)))/det
         else
            fact2 = -rinfin
         end if
      end if

1000  call erclos ( 'mshcrossline' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End mshcrossline'
         write(irefwr,2) 'fact1, fact2', fact1, fact2

      end if  ! ( debug )

      end
