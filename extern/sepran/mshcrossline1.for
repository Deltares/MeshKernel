      subroutine mshcrossline1 ( x1, x2, x3, x4, fact1, fact2, eps )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 26-11-2003
!
!   copyright (c) 2000-2003  "Ingenieursbureau SEPRA"
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
!     Find the intersection of the lines x1-x2 and x3-x4
!     The result is the parameter fact, where the intersection point is defined
!     as x1 + fact*(x2-x1)
!     If no intersection point is found fact gets the value -infinity
!     Two-dimenions only
!     First it is checked if the rectangles around both lines intersect
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
      double precision xmin1, xmin2, ymin1, ymin2, xmax1, xmax2,
     +                 ymax1, ymax2

!     xmax1          maximum x-coordinate of first edge
!     xmax2          maximum x-coordinate of second edge
!     xmin1          maximum y-coordinate of first edge
!     xmin2          maximum y-coordinate of second edge
!     ymax1          maximum x-coordinate of third edge
!     ymax2          maximum x-coordinate of fourth edge
!     ymin1          maximum y-coordinate of third edge
!     ymin2          maximum y-coordinate of fourth edge
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     MSHCROSSLINE   Computes the intersection of two lines
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
      fact1 = -rinfin
      fact2 = -rinfin
      xmin1 = min(x1(1),x2(1))
      xmax1 = max(x1(1),x2(1))
      ymin1 = min(x1(2),x2(2))
      ymax1 = max(x1(2),x2(2))
      xmin2 = min(x3(1),x4(1))
      xmax2 = max(x3(1),x4(1))
      ymin2 = min(x3(2),x4(2))
      ymax2 = max(x3(2),x4(2))
      if ( xmax1.lt.xmin2-eps .or. xmin1.gt.xmax2+eps .or.
     +     ymax1.lt.ymin2-eps .or. ymin1.gt.ymax2+eps ) go to 1000
      call mshcrossline ( x1, x2, x3, x4, fact1, fact2, eps )
1000  end

