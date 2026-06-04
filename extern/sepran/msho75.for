      subroutine msho75 ( xi, yi, xj, yj, x1, y1, x2, y2, ih )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  2.0  date 06-11-2007 Extra checks added, "new" accuracy
!        version  1.0  date 27-06-1996
!
!   copyright (c) 1996-2008  "Ingenieursbureau SEPRA"
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
!     Subroutine to check whether line i-j and line 1-2 have a common
!     point.
! **********************************************************************
!
!                       KEYWORDS
!
!     3d
!     mesh_generation
!     tetrahedron
! **********************************************************************
!
!                       MODULES USED
      use mshconstants
      use mshdummymethods
!
      implicit none

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer ih
      double precision x1, x2, xi, xj, y1, y2, yi, yj

!     ih      o    indicator of result
!                     ih = 0 : there is a common point
!                     ih = 1 : there is no common point
!     x1      i    x-coordinate of point 1
!     x2      i    x-coordinate of point 2
!     xi      i    x-coordinate of point i
!     xj      i    x-coordinate of point j
!     y1      i    x-coordinate of point 1
!     y2      i    x-coordinate of point 2
!     yi      i    x-coordinate of point i
!     yj      i    x-coordinate of point j
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision eps, r1, r2, r3, r4, xs, ys,
     +                 xmin, xmax, ymin, ymax,
     +                 xmi , xma , ymi , yma

!     eps     small number
!     r1      help variable for determination (xs,ys)
!     r2      help variable for determination (xs,ys)
!     r3      help variable for determination (xs,ys)
!     r4      help variables for determining (xs,ys)
!     xma     help variable for extremes line 1-2
!     xmax    help variable for extremes line i-j
!     xmi     help variable for extremes line 1-2
!     xmin    help variable for extremes line i-j
!     xs      x-coordinate of intersection of line segments
!     yma     help variable for extremes line 1-2
!     ymax    help variable for extremes line i-j
!     ymi     help variable for extremes line 1-2
!     ymin    help variable for extremes line i-j
!     ys      y-coordinate of intersection of line segments
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS   deconcatenate name from string of calling routines
!     EROPEN   concatenate name to string of calling routines
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
!     If the line through (x1,y1) and (x2,y2) intersects the line
!     through (xi,yi) and (xj,yj) at a point (xs,ys) and that point
!     lies between (x1,y1) and (x2,y2) as well as between (xi,yi) and
!     (xj,yj) than the two line segments have a point in common.
!     PSEUDO CODE
!     Compute coordinates "boxes"
!     Check whether boxes are disjunct and if so
!         No intersection is possible
!     Elseif xi <> xj then
!        If x1 <> x2 then
!           Determine x-coordinate of common point xs
!           If xs lies between xi and xj then
!              If xs lies between x1 and x2 then
!                 There is a common point
!              Endif
!           Endif
!        Else
!           If yi = yj then
!              If xi and xj lie on opposite sides of xs or
!               yi and yj lie on opposite sides of xs then
!                 There is a common point
!              Endif
!           Else
!              Determine y-coordinate of common point ys
!              If ys lies between yi and yj then
!                 If ys lies between y1 and y2 then
!                    There is a common point
!                 Endif
!              Endif
!           Endif
!        Endif
!     Else
!        If x1 <> x2 then
!           If y1 = y2 then
!              If yi and yj lie on opposite sides of y1 or
!               x1 and x2 lie on opposite sides of xi then
!                 There is a common point
!              Endif
!           Else
!              Determine y-coordinate of common point ys
!              If ys lies between yi and yj then
!                 If ys lies between y1 and y2 then
!                    There is a common point
!                 Endif
!              Endif
!           Endif
!        Else
!           If x1 = xi then: a common point
!        Endif
!     Endif
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      character(len=260) localName
      localName = 'msho75'
      call eropen( localName )

      eps    = 10 * epsmac

!       --- Start with "no intersection"

      ih     = 1

!     --- Check boxes to find out whether an intersection is possible:
!         Determine extreme values line i-j:

      xmin  = min( xi, xj )
      xmax  = max( xi, xj )
      ymin  = min( yi, yj )
      ymax  = max( yi, yj )

!     --- Determine extreme values line 1-2:

      xmi   = min( x1, x2 )
      xma   = max( x1, x2 )
      ymi   = min( y1, y2 )
      yma   = max( y1, y2 )

      if ( xmi>xmax .or. xma<xmin .or.
     +     ymi>ymax .or. yma<ymin ) then

!     --- Intersection not possible (no further search!)

         ih = 1

      else if (abs(xi-xj)>(eps * (abs(xi) + abs(xj)))) then

!     --- xi#xj

         if (abs(x1-x2)>(eps * (abs(x1) + abs(x2)))) then

!        --- xi#xj and x1#x2
!            determine value x-coordinate intersection

            r1 = ( y1 * x2 - y2 * x1 ) / ( x2 - x1 )
            r2 = ( yi * xj - yj * xi ) / ( xj - xi )
            r3 = ( yj - yi ) / ( xj - xi )
            r4 = ( y2 - y1 ) / ( x2 - x1 )
            if ( abs(r3-r4)>eps ) then
               xs = ( r1 - r2 ) / ( r3 - r4 )
               if ( ( xi<xs .and. xj>xs ) .or.
     +              ( xj<xs .and. xi>xs ) .or.
     +              ( abs(xi-xs)<eps )         .or.
     +              ( abs(xj-xs)<eps ) ) then
                  if ( ( x1<xs .and. x2>xs ) .or.
     +                 ( x2<xs .and. x1>xs ) .or.
     +                 ( abs(x1-xs)<eps )         .or.
     +                 ( abs(x2-xs)<eps ) ) then
                     ih = 0
                  end if
               end if
            end if
         else

!        --- x1 = x2  line x = constant

            if ( abs(yi-yj)<eps ) then

!           --- x1 = x2  and  yi = yj

               xs = x1
               ys = yi
               if ( .not. ( ( xi<xs .and. xj<xs ) .or.
     +                      ( xi>xs .and. xj>xs ) .or.
     +                      ( y1<ys .and. y2<ys ) .or.
     +                      ( y1>ys .and. y2>ys ) ) ) then
                  ih = 0
               end if
            else

!           --- Determine y-coordinate intersection

               ys = ( (yj - yi) * x1 + yi * xj - yj * xi ) / (xj - xi)
               if ( ( y1<ys .and. y2>ys ) .or.
     +              ( y2<ys .and. y1>ys ) .or.
     +              ( abs(y1-ys)<eps )         .or.
     +              ( abs(y2-ys)<eps ) ) then
                  if ( ( yi<ys .and. yj>ys ) .or.
     +                 ( yj<ys .and. yi>ys ) .or.
     +                 ( abs(yi-ys)<eps )         .or.
     +                 ( abs(yj-ys)<eps ) ) then
                     ih = 0
                  end if
               end if
            end if
         end if

      else

!     --- xi = xj  line x = constant

         if ( abs ( x1 - x2 )>( eps * ( abs(x1) + abs(x2)))) then

!        --- xi = xj and x1#x2

            if ( abs(y1-y2)<eps ) then
               xs = xi
               ys = y1
               if ( .not. (( yi<ys .and. yj<ys ) .or.
     +                     ( yi>ys .and. yj>ys ) .or.
     +                     ( x1<xs .and. x2<xs ) .or.
     +                     ( x1>xs .and. x2>xs ) ) ) then
                  ih = 0
               end if

            else

!           --- determine y-coordinate intersection

               ys = ( (y2 - y1) * xi + y1 * x2 - y2 * x1 ) / (x2 - x1)
               if ( ( yi<ys .and. yj>ys ) .or.
     +              ( yi>ys .and. yj<ys ) .or.
     +              ( abs(yi-ys)<eps )         .or.
     +              ( abs(yj-ys)<eps ) ) then
                  if ( ( y1<ys .and. y2>ys ) .or.
     +                 ( y2<ys .and. y1>ys ) .or.
     +                 ( abs(y1-ys)<eps )         .or.
     +                 ( abs(y2-ys)<eps ) ) then
                     ih = 0
                  end if
               end if
            end if

         else

!        --- x1 = x2  line also x = constant

            if ( abs(x1 - xi)<( eps * ( abs(x1) + abs(xi)))) ih = 0

         end if

      end if

      call erclos( 'msho75' )

      end
