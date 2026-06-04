      subroutine msho13 ( coor, i1, i2, kdrie, kstapl, kstap, xn, yn )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.2  date 14-02-2010 Checks adjusted (crossing and double)
!        version  2.1  date 27-08-2009 Adjust accuracy for double points
!        version  2.0  date 14-02-2005 Update
!        version  1.4  date 21-02-2003 Layout
!        version  1.3  date 21-02-1997 Adjust for "double" points
!        version  1.2  date 16-02-1994 New norms
!        version  1.1  date 10-10-1990 msho23 replaced by mshg03
!        version  1.0  date 14-04-1989
!
!   copyright (c) 1989-2010  "Ingenieursbureau SEPRA"
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
!     Determine whether lines from i1 and i2 to new point (xn,yn) do or
!     do not cross boundary lines
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     intersection
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
!      include 'SPcommon/cmcdpr'
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer i1, i2, kdrie, kstapl(2,*), kstap
      double precision coor(2,*), xn, yn

!     coor      i   coordinate array
!     i1        i   node of basis element
!     i2        i   node of basis element
!     kdrie     o   number of line that lines to new point
!                   crosses.  0 means no crossing !
!     kstap     i   number of elements in kstapl
!     kstapl    i   current boundary elements array
!     xn        i   x-coordinate new point
!     yn        i   y-coordinate new point
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision xi1, yi1, xi2, yi2, xmin, xmax, ymin, ymax,
     +                 x1, y1, x2, y2, xmi, xma, ymi, yma, xh, yh,
     +                 xpn1, xpn2, ypn1, ypn2,
     +                 dis1, dis2, eps, h1, h2
      integer i, ii1, ii2, ipn1, ipn2, ih

!     dis1      helpdistance
!     dis2      helpdistance
!     eps       accuracy parameter
!     h1        helpvariable for distance
!     h2        helpvariable for distance
!     i         loop variable
!     ih        help indicator
!     ii1       first node of boundary line
!     ii2       second node of boundary line
!     ipn1      node number
!     ipn2      node number
!     x1        x-coordinate point ii1
!     x2        x-coordinate point ii2
!     xh        x-coordinate helppoint
!     xi1       x-coordinate first basis point
!     xi2       x-coordinate second basis point
!     xma       extreme x-value quadrangular envelope boundary line
!     xmax      extreme x-value quadrangular envelope basis line
!     xmi       extreme x-value quadrangular envelope boundary line
!     xmin      extreme x-value quadrangular envelope basis line
!     xpn1      x-coordinate point ipn1
!     xpn2      x-coordinate point ipn2
!     y1        y-coordinate point ii1
!     y2        y-coordinate point ii2
!     yh        y-coordinate helppoint
!     yi1       y-coordinate first basis point
!     yi2       y-coordinate second basis point
!     yma       extreme y-value quadrangular envelope boundary line
!     ymax      extreme y-value quadrangular envelope basis line
!     ymi       extreme y-value quadrangular envelope boundary line
!     ymin      extreme y-value quadrangular envelope basis line
!     ypn1      y-coordinate point ipn1
!     ypn2      y-coordinate point ipn2
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS  Resets old name of previous subroutine of higher level
!     EROPEN  Produces concatenated name of local subroutine
!     MSHO75  Check whether two linesegments have an internal point
!             in common
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
!     Compute coordinates line nodes
!     Determine extreme values of line i1-i2 and new point xn,yn
!     Run through all boundary elements ii1-ii2 and check whether
!     "new" lines i1-n and i2-n intersect line ii1-ii2. Check also
!     whether line midpoint-n intersects line ii1-ii2.
!     Return with number of line that is crossing. Check for double
!     crossings. If it seems to be a double line leave the routine with
!     kdrie = -1. In the other case take the "nearest" line kdrie.
!     If result is kdrie = 0 then no intersection has been found
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Check whether new point is allowed
      call eropen( 'msho13' )

      eps = 1d4 * epsmac

      xi1 = coor( 1, i1 )
      yi1 = coor( 2, i1 )
      xi2 = coor( 1, i2 )
      yi2 = coor( 2, i2 )

!     --- Determine extreme values line i1-i2 and new point

      xmin  = min( xi1, xi2, xn )
      xmax  = max( xi1, xi2, xn )
      ymin  = min( yi1, yi2, yn )
      ymax  = max( yi1, yi2, yn )

      kdrie = 0

!     --- Run through all boundary elements and check whether lines intersect
!         each other

      do i = 2, kstap

         ii1 = kstapl( 1, i )
         ii2 = kstapl( 2, i )

         x1  = coor( 1, ii1 )
         y1  = coor( 2, ii1 )

         x2  = coor( 1, ii2 )
         y2  = coor( 2, ii2 )

         xmi = min(x1,x2)
         xma = max(x1,x2)
         ymi = min(y1,y2)
         yma = max(y1,y2)

         if ( xmi<=xmax .and.
     +        xma>=xmin .and.
     +        ymi<=ymax .and.
     +        yma>=ymin ) then
!           --- Potential line i:

!           --- Check for point i1
            if ( i1/=ii1 .and. i1/=ii2 ) then
!           --- Check for distance and coincidence of points

               dis1 = abs(xi1-x1)+abs(yi1-y1)
               dis2 = abs(xi1-x2)+abs(yi1-y2)

               h1   = abs(xi1)+abs(x1)+abs(yi1)+abs(y1)
               h2   = abs(xi1)+abs(x2)+abs(yi1)+abs(y2)

               if ( dis1>eps*h1 .and. dis2>eps*h2 ) then

                  call msho75( xi1, yi1, xn, yn, x1, y1, x2, y2, ih )

                  if ( ih==0 ) then
                     if ( kdrie == 0 ) then
!                    --- First crossing line:
                        kdrie =  i
                     else if ( kdrie>0 .and. kdrie/=i ) then
!                    --- Double crossing:
!                       -Compute distance to midpoint kdrie
                        ipn1 = kstapl( 1, kdrie )
                        ipn2 = kstapl( 2, kdrie )

!                       -Coordinates midpoint piece kdrie:
                        xpn1  = (coor(1,ipn1)+coor(1,ipn2))/2d0
                        ypn1  = (coor(2,ipn1)+coor(2,ipn2))/2d0

!                       -Coordinates midpoint piece i1 - i2:
                        xpn2  = (xi1+xi2)/2d0
                        ypn2  = (yi1+yi2)/2d0

                        dis1 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                         (ypn1-ypn2)*(ypn1-ypn2)

!                       -Coordinates midpoint piece i:
                        xpn1 = ( x1 + x2 ) / 2d0
                        ypn1 = ( y1 + y2 ) / 2d0

                        dis2 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                         (ypn1-ypn2)*(ypn1-ypn2)

                        if ( dis2<(1d0-eps)*dis1 ) then
!                          --- Take new one i.e. i:
                           kdrie = i
                        else if ( dis2<(1d0+eps)*dis1 ) then
!                          --- Problems: return
                           kdrie = -1
                           goto 1000
                        else
!                          --- Use old value do not change kdrie
                        end if

                     end if
                  end if
               end if
            end if
!           --- Check for point i2
            if ( i2/=ii1 .and. i2/=ii2 ) then
!           --- Check for distance and coincidence of points
               dis1 = abs(xi2-x1)+abs(yi2-y1)
               dis2 = abs(xi2-x2)+abs(yi2-y2)
!              --- Helpaccuracy:
               h1   = abs(xi2)+abs(x1)+abs(yi2)+abs(y1)
               h2   = abs(xi2)+abs(x2)+abs(yi2)+abs(y2)
!              --- If no coïncidence:
               if ( dis1>eps*h1 .and. dis2>eps*h2 ) then
                  call msho75( xi2, yi2, xn, yn, x1, y1, x2, y2, ih )

                  if ( ih == 0 ) then
                     if ( kdrie == 0 ) then
!                    --- First crossing line:
                        kdrie =  i
                     else if ( kdrie>0 .and. kdrie/=i ) then
!                    --- Double crossing:
!                       -Compute distance to midpoint kdrie
                        ipn1 = kstapl( 1, kdrie )
                        ipn2 = kstapl( 2, kdrie )

!                       -Coordinates midpoint piece kdrie:
                        xpn1  = (coor(1,ipn1)+coor(1,ipn2))/2d0
                        ypn1  = (coor(2,ipn1)+coor(2,ipn2))/2d0

!                       -Coordinates midpoint piece i1 - i2:
                        xpn2  = (xi1+xi2)/2d0
                        ypn2  = (yi1+yi2)/2d0

                        dis1 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                         (ypn1-ypn2)*(ypn1-ypn2)

!                       -Coordinates midpoint piece i:
                        xpn1 = ( x1 + x2 ) / 2d0
                        ypn1 = ( y1 + y2 ) / 2d0

                        dis2 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                         (ypn1-ypn2)*(ypn1-ypn2)

                        if ( dis2<(1d0-eps)*dis1 ) then
!                          --- Take new one i.e. i:
                           kdrie = i
                        else if ( dis2<(1d0+eps)*dis1 ) then
!                          --- Problems: return
                           kdrie = -1
                           goto 1000
                        else
!                          --- Use old value do not change kdrie
                        end if
                     end if
                  end if
               end if
            end if
!           --- Finally check for middle point
            xh = ( xi1 + xi2 ) / 2d0
            yh = ( yi1 + yi2 ) / 2d0

            xh = ( eps * xn + (1d0-eps) * xh )
            yh = ( eps * yn + (1d0-eps) * yh )

            call msho75(xh,yh,xn,yn,x1,y1,x2,y2,ih)

            if ( ih == 0 ) then
               if ( kdrie == 0 ) then
!              --- First crossing line:
                  kdrie =  i
               else if ( kdrie>0 .and. kdrie/=i ) then
!              --- Double crossing:
!                 -Compute distance to midpoint kdrie
                  ipn1 = kstapl( 1, kdrie )
                  ipn2 = kstapl( 2, kdrie )

!                 -Coordinates midpoint piece kdrie:
                  xpn1  = (coor(1,ipn1)+coor(1,ipn2))/2d0
                  ypn1  = (coor(2,ipn1)+coor(2,ipn2))/2d0

!                 -Coordinates midpoint piece i1 - i2:
                  xpn2  = (xi1+xi2)/2d0
                  ypn2  = (yi1+yi2)/2d0

                  dis1 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                   (ypn1-ypn2)*(ypn1-ypn2)

!                 -Coordinates midpoint piece i:
                  xpn1 = ( x1 + x2 ) / 2d0
                  ypn1 = ( y1 + y2 ) / 2d0

                  dis2 = (xpn1-xpn2)*(xpn1-xpn2)+
     +                   (ypn1-ypn2)*(ypn1-ypn2)

                  if ( dis2<(1d0-eps)*dis1 ) then
!                    --- Take new one i.e. i:
                     kdrie = i
                  else if ( dis2<(1d0+eps)*dis1 ) then
!                    --- Problems: return
                     kdrie = -1
                     goto 1000
                  else
!                    --- Use old value do not change kdrie
                  end if

               end if
            end if
         end if
      end do

!     --- Checking is ready
1000  call erclos( 'msho13' )

      end
