      subroutine msho14 ( coor, jpn, npoint, itri, i1, i2,
     +                    xn, yn, dista )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 14-02-2005 Update
!        version  2.0  date 16-02-1994 New norms
!        version  1.0  date 14-04-1989
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
!     Determine nearest boundary point to proposed new point.
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     neighbour
!     distance
! **********************************************************************
!
!                       MODULES USED
!
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer i1, i2, itri(*), jpn, npoint
      double precision coor(2,*), xn, yn, dista

!     coor        i     coordinate array
!     dista       o     distance new point to boundary point
!     i1          i     first node basis line
!     i2          i     second node basis line
!     itri        i     array with for each point number of times
!                       that point appears in present boundary
!     jpn         o     node number of nearest point
!     npoint      i     actual number of points in mesh
!     xn          i     x-coordinate new point
!     yn          i     y-coordinate new point
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer ih
      double precision xi, yi, dx, dy, dist

!     dist      Euclidian distance
!     dx        first component distance vector
!     dy        second component distance vector
!     ih        loop variable
!     xi        x-coordinate boundary point
!     yi        y-coordinate boundary point
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
!     Run through all points
!     If point in boundary and not one of the basis points then
!        compute distance new point - boundary point
!        if distance smaller then reference value adjust
!           reference value and node number
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      jpn = 0
      do ih = 1 , npoint
         if ( ih/=i1 .and. ih/=i2 .and. itri(ih)/=0 ) then
            xi = coor( 1, ih )
            yi = coor( 2, ih )
            dx = xn - xi
            dy = yn - yi
            dist = sqrt ( dx*dx + dy*dy )
            if ( dist<dista ) then
               jpn   = ih
               dista = dist
            end if
         end if
      end do
      end
