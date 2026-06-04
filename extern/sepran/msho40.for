      subroutine msho40 ( coor, kelem, npoint, nelem, iperm )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.0  date 15-09-2010 Adjust position continue's
!        version  1.1  date 14-12-2009 Check for double (Plaxis) points
!        version  1.0  date 18-03-2009
!
!   copyright (c) 2009-2010  "Ingenieursbureau SEPRA"
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
!     Check whether all barycenters of triangles are only belonging to
!     their "own" triangle and hence that the mesh is not intersecting.
!     Check also that all nodes are never inside another triangle.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
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
!      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer    npoint, nelem, kelem(3,nelem), iperm
      double precision coor(2,npoint)

!     coor        i   coordinate array
!     iperm       o   result , IPERM = 1 not all barycenters have a unique
!                                        triangle
!                              IPERM = 0 no double triangles
!     kelem       i   element array
!     nelem       i   number of elements in KELEM
!     npoint      i   number of points in COOR
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision x1, y1, x2, y2, x3, y3, xm, ym,
     +                 opp, det, dis, eps
      integer i, i1, i2, i3, ielem, j1, j2, j3, jelem

!     det        area of triangle
!     dis        variable for distance
!     eps        accuracy
!     i          loop variable nodes
!     i1         node
!     i2         node
!     i3         node
!     ielem      loop variable
!     j1         node
!     j2         node
!     j3         node
!     jelem      loop variable
!     opp        area of triangle
!     x1         x-coordinate
!     x2         x-coordinate
!     x3         x-coordinate
!     xm         x-coordinate
!     y1         y-coordinate
!     y2         y-coordinate
!     y3         y-coordinate
!     ym         y-coordinate
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
!     Check for all elements whether barycenter is not in one of the
!     other elements
!     If barycenter is inside give message
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      iperm = 0
      eps   = 1d-5

!     --- Run through all triangles:

      do ielem = 1, nelem

!     --- Determine barycenter:

         i1 = kelem(1,ielem)
         i2 = kelem(2,ielem)
         i3 = kelem(3,ielem)

!        --- Coordinates

         x1 = coor(1,i1)
         y1 = coor(2,i1)
         x2 = coor(1,i2)
         y2 = coor(2,i2)
         x3 = coor(1,i3)
         y3 = coor(2,i3)

!        --- Compute barycentre:

         xm = ( x1 + x2 + x3 ) / 3d0
         ym = ( y1 + y2 + y3 ) / 3d0

!        --- Run through all triangles of kelem:

         do jelem = 1, nelem
            if ( jelem/=ielem ) then
               j1 = kelem(1,jelem)
               j2 = kelem(2,jelem)
               j3 = kelem(3,jelem)

!              --- Coordinates element jelem:

               x1 = coor(1,j1)
               y1 = coor(2,j1)
               x2 = coor(1,j2)
               y2 = coor(2,j2)
               x3 = coor(1,j3)
               y3 = coor(2,j3)

!              --- Compute area of this triangle

               det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)

!              --- Check for barycentre:

               opp = (x1*(y2-ym) + x2*(ym-y1) + xm*(y1-y2))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 100
               opp = (x2*(y3-ym) + x3*(ym-y2) + xm*(y2-y3))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 100
               opp = (x3*(y1-ym) + x1*(ym-y3) + xm*(y3-y1))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 100

!              --- Point (xm,ym) is inside: error

               write(irefwr,*) 'Error: barycentre of element',ielem
               write(irefwr,*) 'is also inside triangle',jelem
               write(irefwr,*) 'Essential error!'
               iperm = 1
            end if
100         continue
         end do
      end do

      if ( iperm==1 ) then
         write(irefwr,*) 'Execution stopped. Please inform SEPRA'
!AvD         call instop
      end if

!     --- Next consider all nodes: they should never be inside a triangle:

      do i = 1 , npoint

!     --- Run through all triangles:

         do ielem = 1, nelem

!        --- Determine nodenumbers:

            i1 = kelem(1,ielem)
            i2 = kelem(2,ielem)
            i3 = kelem(3,ielem)
            if ( i/=i1.and.i/=i2.and.i/=i3 ) then

!           --- Check position of node i

               xm = coor(1,i)
               ym = coor(2,i)
               x1 = coor(1,i1)
               y1 = coor(2,i1)
               x2 = coor(1,i2)
               y2 = coor(2,i2)
               x3 = coor(1,i3)
               y3 = coor(2,i3)

!              --- Check for double points:
!                  Extra check whether coordinates of nodes are exactly the
!                  same,
!                  i.e. check whether it concerns Plaxis points

               dis = (x1-xm)*(x1-xm) + (y1-ym)*(y1-ym)
               det = abs(x1)+abs(xm)+abs(y1)+abs(ym)
               if ( dis<eps * det ) goto 200
               dis = (x2-xm)*(x2-xm) + (y2-ym)*(y2-ym)
               det = abs(x2)+abs(xm)+abs(y2)+abs(ym)
               if ( dis<eps * det ) goto 200
               dis = (x3-xm)*(x3-xm) + (y3-ym)*(y3-ym)
               det = abs(x3)+abs(xm)+abs(y3)+abs(ym)
               if ( dis<eps * det ) goto 200

!              --- Points are different: Compute area of this triangle

               det = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)

!              --- Check for node i:

               opp = (x1*(y2-ym) + x2*(ym-y1) + xm*(y1-y2))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 200
               opp = (x2*(y3-ym) + x3*(ym-y2) + xm*(y2-y3))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 200
               opp = (x3*(y1-ym) + x1*(ym-y3) + xm*(y3-y1))/det
               if ( opp<-eps .or. opp>1d0+eps ) goto 200

!              --- Point (xm,ym) is inside this triangle: error

               write(irefwr,*) 'Error: node',i,'inside element',ielem
               write(irefwr,*) 'Essential error!'
               iperm = 1
            end if
200         continue
         end do
      end do

!     --- Problems encountered?

      if ( iperm==1 ) then
         write(irefwr,*) 'Execution stopped'
         write(irefwr,*) 'Please inform SEPRA'
!AvD         call instop
      end if

      end
