      subroutine msho34( kpoint, coor  , ncurvs, curves,
     +                   kbndpt, nbndpt, numcurvboun,
     +                   boundary, nbound, nholes, userco, nuspnt,
     +                   coaval, ius1, ius2 )
! ======================================================================
!
!        programmer    niek praagman
!        version  1.1  date 07-12-2009 Determination istart adjusted
!        version  1.0  date 10-11-2009
!
!   copyright (c) 2009-2009  "Ingenieursbureau SEPRA"
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
!     Subroutine to determine begin- and end user point of curve on
!     which node kpoint is situated.
!     In node kpoint special coarseness is required.
!     Compute new coarsenesses for begin and end point of curve
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
!     mesh_generation
!     2d
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
!AvDtmp      include 'SPcommon/cmcdpr'
!AvDtmp      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kpoint, ncurvs, curves(ncurvs), numcurvboun,
     +        boundary(2,*),nbound, nholes, nbndpt, kbndpt(nbndpt),
     +        nuspnt, ius1, ius2
      double precision coor(2,*), userco(2,nuspnt),
     +        coaval(nuspnt+1)

!     boundary      i   array containing start and end nodes boundary curves
!     coaval        i   array with first the user given unity of coarseness
!                       and from position 2 to nuspnt+1 the coarsenesses
!                       of all prescribed user points
!     coor          i   coordinate array
!     curves        i   array with number of nodes for each internal curve
!     ius1          i   start node curve
!     ius2          i   end node curve
!     kbndpt        i   array containing numbers of all boundary points on
!                       all boundary curves, internal curves included
!     kpoint        i   nodepoint with too large coarseness
!     nbound        i   number of nodes in boundary curves - 1
!     nbndpt        i   length of boundary array kbndpt
!     ncurvs        i   number of internal curves
!     nholes        i   number of holes in boundary
!     numcurvboun   i   number of boundary curves
!     nuspnt        i   number of userpoints given by user
!     userco        i   coordinate array userpoints
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision :: afst, dis, xt, yt, xst, yst, xen, yen, xk, yk

      integer :: i, istart, j, k, n1, n2

      logical :: debug, found

!     afst         local distance between sequential nodes
!     debug        indicator for debugging
!     dis          Euclidian distance
!     found        indicator for finding position
!     i            loop variable
!     istart       helpvariable to find discontinuties
!     j            loop variable
!     k            loop variable
!     n1           help variable to determine ref number of node
!     n2           help variable to determine reference number of node
!     xen          x-coordinate
!     xk           x-coordinate
!     xst          x-coordinate
!     xt           x-coordinate
!     yen          y-coordinate
!     yk           y-coordinate
!     yst          y-coordinate
!     yt           y-coordinate
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
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Set debug and found

      debug = .false.
      found = .false.

!     debug = .false.

!     --- Coarseness problems in node kpoint. Coordinates are:
      xk = coor(1,kpoint)
      yk = coor(2,kpoint)
      if ( debug ) then
         write(irefwr,*) 'Debug info from msho34'
         write(irefwr,*) 'Node considered',kpoint,'has crds',xk,yk
      end if

!     --- Check contents of all submitted curves and nodes
!         and give warning for coarseness value
      istart = 0

      xst = 9999d0
      yst = 9999d0

      xen = 9999d0
      yen = 9999d0

!     --- Determine start and endpoint of the curve on which kpoint is
!         situated:

      if ( debug ) then
         write(irefwr,*) 'Contents array boundary'
         write(irefwr,*) 'Number of curves is numcurvboun',numcurvboun
         do i = 1, numcurvboun
            n1 = boundary(2,i)
            n2 = boundary(2,i+1)
            write(irefwr,*) 'curve',i,'nodepos',n1,'-',n2
         end do
         write(irefwr,*) 'Contents array',kbndpt
         do i=1,nbndpt
            write(irefwr,*) 'pos',i,'value',kbndpt(i)
         end do
      end if

      i = 0
      do while ( .not.found .and. i<numcurvboun )
         i = i + 1
!        --- Find (local) positions of curve i in boundary
         if ( i<numcurvboun ) then
            n1 = boundary(2,i)
            n2 = boundary(2,i+1)
         else if ( i==numcurvboun) then
            n1 = boundary(2,i)
            n2 = nbound + 1 + nholes
         endif

!        --- If istart = 0: set startvalue
         if ( istart==0 ) istart = kbndpt(n1)

         if ( kbndpt(n2 - 1) == istart .and. n2>n1+1 ) then
              n2 = n2 - 1
              istart = 0
         end if

         if ( debug ) then
            write(irefwr,*) 'Curve number',boundary(1,i)
            write(irefwr,*) 'Nodes sequentially',n1,'up to',n2
         end if

         !debug = .false.
         do j = n1, n2
            if ( kpoint == kbndpt(j) ) then
               if ( debug ) write(irefwr,*) 'On boundary-line',i
               xst = coor(1,kbndpt(n1))
               yst = coor(2,kbndpt(n1))
               xen = coor(1,kbndpt(n2))
               yen = coor(2,kbndpt(n2))
               found = .true.
               !write(irefwr,*) 'Found on line',n1,'- ',n2,'ready'
            end if
         end do

      end do

!     --- If not yet found: internal curves?
      i = 0

      do while ( .not. found .and. i<ncurvs )
         if ( debug ) write(irefwr,*) 'Not yet found'
         i = i+1
!     --- Check also internal lines
         if ( debug )
     +      write(irefwr,*) 'Internal',i,'number of nodes',curves(i)
!        --- Find (local) positions of this curve in boundary
         n1 = n2 + 1
         n2 = n1 + curves(i) - 1
         if ( debug ) write(irefwr,*) 'Intern Curve. Nodes',n1,'to',n2
         do j = n1, n2
            if ( kpoint == kbndpt(j) ) then
               if ( debug ) write(irefwr,*) 'Found on internal-line',i
               xst = coor(1,kbndpt(n1))
               yst = coor(2,kbndpt(n1))
               xen = coor(1,kbndpt(n2))
               yen = coor(2,kbndpt(n2))
               found = .true.
            end if
         end do
      end do

      if ( found .and. debug ) then
         write(irefwr,*) 'Node has been found'
         write(irefwr,*) 'start coordinates are',xst,yst
         write(irefwr,*) 'end coordinates are',xen,yen
      else if ( .not. found ) then
         write(irefwr,*) 'not found: error. Please inform SEPRA'
!AvD         call instop
      end if

!     --- Find userpoint for startpoint and endpoint:
      do i=1, nuspnt

         xt = userco(1,i)
         yt = userco(2,i)

         dis = sqrt( (xst-xt)*(xst-xt)+(yst-yt)*(yst-yt) )

         if ( dis < 1d-8 ) ius1 = i

         dis = sqrt( (xen-xt)*(xen-xt)+(yen-yt)*(yen-yt) )

         if ( dis < 1d-8 ) ius2 = i

      end do

      if ( debug ) then
!     --- If wanted give start and endnodenumber
         write(irefwr,*) 'Curve found has nodes'
         write(irefwr,*) 'Startpunt =',ius1,'en Endpunt =',ius2
         if ( ius1 == ius2 ) then
            write(irefwr,*) 'Error no line found in msho34!'
            write(irefwr,*) 'If your input is ok:'
            write(irefwr,*) 'Please inform SEPRA'
!AvD            call instop
         end if
      end if

      end ! msho34
