      subroutine msho36 ( coor, npoint, istep, ncurvs, curves, cocurvs,
     +                    kstapl, kstap, kbndpt, nbndpt, coarsemin,
     +                    extquanodes, coarse )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  1.1  date 15-07-2009 Quadratic internal lines allowed
!        version  1.0  date 19-06-2008
!
!   copyright (c) 2008-2009  "Ingenieursbureau SEPRA"
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
!     Add nodal points of internal lines to coor and adjust array of
!     boundaries kstapl for mesh generation
! **********************************************************************
!
!                       KEYWORDS
!
!     2d
!     mesh
!     mesh_generation
!     plaxis
!     internal_lines
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
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cconst'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(2,*), cocurvs(2,*), coarse(*), coarsemin
      integer npoint, istep, ncurvs, curves(ncurvs),
     +        kstap, kstapl(2,*), nbndpt, kbndpt(*),
     +        extquanodes(*)

!     coarse        i/o   helparray for coarsenesses in usergiven points
!     coor          i/o   coordinate array
!     cocurvs        i    coordinates nodes on internal lines
!     curves         i    array with number of nodes on internal curves
!                         sequentially
!     extquanodes    i    help array for nodes quadratic elements internal
!                         lines
!     istep          i    stepsize do loop (1=linear, 2=quadratic)
!     kbndpt        i/o   array with local positions of boundary points,
!                         points, nodes on internal lines and
!                         prescribed coarseness points
!     kstap         i/o   length of boundary elements array
!     kstapl        i/o   array with nodes of boundary elements
!     nbndpt        i/o   number of elements in array kbndpt
!     ncurvs         i    number of internal curves
!     npoint        i/o   number of points in contour and internal lines
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, j, iext, npn, nnodes, ibound
      double precision dis, dx, dy, eps, x, y
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     dis            distance
!     dx             x-difference
!     dy             y-difference
!     eps            accuracy
!     i              loop variable
!     ibound         length of array extquanodes
!     iext           help loop variable
!     j              loop variable
!     nnodes         helpvariable to store temporarily number of nodes
!     npn            helpvariable to store nodenumber
!     x              x-coordinate
!     y              y-coordinate
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRININ1        print 2d integer array
!     PRINRL1        Print 2d real vector
! **********************************************************************
!
!                       I/O
!
!    None
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    None
! **********************************************************************
!
!                       PSEUDO CODE
!
!    Run through all internal curves
!    For each curve place nodes in coor
!    Check whether this node was allready in one of the preceding
!    curves (if so consider kbndpt and coor)
!    Augment array kstapl and kbndpt
!
! **********************************************************************
!
!                       DATA STATEMENTS
!
!    None
! ======================================================================
!
      call eropen ( 'msho36' )
      debug = .false.

      eps = 1D-6 * coarsemin

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msho36'
         write(irefwr,100) 'kstap', kstap
100      format ( a, 1x, (10i6) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Compute number of nodes on curves:

      nnodes = 0
      ibound = 1

      do i=1, ncurvs

!     --- Run through line i of the internal curves:

         do j= nnodes + 1, nnodes + curves(i) - 1, istep

!        --- Determine coordinates start node

            npoint = npoint + 1

!           --- Set node number in array kbndpt

            nbndpt = nbndpt + 1
            kbndpt(nbndpt) = npoint

!           --- Fill coor

            coor(1,npoint) = cocurvs(1,j)
            coor(2,npoint) = cocurvs(2,j)

            x = cocurvs(1,j)
            y = cocurvs(2,j)

            dx = cocurvs(1,j+1) - cocurvs(1,j)
            dy = cocurvs(2,j+1) - cocurvs(2,j)

            coarse(npoint) = sqrt( dx*dx + dy*dy )

            npn = 0

!           --- Check nodes so far:

            do iext = 1 , npoint-1

               dx = x - coor(1,iext)
               dy = y - coor(2,iext)

               dis = sqrt(dx*dx + dy*dy)

               if ( dis<eps .and. npn==0 ) then

!              --- Old point

                  npn    = iext
                  npoint = npoint - 1
                  kbndpt(nbndpt) = npn

               end if
            end do

!           --- If quadratic add extra node:

            if ( istep==2 ) then

               npoint = npoint + 1

!              --- Set node number in array kbndpt

               nbndpt = nbndpt + 1
               kbndpt(nbndpt) = npoint

!              --- Fill coor

               coor(1,npoint) = cocurvs(1,j+1)
               coor(2,npoint) = cocurvs(2,j+1)

            end if

!           --- Augment and fill kstapl in both directions

            kstap = kstap + 1

            if ( npn>0 ) then

               kstapl(1, kstap) = npn
               kstapl(2, kstap) = npoint + 1

               kstap = kstap + 1
               kstapl(1, kstap) = npoint + 1
               kstapl(2, kstap) = npn

               if ( istep==2 ) then
!              --- Extra for neighbours array

                  extquanodes(ibound+1) = npn
                  extquanodes(ibound+2) = npoint
                  extquanodes(ibound+3) = npoint+1

                  ibound = ibound + 3
               end if

            else

               kstapl(1, kstap) = npoint + 1 - istep
               kstapl(2, kstap) = npoint + 1

               kstap = kstap + 1
               kstapl(1, kstap) = npoint + 1
               kstapl(2, kstap) = npoint + 1 - istep

               if ( istep==2 ) then
!              --- Extra for neighbours array

                  extquanodes(ibound+1) = npoint - 1
                  extquanodes(ibound+2) = npoint
                  extquanodes(ibound+3) = npoint + 1

                  ibound = ibound + 3
               end if

            end if

         end do

!        --- Add last node of line to coor array and fill
!            kbndpt accordingly:

         j = nnodes + curves(i)

         x = cocurvs(1,j)
         y = cocurvs(2,j)

         npn = 0

!        --- Check nodes so far:

         do iext = 1 , npoint-1

            dx = x - coor(1,iext)
            dy = y - coor(2,iext)

            dis = sqrt(dx*dx + dy*dy)

            if ( dis<eps .and. npn==0 ) then

!           --- Old point, already existing:

               npn    = iext

            end if
         end do

         if ( npn==0 ) then

            npoint = npoint + 1

            nbndpt = nbndpt + 1
            kbndpt(nbndpt) = npoint

            coor(1, npoint) = cocurvs(1,j)
            coor(2, npoint) = cocurvs(2,j)

            dx = cocurvs(1,j) - cocurvs(1,j-1)
            dy = cocurvs(2,j) - cocurvs(2,j-1)

            dis = sqrt( dx*dx + dy*dy )

            coarse( npoint ) = dis

         else

!        --- Old point:

            nbndpt = nbndpt + 1
            kbndpt(nbndpt) = npn

!           --- Adjust kstapl

            kstapl(1,kstap  ) = npn
            kstapl(2,kstap-1) = npn

         end if

         nnodes = nnodes + curves(i)

      end do

!     --- Set number of positions used in array extquanodes:

      if ( istep==2 ) extquanodes(1) = ibound

1000  call erclos ( 'msho36' )

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'kstapl at the end of msho36'
         do i=1, kstap
            write(irefwr,*) i,kstapl(1,i),kstapl(2,i)
         end do
         write(irefwr,*) 'kbndpt at end msho36'
         write(irefwr,*) (kbndpt(i),i=1,nbndpt)

         write(irefwr,*) 'coor at the end of msho36'
         do i=1,npoint
            write(irefwr,*) i,coor(1,i),coor(2,i)
         end do

         write(irefwr,*) 'End debug info msho36'

      end if

      end
