      subroutine mshtrans2dsur ( coor, npoint, xmint, ymint, tran,
     +                           ncoar, coar, ncurvs, curves, cocurvs,
     +                           userco, nuspnt )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.1  date 10-11-2009 Add array of userpoints
!        version  1.0  date 19-08-2009
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
!     Transform 2d region to region of unit length in first quadrant
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     transformation
!     map
!     coordinate
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      use msherror
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cconst'
!      include 'SPcommon/cinout'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer npoint, ncoar, ncurvs, nuspnt, curves(ncurvs)
      double precision coor(2,npoint), xmint, ymint, tran, coar(3,*),
     +                 cocurvs(2,*), userco(2,nuspnt)

!     coar          i/o   array containing coordinates and coarseness of
!                         special points to be used in later calculations:
!                         positions of these points are fixed!
!     cocurvs       i/o   array containing the crds of the fixed line nodes
!     coor          i/o   array containing the coordinates of the points
!     curves         i    the curve numbers of the extra internal curves
!     ncoar          i    number of special points in array coar
!     ncurvs         i    number: extra internal curves submitted by user
!     npoint         i    Number of nodal points in boundary
!     nuspnt         i    Number of userpoints
!     tran           o    scaling parameter
!     userco        i/o   array containing the coords of the user points
!     xmint          o    min x for transformation
!     ymint          o    min y for transformation
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, i1, i2, nnodes
      double precision xmax, ymax
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     i              Counting variable
!     i1             Help parameter to store a constant
!     i2             Help parameter to store a constant
!     nnodes         Number of nodes in curve
!     xmax           Maximum x-value
!     ymax           Maximum y-value
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     PRINRL1        Print 2d real vector
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
      call eropen ( 'mshtrans2dsur' )
      debug = .false. !.and. ioutp>=0 !AvD ioutp?
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshtrans2dsur'
  1      format ( a, 1x, (10i6) )
  2      format ( a, 1x, (5d12.4) )

      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000

!     --- Replace all nodes temporarily to the first (positive) quadrant:
!         (Determine translation parameters x- and y-tran)

      xmint = RINFIN
      xmax  = -RINFIN
      ymint = RINFIN
      ymax  = -RINFIN
      do i = 1, npoint
          xmint = min(xmint,coor(1,i))
          xmax  = max(xmax ,coor(1,i))
          ymint = min(ymint,coor(2,i))
          ymax  = max(ymax ,coor(2,i))
      end do

      if ( debug ) then
         write(irefwr,2) 'xmin, xmax', xmint, xmax
         write(irefwr,2) 'ymin, ymax', ymint, ymax
      end if

!     --- Make temporary transformation:

      if ( xmax - xmint>ymax - ymint ) then
         tran = xmax - xmint
      else
         tran = ymax - ymint
      end if

!     --- Scale all submitted points of coor array:

      coor(1,:) = 1d0 + ( coor(1,:) - xmint ) / tran
      coor(2,:) = 1d0 + ( coor(2,:) - ymint ) / tran

!AvD      if ( debug )
!     +   call prinrl1 ( coor, npoint, 2,
!     +                  'New coordinates, after scaling' )

!     --- Scale all user points (userco array):

      userco(1,:) = 1d0 + ( userco(1,:) - xmint ) / tran
      userco(2,:) = 1d0 + ( userco(2,:) - ymint ) / tran

!AvD      if ( debug )
!     +   call prinrl1 ( userco, nuspnt, 2,
!     +                  'New userpoint coordinates, after scaling' )

      if ( ncoar>0 ) then
!     --- ncoar>0, change internal points

         coar(1,1:ncoar) = 1d0 + ( coar(1,1:ncoar) - xmint ) / tran
         coar(2,1:ncoar) = 1d0 + ( coar(2,1:ncoar) - ymint ) / tran
         coar(3,1:ncoar) = coar(3,1:ncoar) / tran

         if ( debug )
     +      call prinrl1 ( coar, ncoar, 3,
     +                     'Internal points after scaling' )

      end if  ! ( ncoar>0 )

      if ( ncurvs>0 ) then

!     --- Scale nodes on submitted internal curves
!         Start computing number of nodes on the curves:

         nnodes = 0

         do i = 1, ncurvs
            if ( debug )
     +         write(irefwr,1) 'internal curve', i

!           --- Run through line i of the internal curves:

            i1 = nnodes + 1
            i2 = nnodes + curves(i)
            cocurvs(1,i1:i2) = 1d0+(cocurvs(1,i1:i2) - xmint) / tran
            cocurvs(2,i1:i2) = 1d0+(cocurvs(2,i1:i2) - ymint) / tran
!AvD            if ( debug )
!     +         call prinrl1 ( cocurvs(1,i1), curves(i), 2,
!     +                        'coordinates' )
            nnodes = i2
         end do  ! i = 1, ncurvs

      end if  ! ( ncurvs>0 )

1000  call erclos ( 'mshtrans2dsur' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End mshtrans2dsur'

      end if  ! ( debug )

      end
