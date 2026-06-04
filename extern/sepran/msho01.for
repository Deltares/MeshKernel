      subroutine msho01( kbound, nbound, istart, ibuur, coarse,
     +                   coor, npoint, coarsemin, coarsemax, coar,
     +                   ncoar )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  4.1  date 29-08-2008 No coarseness 0 allowed
!        version  4.0  date 26-01-2005 Extra parameters coar and ncoar
!        version  3.0  date 08-03-2003 Extra parameters coarsemin, coarsemax
!        version  2.1  date 08-02-2001 Layout
!        version  2.0  date 03-02-1994 New norms
!        version  1.1  date 07-11-1990 New error message
!        version  0.1  date 12-04-1989
!
!   copyright (c) 1989-2008  "Ingenieursbureau SEPRA"
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
!     Subroutine to check whether the boundary-elements of the area
!     considered are correct and to fill in array coarse for each point
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     coarseness
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      
      implicit none
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coarse(*), coor(*), coarsemin, coarsemax,
     +                 coar(3,*)
      integer kbound(*), nbound, istart(*), ibuur(*), npoint, ncoar

!     coar          i/o   array containing coordinates and coarseness of
!                         special points to be used in later calculations:
!                         positions of these points are fixed!
!                         In ths subroutine the coarseness may be changed
!     coarse         o    array for coarsenesses of each point
!     coarsemax      o    Maximum coarseness at the boundary
!     coarsemin      o    Minimum coarseness at the boundary
!     coor           i    coordinate array  ( Standard SEPRAN )
!     ibuur          o    array with neighbours
!     istart         o    array containing the starting addresses of
!                         the neighbours in ibuur
!     kbound         i    array with boundary elements (lines)
!     nbound         i    number of lines in kbound
!     ncoar          i    Number of internal points in surface
!     npoint         i    number of points already created
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision afstan, afst, xp, yp, coarsemean
      integer itotal, i, jknoop, k, i1, i2, no, ni, ntal, ierror

!     afst           distance between two points
!     afstan         sum of distances for point considered
!     coarsemean     Average coarseness at the boundary
!     i              loop variable
!     i1             node number
!     i2             node number
!     ierror         count for errors
!     itotal         local length of istart
!     jknoop         local node number
!     k              loop variable
!     ni             address in neighbour array
!     no             address in neighbour array
!     ntal           number of neighbours of point considered
!     xp             x-coordinate
!     yp             y-coordinate
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS     Deconcatenate name from string of calling subroutines
!     EROPEN     Concatenate name to string of calling subroutines
!     ERREAL     Place real value in error array
!     ERRSUB     Submit error message
!     MSHO02     Place neighbour points in array of neighbours
!     MSHO03     Compute Euclidian distance
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     905    Number of neighbours of a boundary point is other than
!            zero or two !
! **********************************************************************
!
!                       PSEUDO CODE
!
!     Count all neighbours by running through array KBOUND and check
!     the final values
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen( 'msho01' )

!     --- Initialize coarsemin and coarsemax

      coarsemin = rinfin
      coarsemax = -rinfin

!     --- Fill istart

      do i = 1, npoint
         istart(i) = 0
      end do
      do i = 1, nbound
         i1 = kbound(2*i-1)
         i2 = kbound(2*i  )
         istart(i1) = istart(i1) + 1
         istart(i2) = istart(i2) + 1
      end do

!     --- Check and rearrange

      itotal = 0
      do i = 1, npoint
         itotal    =  itotal + istart(i)
         istart(i) =  itotal
      end do

!     --- Clear and fill array ibuur

      itotal = istart(npoint)
      do i = 1, itotal
         ibuur(i) = 0
      end do
      do i = 1, nbound
         i1 = kbound(2*i-1)
         i2 = kbound(2*i  )
         call msho02( istart, ibuur, i1, i2 )
         call msho02( istart, ibuur, i2, i1 )
      end do

!     --- Check ibuur

      ierror = 0
      no     = 0
      do i = 1, npoint
          ni = istart(i)

!         --- ntal is number of neighbours

          ntal   = ni - no

!         --- check number of neighbours

          if ( ntal==0 .or. ntal==2 ) then

!         --- Right number of neighbours !

          else

!         --- Error in number of neighbours

               xp = coor(2*i-1)
               yp = coor(2*i  )
               call erreal( xp, 1 )
               call erreal( yp, 2 )
               call errsub ( 905, 0, 2, 0 )
               ierror = ierror + 1
          end if
          afstan = 0d0
          if ( ntal>0 ) then
               do k = no + 1, no + ntal
                   jknoop = ibuur(k)
                   if ( jknoop>0 ) then
                        call msho03( i, jknoop, coor, afst )
                        afstan = afstan + afst
                   end if
               end do
               coarse(i) = afstan / ntal
               coarsemin = min(coarsemin,coarse(i))
               coarsemax = max(coarsemax,coarse(i))
          else
               coarse(i) = 0
          end if
          no  =  ni
      end do

!     --- Fill coarseness in coar in case this has not been filled

      coarsemean = 0.5d0*(coarsemin+coarsemax)
      do i = 1, ncoar
         if ( coar(3,i)<=1d-6 ) coar(3,i) = coarsemean
      end do  ! i = 1, ncoar

      call erclos( 'msho01' )
      end
