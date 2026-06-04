      subroutine mshcopyboun ( jpnt, nbound, bcord, coor, kbndpt,
     +                         kbound, inside, nbn, fillkbound)
! ======================================================================
!
!        programmer    Guus Segal
!        version  2.0  date 08-02-2001 Extra parameter ihelp
!        version  1.0  date 03-04-1999
!
!   copyright (c) 1999-2001  "Ingenieursbureau SEPRA"
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
!     Fill the arrays kbndpt, kbound and coor with information concerning
!     the nodes on the boundary
!     kbndpt and kbound are only filled if fillkbound is true
! **********************************************************************
!
!                       KEYWORDS
!
!     boundary
!     copy
!     mesh
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
      integer jpnt, nbound, kbndpt(*), kbound(*), inside, nbn!,
!     +        ihelp(0:*)
      double precision bcord(2,*), coor(2,*)
      logical fillkbound

!     bcord          i    array containing the coordinates of the boundaries
!     coor           o    array of length ndim x npoint containing the
!                         co-ordinates of the nodal points with respect to the
!                         global numbering
!                         In this subroutine only the boundary is filled.
!                         It is supposed that the points on the boundary are
!                         also the first points in coor
!     fillkbound     i    if true the array kbndpt and kbound must be filled
!     ihelp               integer work array of length nparts+1, where nparts
!                         is the number of closed boundaries
!                         ihelp(0) = 0
!                         ihelp(j) is the last relative number of part j with
!                         respect to kbndpt
!     inside         o    indicator how many internal regions have to be
!                         considered
!     jpnt           o    last point number used
!     kbndpt         o    array containing the nodal point numbers of the
!                         boundary
!     kbound         o    help array for boundary lines, length 2*npunt
!     nbn            o    number of boundary points
!     nbound         i    number of points in bcord
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision, allocatable, dimension(:)     :: ihelp

      double precision dx, dy, ref, xst, yst, xp, yp, ds
      integer i, j, jst, ifirst, ilast, nparts
      logical finished

!     ds             distance
!     dx             x-distance
!     dy             y-distance
!     finished       Indicates if the loop has been finished (true) or not
!                    (false)
!     i              Counting variable
!     ifirst         First node in part
!     ilast          Last node in part
!     j              Counting variable
!     jst            nodal point number of first point of local boundary
!     nparts         number of closed parts the boundary consists of
!     ref            reference length
!     xp             x-coordinate of point
!     xst            x-coordinate of first point
!     yp             y-coordinate of point
!     yst            y-coordinate of first point
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
!     initialize parameters
!     Compute reference distance
!     Detect the number of parts the boundary consists of and store information
!     in ihelp
!
!     For all points in the boundary do
!     if ( first point on a new part of the boundary ) then
!        raise node number
!        Copy coor
!        Fill kbndpt if necessary
!     else ( next point on the boundary )
!        Fill kbound
!        If ( last point of closed boundary ) then
!           reset j
!           Fill kbndpt if necessary
!           Fill kbound if necessary
!        else ( not last point of closed boundary )
!           Fill kbndpt if necessary
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- initialize parameters

      nbn    = 0
      jpnt   = 0
      inside = -1

!     --- Allocate helparrays:

      ! Worst-case estimate of nr of closed boundary polygons
      nparts = 1+.5*nbound
      allocate( ihelp  ( 0:nparts ) )

!     --- Compute reference distance

      dx = bcord(1,2) - bcord(1,1)
      dy = bcord(2,2) - bcord(2,1)

      ref = 1d-5 * sqrt( dx*dx + dy*dy )
      j = 0

!     --- Detect number of closed parts the boundary consists of and store in
!         ihelp

      ihelp(0) = 0
      finished = .false.
      nparts = 0
      ifirst = 1
      do while ( .not. finished )

!     --- Not all parts have been found
!         Store starting node
!         ilast is last point in closed boundary

         xst = bcord(1,ifirst)
         yst = bcord(2,ifirst)

!        --- Find last point that coincides with starting point

         ilast = nbound
         do i = ifirst+1, nbound

            xp = bcord(1,i)
            yp = bcord(2,i)
            dx = xp - xst
            dy = yp - yst
            ds = sqrt( dx*dx + dy*dy )
            if ( ds.le.ref ) ilast = i

         end do

!        --- Raise nparts

         nparts = nparts+1
         ihelp(nparts) = ilast
         if ( ilast.eq.nbound ) finished = .true.
         ifirst = ilast+1

      end do

!     --- Loop over all parts

      do j = 1, nparts
         ifirst = ihelp(j-1)+1
         ilast  = ihelp(j)
         do i = ifirst, ilast
            xp  =  bcord(1,i)
            yp  =  bcord(2,i)

            if ( i.eq.ihelp(j-1)+1 ) then

!           --- first point on a new part of the boundary

               if ( inside.lt.1 ) inside = inside + 1
               jst  =  jpnt+1
               jpnt = jpnt + 1
               coor(1,jpnt) = xp
               coor(2,jpnt) = yp
               if ( fillkbound ) kbndpt(i) = jpnt

            else

!           --- next point on the boundary

               if ( fillkbound ) then

!              --- fill kbound

                  nbn  =  nbn + 1
                  kbound(2*nbn-1) = jpnt
                  kbound(2*nbn  ) = jpnt + 1

               end if

               if ( i.eq.ihelp(j) ) then

!              --- end point of a closed boundary found

                  if ( fillkbound ) then
                     kbound(2*nbn) = jst
                     kbndpt(i)     = jst
                  end if

               else

!              --- not end point of a closed boundary

                  jpnt = jpnt + 1
                  coor(1,jpnt) = xp
                  coor(2,jpnt) = yp
                  if ( fillkbound ) kbndpt(i) = jpnt

               end if

            end if

         end do

      end do
      
      deallocate( ihelp )


      end
