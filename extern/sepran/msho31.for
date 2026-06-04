      subroutine msho31 ( coor, npoint, kelem, nelem, kbound, nbn,
     +                    extquanodes )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.2  date 29-10-2010 Dynamic allocation helparrays
!        version  3.1  date 15-07-2009 Quadratic internal lines allowed
!        version  3.0  date 15-02-2005 Update
!
!   copyright (c) 1990-2010  "Ingenieursbureau SEPRA"
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
!     Create quadratic triangles in case of an existing
!     linear triangular SEPRAN mesh as given in array KELEM.
! **********************************************************************
!
!                       KEYWORDS
!
!     2d
!     mesh_generation
!     triangle
!     quadratic
! **********************************************************************
!
!                       MODULES USED
!
      use mshdummymethods
      use mshconstants
      use msherror
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cconst'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(*)
      integer npoint, kelem(*), nelem, nbn, kbound(*), extquanodes(*)

!     coor        i/o     coordinates array
!     extquanodes  i      helparray for quadratic internal line elements
!     kbound       i      array containing the endnodes of the boundary
!                         elements (length 2*nbn)
!     kelem       i,o     array containing the triangular linear elements
!                         at input and quadratic at output
!     nbn          i      number of boundary elements
!     nelem        i      number of elements
!     npoint      i,o     at input total number of points in triangular
!                         linear mesh, at output inj quadratic mesh
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer error, i, ip(6), meshtp, nummer, iextra, jstart, jstend,
     +        ii, ibrp, itotal, npnold, ibrlen, npt, nbr, ia, ib,
     +        ileng
      integer, allocatable, dimension(:)              :: istart
      integer, allocatable, dimension(:)              :: ibrpnt
      integer, allocatable, dimension(:)              :: kelemh

!     error          Return value of allocate or deallocate
!     i              loop variable
!     ia             first node
!     ib             second node
!     ibrlen         declared length of neighbour array
!     ibrp           guess of number of neighbours for each point
!     ibrpnt         helparray containing the neighbours of each point
!     iextra         helpvariable for case of more than five neighbours
!     ii             loop variable
!     ileng          length helparray for quadratic case
!     ip             array to store node numbers
!     istart         helparray with starting positions for each node in ibrpnt
!     itotal         total number of nodes
!     jstart         start of do variable for neighbours
!     jstend         end of do variable for checking neighbours
!     kelemh         helparray used if node has more than five neighbours
!     meshtp         element type according to SEPRAN rules
!     nbr            number of neighbours
!     npnold         old number of points (before creation of quadratic
!                    elements)
!     npt            helpvariable for number of positions
!     nummer         helpvariable for new points
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERALLOC        Produce error message in case allocate went wrong
!     ERCLOS         Resets old name of previous subroutine of higher level
!     ERDEALLOC      Produce error message in case deallocate went wrong
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSH401         Run through group of elements to determine neighbours
!     MSH402         Place nodes of one line in neighbour array
!     MSH403         Determine number refinement point in line piece
!     MSH406         Place quadratic triangle in element array
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     1274  Error extra points not recognized
!     1275  Not enough space reserved for neighbours
! **********************************************************************
!
!                       PSEUDO CODE
!
!     Create quadratic elements and needed points in two steps
!     First place the old already existing midpoints in array ibrpnt
!     next create the new points and place in ibrpnt
!     finally determine the new elements
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'msho31' )

!     --- Allocate space for helparrays

      allocate ( istart(npoint+1), ibrpnt((npoint+1)*10),
     +           kelemh(npoint), stat = error )
      if ( error/=0 ) call eralloc ( error, 12*npoint+11, 'istart' )
      if ( ierror/=0 ) go to 1000

      istart = 0     !  Clear array istart
      ibrpnt = 0     !  Clear array ibrpnt

      ibrp = 10

!     --- Set used length kelemh

      kelemh(1) = 0

!     --- Special value for meshtp

      meshtp = 3

!     --- Fill arrays

      call msh401( npoint, kelem, nelem, istart, ibrp, ibrpnt, kelemh,
     +             meshtp, 3    , 1    , 0     , 0    )

!     --- Make shift in array for all neighbour points

      iextra = 0
      nummer = 0

!     --- First skip all empty places

      do i = 1 , npoint
         jstart = ibrp * (i-1) + 1
         jstend = ibrp *  i
         do ii = jstart, jstend
            if ( ibrpnt(ii)/=0 ) then

!           --- Place neighbour point

               nummer         = nummer + 1
               ibrpnt(nummer) = ibrpnt(ii)
            end if
         end do
      end do

!     --- Keep total number

      itotal = nummer

!     --- Number of positions needed in ibrpnt is 2 * nummer and number of
!         positions declared is ibrlen

      ibrlen = ibrp * ( npoint + 1 )

!     --- Check length

      if ( ibrlen<2 * nummer ) then

!     --- Error, not enough positions

           call errint ( 2*nummer, 1 )
           call errint ( ibrlen  , 2 )
           call errsub ( 1275, 2, 0, 0 )
           go to 1000
      end if

!     --- Make shift

      npt    = ibrlen

!     --- Run through all basis nodes

      do i = npoint, 1, -1

!     --- Determine number of neighbours of i

         nbr = istart(i)

!        --- Check for extra points

         if ( nbr>ibrp ) then

!        --- There are extra points

            iextra = 1

!           --- Make reservation of positions

            do ib = 1 , nbr-ibrp
               ibrpnt(npt) = 0
               npt = npt - 1
            end do

!           --- Add to itotal

            itotal = itotal + nbr - ibrp
            nbr = ibrp
         end if

         do ib = 1 , nbr
                ibrpnt( npt ) = ibrpnt( nummer )
                npt = npt - 1
                nummer = nummer - 1
         end do

         if ( npt<nummer ) then

!        --- Error, too few positions

            call errsub ( 1274, 0, 0, 0 )
            go to 1000
         end if
      end do

      if ( ibrlen<2 * itotal ) then

!     --- Error, not enough positions declared

         call errint ( 2*itotal, 1 )
         call errint ( ibrlen  , 2 )
         call errsub ( 1275, 2, 0, 0 )
         go to 1000
      end if

!     --- Finally shift all points

      do i = 1 , itotal
         ibrpnt(i) = ibrpnt( i+npt )
      end do

      nummer = itotal
      jstend = 0

      do i = 1, npoint
         jstart = jstend + 1
         jstend = jstend + istart(i)
         istart(i) = jstart
      end do

      itotal           =          0
      istart(npoint+1) = jstend + 1

!     --- Extra neighbours

      if ( kelemh(1)/=0 ) then
         if ( iextra==0 ) then
            call errsub ( 1274, 0, 0, 0 )
            go to 1000
         end if

!        --- Place extra neighbours (in total more than ibrp/2)

         itotal = kelemh(1)

!        --- Place these extra nodes

         do i = 1 , itotal
            ia = -kelemh( 2 * i     )
            ib =  kelemh( 2 * i + 1 )
            call msh402( npoint, ibrp, ibrpnt, istart, ia, ib,
     +                   kelemh )
         end do
      end if

!     --- Check declared length of array ibrpnt

      if ( 2 * nummer>ibrlen ) then
         call errint ( 2*nummer, 1 )
         call errint ( ibrlen  , 2 )
         call errsub ( 1275, 2, 0, 0 )
         go to 1000
      end if

!     --- Set number of points for copiing later

      itotal = istart(npoint+1)-1

      do i = 2, npoint+1
         istart(i) = 2 * istart(i)-1
      end do

!     --- Make place for new nodes

      do i = itotal, 1, -1
         ibrpnt(2*i  ) = 0
         ibrpnt(2*i-1) = ibrpnt(i)
      end do

!     --- Place new nodal-point numbers
!         first the old numbers have to be placed again (determine the number
!         while realizing that several gaps may exist in the boundary)

      do i = 1 , nbn
         ip(1) = kbound( 2*i - 1 )
         ip(2) = kbound( 2*i     )
         ip(3) = ip(1) + 1
         if ( ip(1)<ip(2) ) then
            ip(1) = kbound( 2*i     )
            ip(2) = kbound( 2*i - 1 )
         end if
         jstart = istart( ip(1)     )
         jstend = istart( ip(1) + 1 ) - 2
         do ii = jstart , jstend , 2
            if ( ibrpnt(ii)==ip(2) ) ibrpnt(ii+1) = ip(3)
         end do
      end do

!     --- Determine length of part extra array for internal
!         quadratic elements to be used:

      ileng = (extquanodes(1) - 1)/3

!     --- Extra loop internal quadratic elements

      do i = 1 , ileng
         ip(1) = extquanodes( 3*i - 1 )
         ip(2) = extquanodes( 3*i + 1 )
         ip(3) = extquanodes( 3*i     )
         if ( ip(1)<ip(2) ) then
            ip(1) = extquanodes( 3*i + 1 )
            ip(2) = extquanodes( 3*i - 1 )
         end if
         jstart = istart( ip(1)     )
         jstend = istart( ip(1) + 1 ) - 2
         do ii = jstart , jstend , 2
            if ( ibrpnt(ii)==ip(2) ) ibrpnt(ii+1) = ip(3)
         end do
      end do

!     --- place new nodal-point numbers on all internal element-boundaries

      npnold = npoint

      do i = 1 , npnold
         jstart = istart(i)
         jstend = istart(i+1) - 2
         do ii = jstart , jstend , 2
            ip(1) = i
            ip(2) = ibrpnt( ii )
            ip(3) = ibrpnt( ii + 1 )
            if ( ip(2)>0 .and. ip(3)==0 ) then

!           --- Place new point and compute coordinates

               npoint         = npoint + 1
               ibrpnt( ii+1 ) = npoint
               coor(2*npoint-1)=(coor(2*i-1)+coor(2*ip(2)-1)) / 2d0
               coor(2*npoint  )=(coor(2*i  )+coor(2*ip(2)  )) / 2d0
            end if
         end do
      end do

!     --- Create the new elements

      do i = nelem , 1 , -1
         ip(1) = kelem( 3*i - 2 )
         ip(3) = kelem( 3*i - 1 )
         ip(5) = kelem( 3*i     )
         call msh403( npnold, ip(1), ip(3), ibrpnt, istart, ip(2) )
         call msh403( npnold, ip(3), ip(5), ibrpnt, istart, ip(4) )
         call msh403( npnold, ip(5), ip(1), ibrpnt, istart, ip(6) )
         call msh406( kelem, i,
     +                ip(1), ip(2), ip(3), ip(4), ip(5), ip(6) )
      end do

!     --- Finally deallocate helparrays

      deallocate( istart, ibrpnt, kelemh, stat = error )
      if ( error/=0 ) call erdealloc ( error, 'istart' )

!     --- Close routine

1000  call erclos ( 'msho31' )

      end
