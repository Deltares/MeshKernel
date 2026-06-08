      subroutine mshoce ( jnew, coor, kmeshc, inpelm, nbound, bcord,
     +                    kbndpt, boundary, numcurvboun, npoint, nelem,
     +                    holeinfo, nholes, ncoar, coar, userpoints,
     +                    isurnr, numextcurves, numnodextcurvs,
     +                    curvenumbers, rinput, nuspnt, ndim ) bind(c)
! ======================================================================
!
!        programmer    Niek Praagman
!        version  9.1  date 02-11-2010 Dynamical allocation workarrays
!        version  9.0  date 02-09-2010 Remove intarmsh
!        version  8.0  date 27-01-2010 Extra parameter intarmsh
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
!     Subroutine to make a mesh of triangles for an area which is
!     given by the line elements of BCORD. Elements are linear or qua-
!     dratic.(INPELM = 3 , 6 or 7)
!     MSHOCE is the driving subroutine.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh
!     mesh_generation
!     surface
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
!      include 'SPcommon/cinout'
!      include 'SPcommon/cgobs'
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cconst'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision coor(*), bcord(*), coar(*), rinput(*)
      integer kmeshc(*), inpelm, nbound, kbndpt(*), npoint, nelem,
     +        numcurvboun, boundary(2,numcurvboun),
     +        nholes, holeinfo(2,0:nholes+1),
     +        ncoar, userpoints(*), isurnr, numextcurves,
     +        numnodextcurvs(*), curvenumbers(*), nuspnt, ndim
      logical jnew

!     bcord          i    array containing the coordinates of the boundaries
!     boundary       i    array of size 2 x numcurvboun containing
!                         information of the boundary of the surface
!                         (1,i) contains the i-th curve number (including sign)
!                         (2,i) contains the starting address in kbndpt
!     coar           i    array containing the coordinates as well as the
!                         local coarseness of internal points to be placed
!                         in the surface
!                         Length 3 x ncoar
!                         (1,i):  x-coordinate of point i
!                         (2,i):  y-coordinate of point i
!                         (3,i):  coarseness of point i
!     coor           o    array containing the coordinates of the nodes
!     curvenumbers   i    Contains the curve numbers of the extra internal
!                         curves
!     holeinfo       o    Integer array of size 2 x (nholes+2)
!                         Contains information about holes in surfaces and where
!                         a new part of the boundary starts
!                         Pos (1,i-1) contains the local sequence number of the
!                               first curve of part i provided with a sign
!                               If the sign is positive the boundary is
!                               created counter clockwise, otherwise clockwise
!                         Pos (2,i-1) contains the local sequence number of the
!                               first node of of part i
!                         Pos (1,nholes+1) contains the number of curves of the
!                         boundary
!     inpelm         i    number of points per element
!     isurnr         i    Surface sequence number
!     jnew           i    indication whether the mesh is new (true) or old
!     kbndpt         o    array containing the nodal point numbers of the
!                         boundary
!     kmeshc         o    array containing the tetrahedrons
!     nbound         i    number of points in bcord
!     ncoar          i    Number of internal points in surface
!     ndim           i    Number of Euclidian dimensions in original prob
!     nelem          o    Number of elements
!     nholes         o    Number of holes in the boundary of a surface
!     npoint         o    Number of nodal points
!     numcurvboun    i    Number of curves in the boundary of the surface
!     numextcurves   i    Number of extra (internal) curves
!     numnodextcurvs i    Contains the number of nodes for each internal
!                         curve in the sequence of these curves
!     nuspnt         i    Number of userpoints prescribed by user
!     rinput         i    Real input array
!     userpoints     i    integer array of length ncoar containing the
!                         user point numbers of the extra points
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer error, npunt, jpnt, nbn, inside, i, nbndpt,
     +        istep, lenextracrvs, nnodes

      integer, allocatable, dimension (:) :: extquanodes
      integer, allocatable, dimension (:) :: kbound
      double precision  timebefore
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     error          Return value of allocate or deallocate
!     extquanodes    ?
!     i              loopvariable
!     inside         indicator how many internal regions have to be considered
!     istep          indicator for step: 1=linear, 2=quadratic
!     jpnt           point number
!     kbound         array for boundary lines, length 2*npunt
!     lenextracrvs   Number of nodes in extra internal curves
!     nbn            number of boundary points
!     nbndpt         help varaible to store the number of boundary points
!     nnodes         Number of nodes
!     npunt          number of nodal points
!     timebefore     Time at entrance of subroutine
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERALLOC        Produce error message in case allocate went wrong
!     ERCLOS         Resets old name of previous subroutine of higher level
!     ERDEALLOC      Produce error message in case deallocate went wrong
!     EROPEN         concatenate names for error messages
!     ERSETTIME      Get the CPU time at start of subroutine
!     MSHCOPYBOUN    Copy bcord into coor
!     MSHD7          Places coorindates of centroids of 7-point triangles
!                    in array coor
!     MSHO2D         driving subroutine to triangulate region
!     MSHO31         place extra nodes for quadratic case
!     MSHO32         determine highest node number for corner
!     MSHREP         driving subroutine for repositioning
!     MSHTRIANCENTR  Create centre point in 6-point triangles
!     PRININ         print 1d integer array
!     PRININ1        print 2d integer array
!     PRINRL1        Print 2d real vector
!     PRINTTIME      Print the CPU time at end of subroutine
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     none
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
      call eropen('mshoce')

      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshoce'
         write(irefwr,1) 'isurnr, nbound', isurnr, nbound
         call prinrl1 ( bcord, nbound, 2, 'bcord' )
         write(irefwr,1) 'ncoar, numextcurves', ncoar, numextcurves
         if ( ncoar>0 ) then

!        --- ncoar>0

            call prinin ( userpoints, ncoar, 'userpoints' )
            call prinrl1 ( coar, ncoar, 2, 'coar' )
         end if  ! ( ncoar>0 )
         if ( numextcurves>0 ) then

!        --- numextcurves>0

            lenextracrvs = 0
            do i = 1, numextcurves
               lenextracrvs = lenextracrvs+numnodextcurvs(i)
            end do  ! i = 1, numextcurves
            call prinin ( numnodextcurvs, numextcurves,
     +                    'numnodextcurvs' )
            call prinin ( curvenumbers, numextcurves, 'curvenumbers' )
            call prinrl1 ( bcord(2*nbound+1), lenextracrvs, 2,
     +                     'coordinates extra curves' )
         end if  ! ( ncoar>0 )
  1      format ( a, 1x, (10i6) )
      end if  ! ( debug )
      if ( ierror/=0 ) go to 1000
      if ( itime>0 ) call ersettime ( timebefore )

!     --- Estimated number of points
      npunt = npoint
!     --- Fill array coor with boundary points
!         Allocate temporary space for help array kbound
      allocate( kbound( 2 * npunt ), stat = error )
      if ( error/=0 ) call eralloc ( error, 2*npunt, 'kbound' )

      call mshcopyboun ( jpnt, nbound, bcord, coor, kbndpt,
     +                   kbound, inside, nbn, .true. )

      if ( ierror/=0 ) go to 1000

      if ( jnew ) then
         npoint = jpnt
         istep = 1

!        --- Check for quadratic

         if ( inpelm==6 .or. inpelm==7 ) then

!        --- Adjust istep, boundary array and eventually extra curves:

            istep = 2

            nbn = nbn / 2

            do i = 1, nbn
               kbound(2*i-1) = kbound(4*i-3)
               kbound(2*i  ) = kbound(4*i  )
            end do

!           --- Extra curves?

            nnodes = 1

            if ( numextcurves>0 ) then

!           --- Determine the length of extra array for curves:

               nnodes = 0

               do i=1, numextcurves

!              --- Run through line i of the internal curves:

                  nnodes = nnodes + numnodextcurvs(i)
               end do

!              --- Extra array needed to place temporarily the quadratic
!                  elements of the internal lines:
!                  Length needed: 2 * nnodes

            end if

            allocate( extquanodes( 2 * nnodes ), stat = ierror)
            if ( error/=0 ) call eralloc (error,2*nnodes,'extquanodes')
            extquanodes(1) = 0

            if ( debug ) then
               write(irefwr,*) 'extra points = ',nnodes
               write(irefwr,*) 'numextcurves = ',numextcurves
            end if
         else
            allocate( extquanodes( 1 ), stat = ierror)
            extquanodes(1) = 0
         end if

!        --- Set initial number of points

         nbndpt = nbound

         call msho2d ( coor, npoint, kbound, nbn, kmeshc,
     +                 nelem, boundary, numcurvboun, npunt,
     +                 inside, holeinfo, nholes, .true.,
     +                 kbndpt, nbndpt, coar, ncoar,
     +                 userpoints, bcord(2*nbound+1), numextcurves,
     +                 numnodextcurvs, curvenumbers, istep,
     +                 extquanodes, isurnr, rinput,
     +                 nuspnt, ndim )

         if ( debug ) then
            write(irefwr,*) 'Contents array kbndpt after calling msho2d'
            write(irefwr,*) (kbndpt(i),i=1,nbndpt)
            write(irefwr,*) 'First part coor after calling msho2d'
            do i=1, nbndpt
               write(irefwr,*) i,coor(2*i-1),coor(2*i)
            end do
            write(irefwr,*) 'end first part coordinates array coor'
            write(irefwr,*) 'number of extra curves = ',numextcurves
            write(irefwr,*) 'number of extra nodes on curves'
            write(irefwr,*) (numnodextcurvs(i),i=1,numextcurves)
            write(irefwr,*) 'it concerns curves with numbers:'
            write(irefwr,*) (curvenumbers(i),i=1,numextcurves)
         end if

         if ( ierror/=0 ) go to 1000

!        --- Special treatment if elements are quadratic

         if ( inpelm==6 .or. inpelm==7 ) then

!        --- inpelm = 6 or 7, extra nodes in elements

            call msho31 ( coor, npoint, kmeshc, nelem, kbound, nbn,
     +                    extquanodes)
            if ( inpelm==7 ) then

!           --- inpelm = 7, also fill centroid

                !AvD: We only do triangles (inpelm=3), so disable calls to mshtriancentr and mshd7 (routines are missing)
                !AvD: call mshtriancentr ( kmeshc, npoint, nelem, inpelm )
                if ( ierror/=0 ) go to 1000
                !AvD: call mshd7 ( kmeshc, coor, nelem, 2 )
            end if
            if ( ierror/=0 ) go to 1000
         end if

         deallocate(  extquanodes, stat = error )
         if ( error/=0 ) call erdealloc ( error, 'extquanodes' )

      else

!     --- jnew = .FALSE. hence only repositioning

         nbndpt = jpnt

         !AvD: Disable call mshrep (routine is missing): only used when generating
         !     a mesh by REpositioning (given current mesh with only changed boundary
         !     points, move inner points).
         !AvD: call mshrep ( coor, npoint, kmeshc, nelem, inpelm, nbndpt )

!        --- MSHREP uses ordering of GENERAL hence make a renumbering if
!            inpelm>3

         if ( inpelm==6 .or. inpelm==7 ) then

!        --- inpelm = 6 or 7, extra nodes in elements

            call msho32(kmeshc,nelem,inpelm,npunt)
            if ( ierror/=0 ) go to 1000

            if ( npunt>jpnt ) then
                 npoint = npunt
            else
                 npoint = jpnt
            end if

!           --- Adjust boundary points value and boundary node numbers

            nbn = nbn / 2

!           --- Run through all pieces

            do i = 1 , nbn
                kbound(2*i-1) = kbound(4*i-3)
                kbound(2*i  ) = kbound(4*i  )
            end do

!           --- Extra internal curves?

            if ( numextcurves>0 ) then

!           --- Determine the length of extra array for curves:

               nnodes = 0

!              --- Run through internal curves

               do i=1, numextcurves

!              --- Run through line i of the internal curves:

                  nnodes = nnodes + numnodextcurvs(i)
               end do

!              --- Extra array needed to place temporarily the quadratic
!                  elements of the internal lines:
!                  Length needed: 2 * nnodes

               allocate( extquanodes( 2 * nnodes ), stat = ierror)
               if (error/=0) call eralloc (error,2*nnodes,'extquanodes')
            else

!           --- No extra nodes internal curves:

               allocate( extquanodes( 1 ), stat = ierror)
               extquanodes(1) = 0
            end if

            call msho31( coor, npoint, kmeshc, nelem, kbound, nbn,
     +                   extquanodes)

            if ( inpelm==7 ) then

!           --- inpelm = 7, also fill centroid

                !AvD: We only do triangles (inpelm=3), so disable calls to mshtriancentr and mshd7 (routines are missing)
                !AvD: call mshtriancentr ( kmeshc, npoint, nelem, inpelm )
                if ( ierror/=0 ) go to 1000
                !AvD: call mshd7 ( kmeshc, coor, nelem, 2 )
            end if
            if ( ierror/=0 ) go to 1000

            deallocate(  extquanodes, stat = error )
            if ( error/=0 ) call erdealloc ( error, 'extquanodes' )
         end if

      end if

      if ( igobs==0 .and. inpelm>4 ) then

!     --- Quadratic elements, copy boundary again in order to get a really
!         curved boundary

         call mshcopyboun ( jpnt, nbound, bcord, coor, kbndpt,
     +                      kbound, inside, nbn, .false. )

      end if

      deallocate( kbound, stat = error )
      if ( error/=0 ) call erdealloc ( error, 'kbound' )
      if ( itime==1 ) call printtime('TRIANGLE', timebefore)
1000  call erclos ( 'mshoce' )
      if ( debug ) then
         write(irefwr,*) 'Debug information end mshoce'

!        --- Debug information

         write(irefwr,1) 'nholes, npoint, ncoar, nelem, nbound',
     +                    nholes, npoint, ncoar, nelem, nbound
         call prinin ( kbndpt, nbound, 'kbndpt' )
         if ( nholes>0 ) then
            call prinin1 ( holeinfo(1,0), nholes+2, 2, 'holeinfo' )
         end if  ! ( nholes>0 )
         call prinin1 ( kmeshc, nelem, inpelm, 'kmeshc' )

         write(irefwr,*) 'End msho2d with ',npoint,' nodes'
         do i=1, npoint
            write(irefwr,*) i,'x-y',coor(2*i-1),coor(2*i)
         end do
         write(irefwr,*) 'End mshoce'
      end if
      end
