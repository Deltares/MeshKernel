      subroutine msho2d ( coor, npoint, kbound, nbound, kmeshc, nelem,
     +                    boundary, numcrvboun, npunt, inside, holeinfo,
     +                    nholes, reposition, kbndpt, nbndpt, coar,
     +                    ncoar, userpoints, cocurv, ncurvs, curves,
     +                    crvnrs, istep, extquanodes, isurnr, rinput,
     +                    nuspnt, ndim )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  8.2  date 07-01-2011 Improvement tuning nx, ny
!        version  8.1  date 29-10-2010 Dynamical allocation helparrays
!        version  8.0  date 02-09-2010 Remove intarmsh
!        version  7.4  date 11-03-2010 Use SEPRAN estimate for error 901
!
!   copyright (c) 1989-2011  "Ingenieursbureau SEPRA"
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
!     Subroutine to make a mesh of triangles for a polygon which is
!     given by the lines of KBOUND.
!     Besides the boundary internal nodes and internal curves are also
!     allowed to prescribe internally, i.e. within the domain considered
!     coarseness and fixed positions.
!     If the user prescribes coarsenesses causing great size-differences
!     in subsequent elements then give a warning and an advice what
!     coarseness to use in the submitted userpoints
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
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
!      include 'SPcommon/cmcdpr'
!      include 'SPcommon/cconst'
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cplot'
!      include 'SPcommon/cmesh'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer npoint, npunt, kbound(*), nbound, kmeshc(3,*), nelem,
     +        numcrvboun, boundary(2,*), istep, inside, ncoar,
     +        nuspnt, ndim, nholes, holeinfo(2,0:nholes+1),
     +        userpoints(*), isurnr, nbndpt, ncurvs, curves(ncurvs),
     +        crvnrs(ncurvs), kbndpt(*), extquanodes(*)
      double precision coor(2,*), coar(3,*), cocurv(2,*), rinput(*)
      logical reposition

!     boundary       i    array of size 2 x numcrvboun containing
!                         information of the boundary of the surface
!                         (1,i) contains the i-th curve number (including sign)
!                         (2,i) contains the starting address in kbndpt
!     coar           i    array containing coordinates and coarseness of
!                         special points to be used in later calculations:
!                         positions of these points are fixed!
!     cocurv         i    array containing the crds of the fixed line nodes
!     coor          i/o   array containing the coordinates of the points
!     crvnrs         i    array with the the (usergiven)numbers of the internal
!                         curves
!     curves         i    the curvnumbers of the extra internal curves
!     extquanodes    i    helparray to store quadratic elements of internal
!                         lines
!     holeinfo       i    Integer array of size 2 x (nholes+2)
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
!     inside         i    indicator whether there are inside areas or not
!     istep          i    indicator lineair (=1) or quadratic (=2)
!     isurnr         i    Surface sequence number
!     kbndpt         i    array containing the nodal point numbers of the
!                         boundary
!     kbound         i    array containing the boundary lines
!     kmeshc         o    array containing the linear triangles created
!     nbndpt         i    Number of nodal points in boundary, length of array
!                         kbndpt
!     nbound         i    number of boundary elements
!     ncoar          i    number of special points in array coar
!     ncurvs         i    number: extra internal curves submitted by user
!     ndim           i    dimension number of original problem
!     nelem          o    number of triangles in kmeshc
!     nholes         o    Number of holes in the boundary of a surface
!     npoint         i    number of points already created
!                    o    number of points with new points added
!     npunt          i    estimation of number of points in final mesh
!     numcrvboun     i    ?
!     nuspnt         i    Number of user points
!     reposition     i    if true, repositioning may be carried out
!     rinput         i    real input array standard sepran for Probdf
!     userpoints     i    integer array of length ncoar containing the
!                         user point numbers of the extra points
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision eps, dismin, dismax, dism,
     +                 xminloc, xmaxloc, yminloc, ymaxloc,
     +                 xstart, ystart, angle, factor, maxratio,
     +                 surf, surf1, surf2, dx, dy, disold, verh, 
     +                 coa, dist, dista, dist1, adis, bdis, cdis,
     +                 angle1, angle2, valare, xn, yn, xm, ym, e1, e2,
     +                 xnx, yny, x, y, coarsemin, coarsemax, coarmin,
     +                 angleh, tran, xmint, ymint, ratiop
      integer i, nx, ny, ncube, nkmeshc, nrepos, nherha, iperm, kstap,
     +        ih, ipoint, nochan, i1, i2, i3, i4, iex1, iex2, j1, j2,
     +        j3, j4, ipotn1, ipotn2, ja, no, kdrie, icheck, jpn,
     +        ii1, ii2, inear, n1, n2, nc, itot, leng, jelem, nelemi,
     +        ielem, icase, ma, ichan, ifill, nipnt, icount,
     +        lenbnd, imultiply, jcoars, npndef, npntmp,
     +        neldef, neltmp, icancel, idebug
      double precision, allocatable, dimension(:)     :: chelp
      double precision, allocatable, dimension(:)     :: userco
      double precision, allocatable, dimension(:)     :: coaval
      double precision, allocatable, dimension(:,:)   :: coortmp
      double precision, allocatable, dimension(:)     :: cube
      integer, allocatable, dimension(:)              :: icube
      integer, allocatable, dimension(:)              :: itri
      integer, allocatable, dimension(:)              :: jcube
      integer, allocatable, dimension(:,:)            :: meshtmp
      integer, allocatable, dimension(:)              :: kinbnd
      integer, allocatable, dimension(:)              :: kstapl
      logical debug

!     adis           distance of two nodes
!     angle          reference angle value
!     angle1         boundary angle
!     angle2         boundary angle
!     angleh         helpvariable
!     bdis           distance of two nodes
!     cdis           distance of two nodes
!     chelp          help array for coarsenesses of points, and
!                    reference values for areas
!     coa            coarseness
!     coarmin        min coarseness in special points
!     coarsemax      Maximum coarseness at the boundary
!     coarsemin      Minimum coarseness at the boundary
!     coaval         local array to store the coarsenesses of userpoints
!     coortmp        helparray to store coordinates of "best" mesh
!     cube           local array of size nx * ny to store a rectangular 
!                    grid of containing the required coarsenesses  
!     debug          If true debug statements are carried out otherwise
!                    they are not
!     dism           mean value distances to be used for coarseness frame
!     dismax         greatest distance in starting polygon
!     dismin         smallest distance in starting polygon
!     disold         old distance
!     dist           distance
!     dist1          distance
!     dista          distance
!     dx             distance
!     dy             distance
!     e1             first component inside pointing normal
!     e2             second component inside pointing normal
!     eps            accuracy
!     factor         reference factor value
!     i              loop variable
!     i1             node number
!     i2             node number
!     i3             node number
!     i4             node number
!     icancel        counter
!     icase          case indicator
!     ichan          helpvariable for number of changings
!     icheck         indicator for crossing
!     icount         Counter
!     icube          array containing for each node its cube number
!     idebug         Checking is performed if icount>idebug
!     ielem          loop variable for elements
!     iex1           node to be considered, neighbour
!     iex2           node to be considered, neighbour
!     ifill          indicator whether diagonals of quadrilaterals
!                    have been changed (1) or not (0)
!     ih             loop variable
!     ii1            node to be considered
!     ii2            node to be considered
!     imultiply      Multiplication factor
!                    In this program this factor is arbitrarely set to
!                    10 log(coarsemax/coarsemin)
!                    It is used to define the number of searches for
!                    "acceptable" elements
!     inear          indicator how near a boundary point is
!     iperm          indicator whether line is permitted or not
!     ipoint         node number
!     ipotn1         potential point
!     ipotn2         potential point
!     itot           help for summing
!     itri           array containing information whether a nodal point np is
!                    in the actual boundary (itri(np)>0) or not (=0)
!     j1             node number
!     j2             node number
!     j3             node number
!     j4             node number
!     ja             help indicator
!     jcoars         indicator whether coarseness is given (1) or not (0)
!     jcube          array indicating whether a cube is outside the volume (0),
!                    is a boundary cube ( 1 ) or is inside ( 2 )
!     jelem          element number
!     jpn            node number
!     kdrie          number of boundary element for crossing
!     kinbnd         array with extra internal boundary pieces
!     kstap          actual length of array KSTAPL
!     kstapl         array containing all actual boundary elements
!     lenbnd         length axtra internal curvs part in kinbnd
!     leng           actual length of repos array
!     ma             maximum value
!     maxratio       maximum ratio allowed in two subsequent elements
!     meshtmp        helparray to store temporarily elements of "best" mesh
!     n1             actual number of blocks in x-direction
!     n2             actual number of blocks in y-direction
!     nc             cube number
!     ncube          number of reference blocks
!     neldef         defined length of helpelementsarray
!     nelemi         Number of internal cells (elements)
!                    Hence, number of boundary cells = nelem - nelemi
!     neltmp         helpvariable allocation array of elements
!     nherha         count of redivisions performed
!     nipnt          number of initial points
!     nkmeshc        estimated number of elements in MESH
!     no             help indicator
!     nochan         change indicator
!     npndef         defined length of defined help pointsarray
!     npntmp         helpvariable allocation array of points
!     nrepos         estimated length of repositioning array
!     nx             number of reference blocks in x-direction
!     ny             number of reference blocks in y-direction
!     ratiop         min ratio of temporary mesh
!     surf           area value
!     surf1          area value
!     surf2          area value
!     tran           scaling parameter
!     userco         local helparray to store coordinates of all userpoints         
!     valare         reference areavalue
!     verh           ratio
!     x              x-coordinate
!     xm             x-coordinate help node on boundary
!     xmaxloc        extreme (max) x-value for the whole area
!     xminloc        extreme (min) x-value for the whole area
!     xmint          min x for transformation
!     xn             x-coordinate new point
!     xnx            x-coordinate second new point
!     xstart         x-reference coordinate
!     y              y-coordinate
!     ym             y-coordinate help node on boundary
!     ymaxloc        extreme (max) y-value for the whole area
!     yminloc        extreme (min) y-value for the whole area
!     ymint          min y for transformation
!     yn             y-coordinate new point
!     yny            y-coordinate second new point
!     ystart         y-reference coordinate
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSHCHKSTAPL    Check the contents of array kstapl
!     MSHCURVINTERS2 Check if edges in boundary of surface do not intersect
!     MSHO01         Check number of neighbours for all boundary points and
!                    compute coarseness for all points
!     MSHO03         Compute distance between two points
!     MSHO04         Determine extreme coordinate values
!     MSHO05         Determine min and max distance for total area
!     MSHO06         Determine coarseness for each triangle
!     MSHO07         Check whether point belongs to given area
!     MSHO08         Check whether all elements are given right and adjust
!                    boundary elements array
!     MSHO09         Compute midpoint line and the normal vector
!     MSHO10         Fill element in in array and adjust boundary array
!     MSHO11         compute angle
!     MSHO12         Compute neighbour nodes and angles
!     MSHO13         Check for common points with other lines
!     MSHO14         Check whether new point is not too near to
!                    old point
!     MSHO15         Fill element with new point in in array
!     MSHO16         Compute number of neighbours of points
!     MSHO18         Perform repositioning
!     MSHO19         Compute area of triangle
!     MSHO20         Change ordering of boundary
!     MSHO21         Compute reference value for area
!     MSHO22         Check form of triangle
!     MSHO24         Check whether line has point in common with boundary
!     MSHO25         Place smallest element at top of boundary array
!     MSHO26         Determine triangle with given side
!     MSHO27         Determine best diagonal in quadrilateral
!     MSHO28         Check whether boundary points are inside triangle
!     MSHO29         Remove internal points with three or four neighbours
!     MSHO30         Make boundar array for all elements surrounding
!                    element indicated
!     MSHO33         Determine the ratio (2*Rout/Rin) for all triangles of mesh
!     MSHO35         Place special nodes and nodes and elements around
!                    these nodes
!     MSHO36         Add nodal points of internal lines to coor and adjust
!                    array of boundaries kstapl for mesh generation
!     MSHO38         check the coarseness in the total domain considered
!                    in relation to smoothness requirements
!     MSHO40         Check whether all barycenters of triangles are only
!                    belonging to their "own" triangle
!     MSHO42         Check whether line i - j is a piece in one of the internal
!                    curves
!     MSHTRANS2DSUR  Transform 2d region to region of unit length in first
!                    quadrant
!     PLOTBOUNINTER  Plot actual boundary during the process of creation
!                    of a mesh
!     PLOTMESHINTER  Plot actual mesh during the process of creation of a mesh
!     PRININ         print 1d integer array
!     PRININ1        print 2d integer array
!     PRINRL1        Print 2d real vector
! **********************************************************************
!
!                       I/O
!
!     none
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     900   Estimated number of points too small
!     901   Estimated length of element array too small.
!     902   No convergence in MESH generation.
!     903   Either connections array not correct or repositioning
!           array too small
! **********************************************************************
!
!                       PSEUDO CODE
!
!     MSHO2D is meant to divide an arbitrary polygon in triangles
!     The method used is the so-called advancing boundary method
!     Consider all parts that together constitute the boundary of the area
!     that has to be triangulated.
!     A global coarseness is computed over the whole area using a Laplacian
!     method.
!     If special points (array coar) are submitted the division is started by
!     creating special elements around these points.
!     Next the smallest part (segment) of the boundary is determined.
!     First it is considered whether there can be created a triangle with one of
!     the two neighbouring segments
!          If so the new boundary is defined and the
!             process starts again.
!          If not the process goes on
!     A new point is defined using the inside pointing normal. It is con-
!     sidered whether this new point is really inside and far enough away
!     from the other boundary lines.
!          If not so it is considered whether a triangle can be formed using
!          one of the already created points.
!                 If so the boundary is redefined and the process starts
!                 again.
!                 If not the actual considered part is placed at the back
!                 of the boundary array and another part is considered first.
!          If so the new boundary is defined and the process starts.
!     The process ends when there are no parts in the boundary, i.e. the area
!     has been divided totally in triangles or when it is not possible to cre-
!     ate new triangles. In that case a message is given.
!     If the process has terminated normally, a repositioning is performed.
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen( 'msho2d' )
      debug = .false.

!     --- debug = .false.

      if ( debug ) write(irefwr,*) 'At start of msho2d, npunt =',npunt
      idebug = 10000000
      icount = 0

!     --- Copy estimate for number of elements:

      !nkmeshc = nelem
      nkmeshc = 2000000 ! nelem is uninitialised

!     --- Initialize

      neltmp = 0
      npntmp = 0

!     --- Allocate helparrays:

      allocate( chelp  ( npunt ) )
      allocate( icube  ( npunt ) )
      allocate( itri   ( npunt ) )
      allocate( kstapl ( 10 * npunt + 20 ) )

!     --- Fill arrays userco and coaval:

      allocate( userco( 2 * nuspnt ) )
      allocate( coaval( 2 + nuspnt ) )

!     --- Coordinates:

      do i=1, nuspnt
         userco(2*i-1) = rinput(1 + ndim*i)
         userco(2*i  ) = rinput(2 + ndim*i)
      end do

!     --- Coarseness:

      ! AvD: enable this to use coarseness points and coarseness values in rinput array (if nuspnt>0)
      jcoars = 0 !intarmsh(44)

      if ( jcoars==1 ) then

!     --- Coarseness is given:

         do i=1, nuspnt+2
            coaval(  i  ) = rinput(2 + ndim*nuspnt + i)
         end do
      else

!     --- No coarseness given, hence no actions:

         coaval( 1 ) = 1d0
      end if

      if ( debug ) then

!     --- Give coordinates and coarsenesses submitted userpoints:

         write(irefwr,*) 'in msho2d, nuspnt = ',nuspnt
         write(irefwr,*) 'coordinates and coarseness:'
         write(irefwr,*) 'overall coarseness',coaval(1)
         write(irefwr,*) 'maxratio          ',coaval(nuspnt + 2)
         do i=1, nuspnt
            write(irefwr,*) i,userco(2*i-1),userco(2*i),coaval(i+1)
         end do
      end if

!     --- Replace all nodes temporarily to the first (positive) quadrant:
!         (Determine translation parameters x- and y-tran)

      call mshtrans2dsur( coor, npoint, xmint, ymint, tran,
     +                    ncoar, coar, ncurvs, curves, cocurv,
     +                    userco, nuspnt )

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msho2d'
         write(irefwr,*) 'Nodenumbers and coordinates'
         do i=1, npoint
            write(irefwr,*) i,'x = ',coor(1,i),'y = ',coor(2,i)
         end do
         write(irefwr,*) 'Boundary pieces'
         write(irefwr,*) 'End boundary pieces'
      end if
      if ( ierror/=0 ) goto 1000

!     --- Set accuracy

      eps = 10 * epsmac
      xn  = 0d0
      yn  = 0d0

!     --- Check whether each line in KBOUND is connected
!         and determine for each point a coarseness-value
!         (this check concerns the outer boundary and eventually holes
!          as given bij bcord. Internal lines are treated later)

      if ( debug ) then
         write(irefwr,*) 'MSHO2D array coar contains',ncoar,'points'
         write(irefwr,*) 'these points are internal in surf',isurnr
         do i=1, ncoar
            write(irefwr,*) i,'crds',( coar(i1,i),i1=1,3)
         end do
         write(irefwr,*) 'end points array coar for surf',isurnr
      end if

!     --- Compute the coarsenesses:

      call msho01 ( kbound, nbound, itri, kstapl, chelp, coor,
     +              npoint, coarsemin, coarsemax, coar, ncoar )
      if ( ierror/=0 ) goto 1000

!     --- Check that edges in array kbound do not intersect each other

      call mshcurvinters2 ( kbound, nbound, coor, isurnr )
      if ( ierror/=0 ) goto 1000

!     --- Place boundary elements in kstapl:

      kstap = nbound

!     --- Run through all pieces:

      do i = 1, 2*nbound
         kstapl(i)    = kbound(i)
      end do

!     --- Set number of boundary points before internal lines are added

      nipnt = npoint

!     --- Eventually print contents kstapl:

      if ( debug ) then
         write(irefwr,*) 'Array kstapl'
         do i = 1, kstap
            write(irefwr,*) 'Piece',i,':',kstapl(2*i-1),'-',kstapl(2*i)
         end do
      end if

!     --- Internal lines?

      if ( ncurvs>0 ) then

!     --- Place points of extra internal lines in kstapl and check
!         the boundaries again (was already checked):

         call msho36 ( coor, npoint, istep, ncurvs, curves, cocurv,
     +                 kstapl, kstap, kbndpt, nbndpt, coarsemin,
     +                 extquanodes, chelp )
         if ( debug ) then

!        --- Give computed coarsenesses extra nodes due to internal
!            lines

            do i=nipnt+1, npoint
               write(irefwr,*) 'Node',i,'coarseness',chelp(i)
            end do
         end if

!        --- Check for intersections:

         call mshcurvinters2 ( kstapl, kstap, coor, isurnr )
         if ( ierror/=0 ) then
            write(irefwr,*)
     +      'Intersection(s) found due to internal lines.
     +       Error! Please check your input'
            goto 1000
         end if

!        --- Eventually print kstapl again:

         if ( debug ) then
            write(irefwr,*) 'Start MSHO2D with ',npoint,' points'
            do i=1, npoint
               write(irefwr,*) i,'crd-s',(coor(ierror,i),ierror=1,2)
            end do
            write(irefwr,*) 'Check! ok? Check kstapl yourself '
            write(irefwr,*) 'Array kstapl, kstap = ',kstap
            do i=1, kstap
               write(irefwr,*) i,kstapl(2*i-1),kstapl(2*i)
            end do
         end if
         ierror = 0
      end if

!     --- Place all boundary pieces in new array kinbnd,
!         determine length (use all kstapl pieces) and allocate:

      lenbnd = 2 * kstap

!     --- Kinbnd has to be kept for later use, allocate:

      allocate ( kinbnd ( lenbnd ) )

!     --- Fill array:

      do i = 1, lenbnd
         kinbnd(i) = kstapl(i)
      end do

!     --- Eventually print array kinbnd:

      if ( debug ) then
         write(irefwr,*) 'array kinbnd'
         write(irefwr,*) (kinbnd(i),i=1,lenbnd)
         write(irefwr,*) 'end array kinbnd'
      end if

!     --- The minimum and maximum coarseness at the boundary is used to define
!         the multiplication factor

      imultiply = max(1,nint(log10(coarsemax/coarsemin)))
      if ( debug )
     +    write(irefwr,*) 'coarsemin, coarsemax, imultiply ',
     +                     coarsemin, coarsemax, imultiply
      if ( ierror/=0 ) goto 1000

!     --- Determine extreme coordinate-values for this surface

      call msho04 ( coor, npoint, xminloc, xmaxloc, yminloc, ymaxloc )
      if ( ierror/=0 ) goto 1000

!     --- Determine max and min distance for total area
!         ( First set nodes internal curves etc to zero)

      call msho05 ( chelp, nipnt, dismin, dismax )

      if ( ierror/=0 ) goto 1000

!     --- Determine global frame
!         Determine number of reference quadrilaterals for this frame
!         Use mean value for coarse frame:

      dism = (dismax + 2 * dismin) / 3d0

!     --- Compute number of "boxes":

      nx = int( (xmaxloc-xminloc) / (0.9998 * dism) ) + 1
      ny = int( (ymaxloc-yminloc) / (0.9998 * dism) ) + 1

!     --- Number of quadrilaterals

      ncube = nx * ny
   
      if ( debug ) 
     +     write(irefwr,*) 'nx, ny, ncube, dism',nx,ny,ncube,dism

!     --- Allocate space for arrays jcube and cube:

      allocate( jcube( ncube ) )
      allocate(  cube( ncube ) )

!     --- Determine for each centre a dis value

      xstart = ( xmaxloc + xminloc - (nx)*dism ) / 2d0
      ystart = ( ymaxloc + yminloc - (ny)*dism ) / 2d0

!     --- Place startvalues a little "upwards":

      xstart = xstart + dism/12d0
      ystart = ystart + dism/12d0

      if ( debug ) then
         write(irefwr,*) 'MSHO2D:'
         write(irefwr,*) 'xminloc =',xminloc,' xmaxloc =',xmaxloc
         write(irefwr,*) 'yminloc =',yminloc,' ymaxloc =',ymaxloc
         write(irefwr,*) 'x en y start in cube 1',xstart,ystart
         write(irefwr,*) 'x last cube is ',xstart+nx*dism
         write(irefwr,*) 'y last cube is ',ystart+ny*dism
         write(irefwr,*) 'dismin = ',dismin,' dismax = ',dismax
      end if

!     --- Determine coarseness for each quadrilateral

      if ( debug )
     +   write(irefwr,*) 'Determine coarseness for surf ',isurnr

      call msho06 ( npoint, coor, dism, xstart, ystart, nx, ny,
     +              icube, chelp, cube, jcube, kbound, nbound,
     +              coar, ncoar, ncurvs, curves, cocurv )

!     --- Check all mutual coarsenesses: userpoint, internal_points and
!         internal_curves

      if ( jcoars==1 ) then
         maxratio = coaval( nuspnt + 2 )
      else
         maxratio = 0
      end if

!     --- Check coarseness ratio's only in 2D problems and only if user
!         has submitted a value (i.e. maxratio is not zero!)

      if ( ndim==2 .and. maxratio>1d0 ) then

!     --- Check submitted coarsenesses of user:

         call msho38 ( npoint, coor  , dism   , xstart, ystart,
     +                 kbndpt, nbndpt, numcrvboun     , boundary,
     +                 nbound, nholes, nx     , ny    , icube ,
     +                 chelp , cube  , jcube  , kstapl, kstap ,
     +                 ncurvs, curves, cocurv , isurnr, userco,
     +                 nuspnt, coaval, tran )
      end if

      if ( ierror/=0 ) goto 1000

!     --- Compute reference area value for each quadrilateral

      deallocate( chelp )

!     --- Next allocate and use array chelp for area of blocks: 

      allocate( chelp  ( ncube ) )

      call msho21 ( cube, ncube, chelp )

      if ( ierror/=0 ) goto 1000

!     --- Set number of old (fixed!) points
!         (Now nodes on internal lines are not excluded)

      nipnt = npoint

!     --- Repositioning

      nrepos = 10 * npunt + 20

!     --- Clear array itri

      do i = 1, npunt
         itri(i) = 0
      end do

!     --- Fill array itri

      do i = 1, 2*kstap
         ipoint       = kstapl(i)
         itri(ipoint) = itri(ipoint) + 1
      end do

!     --- Eventually check array kstapl:

      nelem  = 0
      nelemi = 0

!     --- First step: check whether there are special points that should be
!         treated first. If so fill elements.

      coarmin = coarsemin

      if ( ncoar>0 ) then

!     --- Place special nodes and fill elements:

         call msho35 ( npoint, coor, xstart, ystart, dism,
     +                 coar, ncoar, icube, nx, kmeshc, nelem,
     +                 kstapl, kstap, itri, isurnr, userpoints,
     +                 kbndpt, nbndpt )

!        --- Adjust number of initial ("fixed") points
!            (Make choice: only the point itself, or also
!             the surrounding points)
!            --- In case of also the surrounding points
!            nipnt  = nipnt + 7 * ncoar
!            --- Choice: only the point itself:

         nipnt  = nipnt + 7 * ncoar
         nelemi = nelem

         do i=1, ncoar
            if ( coar(3,i)<coarmin ) coarmin = coar(3,i)
         end do

      end if

      if ( ierror/=0 ) goto 1000

!     --- Check whether all elements are given right and adjust array kstapl

      call msho08 ( kstapl, kstap, coor, xstart, dismin, holeinfo,
     +              nholes, .false. )

      if ( ierror/=0 ) goto 1000

!     --- Part for division
!         nochan : parameter for counting of changings
!         nherha : parameter for redivisions

      nochan = 0
      nherha = 0
      iperm  = 0

!     --- New value of ichan, used to be ichan = kstap-25
!         Makes only sense for small meshes, i.e kstap small
!         By trial and error we found this to be a working value
!         In case of a mesh with a very large ratio of the elements
!         it was necessary to use a higher minimum value for ichan than in the
!         case of smooth subdivisions of the boundary

      if ( debug ) then
         write(irefwr,*) 'Nodes of boundary as submitted by user'
         do i=1, npoint
            write(irefwr,*) 'Nodenmbr',i,'Crds',coor(1,i),coor(2,i)
         end do
         write(irefwr,*) 'end of given nodes'
      end if

300   ichan = min(imultiply*25, kstap/2+1)

      if ( npoint>npunt ) then
         if ( debug ) write(irefwr,*) 'Error while kstap = ',kstap
         call errint ( npunt , 1 )
         call errint ( npoint, 2 )
         call errsub ( 900, 2, 0, 0 )
      end if

      if ( nelem>nkmeshc ) then
        write(irefwr,*) 'Estimated length of element array too small'
        write(irefwr,*) 'Number of elements found is more than the'
        write(irefwr,*) 'used estimate',nkmeshc,'in SEPRAN.'
        write(irefwr,*) 'Most probable cause: too much difference'
        write(irefwr,*) 'in the coarseness along the boundary.'
        write(irefwr,*) 'Try to divide your original region in'
        write(irefwr,*) 'smaller sections with less coarseness variety'
        write(irefwr,*) 'or check and change the local coarsenesses in'
        write(irefwr,*) 'the userpoints of the original region.'
        write(irefwr,*) 'If nothing helps: Please inform SEPRA'
         call errsub ( 901, 0, 0, 0 )
      end if

      if ( ierror/=0 ) goto 1000

!     --- Here follow helpstatements for debugging, normally these statements
!         are skipped, using debug = .false.:

      debug = .false.

!     --- debug = .false.

      if ( debug .and. nochan==0 .and. nelem>0 ) then
           i1 = kmeshc(1,nelem)
           i2 = kmeshc(2,nelem)
           i3 = kmeshc(3,nelem)
           call msho19( coor, i1, i2, i3, surf1 )
           write(irefwr,*) 'element',nelem,' - ',i1,i2,i3
           if ( 0>1 ) then
              do i=1, kstap
                 i1 = kstapl(2*i-1)
                 i2 = kstapl(2*i)
                 write(irefwr,*) 'Stuk',i1,i2
              end do
              write(irefwr,*) 'coordinates'
              do i=1,kstap
                 i1 = kstapl(2*i-1)
                 xm = coor(1,i1)
                 ym = coor(2,i1)
                 write(irefwr,*) i1,xm,ym
              end do
              write(irefwr,*) 'End of boundary'
           end if
      end if

      if ( nochan==0 .and. nelem==-3 ) then
           write(irefwr,*) 'elem',nelem,'nodes',(kmeshc(i,nelem),i=1,3)
           i1 = kmeshc(1,nelem)
           i2 = kmeshc(2,nelem)
           i3 = kmeshc(3,nelem)
           call msho19( coor, i1, i2, i3, surf1 )
           write(irefwr,*) 'surface of element = ',surf1
           if ( surf1<10*epsmac ) then
              write(irefwr,*) 'Results, kstap =',kstap
              do i=1, npoint
                 xm = coor(1,i)
                 ym = coor(2,i)
                 write(irefwr,'(i6,(2f10.5))' ) i,xm,ym
              end do
              write(irefwr,*) 'Elements'
              do i=1,nelem
                 i1 = kmeshc(1,i)
                 i2 = kmeshc(2,i)
                 i3 = kmeshc(3,i)
                 write(irefwr,'(4i5)') i,i1,i2,i3
              end do
              write(irefwr,*) 'the end '
              go to 950
           end if
      end if
      debug = .false.

      iperm = 0

!     --- Statements to indicate at which stage the subdivision is

      if ( debug .and. (nelem-1000*(nelem/1000))==0.and.nochan==0 ) then
         write(irefwr,*) 'nelem = ',nelem
         write(irefwr,*) 'kstap = ',kstap
      end if

      if ( debug .and. nelem>100 ) then

!     --- Check kstapl

         write(irefwr,*) 'Check kstapl, nelem =',nelem,'kstap',kstap
         call mshcurvinters2 ( kstapl, kstap, coor, isurnr )
         if ( ierror/=0 ) then
            write(irefwr,*)
     +      'Intersection, error! nelem = ',nelem
            goto 1000
         end if
      end if

      if ( debug .and. (nelem-1000*(nelem/1000))==0.and.nochan==0) then
         write(irefwr,*) 'Control: nelem = ',nelem

!        --- Check whether all elements are given right and adjust array kstapl

         write(irefwr,*) 'Check via msho08 '
         call msho08 ( kstapl, kstap, coor, xstart, dismin, holeinfo,
     +                 nholes, .true. )
         if ( ierror/=0 ) goto 1000
         write(irefwr,*) 'No error during check for nelem = ',nelem
         call msho40 ( coor, kmeshc, npoint, nelem, iperm )
      end if

      if ( debug .and. nochan>0 ) then
         write(irefwr,*) 'nelem = ',nelem,' nochan = ',nochan
         write(irefwr,*) 'no element for',i1,'- ',i2
      end if
      if ( debug .and. isurnr==-15 ) then
         write(irefwr,*) 'surface number is',isurnr
         write(irefwr,*) 'nelem = ',nelem,'kstap = ',kstap
         if ( nochan==0 .and. nelem>0 ) then
            write(irefwr,*) (kmeshc(i,nelem),i=1,3)
         end if
      end if

      if ( debug .and. nochan==0 .and. nelem>0 ) then

!     --- Print all elements one after one

         write(irefwr,*) 'Element',nelem,(kmeshc(i,nelem),i=1,3)
      end if

      debug = .false.

!     --- Determine value for first group of boundary pieces to be considered:

      if ( nochan<kstap-1 ) then
         ichan = min( 25, (kstap - nochan - 1) )
      else
         ichan = 0
      end if

      icount = icount+1

      if ( debug .and. icount>=idebug ) then
         write(irefwr,*) 'extra check, kstap, npoint, icount ',
     +                    kstap, npoint , icount

         call msho08 ( kstapl, kstap, coor, xstart, dismin,
     +                 holeinfo, nholes, .true. )

!        --- The next statements can be used if the process "hangs"

         write(irefwr,*) 'kstap, nochan,nelem ', kstap, nochan,nelem
      end if

      if ( debug .and. nelem==-20000 ) then
         write(irefwr,*) 'nelem = ',nelem,' checks'
         call mshchkstapl ( kstapl, kstap, 'geval nelem' )
         write(irefwr,*) 'na mshchkstapl ierror =',ierror
         call mshcurvinters2 ( kstapl, kstap, coor, isurnr )
         write(irefwr,*) 'na mshcurvinters2, ierror =',ierror
         write(irefwr,*) 'nog ok?'
      end if

      if ( debug ) then

!     --- Debug information

      end if

      debug = .false.

!     --- Check for smallest element:

      if ( ichan>0 ) call msho25 ( kstapl, ichan, coor )

      angle  =  -0.1d0

!     --- Remark: factor determines the number of elements
!         Decreasing factor to for example 0.8 gives a much smoother mesh,
!         but also a lot of extra points and elements

      factor =   1d0

!     --- Change values in case of no convergence

      if ( nochan>10 .or. nochan>kstap ) then
         angle  =  -0.7d0
         factor =   0.6d0
      end if

      if ( nochan >2 * kstap ) then
         angle  = -0.95d0
         factor =  0.50d0
      end if

!     --- If no convergence possible stop

      if ( nochan>3*kstap ) then
         !write(irefwr,*) 'nochan =',nochan,' and kstap =',kstap

!        --- nochan too large

         if ( debug .or. 1>0 ) then
            write(irefwr,*) 'Error 902, see information'

            do i=1, kstap
               write(irefwr,*) kstapl(2*i-1),kstapl(2*i)
            end do
            do i=1, kstap
               ja = kstapl(2*i-1)
               write(*,*) 'P',ja,'=',coor(1,ja),coor(2,ja)
            end do
            write(irefwr,*) 'npoint = ',npoint
            write(irefwr,*) 'Enough info?'
         end if

         call errint ( nochan, 1 )
         call errint ( 3*kstap, 2 )
         call errint ( kstap, 3 )
         call errsub ( 902, 3, 0, 0 )

      end if
      if ( ierror/=0 ) goto 1000

!     --- Check for ready

      if ( kstap==0 ) goto 500

!     --- Transform in triangles, run through kstapl

      i1 = kstapl(1)
      i2 = kstapl(2)

!     --- Change array itri

      itri(i1) = itri(i1) - 1
      itri(i2) = itri(i2) - 1

!     --- First check whether triangle can be formed with other
!         boundary lines, compute nodes and angles

      call msho12 ( coor, kstapl, kstap, i1, i2, iex1, iex2,
     +              angle1, angle2 )

!     --- In emergency case:

      if ( nochan>kstap .and. 0>1 ) then

!     --- Try shortcut, compute surfaces:

         call msho19( coor, iex1, i1, i2, surf1 )
         call msho19( coor, i1, i2, iex2, surf2 )

         write(irefwr,*) 'Shortcut, Emergency used:'
         if ( surf1>0 .and. surf1>surf2 ) then
            call msho10( kmeshc, nelem, i1, i2, iex1,
     +                   kstapl, kstap, itri )
            write(irefwr,*) i1,i2,iex1
            nochan = 0

            goto 300

         else if ( surf2>0 .and. surf2>surf1 ) then

            call msho10( kmeshc, nelem, i1, i2, iex2,
     +                   kstapl, kstap, itri )
            write(irefwr,*) i1,i2,iex2
            nochan = 0

            goto 300

         end if
      end if

!     --- Check for errors

      if ( ierror/=0 ) goto 1000

      ja = 0

      if ( nochan>kstap .and. kstap>4 ) then

!     --- Combination of difficult triangles ?

         call msho12( coor, kstapl, kstap, iex1, i1, j1, j2, e1, e2 )
         if ( j1==iex2 ) ja = 1
      end if
      if ( ierror/=0 ) goto 1000

      if ( kstap==4 .or. ja==1 ) then

!     --- Special situation

         call msho19( coor, iex1, i1, i2, surf1 )
         call msho19( coor, i2, iex2, iex1, surf2 )

         if ( surf1<0 .or. surf2<0 ) then

!        --- Difficult situation, first next triangle

            call msho19( coor, i1, i2, iex2, surf1 )
            if ( surf1>0 ) then
               call msho10( kmeshc, nelem, i1, i2, iex2, kstapl,
     +                      kstap, itri )
               nochan = 0
               goto 300
            end if
         end if
      end if
      if ( ierror/=0 ) goto 1000

      if ( iex1==iex2 ) then

!     --- Consider triangle i1 - i2 - iex1 = iex2, compute area

         call msho19( coor, i1, i2, iex1, surf )

         iperm = 0
         if ( inside==1 )
     +        call msho28( coor, i1, i2, iex1, xn, yn, npoint,
     +                     itri, iperm )

         if ( surf>=0 .and. iperm==0 ) then

            call msho10( kmeshc, nelem, i1, i2, iex1,
     +                   kstapl, kstap, itri )
            nochan = 0

            goto 300

         end if

         if ( surf<=0 .and. kstap==3 ) then

!        --- Error: may be that later on repositioning solves this problem:

            write(irefwr,*) 'Emergency: Temporary bad element created'
            write(irefwr,*) 'element has nodes'
            write(irefwr,*) i1,' - ',i2,' - ',iex1
            call msho10( kmeshc, nelem, i1, i2, iex1,
     +                   kstapl, kstap, itri )
            nochan = 0

            goto 300

         end if

      end if

!     --- Check for angles not greater than 180 degrees and determine
!         which points are potential for new triangles

      call msho19( coor, i1, i2, iex1, surf1 )
      if ( ierror/=0 ) goto 1000

      if ( surf1<eps ) then
         angle1 = -1.0d0
         ipotn1 = 1
      else
         ipotn1 = 0
      end if

      call msho19( coor, i1, i2, iex2, surf2 )
      if ( ierror/=0 ) goto 1000

      if ( surf2<eps ) then
         angle2 = -1.0d0
         ipotn2 = 1
      else
         ipotn2 = 0
      end if

      if ( angle1>angle .and. angle2>=angle1 ) then

         call msho11 ( i1, iex1, iex2, coor, angleh )

         if ( angleh<-0.5 ) then

!        --- Try i1 - i2 - iex1

            ja = 1

!           --- Adjust computation (surface is always computed)

            call msho19 ( coor, i2, iex2, iex1, surf )
            if ( surf<100 * eps ) ja = 0

!           --- Check if none of the boundary nodes is in the new triangle:

            iperm = 0
            if ( ja==1 )
     +           call msho28 ( coor, i1, i2, iex1, xn, yn, npoint,
     +                      itri, iperm )

            if ( ja==1 .and. iperm==0 ) then

!           --- Check whether i2 - iex1 is permitted

               call msho24 ( kstapl, kstap, coor, i2, iex1, icheck )

               if ( icheck==0 ) then

                  call msho10( kmeshc, nelem, i1, i2, iex1, kstapl,
     +                         kstap, itri )
                  nochan = 0

                  goto 300

               end if

            else

               surf1 = rinfin

            end if

         else   ! (angleh>angle )

!        --- Try i1 - i2 - iex2

            ja = 1

!           --- Adjust computation (surface is always computed)

            call msho19 ( coor, i1, i2, iex2, surf )
            if ( surf<100 * eps ) ja = 0

!           --- Check if none of the boundary nodes in the new triangle:

            iperm = 0
            if ( ja==1 )
     +           call msho28 ( coor, i1, i2, iex2, xn, yn, npoint,
     +                      itri, iperm )

            if ( ja==1 .and. iperm==0 ) then

!           --- Check whether i1 - iex2 is permitted

               call msho24 ( kstapl, kstap, coor, i1, iex2, icheck )

               if ( icheck==0 ) then

                  call msho10( kmeshc, nelem, i1, i2, iex2, kstapl,
     +                         kstap, itri )
                  nochan = 0

                  goto 300

               end if

            else

               surf2 = rinfin

            end if

         end if

      end if

      if ( ierror/=0 ) goto 1000

      if ( surf1<eps ) surf1 = rinfin

      if ( angle2>angle .and. angle1>=angle2 ) then

         call msho11 ( i2, iex2, iex1, coor, angleh )

         if ( angleh<-0.5 ) then

!        --- Try i1 - i2 - iex2

            ja = 1

!           --- Adjust computation (surface is always computed)

            call msho19 ( coor, i1, iex2, iex1, surf )
            if ( surf<100 * eps ) ja = 0

!           --- Check if none of the boundary nodes is in the new triangle:

            iperm = 0
            if ( ja==1 )
     +           call msho28 ( coor, i1, i2, iex2, xn, yn, npoint,
     +                      itri, iperm )

            if ( ja==1 .and. iperm==0 ) then

!           --- Check whether i1 - iex2 is permitted

               call msho24 ( kstapl, kstap, coor, i1, iex2, icheck )

               if ( icheck==0 ) then
                  call msho10( kmeshc, nelem, i1, i2, iex2, kstapl,
     +                         kstap, itri )
                  nochan = 0
                  goto 300
               end if
            else
               surf2 = rinfin
            end if
         else   ! (angleh>angle )

!        --- Try i1 - i2 - iex1

            ja = 1

!           --- Adjust computation (surface is always computed)

            call msho19 ( coor, i1, i2, iex1, surf )
            if ( surf<100 * eps ) ja = 0

!           --- Check if none of the boundary nodes is in the new triangle:

            iperm = 0
            if ( ja==1 )
     +           call msho28 ( coor, i1, i2, iex1, xn, yn, npoint,
     +                      itri, iperm )
            if ( ja==1 .and. iperm==0 ) then

!           --- Check whether i2 - iex1 is permitted

               call msho24 ( kstapl, kstap, coor, i2, iex1, icheck )

               if ( icheck==0 ) then
                  call msho10( kmeshc, nelem, i1, i2, iex1, kstapl,
     +                         kstap, itri )
                  nochan = 0
                  goto 300
               end if
            else
               surf1 = rinfin
            end if
         end if
      end if

      if ( ierror/=0 ) goto 1000

      if ( surf2<eps ) surf2 = rinfin

!     --- No special reason for triangle with boundary lines ;
!         Determine normal and position of preferable point (xm,ym) for
!         line i1 - i2
!         Check coarseness of "old" line

      dx = coor( 1,i2 ) - coor( 1,i1 )
      dy = coor( 2,i2 ) - coor( 2,i1 )

      disold = sqrt( dx*dx + dy*dy )

!     --- Required coarseness in this cube:

      coa    = ( cube(icube(i1)) + cube(icube(i2)) ) / 2d0
      disold = ( disold + 2 * coa ) / 3d0

!     --- Compute middle point of line:

      call msho09( coor, i1, i2, e1, e2, xm, ym )

!     --- First estimate of point

      xn = xm + 0.95d0 * disold * e1
      yn = ym + 0.95d0 * disold * e2

!     --- Determine helpvalue for cube estimate:

      n1 = int ( (xn-xstart) / dism )
      if ( n1 < 0 ) then
           n1 = 0 
      else if ( n1 > nx ) then
           n1 = nx
      end if
      n2 = int ( (yn-ystart) / dism )
      if ( n2 < 0 ) then
           n2 = 0
      else if ( n2 > ny - 1 ) then
           n2 = ny - 1
      end if 
      nc = 1 + n1 + n2*nx

!     --- Find helpcoarseness

      cdis = cube(nc)

!     --- New guess using both coarsenesses:
      coa = (disold + 2 * cdis)/3d0

      xn = xm + coa * e1
      yn = ym + coa * e2

!     --- Check whether new point belongs to volume

      ja = 0

      call msho07( xn, yn, xstart, coor, kstapl, kstap, ja )

      if ( ierror/=0 ) goto 1000

!     --- Determine extra values for case that new point has to be placed

      if ( cdis > 1.2 * disold ) then
         xnx = xm + coa * e1
         yny = ym + coa * e2
      else if ( cdis < 0.8 * disold ) then
         xnx = xm + 0.7 * coa * e1
         yny = ym + 0.7 * coa * e2
      else
         xnx = xm + 0.85 * coa * e1
         yny = ym + 0.85 * coa * e2
      end if

!     --- Check whether new point is allowed or is there a common
!         point with another line

      call msho13( coor, i1, i2, kdrie, kstapl, kstap, xn, yn )

!     --- Check for double crossinglines:

      if ( kdrie==-1 ) then

!     --- Double crossings: change array kstapl

         if ( debug ) then
            write(irefwr,*) 'Check msho13: Double crossing: skip'
         end if
         call msho20(kstapl,kstap,itri)
         nochan = nochan + 1
         goto 300
      end if

      if ( ierror/=0 ) goto 1000

!     --- Check for impossible situation

      if ( ja==0 .and. kdrie==0 ) then

!     --- Rearrange, impossible

         xn = xn - dismin * 1d-3 + dismin * 1d-4
         yn = yn + dismin * 1d-3 + dismin * 1d-4

         ja = 0
         call msho07( xn, yn, xstart, coor, kstapl, kstap, ja )

         call msho13( coor, i1, i2, kdrie, kstapl, kstap, xn, yn )

         if ( ja==0 .and. kdrie==0 .or. kdrie==-1 ) then

!        --- Really impossible, change

            if ( debug ) then
               write(irefwr,*) 'Check msho13'
               write(irefwr,*) 'Really impossible situation, skip'
               write(irefwr,*) 'If your input is OK inform SEPRA'
            end if
            call msho20( kstapl, kstap, itri )
            nochan = nochan + 1
            goto 300
         end if
      end if
      if ( ierror/=0 ) goto 1000

420   if ( kdrie>0 ) then

!     --- Common point with line kdrie

         ii1 = kstapl( 2*kdrie - 1 )
         ii2 = kstapl( 2*kdrie     )

!        --- Triangle with one of these (old) points
!            ( Check whether points are allowed )

         jpn = 0

         if ( ii1==i2 .and. surf2<rinfin/2. ) then

            call msho19( coor, i1, iex1, ii2, surf )
            if ( surf<-eps .or. ipotn1==1 ) then
               call msho24(kstapl,kstap,coor,i1,ii2,icheck)
               call msho19(coor, i1, i2, ii2, surf)

               if ( icheck==0 .and. surf>eps  ) then
                  call msho10(kmeshc,nelem,i1,i2,ii2,kstapl,
     +                        kstap,itri)
                  nochan = 0
                  goto 300
               else if ( icheck>0 ) then
                  kdrie = icheck
                  goto 420
               end if

            end if

         end if

         if ( ii2==i1 .and. surf1<rinfin/2. ) then

            call msho19(coor,i2,ii1,iex2,surf)
            if ( surf<-eps .or. ipotn2==1 ) then
               call msho24(kstapl,kstap,coor,i2,ii1,icheck)
               call msho19(coor, i1, i2, ii1, surf)

               if ( icheck==0 .and. surf>eps ) then
                  call msho10(kmeshc,nelem,i1,i2,ii1,kstapl,
     +                        kstap,itri)
                  nochan = 0
                  goto 300
               else if ( icheck>0 ) then
                  kdrie = icheck
                  goto 420
               end if

            end if

         end if

         if ( ierror/=0 ) then
            if ( debug ) write(irefwr,*) 'Halfway ierror =',ierror
            goto 1000
         end if

         if ( i2/=ii1 .and. i1/=ii2 ) then

!        --- Point ii1 ?

            dist1 = rinfin

            dx = xnx - coor( 1,ii1 )
            dy = yny - coor( 2,ii1 )

            dist1 = sqrt( dx*dx + dy*dy )
            jpn   = ii1

!           --- Point ii2 ?

            dx = xnx - coor( 1,ii2 )
            dy = yny - coor( 2,ii2 )

            dist  = sqrt( dx*dx + dy*dy )

            if ( dist<dist1 ) then
               dist1 = dist
               jpn   = ii2
            end if

!           --- Point jpn allowed, else place line at bottom of kstapl

            no = 0

            call msho24(kstapl,kstap,coor,i1,jpn,icheck)
            if ( icheck/=0 ) no = 1
            call msho24(kstapl,kstap,coor,i2,jpn,icheck)
            if ( icheck/=0 ) no = 1

            if ( no==1 ) jpn = 0

            if ( jpn/=0 .and. jpn/=iex1 .and.
     +           ipotn1==0 .and. nochan<kstap ) then

               call msho19(coor,iex1,i1,jpn,surf)

               valare = chelp(icube(i1))

               if ( surf<0.15d0 * valare ) then
                  call msho20(kstapl,kstap,itri)
                  nochan = nochan + 1
                  goto 300
               end if

            end if

            if ( jpn/=0 .and. jpn/=iex2 .and.
     +           ipotn2==0 .and. nochan<kstap ) then

               call msho19(coor,i2,iex2,jpn,surf)

               valare = chelp(icube(i2))

               if ( surf<0.15d0 * valare ) then
                  call msho20(kstapl,kstap,itri)
                  nochan = nochan + 1
                  goto 300
               end if

            end if

         end if

         if ( jpn>0 ) then
            call msho19(coor,i1,i2,jpn,surf)
            valare = ( chelp(icube(i1))+chelp(icube(i2)) )/2d0
            if ( nochan<kstap .and.
     +           surf<0.01d0 * valare .or.
     +           surf<10 * eps ) jpn = 0
            if ( jpn==iex1.and.surf1>rinfin/2.) jpn = 0
            if ( jpn==iex2.and.surf2>rinfin/2.) jpn = 0
         end if

         iperm = 0
         if ( jpn>0 )
     +      call msho28(coor,i1,i2,jpn,xn,yn,npoint,itri,iperm)
         if ( iperm==1 ) jpn = 0

         if ( jpn>0 ) then
            call msho10(kmeshc,nelem,i1,i2,jpn,kstapl,
     +                  kstap,itri)
            nochan = 0

         else

!        --- Change array kstapl

            call msho20(kstapl,kstap,itri)
            nochan = nochan + 1
         end if
      else

!     --- New point not to near to old point ? Run through array itri

         jpn   = 0
         dista = 2. * coa

         call msho14(coor,jpn,npoint,itri,i1,i2,xnx,yny,dista)

         if ( dista>0.55d0 * coa ) jpn = 0

         iperm = 0
         if ( jpn/=0 .and. inside==1 )
     +        call msho28(coor,i1,i2,jpn,xn,yn,npoint,itri,iperm)
         if ( iperm==1 ) then
              call msho20(kstapl,kstap,itri)
              nochan = nochan + 1
              goto 300
         end if

         if ( jpn/=0 .and. jpn/=iex1 .and.
     +        ipotn1==0 ) then

            call msho19(coor,iex1,i1,jpn,surf)

            valare = chelp(icube(i1))

            if ( surf<0.15d0 * valare ) then
               call msho20(kstapl,kstap,itri)
               nochan = nochan + 1
               goto 300
            end if

         end if

         if ( jpn/=0 .and. jpn/=iex2 .and.
     +        ipotn2==0 ) then

            call msho19(coor,i2,iex2,jpn,surf)

            valare = chelp(icube(i2))

            if ( surf<0.15d0 * valare ) then
               call msho20(kstapl,kstap,itri)
               nochan = nochan + 1
               goto 300
            end if

         end if

         inear = 0
         no    = 0

         if ( dista<0.55d0*coa ) then
            inear = 1
            call msho24(kstapl,kstap,coor,i1,jpn,icheck)
            if ( icheck/=0 ) no = 1
            call msho24(kstapl,kstap,coor,i2,jpn,icheck)
            if ( icheck/=0 ) no = 1
         end if

         if ( jpn==0 ) then
            ja = 0
         else
            call msho19(coor,i1,i2,jpn,surf)
            valare = ( chelp(icube(i1))+chelp(icube(i2)) )/2.
            if ( nochan<kstap .and.
     +           surf<0.02d0 * valare .or.
     +           surf<10 * eps ) then
               ja = 0
            else
               ja = 1
            end if
         end if

         if ( no==1 ) ja = 0

!        --- If ja = 0 then not allowed

         if ( ja==0 ) dista = rinfin

         if ( dista<0.55d0*coa ) then

!        --- New triangle is i1,i2 with old point jpn

            call msho10(kmeshc,nelem,i1,i2,jpn,kstapl,kstap,
     +                  itri)
            nochan = 0

         else if ( inear==0 ) then

!        --- New triangle is i1,i2,npoint+1 ; fill arrays coor,
!            kmeshc and adjust kstapl and itri

            call msho15(coor,npoint,kmeshc,nelem,kstapl,kstap,
     +                  itri,i1,i2,xnx,yny)
            nochan = 0

!           --- Determine cube for this point

            n1 = int ( (xnx-xstart) / dism )
            n2 = int ( (yny-ystart) / dism )

!     --- Determine helpvalue for cube estimate:

            n1 = int ( (xnx-xstart) / dism )
            if ( n1 < 0 ) then
                 n1 = 0 
            else if ( n1 > nx ) then
                 n1 = nx
            end if
            n2 = int ( (yny-ystart) / dism )
            if ( n2 < 0 ) then
                 n2 = 0
            else if ( n2 > ny - 1 ) then
                 n2 = ny - 1
            end if 

            nc = 1 + n1 + n2*nx

            if ( nc>(nx*ny) ) then
               write(irefwr,*) 'Error in msho2d'
               write(irefwr,*) 'Please inform SEPRA'
               write(irefwr,*) 'Execution NOT terminated' !AvD: do not stop
!AvD               call instop
            end if

            icube(npoint) = nc

         else

!        --- No solution found: change array kstapl

            call msho20(kstapl,kstap,itri)
            nochan = nochan + 1
         end if

      end if

      if ( ierror/=0 ) then
         if ( debug ) write(irefwr,*) 'Error: End of normal loop'
         goto 1000
      end if

!     --- Array kstapl empty ?

      if ( kstap>0 ) goto 300

!     --- Check itri (all entries of this array should by now be zero)

500   itot = 0

      if ( debug ) then
         write(irefwr,*) 'Temporary output after label 500'
         write(irefwr,*) 'nelem = ',nelem,'kstap = ',kstap
         if ( nelem>0 ) then
            write(irefwr,*) 'Listing of',nelem,'elements'
            do i=1, nelem
               write(irefwr,*) i,(kmeshc(i1,i),i1=1,3)
            end do
            write(irefwr,*) 'Output meshgeneration nodes'
            write(irefwr,*) 'End msho2d with ',npoint,' nodes'
            do i=1, npoint
               write(irefwr,*) i,(coor(ja,i),ja=1,2)
            end do

            write(irefwr,*) 'Array kbndpt: ( nbndpt = ',nbndpt,')'
            write(irefwr,*) (kbndpt(i),i=1,nbndpt)
            write(irefwr,*) 'Contents OK?'
         end if
      end if

      do i = 1, npoint
         itot = itot + abs(itri(i))
      end do

      if ( itot/=0 ) then
         write(irefwr,*) 'Error in connections array'
         call errsub ( 903, 0, 0, 0 )
      end if

      if ( npoint>npunt ) then
         call errint ( npoint, 1 )
         call errint ( npunt, 2 )
         call errsub ( 900, 2, 0, 0 )
      end if

      if ( nelem>nkmeshc ) then
         call errsub ( 901, 0, 0, 0 )
      end if
      if ( ierror/=0 ) goto 1000

!     --- First solution found, save before extra actions:

      if ( nherha==0 ) then

!     --- Allocate space for helparrays coortmp and meshtmp:

!         --- Length coortmp is npoint:

         npndef = npoint + int(5 + npoint/10)
         npntmp = npoint

!        --- Length meshtmp is nelem:

         neldef = nelem + int(5 + nelem/10)
         neltmp = nelem

         allocate( coortmp( 2, npndef ) )
         allocate( meshtmp( 3, neldef ) )

!        --- Fill helparrays temporary:
!            --- Coordinates:

         do i = 1, npoint
             coortmp(1,i) = coor(1,i)
             coortmp(2,i) = coor(2,i)
         end do

!        --- Elements:
         do i = 1, nelem 
             meshtmp(1,i) = kmeshc(1,i)
             meshtmp(2,i) = kmeshc(2,i)
             meshtmp(3,i) = kmeshc(3,i)
         end do

!        --- Determine characteristic ratiovalue for comparison reasons in
!            the next steps:

         ratiop = rinfin

         do ielem = 1, nelem
            i1 = kmeshc (1,ielem)
            i2 = kmeshc (2,ielem)
            i3 = kmeshc (3,ielem)
            call msho33 ( coor, verh, i1, i2, i3 )
            call msho19 ( coor, i1, i2, i3, surf )
            if ( surf<0 ) then
               verh = -verh
               if ( debug ) then
                  write(irefwr,*) 'nherha = ',nherha
                  write(irefwr,*) 'verh   = ',verh
                  write(irefwr,*) 'Triangle',ielem,' - ',i1,i2,i3
               end if
            end if
            ratiop = min( verh, ratiop )

         end do  ! ielem = 1, nelem

         if ( debug ) then
            write(irefwr,*) 'End first subdivision, ratio = ',ratiop
         end if
      end if

!     --- Determine number of neighbours of points

540   leng = nrepos

      do ih = 1, 3

         call msho16(kmeshc,nelem,npoint,nipnt,itri,kstapl,leng)

!        --- Check length used

         if ( leng>nrepos ) then
            call errsub ( 903, 0, 0, 0 )
            goto 1000
         end if

!        --- Remove internal points with three or four neighbours

         icancel = 1

         do while ( icancel>0 )

            icancel = 0

            call msho29( kmeshc, nelem, npoint, itri, kstapl, nipnt,
     +                   coor, icancel )

            if ( ierror/=0 ) then
               if ( debug ) write(irefwr,*) 'msho29 ierror =',ierror
               goto 1000
            end if

!           --- Rearrange arrays COOR and KELEM

            do i = 1, npoint
               itri(i) = 0
            end do

            do i = 1, nelem
               i1       = kmeshc( 1,i )
               i2       = kmeshc( 2,i )
               i3       = kmeshc( 3,i )

               itri(i1) = 1
               itri(i2) = 1
               itri(i3) = 1
            end do

            itot = 0

            do i = 1, npoint
               x = coor( 1,i )
               y = coor( 2,i )

               if ( itri(i)==1 .or. i<=nipnt ) then
                  itot = itot + 1
                  coor( 1,itot ) = x
                  coor( 2,itot ) = y
                  itri(i)        = itot
               end if

            end do

            ja     = npoint - itot
            npoint = itot

!           --- Renumber KELEM

            do i = 1, nelem

               i1 = kmeshc( 1,i )
               i2 = kmeshc( 2,i )
               i3 = kmeshc( 3,i )

               kmeshc( 1,i ) = itri(i1)
               kmeshc( 2,i ) = itri(i2)
               kmeshc( 3,i ) = itri(i3)
            end do

         end do ! while loop

      end do

!     --- Temporary for checking:
!          if ( kstap==0 ) goto 900

!         --- Check angles

      ifill = 0

      do ielem = nelemi+1, nelem

         i1 = kmeshc( 1,ielem )
         i2 = kmeshc( 2,ielem )
         i3 = kmeshc( 3,ielem )

!        --- Check all inside "diagonals":
!            Start computing the angles of each triangle

            call msho03(i1,i2,coor,adis)
            call msho03(i2,i3,coor,bdis)
            call msho03(i3,i1,coor,cdis)

!           --- Check angles

            call msho22( adis, bdis, cdis, icase )

            if ( icase>0 ) then

!           --- Triangle ielem has an angle that is larger than 90 degrees,
!               try to change the diagonal

               ma = 1

               jelem = ielem

               if ( icase==1 ) then

!              --- common line is i2 - i3

                  call msho26(kmeshc,nelem,i2,i3,jelem,i4)

!                 --- Check whether i2-i3 is internal boundary:

                  ja = 0
                  call msho42( kinbnd,lenbnd, i2, i3, ja )
                  if ( ja==1 ) i4 = 0
                  if ( i4>0 .and. jelem/=ielem ) then

!                 --- consider i1 - i2 - i3 - i4

                     j1 = i1
                     call msho27(coor,i1,i2,i3,i4)
                     if ( j1/=i1 ) ifill = 1

                     j1 = i1
                     j2 = i2
                     j3 = i3
                     j4 = i4
                  else
                     ma = 0
                  end if
               else if ( icase==2 ) then

!              --- common line is i3 - i1

                  call msho26(kmeshc,nelem,i3,i1,jelem,i4)

!                 --- Check whether i3-i1 is internal boundary:

                  ja = 0
                  call msho42( kinbnd, lenbnd, i3, i1, ja )
                  if ( ja==1 ) i4 = 0
                  if ( i4>0 .and. jelem/=ielem ) then

!                 --- consider i2 - i3 - i1 - i4

                     j2 = i2
                     call msho27(coor,i2,i3,i1,i4)
                     if ( j2/=i2 ) ifill = 1

                     j1 = i2
                     j2 = i3
                     j3 = i1
                     j4 = i4
                  else
                     ma = 0
                  end if
               else if ( icase==3 ) then

!              --- common line is i1 - i2

                  call msho26(kmeshc,nelem,i1,i2,jelem,i4)

!                 --- Check whether i1-i2 is internal boundary:

                  ja = 0
                  call msho42( kinbnd, lenbnd, i1, i2, ja )
                  if ( ja==1 ) i4 = 0
                  if ( i4>0 .and. jelem/=ielem ) then

!                 --- consider i3 - i1 - i2 - i4

                     j3 = i3
                     call msho27(coor,i3,i1,i2,i4)
                     if ( j3/=i3 ) ifill = 1

                     j1 = i3
                     j2 = i1
                     j3 = i2
                     j4 = i4
                  else
                     ma = 0
                  end if
               end if

!              --- Refill element array

               if ( ma>0 ) then
                  kmeshc( 1,ielem ) = j1
                  kmeshc( 2,ielem ) = j2
                  kmeshc( 3,ielem ) = j3

                  kmeshc( 1,jelem ) = j2
                  kmeshc( 2,jelem ) = j4
                  kmeshc( 3,jelem ) = j3
               end if
            end if
      end do

      if ( ierror/=0 ) then
         write(irefwr,*) 'Just past angle checking, ierror =',ierror
         goto 1000
      end if

!     --- Check this mesh (all diagonals should are optimal now)

      dismin = rinfin

      do ielem = 1, nelem
         i1 = kmeshc (1,ielem)
         i2 = kmeshc (2,ielem)
         i3 = kmeshc (3,ielem)
         call msho33 ( coor, verh, i1, i2, i3 )
         call msho19 ( coor, i1, i2, i3, surf )
         if ( surf<0 ) then
            verh = -verh
            if ( debug ) then
               write(irefwr,*) 'nherha = ',nherha
               write(irefwr,*) 'verh   = ',verh
               write(irefwr,*) 'Triangle',ielem,' - ',i1,i2,i3
            end if
         end if
         dismin = min( verh, dismin )
      end do  ! ielem = 1, nelem

!     --- Check whether this mesh is better:

      if ( dismin>ratiop ) then

!     --- Better mesh:

         if ( debug ) then
            write(irefwr,*) 'ratio new mesh (diagonals) ',dismin
            write(irefwr,*) 'ratio',dismin,'>old ratio',ratiop
            write(irefwr,*) 'safe the better one'
         end if

!        --- Use "best" mesh:
!            Set equal to new number of points:

         npntmp = npoint

         if ( npntmp>npndef ) then
            write(irefwr,*) 'Error in allocation temporary pointsarray'
            write(irefwr,*) 'Please inform SEPRA'
            write(irefwr,*) 'Execution stopped'
            go to 1000
         end if

!        --- Coordinates:

         do i = 1, npoint
             coortmp(1,i) = coor(1,i)
             coortmp(2,i) = coor(2,i)
         end do

!        --- Set number of elements:

         neltmp = nelem
         if ( neltmp>neldef ) then
            write(irefwr,*) 'Error in allocation temporary elemarray'
            write(irefwr,*) 'Please inform SEPRA'
            write(irefwr,*) 'Execution stopped'
            go to 1000
         end if

!        --- Elements:

         do i = 1, nelem
            meshtmp(1,i) = kmeshc(1,i)
            meshtmp(2,i) = kmeshc(2,i)
            meshtmp(3,i) = kmeshc(3,i)
         end do 

         ratiop = dismin
      end if

      if ( debug .and. icount>=idebug ) then
         write(irefwr,*) 'extra check, kstap, npoint, icount ',
     +                    kstap, npoint , icount
         call msho08 ( kstapl, kstap, coor, xstart, dismin,
     +                 holeinfo, nholes, .true. )

!        --- The next statements can be used if the process "hangs"

         write(irefwr,*) 'kstap, nochan,nelem ', kstap, nochan,nelem
      end if

!     --- If convergence seems to be difficult and results are "acceptable"
!         return

      if ( reposition .and. ( nherha>10 .and. ifill==0 ) .or.
     +     nherha>20 ) goto 900

!     --- Repositioning

      nherha = nherha + 1
      leng   = nrepos

      call msho16(kmeshc,nelem,npoint,nipnt,itri,kstapl,leng)

      if ( leng>nrepos ) then
         write(irefwr,*) 'Just after msho16',leng,'>',nrepos
         call errsub ( 903, 0, 0, 0 )
         goto 1000
      end if

      coa = dismin

      ! call msho18(kmeshc,nelem,coor,npoint,nipnt,itri,kstapl,coa)  ! Perform a Laplacian repositioning of the nodes

      if ( debug .and. icount>=idebug ) then
         jtimes = 1
         jtimes = 3
         jtimes = 0
         call msho08 ( kstapl, kstap, coor, xstart, dismin,
     +                 holeinfo, nholes, .true. )

!        --- The next statements can be used if the process "hangs"

         write(irefwr,*) 'kstap, nochan,nelem ', kstap, nochan,nelem
      end if

!     --- Check result of Repositioning

      do i = nelemi+1, nelem
         i1 = kmeshc( 1,i )
         i2 = kmeshc( 2,i )
         i3 = kmeshc( 3,i )

!        --- Check only triangles that are fully inside:

         if ( i1>nipnt .and. i2>nipnt .and. i3>nipnt ) then

            call msho19(coor,i1,i2,i3,surf)

            valare = chelp(icube(i1))+chelp(icube(i2))+chelp(icube(i3))
            valare = valare / 3d0

            if ( surf<0d0 .or.
     +           surf<0.1d0 * valare .and. nherha<50 .or.
     +           surf>1.0d1 * valare .and. nherha<50 ) then

!           --- Redivision for this triangle

               ielem = i

               call msho30( kmeshc, nelem, ielem, kstapl, kstap, npoint,
     +                      itri  , nelemi )

!              --- Make a new division for this array kstapl

               if ( icount>=idebug ) then

                  write(irefwr,*) 'nherha, leng ', nherha, leng
                  call msho08 ( kstapl, kstap, coor, xstart, dismin,
     +                          holeinfo, nholes, .true. )

!                 --- The next statements can be used if the process "hangs"

                  write(irefwr,*) 'kstap, nochan,nelem ',
     +                             kstap, nochan,nelem
               end if

               goto 300

            end if

         end if

      end do

      if ( ierror/=0 ) then
         if ( debug )
     +      write(irefwr,*) 'Just before redivision, ierror =',ierror
         goto 1000
      end if

!     --- No redivision needed, check again for points with four neighbours

      ja = 0

      do i = nipnt+1,npoint
         no = itri(i) - itri(i-1)
         if ( no<5 ) ja = 1
      end do

      if ( debug ) write(irefwr,*) 'At the end of loop nherha = ',nherha

      if ( nherha<5 .or. ja==1 .or. ifill==1 ) goto 540

900   continue

!     --- Ready, no error, but check before leaving whether the temporary
!         mesh as stored, is better by comparison of the ratio-values:

      dismin = rinfin

      do ielem = 1, nelem

         i1 = kmeshc(1,ielem)
         i2 = kmeshc(2,ielem)
         i3 = kmeshc(3,ielem)
         call msho33 ( coor, verh, i1, i2, i3 )
         call msho19 ( coor, i1, i2, i3, surf )

         if ( surf<0 ) verh = -verh

         dismin = min(verh,dismin)

      end do  ! ielem = 1, nelem

!     --- Compare

      if ( debug ) then
         write(irefwr,*) 'Check of last mesh'
         write(irefwr,*) 'ratiop = ',ratiop
         write(irefwr,*) 'dismin = ',dismin
         write(irefwr,*) 'nherha = ',nherha
      end if

      if ( dismin<ratiop - 100 * eps ) then
         if ( debug ) then
            write(irefwr,*) 'Repositioning problems: use old mesh'
            write(irefwr,*) 'ratio old (saved) mesh is ',ratiop
            write(irefwr,*) 'last mesh, nherha=',nherha,'ratio =',dismin
         end if

!        --- Use "best" mesh:
!            Set number of points:

         npoint = npntmp

!        --- Coordinates:
         do i = 1, npoint
             coor(1,i) = coortmp(1,i)
             coor(2,i) = coortmp(2,i)
         end do

!        --- Set number of elements:

         nelem  = neltmp

!        --- Elements:
         do i = 1, nelem
             kmeshc(1,i) = meshtmp(1,i)
             kmeshc(2,i) = meshtmp(2,i)
             kmeshc(3,i) = meshtmp(3,i)
         end do
      end if

!     --- Make backtransformation for the nodes:

      do i = 1, npoint
          coor(1,i) = (coor(1,i) - 1d0)*tran + xmint
          coor(2,i) = (coor(2,i) - 1d0)*tran + ymint
      end do

      if ( debug ) then
         write(irefwr,*) 'New coor at end MSHO2D'
         do i = 1, npoint
            write(irefwr,*) 'Node',i,' x,y = ',coor(1,i),coor(2,i)
         end do
         write(irefwr,*) 'Results Ok?'
         write(irefwr,*) 'Controle values'
         write(irefwr,*) 'Initial guess npunt  = ',npunt
         write(irefwr,*) 'Finally found npoint = ',npoint
         write(irefwr,*) 'Guess for number of elements  ',nkmeshc
         write(irefwr,*) 'Number of elmnts finally found',nelem
      end if

      if ( debug .and. kstap>0 ) then

!     --- Information statements:

          write(irefwr,*)
     +      'End msho2d: Check barycenters of',nelem,'elements'
          iperm = 0
          call msho40 ( coor, kmeshc, npoint, nelem, iperm )
          write(irefwr,*) 'Ready with check, iperm = ',iperm

!         --- Compute total surface of all elements and check for negative
!             values:

         surf = 0

         do i = 1, nelem
            i1 = kmeshc( 1,i )
            i2 = kmeshc( 2,i )
            i3 = kmeshc( 3,i )

!           --- Check triangle and add surface:

            call msho19(coor,i1,i2,i3,surf1)

            if ( surf1<0 ) then
               write(irefwr,*) 'error,surface',i1,i2,i3,'neg = ',surf1
            end if

            surf = surf + surf1

         end do

         write(irefwr,*) 'Surface (all elements added) is ',surf

         write(irefwr,*) 'kstap = ',kstap

         if ( kstap>0 ) then
            write(irefwr,*) 'kstap =',kstap,' end msho2d'
            write(irefwr,*) (kstapl(i),i=1,2*kstap)
         end if

         write(irefwr,*) 'Output meshgeneration elements and nodes'
         write(irefwr,*) 'End msho2d with ',npoint,' nodes'
         do i=1, npoint
            write(irefwr,*) i,(coor(ja,i),ja=1,2)
         end do

         write(irefwr,*) 'ok? '
         write(irefwr,*) 'End msho2d'
         write(irefwr,*) 'Number of elements = ', nelem
         write(irefwr,*) 'Contents array kbndpt, value nbndpt = ',nbndpt
         write(irefwr,*) (kbndpt(i),i=1,nbndpt)

      end if

!     --- Check all triangles of this surface for the ratio-value:

      dismin = rinfin
      ratiop = 0

      do ielem = 1, nelem

         i1 = kmeshc(1,ielem)
         i2 = kmeshc(2,ielem)
         i3 = kmeshc(3,ielem)
         call msho33 ( coor, verh, i1, i2, i3 )
         ratiop = ratiop + abs(verh)

         dismin = min(verh,dismin)

      end do  ! ielem = 1, nelem

      if ( debug ) then

!     --- Compute mean value:

         ratiop = ratiop/nelem
         write(irefwr,*) 'Smallest ratio in mesh =',dismin
         write(irefwr,*) 'Mean ratio in mesh = ',ratiop
      
      end if

!     --- Check whether values have been changed and act accordingly:

      if ( maxratio>1d0 .and. dismin<0.15d0 ) then

!     --- Suggest new values for the coarsenesses of the userpoints

         write(irefwr,*) 'Bad-shaped elements found in surface',isurnr
         write(irefwr,*) 'smallest ratiovalue = ',dismin
         write(irefwr,*) 'Use eventually the following coarseness'
         write(irefwr,*) 'values for the',nuspnt,'submitted userpoints'
         do i=2, nuspnt+1
            write(irefwr,*) 'node',i-1,'coarseness',coaval(i)
         end do
      end if

!     --- Delete the temporary helparrays

950   deallocate( userco  )
      deallocate( coaval  )
      deallocate( chelp   )
      deallocate( cube    )
      deallocate( icube   )
      deallocate( itri    )
      deallocate( jcube   )
      deallocate( kstapl  )

      if ( neltmp>0 ) then
         deallocate( coortmp )
         deallocate( meshtmp )
      end if

      if ( lenbnd>0 ) deallocate ( kinbnd )

!     --- Close

1000  call erclos( 'msho2d' )

      end
