      subroutine msho38( npoint, coor  , dist  , xstart, ystart,
     +                   kbndpt, nbndpt, numcurvboun   , boundary,
     +                   nbound, nholes, nx    , ny    , icube ,
     +                   coarse, cube  , jcube , kstapl, kstap ,
     +                   ncurvs, curves, cocurv, isurnr, userco,
     +                   nuspnt, coaval, tran   )
! ======================================================================
!
!        programmer    niek praagman
!        version  1.1  date 07-12-2009 cmin and cmax etc adjusted
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
!     Subroutine to check the coarseness in surface isurnr in relation to
!     smoothness requirements. If the coarseness changes too quickly locally
!     the program suggests "better" coarsenessvalues for the userpoints
!     which are stored in array coaval.
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
!      include 'SPcommon/cmcdpr'
!      include 'SPcommon/cmcdpi'
!      include 'SPcommon/cconst'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer    npoint, nx, ny, icube(*), jcube(*), ncurvs,
     +           curves(ncurvs), isurnr, kstap, kstapl(2,kstap),
     +           numcurvboun, boundary(2,*),nbound, nholes,
     +           nbndpt, kbndpt(nbndpt), nuspnt
      double precision coor(2,*), dist, xstart, ystart, coarse(*),
     +                 cube(*)  , cocurv(2,*),
     +                 userco(2,*), coaval(nuspnt+2), tran

!     boundary      i   array of size 2 x numcurvboun containing
!                       information of the boundary of the surface
!                       (1,i) contains the i-th curve number (including sign)
!                       (2,i) contains the starting address in kbndpt
!     coarse        i   array with coarsenesses of each point
!     coaval       i/o  array with coarsenesses user points
!     cocurv        i   array with coordinates internal lines
!     coor          i   coordinate array
!     cube          i   array with coarsenesses of each quadrilateral
!     curves        i   array with number of nodes for each internal curve
!     dist          i   an avarage distance of two neighbour points in
!                       the contour
!     icube         i   array with for each node quadrilateral that node
!                       belongs to
!     isurnr        i   number of surface considered
!     jcube         i   array with indication whether quadrilateral is
!                       completely outside, completely inside or a
!                       boundary quadrilateral
!     kbndpt        i   array with node numbers boundary sequentially
!     kstap         i   number of bounday pieces in array kstapl
!     kstapl        i   bounday array (all linepieces sequentially)
!     nbound        i   number of nodes in boundary curves - 1
!     nbndpt        i   total number of points in array kbndpt
!     ncurvs        i   number of internal curves
!     nholes        i   number of holes in surface
!     npoint        i   number of points in boundaries and internal lines
!     numcurvboun   i   number of curves in the boundary of the surface
!     nuspnt        i   number of userprescribed nodes
!     nx            i   number of quadrilaterals in x-direction
!     ny            i   number of quadrilaterals in y-direction
!                       considered
!     tran          i   scalingsfactor
!     userco        i   array with coords userpoints
!     xstart        i   smallest x-coordinate of enveloping quadrilateral
!     ystart        i   smallest y-coordinate of enveloping quadrilateral
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision af1, af2, afst, cmax, cmin,
     +      coa, cochng, csmall, dif, dx1, dx2, eps,
     +      maxratio, val1, val2, x1, x2, x3, xt, y1, y2, y3, yt

      integer i, i1, i2, i3, i4, iallow, icheck, ilrg, iperm,
     +        isml, istart, ius1, ius2, j, k, kpoint, n1, n2

      integer locbound(numcurvboun+ncurvs)

      logical debug

!     af1          local distance parameter
!     af2          local distance parameter
!     afst         local distance between sequential nodes
!     cmax         largest local coarseness
!     cmin         smallest local coarseness
!     coa          helpvariable to store temporarily coarseness value
!     cochng       helpvariable suggested coarseness adaption
!     csmall       helpvariable to store smallest coarseness value
!     debug        logical indicator for debugging
!     dif          variable to store difference in coarsenessvalues
!     dx1          Euclidian distance
!     dx2          Euclidian distance
!     eps          accuracy parameter
!     i            loop variable
!     i1           node number
!     i2           node number
!     i3           node number
!     i4           node number
!     iallow       indicator whether node is allowed
!     icheck       indicator for visibility
!     ilrg         number of coarsenesscube with largest value
!     iperm        indicator whether there are coarsenessproblems
!     isml         number of coarsenesscube with smallest value
!     istart       helpvariable in search process boundary nodes
!     ius1         first userpoint ("startpoint curve")
!     ius2         second userpoint ("endpoint curve")
!     j            loop variable
!     k            loop variable
!     kpoint       helpnumber for node with problemcoarseness
!     maxratio     max ratio allowed between two successive elements
!     n1           help variable reference cube of node
!     n2           help variable reference cube of node
!     val1         helpvariable coarseness
!     val2         helpvariable coarseness
!     x1           x-coordinate
!     x2           x-coordinate
!     x3           x-coordinate
!     xt           x-coordinate
!     y1           y-coordinate
!     y2           y-coordinate
!     y3           y-coordinate
!     yt           y-coordinate
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     INSTOP    Stop execution
!     MSHO03    Compute Euclidian distance
!     MSHO24    Check whether line does not cross outer or inner boundary
!     MSHO34    Find number of internal curve and coarseness values
!     MSHO39    Compare coarseness and Euclidian distance
!     MSHO41    Determine needed coarseness in given node
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
!    Run through all boundary lines and check the ratio's of subsequent
!    boundary elements
!    Run through all possible lines of boundary points and check (if the
!    points are visible!) the mutual ratio's
!    Give warning if ratio's are wrong and compute acceptable coarseness
!    values using the location of the usersubmitted userpoints
!
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Set debug

      debug = .false.

      if ( debug ) then
         write(irefwr,*) 'Debuginfo from msho38'
         write(irefwr,*) 'npoint = ',npoint,'nbndpt =',nbndpt
!     --- Check contents of submitted curves and nodes
         write(irefwr,*) 'Check submitted boundarynodes in curves:'

         istart = 0

         do i = 1 , numcurvboun
!        --- Find (local) positions of this curve in boundary
            if ( i<numcurvboun ) then
               n1 = boundary(2,i  )
               n2 = boundary(2,i+1)
            else if ( i==numcurvboun) then
               n1 = boundary(2,i)
               n2 = nbound + 1 + nholes
            endif

            if ( istart==0 ) istart = kbndpt(n1)

            if ( kbndpt(n2 - 1) == istart ) then
                 n2 = n2 - 1
                 istart = 0
            end if
            ! Give info about curve, its and nodes and coords
            write(irefwr,*) 'Curve',boundary(1,i)
            write(irefwr,*) 'with Nodes',n1,'-',n2
            do j = n1, n2
               xt = coor(1,kbndpt(j) )
               yt = coor(2,kbndpt(j) )
               write(irefwr,*) 'Node',kbndpt(j),xt,yt
            end do
         end do

         ! Also info with respect to internal curves
         do i = 1 , ncurvs
            write(irefwr,*) 'Curve',i,'number of nodes',curves(i)
!           --- Find (local) positions of this curve in surface
            n1 = n2 + 1
            n2 = n1 + curves(i) - 1
            write(irefwr,*) 'Internal Curve with Nodes',n1,'-',n2
            do j = n1, n2
               xt = coor(1,kbndpt(j))
               yt = coor(2,kbndpt(j))
               write(irefwr,*) 'Node',kbndpt(j),xt,yt
            end do
         end do
      end if

      debug = .false.

!     --- Set tuning factor:
      maxratio = coaval( nuspnt+2 )

!     --- Run through all cubes and determine the smallest and largest
!         coarseness in the region considered

      coa  = rinfin
      isml = 0

      do i = 1, nx * ny
!     --- Check each cube that has been filled:
         if ( cube(i)<coa .and. jcube(i)>0 ) then
            coa  = cube ( i )
            isml = i
         end if
      end do

!     --- Set smallest value:

      csmall = cube(isml)

      coa  = 0d0

      do i = 1, nx * ny
!     --- Check each cube
         if ( cube(i)>coa .and. jcube(i)>0 ) then
            coa  = cube ( i )
            ilrg = i
         end if
      end do

!     For future use???
!      clarge = cube(ilrg)

!     --- If Debug determine coordinates of these cubes:
      if ( debug ) then
         write(irefwr,*) 'nx = ',nx,' en ny = ',ny
         write(irefwr,*) 'stepsize dist = ',dist
         write(irefwr,*) 'xstart = ',xstart,'ystart = ',ystart
         write(irefwr,*) 'xend   = ',xstart+nx*dist,'yend  = ',
     +                    ystart+ny * dist
         write(irefwr,*) 'en de y max = ',ystart+ny*dist

!        --- Compute n1 and n2 for cube with smallest coarseness:
         n2 = int( isml/nx )
         n1 = isml - 1 - n2 * nx

         write(irefwr,*) 'Cube smallest coarseness'
         write(irefwr,*) 'x- and y-coordinates centre of cube:'
         write(irefwr,*) ' x=',xstart+(n1+0.5)*dist
         write(irefwr,*) ' y=',ystart+(n2+0.5)*dist

!        --- Compute n1 and n2 for cube with largest coarseness:
         n2 = int( ilrg/nx )
         n1 = ilrg - 1 - n2 * nx

         write(irefwr,*) 'Cube largest coarseness'
         write(irefwr,*) 'coordinates:'
         write(irefwr,*) ' x=',xstart+(n1+0.5)*dist
         write(irefwr,*) ' y=',ystart+(n2+0.5)*dist

      end if

!     --- Check all local linepieces, running along
!                 the boundary

!     --- Initialize locbound array
      do i = 1, numcurvboun
         locbound(i) = 0
      end do

      istart = 1

!     --- Eventually (not used now) check the boundary of this domain
      if ( istart<0 ) then

         do i = 1, kstap

            i1 = kstapl(1,i)
            i2 = kstapl(2,i)

            if ( i<kstap ) then
               i3 = kstapl(1,i+1)
               i4 = kstapl(2,i+1)
            else
               i3 = kstapl(1,istart)
               i4 = kstapl(2,istart)
            end if

            x1 = coor(1,i1)
            y1 = coor(2,i1)

            x2 = coor(1,i2)
            y2 = coor(2,i2)

            dx1 = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )

            if ( i2==i3 ) then
!           --- Check pieces i and i+1

               x3 = coor(1,i4)
               y3 = coor(2,i4)

            else if ( i2==kstapl(1,istart) ) then
!           --- Jump in boundary:

               if (debug) write(irefwr,*) 'Jump in boundary',i1,i2,i3,i4

               x3 = coor(1,istart+1)
               y3 = coor(2,istart+1)

               istart = i+1

            else
!           --- Unexpected situation: error

               write(irefwr,*) 'Unexpected situation'
               write(irefwr,*) 'If your input is ok'
               write(irefwr,*) 'Please inform SEPRA'
               write(irefwr,*) 'Execution terminated'
!AvD               call instop

            end if

            dx2 = sqrt( (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) )

            if ( dx2/dx1 > maxratio .or. dx1/dx2 > maxratio ) then
!           --- Coarseness problems in node i2:

               write(irefwr,*) 'Problem for node',i2,'with crds',x2,y2
!              --- Find curve number for node i2:
               do j = 1 , numcurvboun
                  write(irefwr,*) j,boundary(1,j),boundary(2,j)
!                 --- Find positions in boundary
                  n1 = boundary(2,j)
                  if ( j==numcurvboun) then
                     n2 = nbndpt
                  else
                     n2 = boundary(2,j+1)
                  end if
                  write(irefwr,*) 'Curve',boundary(1,j), 'Nodes'
                  do k  = n1, n2
                     xt = coor(1,k)
                     yt = coor(2,k)
                     write(irefwr,*) 'Node',kbndpt(k),xt,yt

                     if ( i2==kbndpt(k) ) then
                        write(irefwr,*) 'Node is in curve',boundary(1,j)
                        locbound(j) = locbound(j)+1
                     end if
                  end do
               end do
            end if

         end do

         iperm = 1
         do i = 1, numcurvboun
            if ( locbound(i)>0 ) then
               iperm = 0
               write(irefwr,*)
     +          'Coarseness problems curve C with number',boundary(1,i),
     +          'which is part of the boundary of surface',isurnr,
     +          'coarsenesses start- and endnode differ too much!'
            end if
            locbound(i) = 0
         end do

         if ( iperm==0 ) then
            write(irefwr,*) 'Check and adjust your input for the curves'
         end if

      end if

!     --- Second step: check all internal mutual lines
!     --- Check all nodes via internal "visible" lines:

!     --- Run through all prescribed points and consider mutual distances in relation
!         to their local coarsenesses:
!         (Give warning if points are not in one line)

      iperm = 1

      eps = 1d-9 * coaval(1)

      do i =   1, npoint-1
      do j = i+1, npoint
!        --- Be sure that nodes are really in linear mesh:
         if ( coarse(i)> eps .and. coarse(j)>eps ) then
!     --- Next step: check mutual visibility and
!         if visible then check coarsenesses

         if ( debug ) write(irefwr,*) 'Check Point',i,'and point',j

!        --- Compute Euclidian distance
         call msho03( i, j, coor, afst )

!        --- Check whether coarsenesses and mutual distance are
!            realistic. Find cmax and cmin of these two values:
         if ( coarse(i)>coarse(j) ) then
              cmax    = coarse(i)
              cmin    = coarse(j)
              kpoint  =  i
         else
              cmax    = coarse(j)
              cmin    = coarse(i)
              kpoint  =  j
         end if

         coa = abs( (cmax - cmin)/(cmax+cmin) )

!        --- Check possibilities:
         if ( afst<0.2 * csmall ) then
!        --- Check whether line does not cross outer or inner boundary:
            call msho24( kstapl, kstap, coor, i, j, icheck)

            if ( icheck==0 ) then
!           --- No intermediate points, give warning: too small distance

!           --- Find number of internal curve and values of coarseness
               call msho34( kpoint, coor, ncurvs, curves, kbndpt,
     +                      nbndpt, numcurvboun, boundary, nbound,
     +                      nholes, userco, nuspnt, coaval, ius1, ius2 )
               if ( debug ) then
                  write(irefwr,*) 'msho38: problems in surface',isurnr
                  write(irefwr,*) 'Locally coarseness too small'
                  write(irefwr,*) 'Found value ',afst
                  write(irefwr,*) 'Position of nodes'
                  write(irefwr,*) ' Node ',i,' = ',coor(1,i),coor(2,i)
                  write(irefwr,*) ' Node ',j,' = ',coor(1,j),coor(2,j)
                  write(irefwr,*) 'Adjust coarsenesses of nodes '
                  iperm = 0
               end if
            end if
         else if ( coa > 0.1 ) then
!        --- Only action if there is a real coarseness difference
            iallow =  1
            call msho39( cmax, cmin, afst, maxratio, iallow )
            if ( iallow==0 .and. debug )
     +         write(irefwr,*) 'Routine msho39 iallow=0 use cmax =',cmax
            icheck = 0
            if ( iallow==0 ) then
!           --- Check visibility:
               call msho24( kstapl, kstap, coor, i, j, icheck)
            end if

            if ( iallow==0 .and. icheck==0 ) then
!           --- Point is visible and coarseness is not ok
               iperm = 0

               if ( debug ) then
                  write(irefwr,*) 'Warning coarseness nodes',i,'and',j
                  write(irefwr,*) 'differs too much'
                  write(irefwr,*) 'Positions of the nodes are'
                  write(irefwr,*) ' Node ',i,' = ',coor(1,i),coor(2,i)
                  write(irefwr,*) ' Node ',j,' = ',coor(1,j),coor(2,j)
                  write(irefwr,*) 'Distance = ',afst
                  write(irefwr,*) 'coarseness  first point = ',coarse(i)
                  write(irefwr,*) 'coarseness second point = ',coarse(j)
                  write(irefwr,*) 'Adjust Coarseness'

               end if
!              --- Determine coarseness needed in kpoint:
               call msho41( cmin, cmax, afst, maxratio )
!              --- Find start- and end-point of line:
               call msho34( kpoint, coor, ncurvs, curves, kbndpt,
     +                      nbndpt, numcurvboun, boundary, nbound,
     +                      nholes, userco, nuspnt, coaval, ius1, ius2 )

!              --- Determine new coarseness needed in ius1 and ius2:
               if ( debug ) then
                  write(irefwr,*) 'cords kpoint',(coor(k,kpoint),k=1,2)
                  write(irefwr,*) 'cords ius1  ',(userco(k,ius1),k=1,2)
                  write(irefwr,*) 'cords ius2  ',(userco(k,ius2),k=1,2)
                  write(irefwr,*) 'coarsen i1  ',(coaval(ius1+1))
                  write(irefwr,*) 'coarsen i2  ',(coaval(ius2+1))
                  write(irefwr,*) 'coarsen kp  ',cmax
                  write(irefwr,*) 'maxratio =',coaval(1),'tran=',tran
                  write(irefwr,*) (coaval(k),k=2,nuspnt+1)
               end if
!              --- Compute coarseness values in ius1 and ius2:
!                  (Use transformation coefs)
               val1 = coaval(1)*coaval(ius1+1)/tran
               val2 = coaval(1)*coaval(ius2+1)/tran
!              --- Adjust values in coaval array:
               if ( cmax<val1 .and. cmax<val2 ) then
!              --- Use a constant value along the curve
                  coaval(ius1+1) = cmax * tran / coaval(1)
                  coaval(ius2+1) = cmax * tran / coaval(1)
               else
!              --- Interpolate
                  x1   = userco(1,ius1)
                  y1   = userco(2,ius1)

                  x2   = coor(1,kpoint)
                  y2   = coor(2,kpoint)

                  x3   = userco(1,ius2)
                  y3   = userco(2,ius2)

                  af1  = sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) )
                  af2  = sqrt( (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3) )

!                 --- Compute value in kpoint as result of val and val2:
                  dif  = val1 + af1 * (val2-val1)/(af1+af2)
!                 --- Determine the difference:
                  dif  = dif - cmax

                  if ( (val1-dif)<0.5*cmax ) val1 = dif + 0.5*cmax
                  if ( (val2-dif)<0.5*cmax ) val2 = dif + 0.5*cmax
!                 --- Adjust the new values:
                  coaval(ius1+1) = tran/coaval(1)*(val1-dif)
                  coaval(ius2+1) = tran/coaval(1)*(val2-dif)

                  cochng  =  1

                  if ( debug ) then
                     val1 = coaval(1)*coaval(ius1+1)/tran
                     write(irefwr,*) 'coa-value',ius1,'w=',val1
                     val2 = coaval(1)*coaval(ius2+1)/tran
                     write(irefwr,*) 'coa-value',ius2,'w=',val2
                  end if

               end if
            end if
         end if

         end if ! Both nodes in mesh?

      end do
      end do

!     --- Check whether values have been changed and act accordingly:
      if ( debug .and. cochng==1 ) then
!     --- Give new values of coarsenesses userpoints
         write(irefwr,*) 'Suggested coarseness values userpoints'
         do i=2, nuspnt+1
            write(irefwr,*) 'node',i-1,'coarseness',coaval(i)
         end do
      end if

      end
