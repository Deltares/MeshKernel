      subroutine msho06( npoint, coor  , dist  , xstart, ystart,
     +                   nx, ny, icube , chelp , cube  , jcube ,
     +                   kbound, nbound, coar  , ncoar , ncurvs,
     +                   curves, cocurvs )
! ======================================================================
!
!        programmer    niek praagman
!        version  4.1  date 07-01-2011 Tuning improved again
!        version  4.0  date 10-11-2009 New tuning of several parameters
!        version  3.4  date 30-06-2009 other values for xstart and ystart
!        version  3.3  date 10-03-2009 Use smoothing of and do not enlarge
!                                      coarseness finally
!        version  3.2  date 28-08-2008 Use also nodes of internal curves
!        version  3.1  date 02-06-2008 Smooth interpolation coarseness of
!                                      cubes (quadrilaterals)
!        version  3.0  date 11-01-2005 Extra coarseness points, different
!                                      computation coarsenesses blocks
!        version  2.1  date 20-09-2004 Correction upper bound loop
!        version  2.0  date 03-02-1994 New norms
!        version  1.0  date 13-04-1989
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
!     Subroutine to determine for each standard quadrilateral the coarse-
!     ness and for all the boundary points in which quadrilateral they are
!     situated.
!     Also the type of each quadrilateral is determined.
!     At the end of the routine:
!        0 : no points in common with boundary lines
!        1 : belongs partially to inside part
!        2 : belongs totally to inside part
!     Determine (if specified) special inside coarsness points
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

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer    npoint, nx, ny, icube(npoint), jcube(*), ncoar, ncurvs,
     +           nbound, curves(ncurvs), kbound(2,nbound)
      double precision coor(2,*), dist, xstart, ystart, chelp(*),
     +                 cube(nx*ny)  , coar(3,ncoar), cocurvs(2,*)

!     chelp         i   array with coarsenesses of each point
!     coar          i   array with x,y coordinates points and their
!                       coarseness
!     cocurvs       i   array with coordinates internal lines
!     coor          i   coordinate array
!     cube          i   array with coarsenesses of each quadrilateral
!     curves        i   array with number of nodes for each internal curve
!     dist          i   mean distance of two neighbour points in
!                       the contour
!     icube         i   array with for each node quadrilateral that node
!                       belongs to
!     jcube         i   array with indication whether quadrilateral is
!                       completely outside, completely inside or a
!                       boundary quadrilateral
!     kbound        i   array with boundary pieces
!     nbound        i   number of boundary pieces
!     ncoar         i   number of extra coarseness points
!     ncurvs        i   number of internal curves
!     npoint        i   number of points in line elements of contour
!     nx            i   number of quadrilaterals in x-direction
!     ny            i   number of quadrilaterals in y-direction
!                       considered
!     xstart        i   smallest x-coordinate of enveloping quadrilateral
!     ystart        i   smallest y-coordinate of enveloping quadrilateral
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision afst, eps, xp, yp, coarse, coxmin, coxmax,
     +      coymin, coymax, cgem, alpha, deel, cox, coy, coa,
     +      val1, val2, val3, val4, valmax, valmin
      integer i, i1, i2, j, ik, kub, n1, n1i1, n1i2, n1min, n1max,
     +        n2, n2i1, n2i2, n2min, n2max,
     +        nc, ntal, nrkube, nnodes,
     +        ikxmin, ikxmax, ikymin, ikymax, jblok
      logical debug

!     afst         local distance between sequential nodes
!     alpha        helpvalue for interpolation
!     cgem         helpvalue to store mean value
!     coa          smallest coarseness value point
!     coarse       coarseness of points resp cube
!     cox          interpolation value coarseness x-direction
!     coxmax       ending coarseness x-direction
!     coxmin       starting coarseness x-direction
!     coy          interpolation value coarseness y-direction
!     coymax       ending coarseness y-direction
!     coymin       starting coarseness y-direction
!     debug        indicator for debugging
!     deel         check variable to see in which direction values
!                  have been used
!     eps          accuracy
!     i            loop variable
!     ik           loop variable
!     ikxmax       maximum value where x-line stops
!     ikxmin       minimum value where x-line starts
!     ikymax       maximum value where x-line stops
!     ikymin       see ikxmin, now for y
!     j            loop variable
!     jblok        help value for block
!     kub          help for starting address in loop
!     n1           help variable to determine ref number of node
!     n2           help variable to determine reference number of node
!     nc           reference number
!     nnodes       number of nodes in internal lines
!     nrkube       cube number
!     ntal         number of steps
!     val1         coarseness value in cube 1
!     val2         coarseness value in cube 2
!     val3         coarseness value in cube 3
!     val4         coarseness value in cube 4
!     valmax       max value surrounding coarsenesses
!     valmin       min value surrounding coarsenesses
!     xp           x-coordinate
!     yp           y-coordinate
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
!     Run through all boundary points and determine the number of the
!     quadrilateral point belongs to
!     Run through all quadrilaterals and check for each quadrilateral
!     whether a coarseness-value is given
!     Determine for all empty quadrilaterals which are inside the area
!     also a value for the coarseness by interpolation
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Set debug
      debug = .false.

      eps = 10 * epsmac

      coa = rinfin

!     --- Determine for all boundary points the cube-numbers
      do i = 1, npoint

          xp = coor( 1,i )
          yp = coor( 2,i )

!         --- Find number of quadrilateral

          n1 = int ( (xp-xstart)/dist )
          n2 = int ( (yp-ystart)/dist )

!         --- Determine smallest coarseness value

          if ( chelp(i)>eps .and. chelp(i)<coa )
     +         coa = chelp(i)

!         --- This point belongs to quadrilateral with number

          nc = 1 + n1 + n2*nx
          icube(i) = nc

          if ( debug )
     +       write(irefwr,*) 'Point',i,'x=',xp,'y=',yp,'cube',nc

      end do

      if ( debug ) then
         write(irefwr,*) 'Contents icube due to',npoint,'boundary nodes'
         write(irefwr,*) 'nx =',nx,'ny =',ny,'dist =',dist
         write(irefwr,*) 'Grid contains ',nx*ny,' boxes'
         write(irefwr,*) (icube(i),i=1,npoint)
      end if

!     --- Consider for each quadrilateral whether a coarseness-value is given

      do i = 1, nx * ny
          jcube(i) = 0
          cube (i) = 0
          ntal     = 0
          coarse   = 0d0
          do ik = 1, npoint
             if ( chelp(ik)>eps ) then
                nc = icube(ik)
                if ( nc==i ) then
                     ntal   = ntal + 1
                     coarse = coarse + chelp(ik)
                end if
             end if
          end do
          if ( ntal>0 ) then
               coarse   = coarse / ntal
               cube(i)  = coarse
               jcube(i) = 1
          end if
      end do

      if ( debug ) then
         write(irefwr,*) 'Check: cubes and cubevalues, dist =',dist
         do j=1, ny
            yp = ystart + (j-0.5)*dist
         do i=1, nx
            nc = i + (j-1)*nx
            xp = xstart + (i-0.5) * dist
            write(irefwr,*) 'Cube, i=',i,'j= ',j,'nc =',nc
            write(irefwr,*) 'Center of cube',xp,yp,'coarse',cube(nc)
         end do
         end do
      end if

!     --- Fill empty places along boundary contour:
      do i=1, nbound

         i1 = kbound(1,i)
         i2 = kbound(2,i)

         cgem = ( chelp(i1)+chelp(i2) ) / 2d0

         xp = coor(1,i1)
         yp = coor(2,i1)

         n1i1 = int( (xp-xstart) / dist )
         n2i1 = int( (yp-ystart) / dist )

         xp = coor(1,i2)
         yp = coor(2,i2)

         n1i2 = int( (xp-xstart) / dist )
         n2i2 = int( (yp-ystart) / dist )

         n1min = min( n1i1, n1i2 )
         n1max = max( n1i1, n1i2 )

         n2min = min( n2i1, n2i2 )
         n2max = max( n2i1, n2i2 )

         do n1 = n1min, n1max
         do n2 = n2min, n2max

            nc = 1 + n1 + n2 * nx

            if ( cube(nc)<eps ) then
!           --- Treat as boundary cube:
               cube (nc) = cgem
               jcube(nc) = 1
            end if

         end do
         end do

      end do

!     --- Place extra coarse values for internal lines:
      if ( ncurvs>0 ) then

         nnodes = 0

         do i=1, ncurvs
!        --- Run through line:

            do j= nnodes + 1, nnodes + curves(i) - 1
!           --- Determine start and endpoint

               val1 = cocurvs(1,j+1) - cocurvs(1,j)
               val2 = cocurvs(2,j+1) - cocurvs(2,j)

               afst = sqrt(val1*val1+val2*val2)

               xp = ( cocurvs(1,j+1) + cocurvs(1,j) ) / 2d0
               yp = ( cocurvs(2,j+1) + cocurvs(2,j) ) / 2d0

               n1 = int ( (xp-xstart)/dist )
               n2 = int ( (yp-ystart)/dist )

!              --- This point belongs to quadrilateral with number
               nc = 1 + n1 + n2*nx

!              --- Set values of coarseness and type
               if ( jcube(nc)==0 ) then
                  cube (nc)  = afst
                  jcube(nc)  = 1
               else if ( jcube(nc)>0 ) then
                  cube(nc) = ( cube(nc) + afst ) /2d0
               end if
            end do
            nnodes = nnodes + curves(i)
         end do
      end if

      if ( debug ) then
         write(irefwr,*) 'Second cube-test'
         do i=1, nx * ny
            write(irefwr,*) i,jcube(i),cube(i)
         end do
      end if

!     --- Place extra coarse values for internal points:
      do i = 1, ncoar

         xp = coar( 1, i )
         yp = coar( 2, i )

         n1 = int ( (xp-xstart)/dist )
         n2 = int ( (yp-ystart)/dist )

!        --- This point belongs to quadrilateral with number
         nc = 1 + n1 + n2*nx

!        --- Set values of coarseness and type
         cube (nc)  = coar( 3, i )
         jcube(nc)  = 1

      end do

!     --- All prescribed cubes have now indication: jcube(nc) == 1
      if ( debug ) then
         write(irefwr,*) 'third cube-test, nx=',nx,'ny=',ny
         do i=1, nx * ny
            write(irefwr,*) i,jcube(i),cube(i)
         end do
      end if

!     --- Determine for all empty quadrilaterals which are inside the area
!         also a value for the coarseness
      do i = 1, nx
      do j = 1, ny

         nrkube = i + (j-1)*nx

         if ( jcube(nrkube)==0 .and. cube(nrkube)<eps ) then
!        --- No boundary points in this quadrilateral
!            Determine if possible quadrilaterals and coarseness in i-direction

            kub    = 1 + (j-1)*nx
            coxmin = coa
            ikxmin = 0

            do ik = kub, nrkube-1
               if ( jcube(ik)==1 ) then
                    ikxmin = ik
                    coxmin = cube(ik)
               end if
            end do

            kub    = j*nx
            coxmax = coa
            ikxmax = 0

            do ik = kub, nrkube+1, -1
               if ( jcube(ik)==1 ) then
                    ikxmax = ik
                    coxmax = cube(ik)
               end if
            end do

!           --- Determine quadrilaterals and coarseness in j-direction

            kub    = i
            coymin = coa
            ikymin = 0

            do ik = kub, nrkube-nx, nx
               if ( jcube(ik)==1 ) then
                    ikymin = ik
                    coymin = cube(ik)
               end if
            end do

            kub    = i  + ( ny - 1 ) * nx
            coymax = coa
            ikymax = 0

            do ik = kub, nrkube + nx, -nx
               if ( jcube(ik)==1 ) then
                    ikymax = ik
                    coymax = cube(ik)
               end if
            end do

!           --- Determine weighted coarseness for this quadrilateral

            deel  = 0d0
            cox   = 0d0
            coy   = 0d0
            jblok = ikxmax * ikxmin

            if ( jblok/=0 ) then
!           --- Make choice (exponential or linear)
               if ( jblok>0 ) then
!              --- Exponential:
                  alpha = (coxmax/coxmin)**(1d0/(ikxmax-ikxmin))
                  cox   =  coxmin*(alpha **(1d0*(nrkube-ikxmin)))
               else
!              --- Linear:
                  alpha = (nrkube - ikxmin)*(1D0/(ikxmax-ikxmin))
                  cox   =  coxmin + alpha * (coxmax-coxmin)
               end if
               deel  =  deel + 1d0
               cube(nrkube) = cox
            end if

            jblok = ikymax * ikymin

            if ( jblok/=0 ) then
!           --- Make choice (exponential or linear)
               if ( jblok>0 ) then
!              --- Exponential:
                alpha = (coymax/coymin)**(1d0/(1d0*(ikymax-ikymin)/nx))
                coy   =  coymin*(alpha **(1d0*(nrkube-ikymin)/nx))
               else
!              --- Linear:
                alpha = (nrkube - ikymin)*(1D0/(ikymax-ikymin))
                coy   =  coymin + alpha * (coymax-coymin)
               end if
                deel  =  deel + 1d0
                cube(nrkube) = (cube(nrkube) + coy) / deel
            end if

!           --- Set type of this node:
            if ( deel>0.5d0 ) jcube(nrkube) = 2

         end if
      end do
      end do

      if ( debug ) then
         do i = 1, nx
         do j = 1, ny
            nrkube = i + (j-1)*nx
            write(irefwr,*) 'Cube ',nrkube,' i= ',i,'j= ',j
            write(irefwr,*) 'type',jcube(nrkube),'value ',cube(nrkube)
         end do
         end do
      end if

!     --- Fill all quadrilaterals with a minimum value if no value
!         was found
      valmin = rinfin
      valmax = 0d0

!     --- Determine smallest and largest cubecoarseness value
      do i = 1, nx * ny
         if ( cube(i)>eps .and. cube(i)<valmin ) valmin = cube(i)
         if ( cube(i)>eps .and. cube(i)>valmax ) valmax = cube(i)
      end do

      coa = ( valmax + 3 * valmin ) / 4d0
!     --- Set value coa in all "empty" cubes
      do i = 1, nx * ny
         if ( cube(i)<valmin ) then
!        --- Use mix between smallest and largest coarseness value
            if ( debug ) write(irefwr,*) 'empty cube nr =',i
            cube(i) = coa
         end if
      end do

      if ( debug ) then
         write(irefwr,*) 'Fourth cube-test'
         do i=1, nx * ny
            write(irefwr,*) i,jcube(i),cube(i)
         end do
      end if

!     --- Smooth internal, not extremal, values
      if ( nx>2 .and. ny>2 .and. 0>1 ) then

         do ik = 1, 5

!        --- 5 times smoothing should be enough to get "smooth" values

            do i = nx-1, 2, -1
            do j = ny-1, 2, -1

               nrkube = i + (j-1)*nx

               if ( jcube(nrkube)==2 ) then

!              --- Internal cube:

                  val1 = cube(nrkube-nx)
                  val2 = cube(nrkube- 1)
                  val3 = cube(nrkube+ 1)
                  val4 = cube(nrkube+nx)

                  valmin = min( val1, val2, val3, val4 )
                  valmax = max( val1, val2, val3, val4 )

                  if (valmin<cube(nrkube).and.valmax>cube(nrkube)) then

!                 --- Adjust value:

                     cube(nrkube)= (2*valmin+val1+val2+val3+val4)/6d0
                  end if

               end if

            end do
            end do

         end do

      end if

      if ( debug ) then
         write(irefwr,*) 'Fifth cube-test'
         do i=1, nx * ny
            write(irefwr,*) i,jcube(i),cube(i)
         end do
      end if

!     --- Adjust values (not too large differences in neighbouring
!         cubes)

      ik = 1

      do while ( ik==1 .and. nx>2 .and. ny>2 )

         ik = 0

         do i = nx-1, 2, -1
         do j = ny-1, 2, -1

            nrkube = i + (j-1)*nx

            if ( jcube(nrkube)==2 ) then

!           --- Internal cube:

               val1 = cube(nrkube-nx)
               val2 = cube(nrkube- 1)
               val3 = cube(nrkube+ 1)
               val4 = cube(nrkube+nx)

               valmin = min( val1, val2, val3, val4 )

               if ( cube(nrkube)>2.75d0*valmin ) then

!              --- Adjust value:

                  cube(nrkube)=2.75d0*valmin
                  ik = 1
               end if

            end if

         end do
         end do

      end do

!     --- Check in case of debugging

      if ( debug ) then
         write(irefwr,*) 'Values nx = ',nx,' ny = ',ny
         write(irefwr,*) 'dist = ',dist
         do j = 1, ny
         do i = 1, nx

            nrkube = i + (j-1)*nx
            write(irefwr,*)
     +       nrkube,'type',jcube(nrkube),'value ',cube(nrkube)
         end do
         end do
         write(irefwr,*) 'End routine msho06'
      end if

      end
