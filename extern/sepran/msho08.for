      subroutine msho08 ( kstapl, kstap, coor, xstart, dismin,
     +                    holeinfo, nholes, check )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  5.2  date 24-03-2009 Extra debug facilities
!        version  5.1  date 25-11-2008 Check for double lines
!                                      and counterclockwise
!        version  5.0  date 17-02-2005 Accuracy improvement
!        version  4.0  date 24-02-2003 Extra parameter check
!        version  3.0  date 05-04-2002 Extra parameters
!        version  2.0  date 14-02-1994 New norms
!        version  1.0  date 12-04-1989
!
!   copyright (c) 1989-2009  "Ingenieursbureau SEPRA"
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
!     Check whether all boundary elements are given in the right direction.
!     If not then reverse direction of boundary element.
!     If element is double determine direction by sequencing
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     check
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kstap, kstapl(*), nholes, holeinfo(2,0:nholes+1)
      double precision coor(2,*), xstart, dismin
      logical check

!     check          i    Indicates if this is the first call of msho08, where
!                         possibly the boundary must be reversed (false) or
!                         that the call is just for checking (true)
!     coor           i    coordinate array
!     dismin         i    smallest distance in boundary polygon
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
!     kstap          i    number of elements in array kstapl
!     kstapl        i/o   array with boundary elements
!     nholes         i    Number of holes in the boundary of a surface
!     xstart         i    smallest x-value of reference quadrangles
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, i1, i2, ichange, ih1, ih2, j, j1, j2, ja, icount,
     +        iloop
      double precision eps, e1, e2, xn, yn, xm, ym
      logical debug

!     debug          parameter to indicate whether debug actions need
!                    to be taken
!     e1             first component of unity vector
!     e2             second component of unity vector
!     eps            accuracy parameter
!     i              loop variable
!     i1             first node of boundary line
!     i2             second node of boundary line
!     ichange        indicator whether changings have been made
!     icount         Counter
!     ih1            helpvariable for first node
!     ih2            helpvariable for second node
!     iloop          helpvariable to count
!     j              loop variable
!     j1             helpnode
!     j2             helpnode
!     ja             help indicator
!     xm             x-coordinate midpoint
!     xn             x-coordinate extra point
!     ym             y-coordinate midpoint
!     yn             y-coordinate extra point
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSHCHKSTAPL    Check the contents of array kstapl
!     MSHO07         Check whether point belongs to given area
!     MSHO09         Compute midpoint line and the normal vector
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    1274   internal error
!    2434   Internal error in triangle
! **********************************************************************
!
!                       PSEUDO CODE
!
!     The check is to consider the nodes i1 and i2 of the boundary element
!     together with an extra point which should be in the interior of the
!     area given by KBOUND. Here it is assumed that the area is closed.
!     Set accuracy (here relative to 1)
!     Run through all boundary elements
!     Compute midpoint and "inside" pointing unity normal
!     Compute point that should be inside
!     Check whether point belongs to area
!     If point does not belong to area reverse direction
!     of boundary element
!     Check whether all pieces are simple. If pieces are double check
!     the correct direction and adjust eventually.
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen( 'msho08' )

      debug = .false.

!     --- Set accuracy

      eps = 1d-5

      icount = 0

      do i = 1, kstap

         i1 = kstapl( 2*i - 1 )
         i2 = kstapl( 2*i     )

         kstapl(2*kstap+i) = 0

!        --- Compute midpoint line and normal vector ("counterclockwise")

         call msho09 ( coor, i1, i2, e1, e2, xm, ym )

!        --- Determine extra point for checking

         xn = xm + eps * e1 * dismin
         yn = ym + eps * e2 * dismin

!        --- Check whether point belongs to area

         ja = 0

         call msho07( xn, yn, xstart, coor, kstapl, kstap, ja )

!        --- If not inside direction is reversed

         if ( ja==0 ) then

!        --- ja = 0, reverse or give error message

            if ( check ) then

!           --- check is true, error

               call errint ( i1, 1 )
               call errint ( i2, 2 )
               call errsub ( 2434, 2, 0, 0 )

            else

!           --- check is false, start of procedure

               if ( debug)
     +            write(irefwr,*) 'Interchange',i1,'and',i2
               kstapl( 2*i-1 ) = i2
               kstapl( 2*i   ) = i1

            end if

         end if

         if ( holeinfo(2,icount)==i1 ) then

!        --- New part of boundary found

            if ( ja==0 ) holeinfo(1,icount) = -holeinfo(1,icount)
            if ( icount<nholes) icount = icount+1

         end if

         if ( ja==1 ) then

!        --- Check whether linepiece is "double", using point on
!            other side of line:

            xn = xm - eps * e1 * dismin
            yn = ym - eps * e2 * dismin

!           --- Check whether other point belongs to area and if so
!               indicate that situation: by value 1!

            ja = 0

            call msho07( xn, yn, xstart, coor, kstapl, kstap, ja )

            if ( ja==1 ) kstapl(2*kstap+i) = 1

         end if

      end do

!     --- Debug: Check new kstapl (double pieces have value 1)
!         Double pieces are internal lines or parts of a crack

      if ( debug ) then
         do i = 1, kstap
            write(irefwr,*)
     +           'msho08: piece ',i,' value = ',kstapl(2*kstap+i)
         end do
      end if

!     --- Check all linepieces: all pieces should be connected!

      iloop   = 0
      ichange = 1
      ih1     = 0
      ih2     = 0

      do while ( ichange==1 )

         ichange = 0

         if ( debug ) then
            do i=1, kstap
               if ( kstapl(2*kstap+i)==1 ) then
                  write(irefwr,*) 'elmnt',i,'not yet definitive'
               end if
            end do
         end if

!        --- Run through all "correct" pieces and find element that has no
!            connection:

         do i = 1, kstap-1

            if ( kstapl(2*kstap+i)==0 ) then

!           --- Check this correct piece i:

               i1 = kstapl(2*i-1)
               i2 = kstapl(2*i  )

               ih1 = 1
               ih2 = 1

!              --- Check this correct element for connections:

               do j = 1, kstap
                  if ( kstapl(2*kstap+j)==0 .and. j/=i ) then

                     j1 = kstapl(2*j-1)
                     j2 = kstapl(2*j  )

                     if ( j1==i2 ) ih2=0
                     if ( j2==i1 ) ih1=0

                  end if
               end do

               if ( ih1+ih2>0 ) then

!              --- Piece has a problem with connections:

                  if ( debug ) then
                     write(irefwr,*) 'Piece',i
                     if ( ih1>0 ) write(*,*) 'No connection',i1
                     if ( ih2>0 ) write(*,*) 'No connection',i2
                  end if

!                 --- Run through all "doubles" and make connection

                  do j = 1 , kstap

!                 --- Check for all "double" pieces:

                     if ( kstapl(2*kstap+j)==1 ) then

!                    --- Candidate

                        j1 = kstapl(2*j-1)
                        j2 = kstapl(2*j  )

                        if ( j1==i2 .and. ih2>0 .or.
     +                       j2==i1 .and. ih1>0 ) then

!                       --- "Double" piece is connected

                           kstapl(2 * kstap+j) = 0
                           ih1 = 0
                           ih2 = 0

                        else if ( j1==i1 .and. ih1>0 .or.
     +                            j2==i2 .and. ih2>0 ) then

!                       --- "Double" piece is connected, but direction has
!                            to be changed:

                           kstapl(2 * kstap+j) = 0
                           kstapl(2 * j - 1  ) = j2
                           kstapl(2 * j      ) = j1

                           ichange = 1

                           ih1 = 0
                           ih2 = 0

                        end if

                     end if

                  end do

               end if

            end if

         end do

         if ( debug .and. ih1+ih2>0 ) then
            write(irefwr,*) 'element',i1,i2,' has no connection '
         end if

         iloop = iloop + 1

         if ( iloop>1000 ) then
            call errsub (1274, 0, 0, 0 )
            go to 1000
         end if

      end do

      call mshchkstapl( kstapl, kstap, 'Final check in msho08')

1000  call erclos( 'msho08' )

      if ( debug ) then
         write(irefwr,*) 'kstapl at the end of msho08'
         do i=1, kstap
            write(irefwr,*) i,kstapl(2*i-1),kstapl(2*i)
         end do
      end if

      end
