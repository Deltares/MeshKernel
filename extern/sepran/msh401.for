      subroutine msh401 ( npoint, kmeshc , nelem , istart, ibrp,
     +                    ibrpnt, kelemh, ishape, inpelm, niedge,
     +                    nisurf, nivolm )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  5.7  date 18-11-2000 Include connection elements
!        version  5.6  date 18-04-1997 Extension with interface elements
!        version  5.5  date 05-04-1995 add niedge, nisurf and nivolm to indi-
!                                      cate whether lines, surfaces and volu-
!                                      mes are included or not)
!        version  5.4  date 15-02-1995 add npoint for call of MSH402
!        version  5.3  date 22-10-1993 extra parameter inpelm for ishape=-199
!        version  5.2  date 02-02-1993 new requirements layout etc
!        version  5.1  date 22-08-1991 Zdenek: no dimension, cdc, ce
!        version  5.0  date 10-09-1990 variable number IBRP of neighbours
!                                      Deletion of inpelm
!        version  4.0  date 14-03-1989 2D and 3D elements
!
!   copyright (c) 1986-2000  "Ingenieursbureau SEPRA"
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
!     Fill neighbour arrays (istart and ibrpnt) for refinement
!     In istart the addresses where the neighbours of point i are
!     stored in array ibrpnt are stored. In this routine it is assumed
!     that each point has ibrp neighbours. Neighbours of point i are
!     stored from position istart(i) to istart(i+1)-1 . In the final
!     form only the neighbours with nodal numbers smaller than point i
!     are stored , together with the new points . If a point has more
!     than ibrp neighbours, the extra neighbours are temporarily stored
!     in array kelemh. Arrays istart and ibrpnt are adjusted according-
!     ly in subroutine msh400.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     refine
!     neighbour
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
!                       INPUT / OUTPUT PARAMETERS
!
      integer npoint, kmeshc(*), istart(*), ibrpnt(*), kelemh(*), nelem,
     +        ibrp, ishape, inpelm, niedge, nisurf, nivolm

!     ibrp        i    number of assumed neighbours
!     ibrpnt      o    neighbour array
!     inpelm      i    number of nodes in element ishape
!     ishape      i    type of elements to be considered (see SEPRAN
!                      PROGRAMMERS GUIDE)
!     istart      o    array containing the start positions of
!                      array ibrpnt
!     kelemh      o    helparray for extra neighbours
!     kmeshc      i    array containing the elements
!     nelem       i    number of elements
!     niedge      i    number of nodes to be placed on edges
!     nisurf      i    number of nodes to be placed on faces
!     nivolm      i    number of nodes to be placed in volumes
!     npoint      i    number of nodes in mesh originally
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i(27), ie, ielem, it, j, jelem
      

!     i        node numbers of element
!     ie       first node of diagonal
!     ielem    loop variable elements
!     it       second node of diagonal
!     j        loop variable
!     jelem    helppointer
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS   Deconcatenate name of subroutine from string
!     EROPEN   Concatenate name to string of calling subroutines
!     ERRINT   Fills integer in error message
!     ERRSUB   Produces error message
!     INSTOP   Stop SEPRAN execution
!     MSH402   Place neighbour in neighbour array
!     MSH416   Determine correct order in nodenumbers
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     531  REFINE or TRANSF not yet available for this type of element
! **********************************************************************
!
!                       PSEUDO CODE
!
!     PSEUDO CODE
!     Run through all elements of this group
!      Fill node array for each element
!      Place lines, surfaces and/or volumes in neighbour array
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      character(len=260) localName
      localName = 'msh401'
      call eropen( localName )

!     --- Run through all elements of this group

      if ( ishape.eq.-9 ) then

!     --- Special line element containing all points of one side

         do j = 1 , inpelm-1

            i(1) = kmeshc( j   )
            i(2) = kmeshc( j+1 )

            call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                   kelemh )

         end do

      else if ( ishape.eq.1 ) then

!     --- inpelm = 2

         do ielem = 1, nelem

            jelem = 2 * ( ielem - 1 )

            do j = 1 , 2

                i(j) = kmeshc( jelem + j )

            end do

            call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                   kelemh )

         end do

      else if ( ishape.eq.2 ) then

!     --- inpelm = 3

         do ielem = 1, nelem

            jelem = 3 * ( ielem - 1 )

            do j = 1 , 3

                i(j) = kmeshc( jelem + j )

            end do

            call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                   kelemh )

         end do

      else if ( ishape.eq.3 ) then

!     --- inpelm = 3

         do ielem = 1, nelem

            jelem = 3 * ( ielem - 1 )

            do j = 1 , 3

                i(j) = kmeshc( jelem + j )

            end do

            if ( niedge.gt.0 ) then

               call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(1),
     +                      kelemh )

            end if

!           --- Place surface later in array ibrpnt (see MSH443) !

         end do

      else if ( ishape.eq.4 ) then

!     --- inpelm = 6

         do ielem = 1, nelem

            jelem = 6 * ( ielem - 1 )

            do j = 1 , 6

                i(j) = kmeshc( jelem + j )

            end do

            call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(5),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(6),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(1),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(6),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(4),
     +                   kelemh )
            call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(6),
     +                   kelemh )

         end do

      else if ( ishape.eq.5 ) then

!     --- inpelm = 4

         do ielem = 1, nelem

            jelem = 4 * ( ielem - 1 )

            do j = 1 , 4

                i(j) = kmeshc( jelem + j )

            end do

            if ( niedge.gt.0 ) then

               call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(1),
     +                      kelemh )

            end if

            if ( nisurf.gt. 0 .or.
     +           nisurf.eq.-1 .and. niedge.gt.0 ) then

!           --- One extra for new middlepoint

               call msh416( i(1), i(2), i(3), i(4), ie, it )

               it = it * npoint

               call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                      kelemh )

            end if

         end do

      else if ( ishape.eq.6 ) then

!     --- inpelm = 9

         do ielem = 1, nelem

             jelem = 9 * ( ielem - 1 )

             do j = 1 , 9

                 i(j) = kmeshc( jelem + j )

             end do

             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(5),
     +                    kelemh )

             call msh416( i(1), i(2), i(9), i(8), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )

             call msh416( i(2), i(3), i(4), i(9), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )

             call msh416( i(8), i(9), i(6), i(7), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )

             call msh416( i(9), i(4), i(5), i(6), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )

         end do

      else if ( ishape.eq.7 ) then

!     --- inpelm = 7

         do ielem = 1, nelem

             jelem = 7 * ( ielem - 1 )

             do j = 1 , 7

                 i(j) = kmeshc( jelem + j )

             end do

             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(5),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(1),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(2),
     +                    kelemh )

         end do

      else if ( ishape.eq.11 ) then

!     --- inpelm = 4

         do ielem = 1, nelem

             jelem = 4 * ( ielem - 1 )

             do j = 1 , 4

                 i(j) = kmeshc( jelem + j )

             end do

             if ( niedge.gt.0 ) then

                call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(3),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(4),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(4),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                       kelemh )

             end if

!            --- Place the four surfaces later in the neighbour arrays !
!                (Use routine MSH443)

         end do

      else if ( ishape.eq.12 ) then

!     --- inpelm = 10

         do ielem = 1, nelem

             jelem = 10 * ( ielem - 1 )

             do j = 1 , 10

                 i(j) = kmeshc( jelem + j )

             end do

             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(5),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(7),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(10),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(10),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(9), i(10),
     +                    kelemh )

         end do

      else if ( ishape.eq.13 ) then

!     --- inpelm = 8

         do ielem = 1, nelem

             jelem = 8 * ( ielem - 1 )

             do j = 1 , 8

                 i(j) = kmeshc( jelem + j )

             end do

             if ( niedge.gt.0 ) then

                call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(4),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(1),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(5),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(6),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(7),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(8),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(6),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(7),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(8),
     +                       kelemh )
                call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(5),
     +                       kelemh )

             end if

             if ( nisurf.gt. 0 .or.
     +            nisurf.eq.-1 .and. niedge.gt.0 ) then

!            --- Store the six faces

                call msh416( i(1), i(2), i(3), i(4), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

                call msh416( i(1), i(2), i(6), i(5), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

                call msh416( i(2), i(3), i(7), i(6), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

                call msh416( i(3), i(4), i(8), i(7), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

                call msh416( i(1), i(4), i(8), i(5), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

                call msh416( i(5), i(6), i(7), i(8), ie, it )
                it = it * npoint
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

             end if

             if ( nivolm.gt. 0 .or.
     +            nivolm.eq.-1 .and. niedge.gt.0 ) then

                ie = min( i(1), i(7) )
                it = max( i(1), i(7) )

                it = npoint * ( npoint + 1 ) + it
                call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                       kelemh )

             end if

         end do

      else if ( ishape.eq.14 ) then

!     --- inpelm = 27

         do ielem = 1, nelem

             jelem = 27 * ( ielem - 1 )

             do j = 1 , 27

                 i(j) = kmeshc( jelem + j )

             end do

!            --- determine new points in "first plane"

             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(2),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(3),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(4),
     +                    kelemh )
             call msh416( i(1), i(2), i(5), i(4), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(5),
     +                    kelemh )
             call msh416( i(2), i(3), i(6), i(5), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(5),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(6),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(7),
     +                    kelemh )
             call msh416( i(4), i(5), i(8), i(7), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(8),
     +                    kelemh )
             call msh416( i(5), i(6), i(9), i(8), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(9),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(8),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(9),
     +                    kelemh )

!            --- "second plane"

             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(10),
     +                    kelemh )
             call msh416( i(1), i(2), i(11), i(10), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(11),
     +                    kelemh )
             call msh416( i(2), i(3), i(12), i(11), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(12),
     +                    kelemh )
             call msh416( i(1), i(4), i(13), i(10), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(3), i(14),
     +                    kelemh )
             call msh416( i(3), i(6), i(15), i(12), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(13),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(4), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(5), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(6), i(15),
     +                    kelemh )
             call msh416( i(4), i(7), i(16), i(13), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(9), i(14),
     +                    kelemh )
             call msh416( i(6), i(9), i(18), i(15), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(7), i(16),
     +                    kelemh )
             call msh416( i(7), i(8), i(17), i(16), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(8), i(17),
     +                    kelemh )
             call msh416( i(8), i(9), i(18), i(17), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(9), i(18),
     +                    kelemh )

!            --- "third plane"

             call msh402( npoint, ibrp, ibrpnt, istart, i(10), i(11),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(11), i(12),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(10), i(13),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(10), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(11), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(12), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(12), i(15),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(13), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(14), i(15),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(13), i(16),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(14), i(16),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(14), i(17),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(14), i(18),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(15), i(18),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(16), i(17),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(17), i(18),
     +                    kelemh )

!            --- "fourth plane"

             call msh402( npoint, ibrp, ibrpnt, istart, i(10), i(19),
     +                    kelemh )
             call msh416( i(10), i(11), i(20), i(19), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(11), i(20),
     +                    kelemh )
             call msh416( i(11), i(12), i(21), i(20), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(12), i(21),
     +                    kelemh )
             call msh416( i(10), i(13), i(22), i(19), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(19), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(20), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(21), i(14),
     +                    kelemh )
             call msh416( i(12), i(15), i(24), i(21), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(22), i(13),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(22), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(23), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(24), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(24), i(15),
     +                    kelemh )
             call msh416( i(13), i(16), i(25), i(22), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(25), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(26), i(14),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(27), i(14),
     +                    kelemh )
             call msh416( i(15), i(18), i(27), i(24), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(25), i(16),
     +                    kelemh )
             call msh416( i(16), i(17), i(26), i(25), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(26), i(17),
     +                    kelemh )
             call msh416( i(17), i(18), i(27), i(26), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(27), i(18),
     +                    kelemh )

!            --- "fifth plane"

             call msh402( npoint, ibrp, ibrpnt, istart, i(19), i(20),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(20), i(21),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(19), i(22),
     +                    kelemh )
             call msh416( i(19), i(20), i(23), i(22), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(20), i(23),
     +                    kelemh )
             call msh416( i(20), i(21), i(24), i(23), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(21), i(24),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(22), i(23),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(23), i(24),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(22), i(25),
     +                    kelemh )
             call msh416( i(22), i(23), i(26), i(25), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(23), i(26),
     +                    kelemh )
             call msh416( i(23), i(24), i(27), i(26), ie, it )
             it = it * npoint
             call msh402( npoint, ibrp, ibrpnt, istart, ie, it,
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(24), i(27),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(25), i(26),
     +                    kelemh )
             call msh402( npoint, ibrp, ibrpnt, istart, i(26), i(27),
     +                    kelemh )

         end do

      else if ( ishape.eq.50 ) then

!     --- inpelm = 50

         do ielem = 1, nelem

            jelem = 4 * ( ielem - 1 )

            do j = 1, 4

                i(j) = kmeshc( jelem + j )

            end do

            if ( niedge.gt.0 ) then

               call msh402( npoint, ibrp, ibrpnt, istart, i(1), i(3),
     +                      kelemh )
               call msh402( npoint, ibrp, ibrpnt, istart, i(2), i(4),
     +                      kelemh )

            end if

!           --- Place surface later in array ibrpnt (see MSH443) !

         end do

      else

!     --- not yet available

         call errint ( ishape, 1 )
         call errsub ( 531, 1, 0, 0 )
         go to 1000

      end if

1000  call erclos('msh401')

      end

