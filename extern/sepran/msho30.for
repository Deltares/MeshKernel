      subroutine msho30 ( kelem, nelem, ielem, kstapl, kstap, npoint,
     +                    itri , nelmfix)
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.3  date 09-04-2009 use nelmfix i.s.o. nelemi
!        version  3.2  date 18-03-2005 Extra parameter nelemi
!        version  3.1  date 15-02-2005 Update
!        version  3.0  date 11-01-2005 Remove nipnt
!        version  2.1  date 14-08-2002 Consider Plaxis boundaries
!        version  2.0  date 28-01-1994 Complete new version
!        version  1.1  date 19-02-1993 New layout
!        version  1.0  date 13-06-1989
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
!     Subroutine to refill arrays KSTAPL and ITRI after a Repositioning
!     error
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
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kelem(3,*), nelem, ielem, kstapl(2,*), kstap,
     +        npoint, itri(*), nelmfix

!     ielem          i    element to be cancelled
!     itri           o    indicator whether point i is in boundary or not
!     kelem         i/o   New KMESH part c with respect to the surface elements
!                         All elements have the same number of nodes
!     kstap          o    number of line segments in KSTAP
!     kstapl         o    new array of boundary segments
!     nelem         i/o   Number of elements
!     nelmfix        i    Number of elements that is fixed
!     npoint         i    Number of nodal points
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, ia, ib, ic, nextra, i1, i2, i3, j1, j2, jp, jtal,
     +        iaan, is

!     i              counting variable
!     i1             help node number
!     i2             help node number
!     i3             help node number
!     ia             node number of triangle to be considered
!     iaan           indicator whether segment has been cancelled or not
!     ib             node number of triangle to be considered
!     ic             node number of triangle to be considered
!     is             counting variable
!     j1             segment node
!     j2             segment node
!     jp             indicator whether element should be reconsidered
!                    or not
!     jtal           temporary number of segments in kstapl
!     nextra         number of elements which are not cancelled
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
! **********************************************************************
!
!                       I/O
!
!     none
! **********************************************************************
!
!                       ERROR MESSAGES
!
! **********************************************************************
!
!                       PSEUDO CODE
!
!     determine the three node numbers of the triangle to be cancelled
!     run through all elements
!         if element has at least one node in common with triangle
!         cancel element and fill the three segments in kstapl
!     refill itri
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Determine the three node numbers of the triangle to be skipped

      ia = kelem( 1, ielem )
      ib = kelem( 2, ielem )
      ic = kelem( 3, ielem )

!     --- Set parameters KSTAP and NEXTRA

      kstap  = 0
      nextra = nelmfix

!     --- Run through all elements that have at least one point in common
!         with the triangle

      do i = nelmfix + 1 , nelem

          i1  =  kelem( 1, i )
          i2  =  kelem( 2, i )
          i3  =  kelem( 3, i )

          jp  = 0

          if ( i1==ia .or. i1==ib .or. i1==ic ) jp = 1
          if ( i2==ia .or. i2==ib .or. i2==ic ) jp = 1
          if ( i3==ia .or. i3==ib .or. i3==ic ) jp = 1

          if ( jp==1 ) then

!         --- add three line segments to kstapl
!             start with i1 - i2

             if ( i1<=-1 .and. i2<=-1 ) then

                kstap = kstap + 1
                kstapl( 1, kstap ) = i1
                kstapl( 2, kstap ) = i2

             else

!            --- Check  whether i1 - i2 already belongs to kstapl

                jtal = 0
                iaan = 0

                do is = 1, kstap
                   j1 = kstapl( 1, is )
                   j2 = kstapl( 2, is )
                   if ( j1==i2 .and. j2==i1 ) then

!                  --- do nothing, i.e. skip segment

                      iaan = 1
                   else
                      jtal = jtal + 1
                      kstapl( 1, jtal ) = j1
                      kstapl( 2, jtal ) = j2
                   end if
                end do
                if ( iaan==0 ) then

!               --- add segment

                   jtal = jtal + 1
                   kstapl( 1, jtal ) = i1
                   kstapl( 2, jtal ) = i2
                end if
                kstap = jtal
             end if

!            --- next check whether i2 - i3 already belongs to kstapl

             jtal = 0
             iaan = 0

             do is = 1, kstap
                j1 = kstapl( 1, is )
                j2 = kstapl( 2, is )
                if ( j1==i3 .and. j2==i2 ) then

!               --- do nothing, i.e. skip segment

                   iaan = 1
                else
                   jtal = jtal + 1
                   kstapl( 1, jtal ) = j1
                   kstapl( 2, jtal ) = j2
                end if
             end do
             if ( iaan==0 ) then

!            --- add segment

                jtal = jtal + 1
                kstapl( 1, jtal ) = i2
                kstapl( 2, jtal ) = i3
             end if
             kstap = jtal

!            --- finally check whether i3 - i1 already belongs to kstapl

             jtal = 0
             iaan = 0

             do is = 1, kstap
                j1 = kstapl( 1, is )
                j2 = kstapl( 2, is )
                if ( j1==i1 .and. j2==i3 ) then

!               --- do nothing, i.e. skip segment

                   iaan = 1
                else
                   jtal = jtal + 1
                   kstapl( 1, jtal ) = j1
                   kstapl( 2, jtal ) = j2
                end if
             end do
             if ( iaan==0 ) then

!            --- add segment

                jtal = jtal + 1
                kstapl( 1, jtal ) = i3
                kstapl( 2, jtal ) = i1
             end if
             kstap = jtal

          else

!         --- triangle is outside interesting area

             nextra = nextra + 1
             kelem( 1, nextra ) = i1
             kelem( 2, nextra ) = i2
             kelem( 3, nextra ) = i3
          end if
      end do
      nelem = nextra

!     --- Refill array ITRI

      do i = 1 , npoint
          itri(i) = 0
      end do

      do i = 1 , kstap

          i1 = kstapl( 1, i )
          i2 = kstapl( 2, i )

          itri(i1) = itri(i1) + 1
          itri(i2) = itri(i2) + 1

      end do

      end
