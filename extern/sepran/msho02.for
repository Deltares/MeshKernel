      subroutine msho02( istart, ibuur, i, j)
! ======================================================================
!
!         programmer   Niek Praagman
!         version 2.1  date  17-08-1997  Check coincidence of double
!                                        Plaxis lines
!         version 2.0  date  03-02-1994 New norms
!         version 1.0  date  12-04-1989
!
!
!
!   copyright (c) 1989-1997  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard or soft copy granted only by written license obtained
!   from "ingenieursbureau SEPRA".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system (e.g. , in memory, disk or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
!
! **********************************************************************
!
!                       DESCRIPTION
!
!     Subroutine to fill array of neighbour points
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     neighbour
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      use mshdummymethods
      
      implicit none
      integer istart(*), ibuur(*), i, j

!     i         i   first node of segment to be considered
!     ibuur     o   neighbour array
!     istart    i   array with startaddresses for neighbours of
!                   each point
!     j         i   second node of line to be considered
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer ist, kn, ien, ib

!     ib         node number of neighbour
!     ien        end address of neighbours in array IBUUR
!     ist        starting address of neighbours in array IBUUR
!     kn         loop variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS     Deconcatenate name from string of calling subroutines
!     EROPEN     Concatenate name to string of calling subroutines
!     ERRSUB     Submit error message
!     INSTOP     Stop execution of SEPRAN
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     906    Double line in the boundary
! **********************************************************************
!
!                       PSEUDO CODE
!
!     find starting position for neighbours of point i
!
!     run through all points; if node=0 then place new point, else
!     if node = j give error message
!
! **********************************************************************
!
!                       MODULES USED
!
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen( 'msho02' )

!     --- Consider line i - j

      if ( i.eq.1 ) then
         ist = 1
      else
         ist = istart(i-1) + 1
      end if

      ien = istart(i)

      do kn = ist, ien

         ib = ibuur(kn)

         if ( ib.eq.0 ) then

            ibuur(kn) = j

!           --- Ready, point has been placed

            goto 1000

         else if ( ib.eq.j ) then

!        --- Double line found, do nothing more

            go to 1000

         end if

      end do

1000  call erclos( 'msho02' )

      end
