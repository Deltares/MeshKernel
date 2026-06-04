      subroutine msh402( npoint, ibrp, ibrpnt, istart, i1, i2, kelemh )
! ======================================================================
!
!         programmer   Niek Praagman
!         version 2.5  date 14-02-1995 use alternative of neighbours
!                                      and npoint
!         version 2.4  date 03-03-1993 layout adjustments, new norms
!         version 2.3  date 18-08-1992 iabs==>abs
!         version 2.2  date 22-08-1991 Zdenek: no dimension, cdc, ce
!         version 2.1  date 10-09-1990 variable length for IBRPNT
!         version 2.0  date 01-03-1989
!
!
!
!   copyright (c) 1986-1995  "Ingenieursbureau SEPRA"
!   permission to copy or distribute this software or documentation
!   in hard or soft copy granted only by written license obtained
!   from "Ingenieursbureau SEPRA".
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
!     Fill the nodal point numbers in the right positions of array ibrpnt
!
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     neighbour
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      use mshdummymethods
      
      implicit none
      integer npoint, ibrpnt(*), i1, i2, ibrp, istart(*), kelemh(*)

!     i1          i    first node number to be considered
!     i2          i    second node number to be considred
!     ibrp        i    number of places available in IBRPNT for each
!                      point
!     ibrpnt     i,o   array containing sequentially the neighbour-
!                      points
!     istart      i    array of startposition for each nodal point
!                      in ibrpnt
!     kelemh      o    array to store extra neighbours if a point has
!                      more than ibrp neighbours
!     npoint      i    number of nodes before operation
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer idif, ih, ihigh, il, ilow, j, ja, jstart, jstend

!     idif    difference
!     ih      number of neighbour
!     ihigh   highest node number of i1,i2
!     il      loop variable
!     ilow    lowest node number of i1,i2
!     j       loop variable
!     ja      local indicator which storage system is used
!     jstart  startposition in ibrpnt
!     jstend  endposition in ibrpnt
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS  Deconcatenate name from string of calling subroutines
!     EROPEN  Concatenate name to string of calling subroutines
!     ERRSUB  Produce error message
!     INSTOP  Stop SEPRAN execution
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     1274    No place left for neighbour
! **********************************************************************
!
!                       PSEUDO CODE
!
!     Check whether special situation (i1 < 0)
!
!     If special situation ja = 1 else ja = 0
!
!     Determine lowest and highest node number if both are smaller
!     than npoint else change roles
!
!     Store : ilow is neighbour of ihigh
!
!     Determine places in array with neighbours
!
!     If ja=1 than start situation else rearranged situation
!
!     Check whether new point
!
!        If new than place number in array of neighbours
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
      character(len=260) localName
      localName = 'msh402'
      call eropen( localName )

      if ( i1.lt.0 ) then

           ja =   0
           i1 = -i1

      else

           ja =   1

      end if

!     --- Determine ilow and ihigh

      if ( i1.gt.i2 ) then

!     --- Check for npoint

           if ( i1.gt.npoint ) then

                ilow   = i1
                ihigh  = i2

           else

                ilow   = i2
                ihigh  = i1

           end if

      else

!     --- Check for npoint

           if ( i2.gt.npoint ) then

                ilow   = i2
                ihigh  = i1

           else

                ilow   =  i1
                ihigh  =  i2

           end if

      end if

!     --- ilow is neighbour of ihigh

      if ( ja.eq.1 ) then

           jstart = ibrp * ( ihigh - 1 ) + 1
           jstend = jstart + ibrp  - 1

      else

           jstart = istart( ihigh     )
           jstend = istart( ihigh + 1 ) - 1

      end if

!     --- Run through available positions

100   ih = ibrpnt( jstart )

      if ( ih.eq.0 ) then

           ibrpnt( jstart ) = ilow
           istart( ihigh  ) = istart( ihigh ) + ja

!          --- Ready

           goto 1000

      else if ( ih.eq.ilow ) then

!     --- Line already stored

           goto 1000

      end if

      jstart = jstart + 1

      if ( jstart.gt.jstend ) then

           if ( ja.eq.0 ) then

              call errsub ( 1274, 0, 0, 0 )

           end if

           j = kelemh(1)

!          --- Check whether this line is already in kelemh

           do il = 1 , j

                idif = abs ( kelemh(2*il  ) - ihigh ) +
     +                 abs ( kelemh(2*il+1) - ilow  )

                if ( idif.eq.0 ) goto 1000

           end do

!          --- Place in kelemh

           kelemh( 2 * j + 2 )  = ihigh
           kelemh( 2 * j + 3 )  = ilow
           kelemh( 1         )  = j + 1
           istart( ihigh     )  = istart( ihigh ) + 1

!          --- Ready

           goto 1000

      end if

      goto 100

1000  call erclos('msh402')

      end

