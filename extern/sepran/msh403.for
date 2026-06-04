      subroutine msh403( npoint, i1, i3, ibrpnt, istart, i2 )
! ======================================================================
!
!         programmer   niek praagman
!         version 2.1  date  22-07-1998 Layout
!         version 2.0  date  14-02-1995 npoint added; other storage of
!                                       neighbour nodes
!         version 1.3  date  08-04-1994 initialize i2
!         version 1.2  date  03-03-1993 adjustment layout and error 949
!         version 1.1  date  22-08-1991 Zdenek: no dimension, cdc, ce
!         version 1.0  date  01-02-1986
!
!
!
!   copyright (c) 1986-1998  "Ingenieursbureau SEPRA"
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
!     Determine the nodal point number i2 of a new point between the nodes
!     i1 and i3 as stored in array ibrpnt
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
      integer i1, i2, i3, ibrpnt(*), istart(*), npoint

!     i1         i    first node number to be considered
!     i2         o    node number of intermediate point
!     i3         i    last node number of line to be considered
!     ibrpnt     i    array containing the neighbours of the
!                     nodes sequentially
!     istart     i    array containing the startaddresses for
!                     each original nodal point
!     npoint     i    number of nodes in original mesh
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer ihelp, ihigh, ilow, jstart, jstend

!     ihelp    place to store node number temporarily
!     ihigh    highest node number
!     ilow     lowest node number
!     jstart   start position neighbour points
!     jstend   last position neighbour points
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS  Deconcatenate name from string of calling subroutines
!     EROPEN  Concatenate name to string of calling subroutines
!     ERRSUB  Produce error message
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     949    No intermediate neighbour found
! **********************************************************************
!
!                       PSEUDO CODE
!
!     Determine position taking the value of npoint into account
!
!     Check new neighbour number
!
!     If no point is found, error !
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
      localName = 'msh403'
      call eropen( localName )

      i2 = 0

!     --- determine place

      if ( i1.gt.i3 ) then

!     --- Check for npoint

         if ( i1.gt.npoint ) then

            ilow   =  i1
            ihigh  =  i3

         else

            ilow   =  i3
            ihigh  =  i1

         end if

      else

!     --- Check for npoint

         if ( i3.gt.npoint ) then

            ilow   =  i3
            ihigh  =  i1

         else

            ilow   =   i1
            ihigh  =   i3

         end if

      end if

      jstart      =   istart( ihigh     )
      jstend      =   istart( ihigh + 1 ) - 1

100   ihelp   =   ibrpnt(jstart)

      if ( ihelp.eq.ilow ) then

!     --- Point found

         i2 = ibrpnt( jstart + 1 )

         goto 1000

      end if

      jstart  =  jstart + 2

      if ( jstart.le.jstend ) goto 100

!     --- Error, no point number found

      call errsub ( 949, 0, 0, 0 )

1000  call erclos('msh403')

      end

