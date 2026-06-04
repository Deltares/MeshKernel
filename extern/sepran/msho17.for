      subroutine msho17( ibuurp, iwork, ih, ij )

! =======================================================================
!
!
!         programmer   niek praagman
!
!         version 2.0    date  21-02-1994 New norms
!         version 1.0    date  14-04-1989
!
!
!
!   copyright (c) 1989-1994  "ingenieursbureau sepra"
!   permission to copy or distribute this software or documentation
!   in hard or soft copy granted only by written license obtained
!   from "ingenieursbureau sepra".
!   all rights reserved. no part of this publication may be reproduced,
!   stored in a retrieval system (e.g. , in memory, disk or core)
!   or be transmitted by any means, electronic, mechanical, photocopy,
!   recording, or otherwise, without written permission from the
!   publisher.
!
! ***********************************************************************
!
!     DESCRIPTION
!
!     Fill nodal point number ij in right position of array ibuurp
!
! ***********************************************************************
!
!     KEYWORDS
!
!     mesh_generation
!     neighbour
! ***********************************************************************
!
!     INPUT / OUTPUT PARAMETERS
!
      implicit none
      integer ibuurp(*),iwork(*),ih,ij

!     ibuurp  i,o    array containing the nodal point numbers of
!                    the neighbours
!     ih       i     nodal point number to be considered
!     ij       i     nodal point number of neighbour
!     iwork    i     pointer array
!
! ***********************************************************************
!
!     COMMON BLOCKS
!
! ***********************************************************************
!
!     LOCAL PARAMETERS
!
      integer jstart
!
!     jstart    startaddress for present point ih in array IBUURP
!
! ***********************************************************************
!
!     SUBROUTINES CALLED
!
! ***********************************************************************
!
!     ERROR MESSAGES
!
! ***********************************************************************
!
!     PSEUDO CODE
!
!     trivial
!
! =======================================================================

      jstart=0

      if ( ih.ne.1 ) jstart = iwork(ih-1)

100   jstart = jstart + 1

!     If point is already placed then return

      if ( ibuurp ( jstart ).eq.ij )  goto 200

!     If not yet last neighbour, then check next

      if ( ibuurp ( jstart ).ne. 0 )  goto 100

!     Position is empty, place new neighbour

      ibuurp ( jstart ) = ij

200   continue

      end
