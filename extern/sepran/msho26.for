      subroutine msho26( kelem, nelem, i2, i3, jelem, i4 )

! =======================================================================
!
!
!         programmer   niek praagman
!
!         version 3.0  date  02-02-2010 Redesign
!         version 2.0  date  22-02-1994 New norms
!         version 1.0  date  17-04-1989
!
!
!
!   copyright (c) 1989-2010  "ingenieursbureau sepra"
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
!     Determine triangle where i2 - i3 is side of
!
! ***********************************************************************
!
!     KEYWORDS
!
!     mesh_generation
!     triangle
!
! ***********************************************************************
!
!     INPUT / OUTPUT  PARAMETERS
!
      integer kelem(*), nelem, i2, i3, jelem, i4

!     i2       i     node number
!     i3       i     node number
!     i4       o     extra node to form triangle
!     jelem    o     number of element found
!     kelem    i     element array
!     nelem    i     number of elements in kelem
!
! ***********************************************************************
!
!     COMMON BLOCKS
!
! ***********************************************************************
!
!     LOCAL PARAMETERS
!
      integer ielem, j, j1, j2, j3
!
!     j         loop variable
!     j1        local element node
!     j2        local element node
!     j3        local element node
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
!     METHOD
!
! =======================================================================

      ielem = jelem

      jelem = 0

      do 200 j = 1, nelem

          if ( j /= ielem ) then
!         --- Find triangle j with the same side i2 - i3 that is not equal to
!             element ielem:

             j1 = kelem(3*j-2)
             j2 = kelem(3*j-1)
             j3 = kelem(3*j  )

             if ( j1.eq.i3 .and. j2.eq.i2 ) then
                  jelem = j
                  i4    = j3
                  goto 500
             else if ( j2.eq.i3 .and. j3.eq.i2 ) then
                  jelem = j
                  i4    = j1
                  goto 500
             else if ( j3.eq.i3 .and. j1.eq.i2 ) then
                  jelem = j
                  i4    = j2
                  goto 500
             end if

          end if

200   continue

      i4 = 0

500   continue

      end
