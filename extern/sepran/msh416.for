      subroutine msh416( i1, i2, i3, i4, ie, it )
! =======================================================================
!
!
!         programmer   niek praagman
!         version 1.3  date  06-04-1995 ie is min, it is max
!         version 1.2  date  04-02-1993 New norms
!         version 1.1  date  23-02-1989
!
!
!
!   copyright (c) 1989 - 1995  "ingenieursbureau sepra"
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
!     Determine begin and endpoint of (main-)diagonal in quadrilateral
!
! ***********************************************************************
!
!     KEYWORDS
!
!     quadrilateral
!     diagonal
!
! ***********************************************************************
!
!     INPUT / OUTPUT PARAMETERS
!
      implicit none

      integer i1, i2, i3, i4, ie, it

!     i1      first node of quadrilateral
!     i2      second
!     i3      third
!     i4      fourth
!     ie      first node of chosen diagonal
!     it      second node of chosen diagonal
!
! ***********************************************************************
!
!     COMMON BLOCKS
!
! ***********************************************************************
!
!     LOCAL PARAMETERS
!
      integer im

!     im      smallest node value of quadrilateral
!
! ***********************************************************************
!
!     SUBROUTINES CALLED
!
!
! ***********************************************************************
!
!     ERROR MESSAGES
!
! ***********************************************************************
!
!     METHOD
!
!     PSEUDO CODE
!
! =======================================================================

      im = min (i1,i2,i3,i4)

      if ( im.eq.i1 .or. im.eq.i3 ) then

           ie = min( i1, i3 )
           it = max( i1, i3 )

      else

           ie = min( i2, i4 )
           it = max( i2, i4 )

      end if

      end
