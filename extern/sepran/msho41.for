      subroutine msho41( cmin, cmax, dist, maxratio )
! ======================================================================
!
!        programmer    niek praagman
!        version  1.0  date 10-11-2009
!
!   copyright (c) 2009-2009  "Ingenieursbureau SEPRA"
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
!     Compute the max allowed coarseness cmax in point B that is in
!     accordance with given coarseness cmin in a point A. The
!     Euclidian distance from A to B is dist and the max allowed ratio
!     between two linepieces is maxratio.
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
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision cmax, cmin, dist, maxratio

!     cmax        i/o   i: given coarseness in point B
!                       o: max allowed coarseness in point B
!     cmin         i    coarseness point A
!     dist         i    Euclidian distance A to B
!     maxratio     i    max allowed multiplication in two neighbouring
!                       elements (i.e. linepieces)
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision afst, som

!     afst         local distance between sequential nodes
!     som          temporary distance
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
!     Trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Check possibilities:

      if ( dist < cmin ) then
!     --- Length dist too small:
         cmax = cmin
      else if ( dist>cmin ) then
!     --- Check whether there is enough space from small to large:
         afst = 0.65 * cmin
         som  = afst

         do while ( som<dist )

            afst = afst * maxratio
            som  = som + afst

         end do

!        --- Set value of cmax to last distance:

         if ( afst < cmax ) cmax = afst

      end if

      end
