      subroutine msho39( cmax, cmin, dist, maxratio, iallow )
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
!     Subroutine to check whether the given coarseness in two points is
!     realistic with respect to their Euclidian distance
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
      use mshconstants
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer    iallow
      double precision cmax, cmin, dist, maxratio

!     cmax       i,o  largest coarsenes allowed at "end" of line
!     cmin        i   smallest coarseness at begin of line
!     dist        i   Euclidian distance between start and endpoint curve
!     maxratio    i   max ratio allowed between two successive elements
!     iallow      o   indicator whether present coarsenessesare
!                     0: not allowed, hence cmax is adjusted
!                     1: allowed, no changes in cmax
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision afst, som

!     afst         distance between sequential nodes
!     som          helpvariable addedlength of all elements in loop
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
!       Trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
!     --- Check possibilities:

      iallow = 1

      if ( dist < cmin ) then
!     --- Length dist too small:
         iallow = 0
!     --- Adjust cmax:
         cmax = dist
      else if ( dist>cmin ) then
!     --- Check whether there is enough space from small to large:
         afst = 0.6 * cmin
         som  = afst

         do while ( som<dist )

            afst = afst * maxratio
            som  = som + afst

         end do

!        --- Check whether the last element causes problems:

         if ( afst<0.5 * cmax ) then
!        --- Value of cmax has to be adjusted:
            iallow=0
            cmax = afst
         end if

      end if

      end
