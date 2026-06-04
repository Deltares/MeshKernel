      subroutine msho05 ( dist, npoint, dismin, dismax )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 11-02-2005 Update min, max determination
!        version  2.0  date 03-02-1994 New norms
!        version  1.1  date 04-01-1990 Quadratic case added
!        version  1.1  date 12-04-1989
!
!   copyright (c) 1989-2005  "Ingenieursbureau SEPRA"
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
!     Subroutine to determine the extreme coarsenesses of the points on
!     the boundary of the area given
! **********************************************************************
!
!                       KEYWORDS
!
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
!      include 'SPcommon/cmcdpr'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      double precision dist(*), dismin, dismax
      integer npoint

!     dismax    o    largest distance
!     dismin    o    smallest distance
!     dist      i    array containing the representative distances
!                    for the surface-points
!     npoint    i    number of points in boundary elements
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision dis, eps
      integer i

!     dis      coarseness of node considered
!     eps      accuracy
!     i        loop variable
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
!     Compute extreme coarseness values
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      eps    = 10 * epsmac

      dismin =  rinfin
      dismax = -rinfin

      do i = 1, npoint

         dis = dist(i)

         if ( dis>eps ) then

              dismin = min( dismin, dis)
              dismax = max( dismax, dis)

         end if

      end do

      end
