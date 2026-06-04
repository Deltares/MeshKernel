      subroutine msho22( a, b, c, icase )
! ======================================================================
!
!        programmer    niek praagman
!        version  2.2  date 01-02-2010 Remove jelem
!        version  2.1  date 11-01-2005 Adjust value jelem
!        version  2.0  date 21-02-1994 New norms
!        version  1.0  date 17-04-1989
!
!   copyright (c) 1989-2010  "Ingenieursbureau SEPRA"
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
!     Compute all angles of triangle with sides a,b,c of element
!     and check whether these angles are too large i.e. > 90 degrees
!     Indicate which side it concerns.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     angle
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
      double precision a, b, c
      integer icase

!     a      i   length of first side of triangle
!     b      i   length of second side of triangle
!     c      i   length of third side of triangle
!     icase  o   indicator if (>0) an angle and which angle is too large
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      double precision angle1, angle2, angle3, aa, bb, cc, q

!     aa      length of side one squared
!     angle1  first angle
!     angle2  second angle
!     angle3  third angle
!     bb      length of side two squared
!     cc      length of side three squared
!     q       constant
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
!     --- Square the three sides

      aa = a * a
      bb = b * b
      cc = c * c

!     --- Use cosine rule to compute the cosine of each of the three angles

      angle1 = ( cc + aa - bb ) / ( 2 * c * a )
      angle2 = ( aa + bb - cc ) / ( 2 * a * b )
      angle3 = ( bb + cc - aa ) / ( 2 * b * c )

      q      = -0.01d0

!     --- Suppose that all angles are < 90 degrees:
      icase  = 0

!     --- Next check whether this is correct:
      if ( angle1.lt.q ) then
           icase = 1
      else if ( angle2.lt.q ) then
           icase = 2
      else if ( angle3.lt.q ) then
           icase = 3
      end if

      end
