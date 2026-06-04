      subroutine msho12 ( coor, kstapl, kstap, i1, i2, iex1, iex2,
     +                    angle1, angle2 )
! ======================================================================
!
!        programmer    Niek Praagman
!        version  3.3  date 14-02-2010 Higher accuracy
!        version  3.2  date 26-01-2005 Layout
!        version  3.1  date 21-02-2003 Layout
!        version  3.0  date 14-02-1994 New norms
!        version  2.0  date 01-02-1990 only real neighbours !
!        version  1.0  date 14-04-1989
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
!     Determine for line i1 - i2 the other lines where i1 resp i2 are
!     nodes off and compute the angles between these lines and the original
!     line
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     2d
!     neighbour
!     angle
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use mshdummymethods
      
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!      include 'SPcommon/cmcdpi'

! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kstap, kstapl(2,kstap), i1, i2, iex1, iex2
      double precision coor(*), angle1, angle2

!     angle1      o  angle between basis line and iex1 line
!     angle2      o  angle between basis line and iex2 line
!     coor        i  coordinate array
!     i1          i  first node segment to be considered
!     i2          i  second node segment to be considered
!     iex1        o  node of line adjacent to i1
!     iex2        o  node of line adjacent to i2
!     kstap       i  number of line-elements in kstapl
!     kstapl      i  boundary array
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer ii1, ii2, i
      double precision angle, surf
      logical debug

!     angle          help variable to determine angle
!     debug          If true debug statements are carried out otherwise
!                    they are not
!     i              loop variable
!     ii1            node number
!     ii2            node number
!     surf           help variable for determination of surface
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     ERRINT         Put integer in error message
!     ERRSUB         Error messages
!     MSHO11         compute angle
!     MSHO19         Compute area of triangle
!     PRININ1        print 2d integer array
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!    2433   Internal error in triangle
! **********************************************************************
!
!                       PSEUDO CODE
!
!     trivial
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'msho12' )

      debug = .false.

      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from msho12'

      end if

      iex1 = 0
      iex2 = 0

!     --- Initialize angle values

      angle1 = -1.2d0
      angle2 = -1.2d0

      do i = 1 , kstap

         ii1 = kstapl(1,i)
         ii2 = kstapl(2,i)

         if ( i1==ii2 ) then

!        --- Compute surface of area

              call msho19( coor, ii1, i1, i2, surf )

              if ( surf<=0 ) then
                   angle = -1.1d0
              else

!             --- Compute angle ii1-i1-i2

                   call msho11( ii1, i1, i2, coor, angle )
              end if
              if ( angle>angle1 ) then
                   iex1   = ii1
                   angle1 = angle
              end if
         end if

         if ( i2==ii1 ) then

!        --- Compute surface of area

              call msho19( coor, i1, i2, ii2, surf )
              if ( surf<=0 ) then
                   angle = -1.1d0
              else

!             --- Compute angle i1-i2-ii2

                   call msho11( i1, i2, ii2, coor, angle )
              end if
              if ( angle>angle2 ) then
                   iex2   = ii2
                   angle2 = angle
              end if
         end if
      end do

      call erclos ( 'msho12' )

      if ( debug ) write(irefwr,*) 'End msho12'

      if ( iex1==0 .or. iex2==0 ) then

!     --- Problem in msho12

         call errint ( i1, 1 )
         call errint ( i2, 2 )
         call errint ( iex1, 3 )
         call errint ( iex2, 4 )
         call errsub ( 2433, 4, 0, 0 )

      end if

      end
