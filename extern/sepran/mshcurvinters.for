      subroutine mshcurvinters ( ncurvs, iwork, xbox, icurvs, curves )
! ======================================================================
!
!        programmer    Guus Segal
!        version  1.0  date 28-11-2005
!
!   copyright (c) 2005-2005  "Ingenieursbureau SEPRA"
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
!     Give an error if two curves intersect
!     At this moment 2d only
! **********************************************************************
!
!                       KEYWORDS
!
!     intersection
!     curve
!     2d
! **********************************************************************
!
!                       MODULES USED
!
      use mshconstants
      use msherror
      
      implicit none
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer ncurvs, iwork(*), icurvs(*)
      double precision xbox(2,ncurvs,2), curves(2,*)

!     curves         i    array containing the coordinates of the curves
!     icurvs         i    array containing the number of points in the curves
!                         accumulated
!     iwork          i    integer work array of length ncurvs
!                         iwork(i) = -1: curve is empty
!                         iwork(i) = 0: curve is single
!                         iwork(i) = 1: curve is compound
!     ncurvs         i    Number of curves in mesh
!     xbox           i    Is used to store the minimum and maximum coordinates
!                         of all curves
!                         xbox(i,j,1)  minimum i^th coordinate of curve j
!                         xbox(i,j,2)  maximum i^th coordinate of curve j
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer icurnr, jcurnr
      logical debug

!     debug          If true debug statements are carried out otherwise
!                    they are not
!     icurnr         Curve sequence number
!     jcurnr         Curve sequence number
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS         Resets old name of previous subroutine of higher level
!     EROPEN         Produces concatenated name of local subroutine
!     MSHCURVINTERS1 Give an error if 2 curves intersect
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
!     For all single curves icurnr do
!        for all single curve with sequence number>icurnr do
!           If ( possible intersection ) then
!              if ( subparts of curve intersect ) then
!                 give error message
! **********************************************************************
!
!                       DATA STATEMENTS
!
! ======================================================================
!
      call eropen ( 'mshcurvinters' )
      debug = .false.
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'Debug information from mshcurvinters'

      end if
      if ( ierror/=0 ) go to 1000

!     --- For all single curves icurnr do

      do icurnr = 1, ncurvs

         if ( iwork(icurnr)==0 ) then

!        --- single non-empty curve found
!            for all single curve with sequence number>icurnr do

            do jcurnr = icurnr+1, ncurvs

               if ( iwork(jcurnr)==0 ) then

!              --- single non-empty curve found
!                  If ( possible intersection ) then
!                  xmini<xmaxj, xminj<xmaxi

                  if ( xbox(1,icurnr,1)>xbox(1,jcurnr,2) ) go to 200
                  if ( xbox(1,jcurnr,1)>xbox(1,icurnr,2) ) go to 200

!                 --- ymini<ymaxj, yminj<ymaxi

                  if ( xbox(2,icurnr,1)>xbox(2,jcurnr,2) ) go to 200
                  if ( xbox(2,jcurnr,1)>xbox(2,icurnr,2) ) go to 200

!                 --- if ( subparts of curve intersect ) then
!                        give error message

                  call mshcurvinters1 ( icurnr, jcurnr, icurvs, curves )

               end if  ! ( iwork(jcurnr)==0 )

200         end do  ! jcurnr = icurnr+1, ncurvs

         end if  ! ( iwork(i)==0 )

      end do  ! icurnr = 1, ncurvs

1000  call erclos ( 'mshcurvinters' )
      if ( debug ) then

!     --- Debug information

         write(irefwr,*) 'End mshcurvinters'

      end if

      end

