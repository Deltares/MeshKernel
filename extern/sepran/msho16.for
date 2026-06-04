      subroutine msho16 ( kelem, nelem, npoint, nipnt, iwork,
     +                    ibuurp, leng )
! ======================================================================
!
!        programmer    niek praagman
!        version  3.0  date 14-02-2005 Update
!        version  2.0  date 21-02-1994 New norms
!        version  1.0  date 14-04-1989
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
!     Fill array iwork which contains the number of neighbour points for
!     each point accumulated. Fill array ibuurp which contains the nodal
!     point numbers of the neighbours for the case of triangles.
! **********************************************************************
!
!                       KEYWORDS
!
!     mesh_generation
!     neighbour
! **********************************************************************
!
!                       MODULES USED
!
      use mshdummymethods
       
      implicit none
! **********************************************************************
!
!                       COMMON BLOCKS
!
!     None
! **********************************************************************
!
!                       INPUT / OUTPUT PARAMETERS
!
      integer kelem(*), nelem, npoint, nipnt, iwork(*), ibuurp(*),
     +        leng

!     ibuurp   o    neighbour nodal points
!     iwork    o    array containing pointers for array ibuurp
!     kelem    i    element array
!     leng     o    needed length array ibuurp
!              i    declared length
!     nelem    i    number of elements in the mesh
!     nipnt    i    number of boundary points
!     npoint   i    number of points in the mesh
! **********************************************************************
!
!                       LOCAL PARAMETERS
!
      integer i, isum, is, i1, i2, i3

!     i        loop variable
!     i1       node number
!     i2       node number
!     i3       node number
!     is       position variable
!     isum     help variable
! **********************************************************************
!
!                       SUBROUTINES CALLED
!
!     ERCLOS   deconcatenate name from calling string
!     EROPEN   concatenate name to calling string
!     ERRSUB   submit error message
!     INSTOP   terminate execution SEPRAN
!     MSHO17   fill neighbours in in neighbour array
! **********************************************************************
!
!                       I/O
!
! **********************************************************************
!
!                       ERROR MESSAGES
!
!     903   Declared length of array IBUURP is not enough
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
      call eropen( 'msho16' )

!     --- Clear array iwork

      do i = 1, npoint
         iwork(i)=0
      end do

!     --- Determine number of triangles point is related to

      do i= 1, nelem * 3
         iwork(kelem(i)) = iwork(kelem(i))+1
      end do

!     --- Add 1 to all boundary points

      do i = 1, nipnt
         iwork(i) = iwork(i) + 1
      end do

!     --- Rearrange and accumulate

      isum = 0
      do i = 1, npoint
         iwork(i) = iwork(i)+isum
         isum     = iwork(i)
      end do
      if ( isum>leng ) then

!     --- Declared length too small

           call errsub ( 903, 0, 0, 0 )
      end if

!     --- Clear ibuurp

      do i = 1, isum
         ibuurp(i) = 0
      end do
      leng = isum

!     --- Fill array ibuurp

      do i = 1 , nelem
         is = 3*(i-1)
         i1 = kelem(is+1)
         i2 = kelem(is+2)
         i3 = kelem(is+3)
         call msho17(ibuurp,iwork,i1,i2)
         call msho17(ibuurp,iwork,i1,i3)
         call msho17(ibuurp,iwork,i2,i1)
         call msho17(ibuurp,iwork,i2,i3)
         call msho17(ibuurp,iwork,i3,i1)
         call msho17(ibuurp,iwork,i3,i2)
      end do
      call erclos( 'msho16' )
      end
