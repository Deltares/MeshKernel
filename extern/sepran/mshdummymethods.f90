module mshdummymethods

implicit none

private

public eropen
public ersettime
public erclos
public errsub
public errint
public errchr
public erreal
public prinrl1
public prinin
public prinin1
public eralloc
public erdealloc
public printtime

!===============================================================================
contains
!
!------------------------------------------------------------------------------
subroutine eropen(aErrorString)
    character(len=*) :: aErrorString
    
end subroutine

subroutine ersettime(aErrorTime)
    doubleprecision :: aErrorTime
    integer i
    i = 1
end subroutine

subroutine erclos(aErrorString)
    character :: aErrorString
    integer i
    i = 1
end subroutine

subroutine errsub(aError1, aError2, aError3, aError4)
use msherror
    integer :: aError1
    integer :: aError2
    integer :: aError3
    integer :: aError4
    integer i
       
    ierror = 1
end subroutine

subroutine errint(aError1, aError2)
    integer :: aError1
    integer :: aError2
    integer i
    i = 1
end subroutine

subroutine errchr(aErrorString, aErrorNumber)
    character :: aErrorString
    integer :: aErrorNumber
    integer i
    i = 1
    
end subroutine

subroutine erreal(aErrorValue, aErrorNumber)
    doubleprecision :: aErrorValue
    integer :: aErrorNumber
    integer i
    i = 1
end subroutine

subroutine prinrl1( rarr, n1, n2, name)
    double precision, intent(in) :: rarr(*)
    integer, intent(in)          :: n1, n2
    character(len=*), intent(in) :: name

    write (*,*) 'prinrl1 not implemented yet'
end subroutine

subroutine prinin( iarr, n1, name)
    integer, intent(in)          :: iarr(*)
    integer, intent(in)          :: n1
    character(len=*), intent(in) :: name

    write (*,*) 'prinin not implemented yet'
end subroutine

subroutine prinin1( iarr, n1, n2, name)
    integer, intent(in)          :: iarr(*)
    integer, intent(in)          :: n1, n2
    character(len=*), intent(in) :: name

    write (*,*) 'prinin1 not implemented yet'
end subroutine

subroutine eralloc ( ierror, numel, name)
    integer, intent(in)          :: ierror
    integer, intent(in)          :: numel
    character(len=*), intent(in) :: name

    write (*,*) 'eralloc not implemented yet'
end subroutine

subroutine erdealloc ( ierror, name)
    integer, intent(in)          :: ierror
    character(len=*), intent(in) :: name

    write (*,*) 'erdealloc not implemented yet'
end subroutine

subroutine printtime(text, time)
    character(len=*), intent(in) :: text
    double precision, intent(in) :: time

    write (*,*) 'printtime not implemented yet'
end subroutine

end module