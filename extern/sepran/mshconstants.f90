module mshconstants

double precision, parameter :: EPSMAC = 1.0d-15
double precision, parameter :: SQREPS = 1.0d-15
double precision, parameter :: RINFIN = 1.0d77
integer, parameter :: IREFWR = 66
integer :: ITIME = 1
integer :: IGOBS = 1
integer :: JTIMES = 1

end module

module msherror

integer :: IERROR = 0

end module
