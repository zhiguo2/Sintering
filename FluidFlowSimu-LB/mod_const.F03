!!!--- Module file ---!!!
!!!-Description:
!!!   Common consts that are used in programs.
!!!-Details
!!!   Digits precision control
!!!   Mathematic constants
!!!   Physics constants
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   Basic module


MODULE const_mod
IMPLICIT NONE

!*** Dimension ***!
INTEGER, PARAMETER     :: nd = 3      ! spatial dimension fixed as 3


!*** Digits precision control ***!
! int kind: short(1), ordinary(4) and long(8)
INTEGER, PARAMETER      :: ST = SELECTED_INT_KIND(r=2)
INTEGER, PARAMETER      :: OD = SELECTED_INT_KIND(r=9)
INTEGER, PARAMETER      :: LG = SELECTED_INT_KIND(r=18)

! real kind: single(4), double(8) and quardra(16) precision
INTEGER, PARAMETER      :: SGL = SELECTED_REAL_KIND(r=37)
INTEGER, PARAMETER      :: DBL = SELECTED_REAL_KIND(r=307)
INTEGER, PARAMETER      :: QUD = SELECTED_REAL_KIND(r=4931)

! minimum real for double precision
REAL(DBL), PARAMETER    :: EPS = EPSILON(1d0)


!*** Mathematics constants ***!
! PI
REAL(DBL), PARAMETER    :: PI = 3.141592653589793238462643383


!*** Physics constants ***!
! RGC - gas constant
REAL(DBL), PARAMETER    :: RGC = 8.314  ! [J/mol.K]


END MODULE const_mod
