function r1mach(i)
!
!*******************************************************************************
!
! R1MACH returns single precision machine constants.
!!
! Assume that single precision numbers are stored with a mantissa of T digits
! in base B, with an exponent whose value must lie between EMIN and EMAX.  Then
! for values of I between 1 and 5, R1MACH will return the following values:
!
! R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
! R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
! R1MACH(3) = B**(-T), the smallest relative spacing.
! R1MACH(4) = B**(1-T), the largest relative spacing.
! R1MACH(5) = log10(B)
!
! To alter this function for a particular environment, the desired set of data
! statements should be activated by removing the C from column 1.
!
! On rare machines a STATIC statement may need to be added.  But probably more
! systems prohibit it that require it.
!
! For IEEE-arithmetic machines (binary standard), the first set of constants
! below should be appropriate.
!
! Where possible, octal or hexadecimal constants have been used to specify the
! constants exactly which has in some cases required the use of EQUIVALENCED
! integer arrays.  If your compiler uses half-word integers by default
! (sometimes called integer*2), you may need to change INTEGER to INTEGER*4 or
! otherwise instruct your compiler to use full-word integers in the next 5
! declarations.
!
integer diver(2)
integer i
integer large(2)
integer log10(2)
real r1mach
integer right(2)
real rmach(5)
integer small(2)
!
equivalence (rmach(1),small(1))
equivalence (rmach(2),large(1))
equivalence (rmach(3),right(1))
equivalence (rmach(4),diver(1))
equivalence (rmach(5),log10(1))
!
! IEEE arithmetic machines, such as the ATT 3B series, Motorola 68000 based
! machines such as the SUN 3 and ATT PC 7300, and 8087 based micros such as
! the IBM PC and ATT 6300.
!
! data small(1) /     8388608 / !disabled in favor of the following set (11/30/06)
! data large(1) /  2139095039 / !disabled in favor of the following set (11/30/06)
! data right(1) /   864026624 / !disabled in favor of the following set (11/30/06)
! data diver(1) /   872415232 / !disabled in favor of the following set (11/30/06)
! data log10(1) /  1050288283 / !disabled in favor of the following set (11/30/06)
!
! IBM PC - Professional FORTRAN and Lahey FORTRAN:  The values below are a little 
! better than the values above.
!
data small(1)/ Z'00800000' /
data large(1)/ Z'7F7FFFFF' /
data right(1)/ Z'33800000' /
data diver(1)/ Z'34000000' /
data log10(1)/ Z'3E9A209A' /
!
return
end function r1mach
