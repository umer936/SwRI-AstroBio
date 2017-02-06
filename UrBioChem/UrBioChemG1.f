      program UrBioChemG1
!
!     UrBioChemG1:  
!     G1 mean gasphase with H, He, C, N, and O only.
!     ******************************************************************
!     Have FChem and Rates been replaced with the new FChem and Rates
!     from FORTRAN Equation Writer Subroutins (EWS.f)?
!     LenTbl and iCRTbl are read from file SpecReac.  The file, SRT, 
!     produced by Eqwrt can be copied into UrBoChemG1 and renamed 
!     SpecReac. 
!     MxCRTbl should be larger than LenTbl.  If MxCRTbl is changed, 
!     its octal equivalent in common/lc4/ in program olload must also 
!     be changed.
!     If dimensions of Sig in lc3 are changed, then its octal eqivalent
!     in lc3 in program olload must also be changed.
!     Input is obtained from hydra disk by using:  switch ComInp Input.
!     ******************************************************************
!     NYDim is the maximum size of the array for the number of Equs.
!     If NYDim is changed, then change the numerical value in the write 
!     statemens fmt = "(a9, 600i9)") accordingly.
!     NReac is the number of chemical reactions in the system.
!     N is the number of species in the system including M.
!     rNuc is the mean radius of the protosun in km. 
!     rh is the heliocentric distance in AU.
!     Temp is the temperature in degree kelvin. 
!     t is the time in sec (1 year = 3.155693e+07 sec).
! **********************************************************************
!     Use tOut = log10(tOut + 0.1
!     and
!     tOut = 10.0**tOut?            ! ?  See UMcse code ?
! **********************************************************************
      parameter (NyDim = 600, NyDim10 = 610, NyDim50 = 650, MxNt = 500, 
     1  MxPlPts = 1000, MxCRTbl = 30000, MxStoich = 18, 
     3  SSAge = 4.568E+09)                                 ! in years.  
      dimension jStoich(MxStoich, NYDim)
      real :: k
      dimension tlog(MxNt)
      dimension YP(NYDim)
      integer, parameter :: MxReac = 8000, LTabl = 30000, LIW = 160000, 
     1  LCW = 120000                                        ! NSym = 7 
      dimension NAt(NYDim), NDegF(NYDim), T3(NYDim), P3(NYDim)
      real :: MxtLim
      real, dimension(MxPlPts) :: PlOut
      real, dimension(MxPlPts, NYDim) :: PlotData
      character*4 Cat, SPD, SeD, SPI, SDI        ! Neded for solar distance scaling
      character*8 Prog
c      character*8 Speci, Symbol(7,MxReac), SpecIn, x, Quelle
      character*8 Speci, SpecIn, x, Quelle
      character(len = 92) :: ColHeads
!
      common/lc1/ iStoch(MxStoich, NYDim), WtM(NYDim)      ! lc1  not needed?
      common/lc2/ A(MxReac), B(MxReac), C(MxReac), TLo(MxReac), 
     1  THi(MxReac), Cat(MxReac), Quelle(MxReac) 
      common/lc4/ iCRTbl(MxCRTbl)
      common/char1/ Speci(NYDim)
      common/char2/ Symbol(7, MxReac)
      common/stats/ avgh, hused, avgord, nqused, imxerr, NStep, nfe, 
     1  nje        !, NtSteps
      common/stcom9/ k(MxReac)
      common/stcom11/ Q(NYDim)
      common/rivst/ tOut, N, t, Y(NYDim)
      common/sdriv/ NState, Ntask, eps, ewt(NYDim), iError, mf, ml, mu, 
     1  MxOrd, hmax, Work(NYDim10, NYDim10), LenW, iWork(NYDim50), NDE, 
     2  LOut, MxStep, NSpeci, mint, miter, impl, NRoot, LeniW, Meth 
      common/gph/ Nt
!
      data secperyr/3.155693e+07/  ! Number of seconds in an astronomical year.  
!
      external FChem
!
      open(unit = 7, file = "SpecReac", status = "old")
      open(unit = 8, file = "Results", status = "replace")
      open(unit = 9, file = "Plot", status = "replace")
      open(unit = 10, file = "Input", status = "replace")
      open(unit = 11, file = "Output", status = "replace")
! ***
      read(unit = 7, fmt = "(5i4, f5.2)") NReac, NSpeci, NDE, NHot, ! SpecReac 
     1  iPrnt, rh
      NAll = NSpeci + NHot
      read(unit = 7, fmt = "(a8, 18i3, f9.4, i3, i2, f7.2, e10.3)") 
     1  (Speci(i), (iStoch(j, i), j = 1, MxStoich), WtM(i), NAt(i), 
     2  NDegF(i), T3(i), P3(i), i = 1, NAll) 
!
      if(NSpeci .gt. NYDim) then 
        write(unit = *, fmt = *) " NSpeci > NYDim.  Increase NYDim."
        stop
      else
      end if
      itime = 1
      t = 0.0
      dtlg = 0.05             ! Logarithmic time step.  
!      dtlg = 1.0             ! Logarithmic time step.  
!      dtlg = 0.50            ! Logarithmic time step.  
!      dtlg = 0.20            ! Logarithmic time step.  
!      dtlg = 0.10            ! Logarithmic time step.  
!      dtlg = 0.05            ! Logarithmic time step.  
      Nt = 0
! ***
      tlg = -1.0
  !    tlg0 = -1.0
      tlg0 = 0.0
! ***
      tlog(1) = tlg0
      write(unit = *, fmt = *) 
     1  "Type in the gas number density [1/cm^3] in this format:  ", 
     1  "1pe8.2,"
!      TotDens = 2.5e+12
!      TotDens = 2.5e+03
      TotDens = 3.70004E+13  ! @ 1 AU, Midplane
!      read(unit = *, fmt = "(1pe8.2)") TotDens
      write(unit = *, fmt = *) 
     1  "Temperature [K] in this format:  1pe8.2, "
!      Temp = 1.00e+02
!      Temp = 60.
      Temp = 353.44          ! @ 1 AU, Midplane
!      read(unit = *, fmt = "(1pe8.2)") Temp
      write(unit = *, fmt = *)
     1  "and total time span, tLim, in sec in this format:  1pe8.2"
!      read(unit = *, fmt = "(1pe8.2)") tLim
      MxtLim = SSAge*secperyr                       ! in sec
!      tLim = 1.0E+17                   ! sec
!      tLim = 1.0E+14                   ! sec
!      tLim = 1.0E+11                   ! sec
!      tLim = 1.0E+10                   ! sec
!      tLim = 1.0E+09                   ! sec
!      tLim = 1.0E+08                   ! sec
!      tLim = 1.0E+07                   ! sec
      tLim = 1.0E+04                   ! sec
!      tLim = 1.0E+01                   ! sec
      if(tLim .gt. MxtLim) then
        write(unit = *, fmt = *) " tLim = ", tLim, " sec,  "
        write(unit = *, fmt = *) " exceeds the age of the Solar System,"
        write(unit = *, fmt = *) " MxtLim = ", MxtLim, " sec." 
        stop
      else
      end if
      write(unit = *, fmt = "(a17, 1pe9.2, a13, e9.2, a10, e9.2, a8, 
     1  e9.2, a5)") "Number density =", TotDens, " cm^(-3), T =", 
     2  Temp, " K, tLim =", tLim/secperyr, " years =", tLim, " sec."
      NtSteps = log10(tLim)/dtlg + int(1./dtlg) + 1
      write(unit = *, fmt = "(a14, 1pe10.3, a9, e9.3, a1)") 
     1  " log10(tLim) =", log10(tLim), ", dtlg = ", dtlg, "."
      if(NtSteps .gt. MxPlPts) then
        write(unit = *, fmt = "(a9, 1pe10.3, a12, e10.3)") 
     1    "NtSteps = ", NtSteps, "> MxPlPts = ", MxPlPts
        write(unit = *, fmt = *) " Increase MxPlPts."
        stop
      else
      end if
      write(unit = *, fmt = "(a18, 1pe10.3, a40, i4, a1)") 
     1  "Total time span =", tLim, " sec.  The total number of time step
     2s = ", NtSteps, "."
!     Set element abundances.  The default is solar abundances.  
!
!     Log Solar photospheric number abundances:  Asplund, M., Grevesse, 
!       N., Sauval, A. J., Scott, P.,  The chemical composition of the 
!       Sun, Annual Rev. Astron. Astrophys. 47, 481-522 (2009).
!
!     H  : 12.00, He : 10.93, Li :  1.05, Be :  1.38, B  :  2.70 
!     C  :  8.43, N  :  7.83, O  :  8.69, F  :  4.56, Ne :  7.93 
!     Na :  6.24, Mg :  7.60, Al :  6.45, Si :  7.51, P  :  5.41 
!     S  :  7.12, Cl :  5.50, Ar :  6.40, K  :  5.03, Ca :  6.34 
!     Sc :  3.15, Ti :  4.95, V  :  3.93, Cr :  5.64, Mn :  5.43 
!     Fe :  7.50, Co :  4.99, Ni :  6.22, Cu :  4.19, Zn :  4.56 
!
      YSum = 0.0
      write(unit = 11, fmt = "(a18, a15, 1pe10.3, a10, a6, e10.3, a2)") 
     1  "Solar Nebula G1:  ", "Total Density =", TotDens, " cm^(-3), ", 
     2  "Temp =", Temp, "K"
      write(unit = 11, fmt = "(a8, i4)")
     1  "NSpeci =", NSpeci
      write(unit = 8, fmt = "(a16)") " Solar Nebula G1"
      do i = 2, NSpeci       ! Calculate abundances, excluding electrons.
        Y(i) = 1.00e-10
        YP(i) = 0.0
        Q(i) = 0.0
      end do
!      write(unit = *, fmt = *) "Relative input densities of elements:"
  !    do j = 1, NSpeci
!        if(Speci(j) == "H ") then
!          write(unit = *, fmt = *) " if ", Speci(j), "== H"
!          Y(j) = 1.00
!          YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!        else
!          if(Speci(j) == "He") then
!            Y(j) = 10.0**(10.93 - 12.00)
!            YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          else
!            if(Speci(j) == "C ") then
!              Y(j) = 10.0**(8.43 - 12.00)
!              YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!            else
!              if(Speci(j) == "N ") then
!                Y(j) = 10.0**(7.83 - 12.00)
!                YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!              else
!                if(Speci(j) == "O ") then
!                  Y(j) = 10.0**(8.69 - 12.00)
!                  YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                else
!                  if(Speci(j) == "F ") then
!                    Y(j) = 10.0**(4.56 - 12.00)
!                    YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                  else
!                    if(Speci(j) == "Na") then
!                      Y(j) = 10.0**(6.24 - 12.00)
!                      YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                    else
!                      if(Speci(j) == "Si") then
!                        Y(j) = 10.0**(7.51 - 12.00)
!                        YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                      else
!                        if(Speci(j) == "P ") then
!                          Y(j) = 10.0**(5.44 - 12.00)
!                          YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                        else
!                          if(Speci(j) == "S ") then
!                            Y(j) = 10.0**(7.12 - 12.00)
!                            YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                          else
!                            if(Speci(j) == "Cl") then
!                              Y(j) = 10.0**(5.50 - 12.00)
!                              YSum = YSum + Y(j)
!          write(unit = *, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!          write(unit = 8, fmt = "(a5, i3, 1x, a2, a5, f10.7, a8, 
!     1 f10.7)") "j = ", j, Speci(j), " Y = ", Y(j), " YSum = ", YSum
!                            else
!                            end if
!                          end if
!                        end if
!                      end if
!                    end if
!                  end if
!                end if
!              end if
!            end if
!          end if
!        end if
! *****
!      end do
! *****
! The folowing are aproximate initil abundnces basd on solar abundces of H, 
! He, C, N, and O.  We assumed that 1/5 of the total C is in each C, CH, CH2, 
! CH3 , and CH4.  Similar for N (1/4), and O (1/3).  The total H is then 
! adjusted to make the sum of all hydrogen bearing species to add up to 
! 1.000E+12.  The electrons are assumed to be arbitrarily 1.000E+11.  H+ is 
! equal to e- to maintain neutrality.
!

#include "Y.h"


      write(unit = *, fmt = "(a12, 1pe10.3, a9)") "M:  YSum = ", YSum, 
     1  " cm^(-3)."

! *********
! ARE THESE ABUNDANCE RATIOS CONSISTENT WITH SOLAR VAULES OF H;He;C;N;O?
! DEVELOP AN INPUT FILE FOR INITIAL ABUNDANCE VALUES.
! *********
! Renormalize
      write(unit = *, fmt = *) "Normalized densities of starting ",
     1  "composition:"
      Y(1) = TotDens*Y(1)/YSum 
      write(unit = *, fmt = "(1x, a8, a4, 1pe10.3, a8)") Speci(1), 
     1  " Y =", Y(1), " cm^(-3)" 
      do j = 2, NDE          ! Excludes electrons and total number M.
        if(Y(j) .gt. 1.01e-10) then
          Y(j) = TotDens*Y(j)/YSum 
          write(unit = *, fmt = "(1x, a8, a4, 1pe10.3, a8)") Speci(j), 
     1      " Y =", Y(j), " cm^(-3)" 
          write(unit = 8, fmt = "(1x, a8, a5, 1pe10.3, a8)") Speci(j), 
     1      " Y = ", Y(j), " cm^(-3)" 
        else
        end if
      end do
      Y(NSpeci) = 0.0
      do i = 2, NSpeci - 1        ! This excludes the density of e and M.
        Y(NSpeci) = Y(NSpeci) + Y(i)
      end do
      write(unit = *, fmt = "(1x, a8, a4, 1pe10.3, a8)") Speci(NSpeci), 
     1  " Y =", Y(NSpeci),  " cm^(-3) X" ! This corresponds to the density of M
      write(unit = 8, fmt = "(1x, a8, a5, 1pe10.3, a8)") Speci(NSpeci), 
     1  " Y = ", Y(NSpeci),  " cm^(-3)"          ! This corresponds to the density of M
      write(unit = 8, fmt = "(1x, a8, i5)") !  a8, a5, e10.3, a8)")  
     1  "NSpeci =", NSpeci !Speci(NSpeci), " Y = ", Y(NSpeci),  " cm^(-3)" !This corresponds to the density of M
      LenW = 2*NYDim*NYDim + (MxOrd + 5)*NYDim + 2*NRoot + 250
      LeniW = NYDim50
!
      write(unit = *, fmt = "(a8, i5, a10, i4, a7, i4, a8, i3, a9, i2, 
     1  a6, f5.2)") " NReac =", NReac, ", NSpeci =", NSpeci, ", NDE =", 
     2  NDE, ", NHot =", NHot, ", iPrnt =", iPrnt, ", rh =", rh 
      if(iPrnt .gt. 1) then
        write(unit = *, fmt = "(a9, i4, a11, i4, a10, i3, a11, i5, a10, 
     1    i4)") " NReac = ", NReac, ", NSpeci = ", NSpeci,  ", 
     2    Nhot = ", NHot, ", iPrnt = ", iPrnt, ", MxReac = ", 
     3    MxReac, ", NYDim = ", NYDim
      else
      end if
      if(NReac .gt. MxReac) then
        write(unit = *, fmt = "(a8, i5, a25, i5)") " NReac = ", NReac, 
     1    " is greater than MxReac =", MxReac
        stop
      else
      end if
      write(unit = *, fmt = "(a7, i4, a9, i4)") " NAll =", 
     1  NAll, ", NYDim =", NYDim
      write(unit = 8, fmt = "(a7, i4, a5, i4, a9, i4)") " NAll = ", 
     1  NAll, ", N =", N, ", NYDim =", NYDim
      if(iPrnt .gt. 2) then
        write(unit = 8, fmt = "(i4, 2x, a8, 18i3, f7.3, f9.4)") 
     1    (i, Speci(i), (iStoch(j, i), j = 1, MxStoich), WtM(i), 
     2    i = 1, NAll)
      else
      end if
      write(unit = 8, fmt = *) " NReac =", NReac
      read(unit = 7, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0,  
     1  0pf7.0, a4, a8)")
     2  ((Symbol(j, i), j = 1, 7), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i), i = 1, NReac) 
      read(unit = 7, fmt = "(i5)") LenTbl
      if(LenTbl .gt. MxCRTbl) then
        write(unit = *, fmt = "(a9, i6, a11, i6)") " LenTbl =", LenTbl, 
     1    ", MxCRTbl =", MxCRTbl  
        write(unit = *, fmt = *) " Increase MxCRTbl to at least ", 
     1    LenTbl
        stop
      else
      end if
      write(unit = 8, fmt = "(a9, i6)") " LenTbl =", LenTbl
      read(unit = 7, fmt = "(16i5)") (iCRTbl(i), i = 1, LenTbl)
!
!     Write rate coefficients, k
      do i = 1, NReac
        k(i) = A(i)*((Temp/300.0)**B(i))*exp(-C(i)/Temp)
        if(Cat(i) .eq. "SPD" .or. Cat(i) .eq. "SeD" .or. Cat(i) .eq. 
     1    "SPI" .or. Cat(i) .eq. "SDI") then
          k(i) = k(i)/rh**2
        else
        end if
      end do
      do i = 1, MxPlPts
        PlOut(i) = 0.0
      end do
  100 tOut = exp(log(10.)*tlg)         ! in sec.  For tlg = 3, tOut = 1000 sec
!  100 tOut = exp(log(10.)*tLim)         ! in sec.  For tlg = 3, tOut = 1000 sec
!  100 tOut = tLim
      PlOut(itime) = tOut
      write(unit = *, fmt = "(a8, i3, a9, 1pe10.3, a5)") " itime =", 
     1  itime, ", PlOut =", PlOut(itime), " sec."
!      write(unit = *, fmt = "(a5, i3, a12, f5.2, a5, 1pe9.2, a8, e9.2, 
!     1  a10, i2, a9, i2, a9, i2, a1)") " Nt =", Nt, ", tlog(Nt) =", 
!     2  tlog(Nt), ", t =", t, ", tOut =", tOut, ", NState =", NState, 
!     3  ", NTask =", NTask, ", NRoot =", NRoot, ","
!      write(unit = *, fmt = "(a6, 1pe9.2, a8, i2, a9, i2, a8, 
!     1  i2, 2(a6, i2), a9, i2, a7, i5)")  " EPS =", EPS, ", mint =", 
!     2  mint, ", miter =", miter, ", impl =", impl, ", ML =", ML, 
!     3  ", MU =", MU, ", MxOrd =", MxOrd, ", NDE =", NDE
!      write(unit = *, fmt = *) " Before call to sdrvb3"
!      write(unit = *, fmt = "(a8, i4, a5, 1pe10.3)") " NYDim =", NYDim, 
!     1  ", t =", t
! xxxxxxxxxx
      write(unit = *, fmt = "(a9, i2, a8, 1pe9.2, 2(a9, i2), a10, i2, 
     1  a8, i2, a9, i2, a8, i2, 2(a6, i2), a9, i2)") 
     2  " NState =", NState, ", tOut =", tOut, ", NTask =", NTask, 
     3  ", NRoot =", NRoot, ", iError =", iError, ", Mint =", Mint, 
     4  ", Miter =", Miter, ", Impl =", Impl, ", Ml =", Ml, ", Mu =", 
     5  Mu, ", MxOrd =", MxOrd
      call sdrvb3(NYDim, t, Y, FChem, NState, tOut, NTask, NRoot, eps, 
     1  ewt, iError, mint, miter, impl, ml, mu, MxOrd, hmax, work, 
     2  LenW, iWork, LeniW, FChem, FChem, NDE, MxStep, FChem)
!      write(unit = *, fmt = *) " After call sdrvb3:  **** t =", t   ,
!     1  " Y(1-3,7,9,43) =", (Y(ii), ii = 1, 3), Y(7), Y(9), Y(43)
!      pause
      write(unit = *, fmt = "(a9, i2, a8, 1pe9.2, a7, e9.2)") 
     1  " NState =", NState, " tOut =", tOut, ", eps =", eps
! *****
      if(NState == 3) then
!        NState = 1
        go to 100
      else
      end if
! *****
!      write(unit = *, fmt = *) " hmax =", hmax, " LenW =", LenW, 
!     1  " iWork(1-20) =", (iWork(ii), ii = 1, 20), " NDE =", NDE, 
!     2  " MxStep =", MxStep
!      pause

! *****
! ksp - check for Nt > 0 before tlog lookup
      if (Nt > 0) then
        if((tlg .lt. tlog(Nt))) then
          go to 200
        else
        end if
      end if
      Nt = Nt + 1
!      tlog(Nt) = tlg + 1.0
      tlog(Nt) = tlg + dtlg
! *****
      write(unit = *, fmt = "(a5, i4, a27, 1pe10.3)") " Nt =", Nt, 
     1  "                tlog(Nt) = ", tlog(Nt)
      write(unit = *, fmt = "(a15, 1pe10.3, a11, e10.3, a20, e10.3)") 
     1  " tlg = ", tlg, " tOut =    ", tOut, " sec, PlOut(itime) = ", 
     2  PlOut(itime)
! *****
  200 continue
      tlg = tlg + dtlg + 1.e-12 
! *****
      write(unit = *, fmt = "(a6, 1pe10.3, a6, i4, a11, i4, a8, e10.3)") 
     1  " tlg =", tlg, ", Nt =", Nt, ", NtSteps =", NtSteps, 
     2  ", tLim =", tLim
! *****
      if(Nt .le. NtSteps) then
        do j = 1, NSpeci - 1
          if(Y(j) .lt. 1.0E-32) then
            Y(j) = 1.0E-32
          else
          end if
          PlotData(itime, j) = Y(j)
        end do
        Y(NSpeci) = 0.0
        do j = 2, NSpeci - 1
          Y(NSpeci) = Y(NSpeci) + Y(j)
        end do
        write(unit = *, fmt = "(a4, 1pe10.3)") " M =", Y(NSpeci)
        PlotData(itime, NSpeci) = Y(NSpeci)
!        pause
!        if(itime .ge. MxPlPts) then
!          write(unit = *, fmt = *) " Increase MxPlPts for PlotData array
!     1"
!          stop
!        else
          itime = itime + 1
! *****
          NState = 1
! *****
          go to 100
!        end if
      else
        if(iPrnt .ge. 0) then
          write(unit = 8, fmt = "(a16)") "Y"
          write(unit = 8, fmt = "(i4, a9, 1pe10.3)") (i, 
     1      Speci(i), Y(i), i = 1, NSpeci)
        else
        end if
        write(unit = *, fmt = "(a7, e10.3)") " t =", t, " Temp =", Temp
        write(8, 280) tOut, Temp, (Speci(i), Y(i), i = 1, NSpeci)
 280    format(" Time =", 1pe10.3, " sec, Temp =", 0pf8.1, " K"/
     1    " Spec. Concentrat.   Spec. Concentrat.   Spec. Concentrat.   
     2  Spec. Concentrat.   Spec. Concentrat."/
     7    (5(1x, a6, 1pe11.4, 2x)))
        write(unit = 9, fmt = "(a9, 600i9)") "   time  ",  ! 600 = NYDim
     1    (j, j = 1, NSpeci)
        write(unit = 9, fmt = "(600(a9))") "   [sec] ",    ! 600 = NYDim
     1    (Speci(i), i = 1, NSpeci)
        t0 = 0.0
        write(unit = *, fmt = "(600(1pe9.2))")             ! Plot, 600 = NYDim
     1    t0, (Y(j), j = 1, NSpeci)                        ! at t = 0
!      pause
! *****
        write(unit = *, fmt = *) " PlOut (Temps) for NtSteps = ", 
     1    NtSteps
        write(unit = *, fmt = "(2(1x, 1pe8.2))") (PlOut(j), j = 1, 
     1    NtSteps)
!        write(unit = 9, fmt = "(a11, 2x, 10(1pe9.2))") "    [cm^-3]", 
!     1    (PlOut(j), j = 1, NtSteps)
!        do i = 1, NSpeci
!          write(unit = 9, fmt = "(i4, 1x, a8, 10(1pe9.2))") i, 
!     1      Speci(i), (PlotData(j, i), j = 1, NtSteps) 
  !      write(unit = 9, fmt = "(600(1pe9.2))") PlOut(Temps), 
  !   1    (PlOut(j), j = 1, NSpeci)
!        do i = 1, MxPlPts
! *****
!        do i = 1, NtSteps
! *****
        do i = 1, NSpeci
          write(unit = 11, fmt = "(a8, i3, a3, 1pe10.3, a3, a8)")
     1      "Y(", i, ") =", PlotData(NtSteps, i), " ! ", Speci(i)
        end do
! *****
        do i = 1, NtSteps + 1
          write(unit = 9, fmt = "(600(1pe9.2))")           ! Plot, 600 = NYDim
     1      PlOut(i), (PlotData(i, j), j = 1, NSpeci)      ! PlOut(i) = time
        end do
        write(unit = *, fmt = *) "      ***** FINISHED *****"
        close(unit = 8)
        close(unit = 9)
        close(unit = 10)
        close(unit = 11)
        stop
      end if
      end program UrBioChemG1
!
!     do 410 i=1, nall
!       ylog(i, nt) = alog10(abs(y(i)))
!410  continue
!
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Adjust ewt downward in BLOCK DATA
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      block data
      parameter (NYDim = 600, NYDim10 = 610, NYDim50 = 650, 
! 8/29/16
     1  MxReac = 8000)
! 8/29/16
      common/sdriv/ NState, Ntask, eps, ewt(NYDim), iError, mf, ml, mu, 
     1  MxOrd, hmax, Work(NYDim10, NYDim10), LenW, iWork(NYDim50), NDE, 
     2  LOut, MxStep, NSpeci, mint, miter, impl, NRoot, LeniW, Meth
! 8/29/16
      common/char2/ Symbol(7, MxReac)
! 8/29/16
!      common/stcom10/ Rat(MxReac), ewt 
      common/stcom11/ Q(NYDim)
! *****
!      data NState/1/, NTask/3/, eps/1.0e-04/, ewt/599*0.001, 1.0e-04/ ! 599 = NYDim - 1
!      data NState/1/, NTask/3/, eps/1.0e-03/, ewt/599*0.001, 1.0e-04/ ! 599 = NYDim - 1
      data NState/1/, NTask/3/, eps/1.0e-02/, ewt/599*0.001, 1.0e-04/ ! 599 = NYDim - 1
!      data NState/1/, NTask/3/, eps/1.0e-02/, ewt/599*0.01, 1.0e-04/ ! 599 = NYDim - 1
!      data NState/1/, NTask/3/, eps/1.0e-02/, ewt/599*0.001, 1.0e-02/ ! 599 = NYDim - 1
!      data NState/1/, NTask/3/, eps/1.0e-01/, ewt/599*0.1, 1.0e-01/ ! 599 = NYDim - 1
!      data NState/1/, NTask/3/, eps/1.0e-01/, ewt/599*0.1, 1.0e-01/ ! 599 = NYDim - 1
!      data ierror/4/, mf/022/, MxOrd/5/, hmax/3.0e+03/, Meth/2/ 
      data ierror/3/, mf/022/, MxOrd/5/, hmax/3.0e+03/, Meth/2/ 
      data LOut/1/, mxStep/3000/, NRoot/0/, mint/2/, miter/2/, impl/0/
      data Q/NYDim*0.0/
      end
!
      subroutine SDrvb3(N, t, Y, F, Nstate, tOut, NTask, NRoot, eps,    sdrvb3 3
     1  ewt, iError, MInt, MIter, Impl, ML, MU, MxOrd, HMax, Work, LenW,sdrvb3 4
     2  iWork, LenIW, Jacobn, FA, NDE, MxStep, G)                       sdrvb3 5
C   Begin prologue sdrvb3                                               sdrvb3 6
C***DATE WRITTEN   790601   (YYMMDD)                                    sdrvb3 7
C***REVISION DATE  841119   (YYMMDD)                                    sdrvb3 8
C***CATEGORY NO.  I1A2,I1A1B                                            sdrvb3 9
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                sdrvb310
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,                     sdrvb311
C             SINGLE PRECISION                                          sdrvb312
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdrvb313
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdrvb314
C***PURPOSE  The function of SDRVB3 is to solve N ordinary differential sdrvb315
C            equations of the form dY(I)/dT = F(Y(I),T), given the      sdrvb316
C            initial conditions Y(I) = YI.  The program has options to  sdrvb317
C            allow the solution of both stiff and non-stiff differentialsdrvb318
C            equations.  Other important options are available.  SDRVB3 sdrvb319
C            uses single precision arithmetic.                          sdrvb320
C***DESCRIPTION                                                         sdrvb321
C                                                                       sdrvb322
C  I.  ABSTRACT  .......................................................sdrvb323
C                                                                       sdrvb324
C    The primary function of SDRVB3 is to solve N ordinary differential sdrvb325
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial      sdrvb326
C    conditions Y(I) = YI.  The program has options to allow the        sdrvb327
C    solution of both stiff and non-stiff differential equations.  In   sdrvb328
C    addition, SDRVB3 may be used to solve:                             sdrvb329
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is sdrvb330
C         a non-singular matrix depending on Y and T.                   sdrvb331
C      2. The hybrid differential/algebraic initial value problem,      sdrvb332
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may sdrvb333
C         depend upon Y and T) some of whose components will be zero    sdrvb334
C         corresponding to those equations which are algebraic rather   sdrvb335
C         than differential.                                            sdrvb336
C    SDRVB3 is to be called once for each output point of T.            sdrvb337
C                                                                       sdrvb338
C  II.  PARAMETERS  ....................................................sdrvb339
C                                                                       sdrvb340
C    The user should use parameter names in the call sequence of SDRVB3 sdrvb341
C    for those quantities whose value may be altered by SDRVB3.  The    sdrvb342
C    parameters in the call sequence are:                               sdrvb343
C                                                                       sdrvb344
C    N      = (Input) The number of dependent functions whose solution  sdrvb345
C             is desired.  N must not be altered during a problem.      sdrvb346
C                                                                       sdrvb347
C    T      = The independent variable.  On input for the first call, T sdrvb348
C             is the initial point.  On output, T is the point at which sdrvb349
C             the solution is given.                                    sdrvb350
C                                                                       sdrvb351
C    Y      = The vector of dependent variables.  Y is used as input on sdrvb352
C             the first call, to set the initial values.  On output, Y  sdrvb353
C             is the computed solution vector.  This array Y is passed  sdrvb354
C             in the call sequence of the user-provided routines F,     sdrvb355
C             JACOBN, FA, USERS, and G.                                 sdrvb356
C                                                                       sdrvb357
C    F      = A subroutine supplied by the user.  The name must be      sdrvb358
C             declared EXTERNAL in the user's calling program.  This    sdrvb359
C             subroutine is of the form:                                sdrvb360
C                   SUBROUTINE F (N, T, Y, YDOT)                        sdrvb361
C                   REAL Y(*), YDOT(*)                                  sdrvb362
C                     .                                                 sdrvb363
C                     .                                                 sdrvb364
C                   YDOT(1) = ...                                       sdrvb365
C                     .                                                 sdrvb366
C                     .                                                 sdrvb367
C                   YDOT(N) = ...                                       sdrvb368
C                   END                                                 sdrvb369
C             This computes YDOT = F(Y,T), the right hand side of the   sdrvb370
C             differential equations.  Here Y is a vector of length at  sdrvb371
C             least N.  The actual length of Y is determined by the     sdrvb372
C             user's declaration in the program which calls SDRVB3.     sdrvb373
C             Thus the dimensioning of Y in F, while required by FORTRANsdrvb374
C             convention, does not actually allocate any storage.  When sdrvb375
C             this subroutine is called, the first N components of Y aresdrvb376
C             intermediate approximations to the solution components.   sdrvb377
C             The user should not alter these values.  Here YDOT is a   sdrvb378
C             vector of length N.  The user should only compute YDOT(I) sdrvb379
C             for I from 1 to N.                                        sdrvb380
C                                                                       sdrvb381
C    NSTATE = An integer describing the status of integration.  The     sdrvb382
C             meaning of NSTATE is as follows:                          sdrvb383
C               1  (Input) Means the first call to the routine.  This   sdrvb384
C                  value must be set by the user.  On all subsequent    sdrvb385
C                  calls the value of NSTATE should be tested by the    sdrvb386
C                  user, but must not be altered.  (As a convenience to sdrvb387
C                  the user who may wish to put out the initial         sdrvb388
C                  conditions, SDRVB3 can be called with NSTATE = 1, andsdrvb389
C                  tOut = t.  In this case the program will return with sdrvb390
C                  NSTATE unchanged, i.e., NSTATE = 1.)                 sdrvb391
C               2  (Output) Means a successful integration.  If a normalsdrvb392
C                  continuation is desired (i.e., a further integration sdrvb393
C                  in the same direction), simply advance TOUT and call sdrvb394
C                  again.  All other parameters are automatically set.  sdrvb395
C               3  (Output)(Unsuccessful) Means the integrator has takensdrvb396
C                  MXSTEP steps without reaching TOUT.  The user can    sdrvb397
C                  continue the integration by simply calling SDRVB3    sdrvb398
C                  again.                                               sdrvb399
C               4  (Output)(Unsuccessful) Means too much accuracy has   sdrvb100
C                  been requested.  EPS has been increased to a value   sdrvb101
C                  the program estimates is appropriate.  The user can  sdrvb102
C                  continue the integration by simply calling SDRVB3    sdrvb103
C                  again.                                               sdrvb104
C               5  (Output) A root was found at a point less than TOUT. sdrvb105
C                  The user can continue the integration toward TOUT by sdrvb106
C                  simply calling SDRVB3 again.                         sdrvb107
C                                                                       sdrvb108
C    TOUT   = (Input) The point at which the solution is desired.  The  sdrvb109
C             position of TOUT relative to T on the first call          sdrvb110
C             determines the direction of integration.                  sdrvb111
C                                                                       sdrvb112
C    NTASK  = (Input) An index specifying the manner of returning the   sdrvb113
C             solution, according to the following:                     sdrvb114
C               NTASK = 1  Means SDRVB3 will integrate past TOUT and    sdrvb115
C                          interpolate the solution.  This is the most  sdrvb116
C                          efficient mode.                              sdrvb117
C               NTASK = 2  Means SDRVB3 will return the solution after  sdrvb118
C                          each internal integration step, or at TOUT,  sdrvb119
C                          whichever comes first.  In the latter case,  sdrvb120
C                          the program integrates exactly to TOUT.      sdrvb121
C               NTASK = 3  Means SDRVB3 will adjust its internal step tosdrvb122
C                          reach TOUT exactly (useful if a singularity  sdrvb123
C                          exists beyond TOUT.)                         sdrvb124
C                                                                       sdrvb125
C    NROOT  = (Input) The number of equations whose roots are desired.  sdrvb126
C             If NROOT is zero, the root search is not active.  This    sdrvb127
C             option is useful for obtaining output at points which are sdrvb128
C             not known in advance, but depend upon the solution, e.g., sdrvb129
C             when some solution component takes on a specified value.  sdrvb130
C             The root search is carried out using the user-written     sdrvb131
C             function G (see description of G below.)  SDRVB3 attempts sdrvb132
C             to find the value of T at which one of the equations      sdrvb133
C             changes sign.  SDRVB3 can find at most one root per       sdrvb134
C             equation per internal integration step, and will then     sdrvb135
C             return the solution either at TOUT or at a root, whicheversdrvb136
C             occurs first in the direction of integration.  The index  sdrvb137
C             of the equation whose root is being reported is stored in sdrvb138
C             the sixth element of IWORK.                               sdrvb139
C             NOTE: NROOT is never altered by this program.             sdrvb140
C                                                                       sdrvb141
C    EPS    = On input, the requested relative accuracy in all solution sdrvb142
C             components.  EPS = 0 is allowed.  On output, the adjusted sdrvb143
C             relative accuracy if the input value was too small.  The  sdrvb144
C             value of EPS should be set as large as is reasonable,     sdrvb145
C             because the amount of work done by SDRVB3 increases as EPSsdrvb146
C             decreases.                                                sdrvb147
C                                                                       sdrvb148
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,        sdrvb149
C             physically meaningful value for the solution.  (Array,    sdrvb150
C             possibly of length one.  See following description of     sdrvb151
C             IERROR.)  Setting EWT smaller than necessary can adverselysdrvb152
C             affect the running time.                                  sdrvb153
C                                                                       sdrvb154
C    IERROR = (Input) Error control indicator.  A value of 3 is         sdrvb155
C             suggested for most problems.  Other choices and detailed  sdrvb156
C             explanations of EWT and IERROR are given below for those  sdrvb157
C             who may need extra flexibility.                           sdrvb158
C                                                                       sdrvb159
C             These last three input quantities EPS, EWT, and IERROR    sdrvb160
C             control the accuracy of the computed solution.  EWT and   sdrvb161
C             IERROR are used internally to compute an array YWT.  One  sdrvb162
C             step error estimates divided by YWT(I) are kept less than sdrvb163
C             EPS in root mean square norm.                             sdrvb164
C                 IERROR (Set by the user) =                            sdrvb165
C                 1  Means YWT(I) = 1. (Absolute error control)         sdrvb166
C                                   EWT is ignored.                     sdrvb167
C                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)sdrvb168
C                                   EWT is ignored.                     sdrvb169
C                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).             sdrvb170
C                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).             sdrvb171
C                    This choice is useful when the solution components sdrvb172
C                    have differing scales.                             sdrvb173
C                 5  Means YWT(I) = EWT(I).                             sdrvb174
C             If IERROR is 3, EWT need only be dimensioned one.         sdrvb175
C             If IERROR is 4 or 5, the user must dimension EWT at least sdrvb176
C             N, and set its values.                                    sdrvb177
C                                                                       sdrvb178
C    MINT   = (Input) The integration method indor.                 sdrvb179
C               MINT = 1  Means the Adams methods, and is used for      sdrvb180
C                         non-stiff problems.                           sdrvb181
C               MINT = 2  Means the stiff methods of Gear (i.e., the    sdrvb182
C                         backward differentiation formulas), and is    sdrvb183
C                         used for stiff problems.                      sdrvb184
C               MINT = 3  Means the program dynamically selects the     sdrvb185
C                         Adams methods when the problem is non-stiff   sdrvb186
C                         and the Gear methods when the problem is      sdrvb187
C                         stiff.  When using the Adams methods, the     sdrvb188
C                         program uses a value of MITER=0; when using   sdrvb189
C                         the Gear methods, the program uses the value  sdrvb190
C                         of MITER provided by the user.  Only a value  sdrvb191
C                         of IMPL = 0 and a value of MITER = 1, 2, 4, orsdrvb192
C                         5 is allowed for this option.  The user may   sdrvb193
C                         not alter the value of MINT or MITER without  sdrvb194
C                         restarting, i.e., setting NSTATE to 1.        sdrvb195
C                                                                       sdrvb196
C    MITER  = (Input) The iteration method indor.                   sdrvb197
C               MITER = 0  Means functional iteration.  This value is   sdrvb198
C                          suggested for non-stiff problems.            sdrvb199
C               MITER = 1  Means chord method with analytic Jacobian.   sdrvb200
C                          In this case, the user supplies subroutine   sdrvb201
C                          JACOBN (see description below).              sdrvb202
C               MITER = 2  Means chord method with Jacobian calculated  sdrvb203
C                          internally by finite differences.            sdrvb204
C               MITER = 3  Means chord method with corrections computed sdrvb205
C                          by the user-written routine named USERS.     sdrvb206
C                          This option allows all matrix algebra and    sdrvb207
C                          storage decisions to be made by the user.    sdrvb208
C                          The routine USERS is called by SDRVB3 when   sdrvb209
C                          certain linear systems must be solved.  The  sdrvb210
C                          user may choose any method to form, store andsdrvb211
C                          solve these systems in order to obtain the   sdrvb212
C                          solution result that is returned to SDRVB3.  sdrvb213
C                          In particular, this allows sparse matrix     sdrvb214
C                          methods to be used.                          sdrvb215
C                          The call sequence for this routine is        sdrvb216
C                                                                       sdrvb217
C                           SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, sdrvb218
C                          8              T, H, EL, IMPL, N, NDE, IFLAG)sdrvb219
C                           REAL Y(*), YH(*), YWT(*),                   sdrvb220
C                          8        SAVE1(*), SAVE2(*), T, H, EL        sdrvb221
C                                                                       sdrvb222
C                          The input variable IFLAG indes what      sdrvb223
C                          action is to be taken.  Subroutine USERS     sdrvb224
C                          should perform the following operations,     sdrvb225
C                          depending on the value of IFLAG and IMPL.    sdrvb226
C                                                                       sdrvb227
C                          IFLAG = 0                                    sdrvb228
C                            IMPL = 0.  USERS is not called.            sdrvb229
C                            IMPL = 1 or 2.  Solve the system           sdrvb230
C                                A*X = SAVE2,                           sdrvb231
C                              returning the result in SAVE2.  The arraysdrvb232
C                              SAVE1 can be used as a work array.       sdrvb233
C                                                                       sdrvb234
C                          IFLAG = 1                                    sdrvb235
C                            IMPL = 0.  Compute, decompose and store thesdrvb236
C                              matrix (I - H*EL*J), where I is the      sdrvb237
C                              identity matrix and J is the Jacobian    sdrvb238
C                              matrix of the right hand side.  The arraysdrvb239
C                              SAVE1 can be used as a work array.       sdrvb240
C                            IMPL = 1 or 2. Compute, decompose and storesdrvb241
C                              the matrix (A - H*EL*J).  The array SAVE1sdrvb242
C                              can be used as a work array.             sdrvb243
C                                                                       sdrvb244
C                          IFLAG = 2                                    sdrvb245
C                            IMPL = 0.   Solve the system               sdrvb246
C                                (I - H*EL*J)*X = H* - YH - SAVE1,      sdrvb247
C                              returning the result in .                sdrvb248
C                            IMPL = 1, or 2.  Solve the system          sdrvb249
C                              (A - H*EL*J)*X = H* - A*(YH + SAVE1)     sdrvb250
C                              returning the result in .                sdrvb251
C                            The array SAVE1 should not be altered.     sdrvb252
C                                                                       sdrvb253
C                          When using a value of MITER = 3, the         sdrvb254
C                          subroutine FA is not required, even if IMPL  sdrvb255
C                          is not 0.  For further information on using  sdrvb256
C                          this option, see section IV-F below.         sdrvb257
C                                                                       sdrvb258
C               MITER = 4  Means the same as MITER = 1 but the A and    sdrvb259
C                          Jacobian matrices are assumed to be banded.  sdrvb260
C               MITER = 5  Means the same as MITER = 2 but the A and    sdrvb261
C                          Jacobian matrices are assumed to be banded.  sdrvb262
C                                                                       sdrvb263
C    IMPL   = (Input) The implicit method indor.                    sdrvb264
C               IMPL = 0 Means solving dY(I)/dT = F(Y(I),T).            sdrvb265
C               IMPL = 1 Means solving A*dY(I)/dT = F(Y(I),T),          sdrvb266
C                        non-singular A (see description of FA below.)  sdrvb267
C                        Only MINT = 1 or 2, and MITER = 1, 2, 3, 4, or sdrvb268
C                        5 are allowed for this option.                 sdrvb269
C               IMPL = 2 Means solving certain systems of hybrid        sdrvb270
C                        differential/algebraic equations (see          sdrvb271
C                        description of FA below.)  Only MINT = 2 and   sdrvb272
C                        MITER = 1, 2, 3, 4, or 5, are allowed for this sdrvb273
C                        option.                                        sdrvb274
C               The value of IMPL must not be changed during a problem. sdrvb275
C                                                                       sdrvb276
C    ML     = (Input) The lower half-bandwidth in the case of a banded  sdrvb277
C             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero    sdrvb278
C             A(R,C).)                                                  sdrvb279
C                                                                       sdrvb280
C    MU     = (Input) The upper half-bandwidth in the case of a banded  sdrvb281
C             A or Jacobian matrix.  (I.e., maximum(C-R).)              sdrvb282
C                                                                       sdrvb283
C    MXORD  = (Input) The maximum order desired. This is .LE. 12 for    sdrvb284
C             the Adams methods and .LE. 5 for the Gear methods.  alsdrvb285
C             value is 12 and 5, respectively.  If MINT is 3, the       sdrvb286
C             maximum order used will be MIN(MXORD, 12) when using the  sdrvb287
C             Adams ods, and MIN(MXORD, 5) when using the Gear      sdrvb288
C             ods.  MXORD must not be altered during a problem.     sdrvb289
C                                                                       sdrvb290
C    HMAX   = (Input) The maximum magnitude of the step size that will  sdrvb291
C             be used for the problem.  This is useful for ensuring thatsdrvb292
C             important details are not missed.  If this is not the     sdrvb293
C             case, a large value, such as the interval length, is      sdrvb294
C             suggested.                                                sdrvb295
C                                                                       sdrvb296
C    WORK                                                               sdrvb297
C    LENW   = (Input)                                                   sdrvb298
C             WORK is an array of LENW real words used                  sdrvb299
C             internally for temporary storage.  The user must allocate sdrvb300
C             space for this array in the calling program by a statementsdrvb301
C             such as                                                   sdrvb302
C                       REAL WORK(...)                                  sdrvb303
C             The following table gives the required minimum value for  sdrvb304
C             the length of WORK, depending on the value of IMPL and    sdrvb305
C             MITER.  LENW should be set to the value used.  The        sdrvb306
C             contents of WORK should not be disturbed between calls to sdrvb307
C             SDRVB3.                                                   sdrvb308
C                                                                       sdrvb309
C      IMPL =   0                   1                   2               sdrvb310
C              ---------------------------------------------------------sdrvb311
C MITER =  0   (MXORD+4)*N +       Not allowed         Not allowed      sdrvb312
C              2*NROOT + 204                                            sdrvb313
C                                                                       sdrvb314
C         1,2  N*N+(MXORD+4)*N     2*N*N+(MXORD+4)*N   N*N+(MXORD+5)*N  sdrvb315
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204  sdrvb316
C                                                                       sdrvb317
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +    sdrvb318
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204    sdrvb319
C                                                                       sdrvb320
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +    sdrvb321
C              (MXORD+5)*N +       (MXORD+6)*N +       (MXORD+6)*N +    sdrvb322
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204    sdrvb323
C              ---------------------------------------------------------sdrvb324
C                                                                       sdrvb325
C    IWORK                                                              sdrvb326
C    LENIW  = (Input)                                                   sdrvb327
C             IWORK is an integer array of length LENIW used internally sdrvb328
C             for temporary storage.  The user must allocate space for  sdrvb329
C             this array in the calling program by a statement such as  sdrvb330
C                       INTEGER IWORK(...)                              sdrvb331
C             The length of IWORK should be at least                    sdrvb332
C               20      if MITER is 0 or 3, or                          sdrvb333
C               N+20    if MITER is 1, 2, 4, or 5, or MINT is 3,        sdrvb334
C             and LENIW should be set to the value used.  The contents  sdrvb335
C             of IWORK should not be disturbed between calls to SDRVB3. sdrvb336
C                                                                       sdrvb337
C    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.    sdrvb338
C             If this is the case, the name must be declared EXTERNAL insdrvb339
C             the user's calling program.  Given a system of N          sdrvb340
C             differential equations, it is meaningful to speak about   sdrvb341
C             the partial derivative of the I-th right hand side with   sdrvb342
C             respect to the J-th dependent variable.  In general there sdrvb343
C             are N*N such quantities.  Often however the equations can sdrvb344
C             be ordered so that the I-th differential equation only    sdrvb345
C             involves dependent variables with index near I, e.g., I+1,sdrvb346
C             I-2.  Such a system is called banded.  If, for all I, the sdrvb347
C             I-th equation depends on at most the variables            sdrvb348
C               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)   sdrvb349
C             then we call ML+MU+1 the bandwith of the system.  In a    sdrvb350
C             banded system many of the partial derivatives above are   sdrvb351
C             automatically zero.  For the cases MITER = 1, 2, 4, and 5,sdrvb352
C             some of these partials are needed.  For the cases         sdrvb353
C             MITER = 2 and 5 the necessary derivatives are             sdrvb354
C             approximated numerically by SDRVB3, and we only ask the   sdrvb355
C             user to tell SDRVB3 the value of ML and MU if the system  sdrvb356
C             is banded.  For the cases MITER = 1 and 4 the user must   sdrvb357
C             derive these partials algebraically and encode them in    sdrvb358
C             subroutine JACOBN.  By computing these derivatives the    sdrvb359
C             user can often save 20-30 per cent of the computing time. sdrvb360
C             Usually, however, the accuracy is not much affected and   sdrvb361
C             most users will probably forego this option.  The optionalsdrvb362
C             user-written subroutine JACOBN has the form:              sdrvb363
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)   sdrvb364
C                   REAL Y(*), DFDY(MATDIM,*)                           sdrvb365
C                     .                                                 sdrvb366
C                     .                                                 sdrvb367
C                     Calculate values of DFDY                          sdrvb368
C                     .                                                 sdrvb369
C                     .                                                 sdrvb370
C                   END                                                 sdrvb371
C             Here Y is a vector of length at least N.  The actual      sdrvb372
C             length of Y is determined by the user's declaration in thesdrvb373
C             program which calls SDRVB3.  Thus the dimensioning of Y insdrvb374
C             JACOBN, while required by FORTRAN convention, does not    sdrvb375
C             actually allocate any storage.  When this subroutine is   sdrvb376
C             called, the first N components of Y are intermediate      sdrvb377
C             approximations to the solution components.  The user      sdrvb378
C             should not alter these values.  If the system is not      sdrvb379
C             banded (MITER=1), the partials of the I-th equation with  sdrvb380
C             respect to the J-th dependent function are to be stored insdrvb381
C             DFDY(I,J).  Thus partials of the I-th equation are stored sdrvb382
C             in the I-th row of DFDY.  If the system is banded         sdrvb383
C             (MITER=4), then the partials of the I-th equation with    sdrvb384
C             respect to Y(J) are to be stored in DFDY(ki,J), where      sdrvb385
C             ki=I-J+MU+1.                                               sdrvb386
C                                                                       sdrvb387
C    FA     = A subroutine supplied by the user if IMPL is 1 or 2, and  sdrvb388
C             MITER is not 3.  If so, the name must be declared EXTERNALsdrvb389
C             in the user's calling program.  This subroutine computes  sdrvb390
C             the array A, where A*dY(I)/dT = F(Y(I),T).                sdrvb391
C             There are two cases:                                      sdrvb392
C                                                                       sdrvb393
C               IMPL=1.                                                 sdrvb394
C               Subroutine FA is of the form:                           sdrvb395
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)     sdrvb396
C                   REAL Y(*), A(MATDIM,*)                              sdrvb397
C                     .                                                 sdrvb398
C                     .                                                 sdrvb399
C                     Calculate ALL values of A                         sdrvb400
C                     .                                                 sdrvb401
C                     .                                                 sdrvb402
C                   END                                                 sdrvb403
C               In this case A is assumed to be a nonsingular matrix,   sdrvb404
C               with the same structure as DFDY (see JACOBN description sdrvb405
C               above).  Programming considerations prevent complete    sdrvb406
C               generality.  If MITER is 1 or 2, A is assumed to be fullsdrvb407
C               and the user must compute and store all values of       sdrvb408
C               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumedsdrvb409
C               to be banded with lower and upper half bandwidth ML and sdrvb410
C               MU.  The left hand side of the I-th equation is a linearsdrvb411
C               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,        sdrvb412
C               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in thesdrvb413
C               I-th equation, the coefficient of dY(J)/dT is to be     sdrvb414
C               stored in A(ki,J), where ki=I-J+MU+1.                     sdrvb415
C               NOTE: The array A will be altered between calls to FA.  sdrvb416
C                                                                       sdrvb417
C               IMPL=2.                                                 sdrvb418
C               Subroutine FA is of the form:                           sdrvb419
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)     sdrvb420
C                   REAL Y(*), A(*)                                     sdrvb421
C                     .                                                 sdrvb422
C                     .                                                 sdrvb423
C                     Calculate non-zero values of A(1),...,A(NDE)      sdrvb424
C                     .                                                 sdrvb425
C                     .                                                 sdrvb426
C                   END                                                 sdrvb427
C               In this case it is assumed that the system is ordered bysdrvb428
C               the user so that the differential equations appear      sdrvb429
C               first, and the algebraic equations appear last.  The    sdrvb430
C               algebraic equations must be written in the form:        sdrvb431
C               0 = F(Y(I),T).  When using this option it is up to the  sdrvb432
C               user to provide initial values for the Y(I) that satisfysdrvb433
C               the algebraic equations as well as possible.  It is     sdrvb434
C               further assumed that A is a vector of length NDE.  All  sdrvb435
C               of the components of A, which may depend on T, Y(I),    sdrvb436
C               etc., must be set by the user to non-zero values.       sdrvb437
C             Here Y is a vector of length at least N.  The actual      sdrvb438
C             length of Y is determined by the user's declaration in thesdrvb439
C             program which calls SDRVB3.  Thus the dimensioning of Y insdrvb440
C             FA, while required by FORTRAN convention, does not        sdrvb441
C             actually allocate any storage.  When this subroutine is   sdrvb442
C             called, the first N components of Y are intermediate      sdrvb443
C             approximations to the solution components.  The user      sdrvb444
C             should not alter these values.  FA is always called       sdrvb445
C             immediately after calling F, with the same values of T    sdrvb446
C             and Y.                                                    sdrvb447
C                                                                       sdrvb448
C    NDE    = (Input) The number of differential equations.  This is    sdrvb449
C             required only for IMPL = 2, with NDE .LT. N.              sdrvb450
C                                                                       sdrvb451
C    MXSTEP = (Input) The maximum number of internal steps allowed on   sdrvb452
C             one call to SDRVB3.                                       sdrvb453
C                                                                       sdrvb454
C    G      = A real FORTRAN function supplied by the user              sdrvb455
C             if NROOT is not 0.  In this case, the name must be        sdrvb456
C             declared EXTERNAL in the user's calling program.  G is    sdrvb457
C             repeatedly called with different values of IROOT to obtainsdrvb458
C             the value of each of the NROOT equations for which a root sdrvb459
C             is desired.  G is of the form:                            sdrvb460
C                   REAL FUNCTION G (N, T, Y, IROOT)                    sdrvb461
C                   REAL Y(*)                                           sdrvb462
C                   GO TO (10, ...), IROOT                              sdrvb463
C              10   G = ...                                             sdrvb464
C                     .                                                 sdrvb465
C                     .                                                 sdrvb466
C                   END                                                 sdrvb467
C             Here, Y is a vector of length at least N, whose first N   sdrvb468
C             components are the solution components at the point T.    sdrvb469
C             The user should not alter these values.  The actual lengthsdrvb470
C             of Y is determined by the user's declaration in the       sdrvb471
C             program which calls SDRVB3.  Thus the dimensioning of Y insdrvb472
C             G, while required by FORTRAN convention, does not actuallysdrvb473
C             allocate any storage.                                     sdrvb474
C                                                                       sdrvb475
C***LONG DESCRIPTION                                                    sdrvb476
C                                                                       sdrvb477
C  III.  OTHER COMMUNICATION TO THE USER  ..............................sdrvb478
C                                                                       sdrvb479
C    A. The solver communicates to the user through the parameters      sdrvb480
C       above.  In addition it writes diagnostic messages through the   sdrvb481
C       standard error handling program XERROR.  That program will      sdrvb482
C       terminate the user's run if it detects a probable problem setup sdrvb483
C       error, e.g., insufficient storage allocated by the user for the sdrvb484
C       WORK array.  Messages are written on the standard error message sdrvb485
C       file.  At installations which have this error handling package  sdrvb486
C       the user should determine the standard error handling file from sdrvb487
C       the local documentation.  Otherwise the short but serviceable   sdrvb488
C       routine, XERROR, available with this package, can be used.  Thatsdrvb489
C       program writes on logical unit 6 to transmit messages.  A       sdrvb490
C       complete description of XERROR is given in the Sandia           sdrvb491
C       Laboratories report SAND78-1189 by R. E. Jones.  Following is a sdrvb492
C       list of possible errors.  Unless otherwise noted, all messages  sdrvb493
C       come from SDRVB3:                                               sdrvb494
C                                                                       sdrvb495
C        No.  Type         Message                                      sdrvb496
C        ---  ----         -------                                      sdrvb497
C         1   Fatal        From SDRVB2: The integration od flag has an  sdrvb498
C                          illegal value.                               sdrvb499
C         2   Warning      The output point is inconsistent with the    sdrvb500
C                          value of NTASK and T.                        sdrvb501
C         3   Warning      Number of steps to reach TOUT exceeds MXSTEP.sdrvb502
C         4   Recoverable  Requested accuracy is too stringent.         sdrvb503
C         5   Warning      Step size is below the roundoff level.       sdrvb504
C         6   Fatal        EPS is less than zero.                       sdrvb505
C         7   Fatal        N is not positive.                           sdrvb506
C         8   Fatal        Insufficient work space provided.            sdrvb507
C         9   Fatal        Improper value for MINT, MITER and/or IMPL.  sdrvb508
C        10   Fatal        The IWORK array is too small.                sdrvb509
C        11   Fatal        The step size has gone to zero.              sdrvb510
C        12   Fatal        Excessive amount of work.                    sdrvb511
C        13   Fatal        For IMPL=1 or 2, the matrix A is singular.   sdrvb512
C        14   Fatal        MXORD is not positive.                       sdrvb513
C        15   Fatal        From SDRVB1: N is greater than 200.          sdrvb514
C        16   Fatal        From SDRVB1: The WORK array is too small.    sdrvb515
C                                                                       sdrvb516
C    B. The first three elements of WORK and the first five elements of sdrvb517
C       IWORK will contain the following statistical data:              sdrvb518
C         AVGH     The average step size used.                          sdrvb519
C         HUSED    The step size last used (successfully).              sdrvb520
C         AVGORD   The average order used.                              sdrvb521
C         IMXERR   The index of the element of the solution vector that sdrvb522
C                  contributed most to the last error test.             sdrvb523
C         NQUSED   The order last used (successfully).                  sdrvb524
C         NSTep    The number of steps taken.                           sdrvb525
C         NFE      The number of evaluations of the right hand side.    sdrvb526
C         NJE      The number of evaluations of the Jacobian matrix.    sdrvb527
C                                                                       sdrvb528
C  IV.  REMARKS  .......................................................sdrvb529
C                                                                       sdrvb530
C    A. Other routines used:                                            sdrvb531
C         SDNTPB, SDZROB, SDSTPB, SDNTLB, SDPSTB, SDCORB, SDCSTB,       sdrvb532
C         SDPSCB, and SDSCLB;                                           sdrvb533
C         SGEFA, SGESL, SGBFA, SGBSL, and SNRM2 (from LINPACK)          sdrvb534
C         R1MACH (from the Bell Laboratories Machine Constants Package) sdrvb535
C         XERROR (from the SLATEC Common Math Library)                  sdrvb536
C       The last seven routines above, not having been written by the   sdrvb537
C       present authors, are not explicitly part of this package.       sdrvb538
C                                                                       sdrvb539
C    B. On any return from SDRVB3 all information necessary to continue sdrvb540
C       the calculation is contained in the call sequence parameters,   sdrvb541
C       including the work arrays.  Thus it is possible to suspend one  sdrvb542
C       problem, integrate another, and then return to the first.       sdrvb543
C                                                                       sdrvb544
C    C. There are user-written routines which are only required by      sdrvb545
C       SDRVB3 when certain parameters are set.  Thus a message warning sdrvb546
C       of unsatisfied externals may be issued during the load or link  sdrvb547
C       phase.  This message should never refer to F.  This message can sdrvb548
C       be ignored if: it refers to JACOBN and MITER is not 1 or 4, or  sdrvb549
C       it refers to FA and IMPL is 0 or MITER is 3, or it refers to    sdrvb550
C       USERS and MITER is not 3, or it refers to G and NROOT is 0.     sdrvb551
C                                                                       sdrvb552
C    D. If this package is to be used in an overlay situation, the user sdrvb553
C       must declare in the primary overlay the variables in the call   sdrvb554
C       sequence to SDRVB3.                                             sdrvb555
C                                                                       sdrvb556
C    E. Changing parameters during an integration.                      sdrvb557
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may  sdrvb558
C       be altered by the user between calls to SDRVB3.  For example, ifsdrvb559
C       too much accuracy has been requested (the program returns with  sdrvb560
C       NSTATE = 4 and an increased value of EPS) the user may wish to  sdrvb561
C       increase EPS further.  In general, prudence is necessary when   sdrvb562
C       making changes in parameters since such changes are not         sdrvb563
C       implemented until the next integration step, which is not       sdrvb564
C       necessarily the next call to SDRVB3.  This can happen if the    sdrvb565
C       program has already integrated to a point which is beyond the   sdrvb566
C       new point TOUT.                                                 sdrvb567
C                                                                       sdrvb568
C    F. As the price for complete control of matrix algebra, the SDRVB3 sdrvb569
C       USERS option puts all responsibility for Jacobian matrix        sdrvb570
C       evaluation on the user.  It is often useful to approximate      sdrvb571
C       numerically all or part of the Jacobian matrix.  However this   sdrvb572
C       must be done carefully.  The FORTRAN sequence below illustrates sdrvb573
C       the od we recommend.  It can be inserted directly into          sdrvb574
C       subroutine USERS to approximate Jacobian elements in rows I1    sdrvb575
C       to I2 and columns J1 to J2.                                     sdrvb576
C              REAL DFDY(N,N), EPSJ, H, R, R0, R1MACH,                  sdrvb577
C             8     SAVE1(N), SAVE2(N), SUM, T, UROUND, Y(N), YJ        sdrvb578
C              SUM = 0.E0                                               sdrvb579
C              DO 10 I = 1,N                                            sdrvb580
C         10     SUM = SUM + SAVE2(I)**2                                sdrvb581
C              UROUND = R1MACH(4)                                       sdrvb582
C              R0 = ABS(H)*SQRT(SUM)*1000.E0*UROUND                     sdrvb583
C              EPSJ = SQRT(UROUND)                                      sdrvb584
C              DO 30 J = J1,J2                                          sdrvb585
C                R = MAX(R0, EPSJ*ABS(Y(J)))                            sdrvb586
C                IF ((Y(J) + R) .EQ. Y(J)) R = UROUND                   sdrvb587
C                YJ = Y(J)                                              sdrvb588
C                Y(J) = Y(J) + R                                        sdrvb589
C                CALL F (N, T, Y, SAVE1)                                sdrvb590
C                Y(J) = YJ                                              sdrvb591
C                DO 20 I = I1,I2                                        sdrvb592
C         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R                  sdrvb593
C         30     CONTINUE                                               sdrvb594
C       Many problems give rise to structured sparse Jacobians, e.g.,   sdrvb595
C       block banded.  It is possible to approximate them with fewer    sdrvb596
C       function evaluations than the above procedure uses; see Curtis, sdrvb597
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,       sdrvb598
C       pp. 117-119.                                                    sdrvb599
C                                                                       sdrvb600
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN       sdrvb601
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.sdrvb602
C***ROUTINES CALLED  SDSTPB,SDNTPB,SDZROB,SGEFA,SGESL,SGBFA,SGBSL,SNRM2,sdrvb603
C                    R1MACH,XERROR                                      sdrvb604
C   End prologue sdrvb3                                                 sdrvb605
      external F, JACOBN, FA, G                                         sdrvb606
      REAL :: AE, BIG, CMPR, EPS, EWT(*), G, GLAST, H, HMAX,            sdrvb607
     8     HSIGN, RE, RELOUT, SIZE, SNRM2, SUM, t, TLAST, tOut,         sdrvb608
     8     TROOT, UROUND, WORK(*), Y(*)                                 sdrvb609
      INTEGER :: IWORK(*)                                               sdrvb610
      LOGICAL :: CONVRG                                                 sdrvb611
      CHARACTER MSG*205                                                 sdrvb612
      DATA IPRNT /0/, RELOUT /.05E0/                                    sdrvb613
      DATA IAVGH, IHUSED, IAVGRD, IEL /1, 2, 3, 4/                      sdrvb614
      DATA IH, IHMAX, IHOLD, IHSIGN, IRC /160, 161, 162, 163, 164/      sdrvb615
      DATA IRMAX, IT, ITOUT, ITQ, ITREND /165, 166, 167, 168, 204/      sdrvb616
      DATA IYH /205/                                                    sdrvb617
      DATA INDMXR, INQUSD, iNStep, INFE, INJE /1, 2, 3, 4, 5/           sdrvb618
      DATA INROOT, ICNVRG, IJROOT, IJTASK, IMNTLD /6, 7, 8, 9, 10/      sdrvb619
      DATA IMTRLD, INQ, INRTLD, INDTRT, INWAIT /11, 12, 13, 14, 15/     sdrvb620
      DATA IMNT, IMTRSV, IMTR, IMXRDS, IMXORD /16, 17, 18, 19, 20/      sdrvb621
      DATA INDPVT /21/                                                  sdrvb622
C***FIRST EXECUTABLE STATEMENT  SDRVB3                                  sdrvb623
C     R1MACH (3) = 1.110223024625157E-16                                sdrvb624
C     UROUND = R1MACH (4)    
! *****
      write(unit = 8, fmt = "(a4, i4, a5, 1pe8.1, a10, i2, a8, e8.1, 
     1  2(a9, i2), a7, e8.1, a8, i2, a9, i2, a8, i2, 2(a6, i2), a7, 
     3  i2, a5, i4)") " N =", N, ", t =", t, ", NState =", NState, ", 
     4  tOut =", tOut, ", NTask =", NTask, ", NRoot =", NRoot, ", 
     4  EPS =", EPS, ", mint =", mint, ", miter =", miter, ", impl =", 
     5  impl, ", ML =", ML, ", MU =", MU, ", MxOrd =", MxOrd, ", 
     6  NDE =", NDE
! *****
!      UROUND = 2.220446049250313E-16   ! 64 bit processor
      uround = 2.384185791e-07         ! 32 bit processor
      IF (NROOT .NE. 0) THEN                                            sdrvb625
C       AE = R1MACH(1)                                                  sdrvb626
        RE = UROUND                                                     sdrvb627
      else  
      END IF                                                            sdrvb628
      IF (EPS .LT. 0.E0) THEN                                           sdrvb629
        WRITE(MSG, '(''SDRVB36FE Illegal input.  EPS,'', E16.8,         sdrvb630
     8  '', is negative.'')') EPS                                       sdrvb631
C       CALL XERROR(MSG, 60, 6, 2)                                      sdrvb632
        RETURN                                                          sdrvb633
      else  
      END IF                                                            sdrvb634
      IF (N .LE. 0) THEN                                                sdrvb635
        WRITE(MSG, '(''SDRVB37FE Illegal input.  Number of equations,'',sdrvb636
     8  I8, '', is not positive.'')') N                                 sdrvb637
C       CALL XERROR(MSG, 72, 7, 2)                                      sdrvb638
        RETURN                                                          sdrvb639
      END IF                                                            sdrvb640
      IF (MXORD .LE. 0) THEN                                            sdrvb641
        WRITE(MSG, '(''SDRVB314FE Illegal input.  Maximum order,'', I8, sdrvb642
     8  '', is not positive.'')') MXORD                                 sdrvb643
C       CALL XERROR(MSG, 67, 14, 2)                                     sdrvb644
        RETURN                                                          sdrvb645
      END IF                                                            sdrvb646
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.        sdrvb647
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))              sdrvb648
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.                      sdrvb649
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.          sdrvb650
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.        sdrvb651
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1)) THEN                           sdrvb652
        WRITE(MSG, '(''SDRVB39FE Illegal input.  Improper value for '', sdrvb653
     8  ''MINT, MITER and/or IMPL.'')')                                 sdrvb654
C       CALL XERROR(MSG, 69, 9, 2)                                      sdrvb655
        RETURN                                                          sdrvb656
      END IF                                                            sdrvb657
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN                          sdrvb658
        LIWCHK = INDPVT - 1                                             sdrvb659
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.    sdrvb660
     8  MITER .EQ. 5) THEN                                              sdrvb661
        LIWCHK = INDPVT + N - 1                                         sdrvb662
      END IF                                                            sdrvb663
      IF (LENIW .LT. LIWCHK) THEN                                       sdrvb664
        WRITE(MSG, '(''SDRVB310FE Illegal input.  Insufficient '',      sdrvb665
     8  ''storage allocated for the IWORK array.  Based on the '')')    sdrvb666
        WRITE(MSG(94:), '(''value of the input parameters involved, '', sdrvb667
     8  ''the required storage is'', I8)') LIWCHK                       sdrvb668
C       CALL XERROR(MSG, 164, 10, 2)                                    sdrvb669
        RETURN                                                          sdrvb670
      else  
      end if                                                            sdrvb671
C                                                Allocate the WORK arraysdrvb672
C                                         IYH is the index of YH in WORKsdrvb673
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN                            sdrvb674
        MAXORD = MIN(MXORD, 12)                                         sdrvb675
      ELSE IF (MINT .EQ. 2) THEN                                        sdrvb676
        MAXORD = MIN(MXORD, 5)                                          sdrvb677
      else 
      end if                                                             sdrvb678
      IDFDY = IYH + (MaxOrd + 1)*N                                      sdrvb679
C                                            IDFDY is the index of DFDY sdrvb680
C                                                                       sdrvb681
      if (MITER .eq. 0 .or. MITER .eq. 3) then                          sdrvb682
        IYWT = IDFDY                                                    sdrvb683
      ELSE IF (MITER .eq. 1 .or. MITER .eq. 2)  then                    sdrvb684
        IYWT = IDFDY + N*N                                              sdrvb685
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN                    sdrvb686
        IYWT = IDFDY + (2*ML + MU + 1)*N                                sdrvb687
      END IF                                                            sdrvb688
! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C                                               IYWT is the index of YWTsdrvb689
      ISAVE1 = IYWT + N                                                 sdrvb690
C                                           ISAVE1 is the index of SAVE1sdrvb691
      ISAVE2 = ISAVE1 + N                                               sdrvb692
C                                           ISAVE2 is the index of SAVE2sdrvb693
      IGNOW = ISAVE2 + N                                                sdrvb694
C                                             IGNOW is the index of GNOWsdrvb695
      ITROOT = IGNOW + NRoot                                            sdrvb696
C                                           ITROOT is the index of TROOTsdrvb697
      IA = ITRoot + NRoot                                               sdrvb698
C                                                   IA is the index of Asdrvb699
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN                           sdrvb700
        LENCHK = IA - 1                                                 sdrvb701
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN sdrvb702
        LENCHK = IA - 1 + N*N                                           sdrvb703
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN sdrvb704
        LENCHK = IA - 1 + (2*ML + MU + 1)*N                             sdrvb705
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN                     sdrvb706
        LENCHK = IA - 1 + N                                             sdrvb707
      end if                                                            sdrvb708
      IF (LENW .LT. LENCHK) THEN                                        sdrvb709
        WRITE(MSG, '(''SDRVB38FE Illegal input.  Insufficient '',       sdrvb710
     8  ''storage allocated for the WORK array.  Based on the '')')     sdrvb711
        WRITE(MSG(92:), '(''value of the input parameters involved, '', sdrvb712
     8  ''the required storage is'', I8)') LENCHK                       sdrvb713
C       CALL XERROR(MSG, 162, 8, 2)                                     sdrvb714
        RETURN                                                          sdrvb715
      end if                                                            sdrvb716
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN                          sdrvb717
        MATDIM = 1                                                      sdrvb718
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN                     sdrvb719
        MATDIM = N                                                      sdrvb720
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN                     sdrvb721
        MATDIM = 2*ML + MU + 1                                          sdrvb722
      end if                                                            sdrvb723
! *****
!      write(unit = 9, fmt = *) (" i = ", i, " Y = ", Y(i), i = 1, 10)
! *****
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN                            sdrvb724
        NDECOM = N                                                      sdrvb725
      ELSE IF (IMPL .EQ. 2) THEN                                        sdrvb726
        NDECOM = NDE                                                    sdrvb727
      end if                                                            sdrvb728
      IF (NSTATE .EQ. 1) THEN                                           sdrvb729
C                                                Initialize parameters  sdrvb730
        if(t .eq. tOut) return                                          sdrvb731
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN                          sdrvb732
          IWORK(IMXORD) = MIN(MXORD, 12)                                sdrvb733
        ELSE IF (MINT .EQ. 2) THEN                                      sdrvb734
          IWORK(IMXORD) = MIN(MXORD, 5)                                 sdrvb735
        END IF                                                          sdrvb736
        IWORK(IMXRDS) = MXORD                                           sdrvb737
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN                          sdrvb738
          IWORK(IMNT) = MINT                                            sdrvb739
          IWORK(IMTR) = MITER                                           sdrvb740
          IWORK(IMNTLD) = MINT                                          sdrvb741
          IWORK(IMTRLD) = MITER                                         sdrvb742
        ELSE IF (MINT .EQ. 3) THEN                                      sdrvb743
          IWORK(IMNT) = 1                                               sdrvb744
          IWORK(IMTR) = 0                                               sdrvb745
          IWORK(IMNTLD) = IWORK(IMNT)                                   sdrvb746
          IWORK(IMTRLD) = IWORK(IMTR)                                   sdrvb747
          IWORK(IMTRSV) = MITER                                         sdrvb748
        end if                                                          sdrvb749
        WORK(IHMAX) = HMAX                                              sdrvb750
        H = (TOUT - t)*(1.0E+00 - 4.0E+00*UROUND)                       sdrvb751
        H = SIGN(MIN(ABS(H), HMAX), H)                                  sdrvb752
! *****
   !     write(unit = *, fmt = *) " IMNT = ", IMNT, " IMTR = ", IMTR, 
   !  1    " IMNTLD = ", IMNTLD, " IMTRLD = ", IMTRLD,  
   !  2    " IWORK (10 - 18) = ", (i, IWORK(i), i = 10, 18)
! *****
        WORK(IH) = H                                                    sdrvb753
        HSIGN = SIGN(1.0E+00, H)                                        sdrvb754
        WORK(IHSIGN) = HSIGN                                            sdrvb755
        IWORK(IJTASK) = 0                                               sdrvb756
        WORK(IAVGH) = 0.0E+00                                           sdrvb757
        WORK(IAVGRD) = 0.0E+00                                          sdrvb758
        IWORK(INQUSD) = 0                                               sdrvb759
        IWORK(iNStep) = 0                                               sdrvb760
        IWORK(INFE) = 0                                                 sdrvb761
        IWORK(INJE) = 0                                                 sdrvb762
        WORK(IT) = t                                                    sdrvb763
        IWORK(ICNVRG) = 0                                               sdrvb764
C                                                 Set initial conditionssdrvb765
        do 30 i = 1, N                                                  sdrvb766
          JYH = I + IYH - 1                                             sdrvb767
 30       WORK(JYH) = Y(I)                                              sdrvb768
! *****
!      write(unit = *, fmt = *) " sdrvb768: ", " JYH =", JYH, " IYH =", 
!     1  IYH,  " Y(i=1-9) =", (Y(i), i = 1, 9)
! *****
        go to 180                                                       sdrvb769
      end if                                                            sdrvb770
C                                             On a continuation, check  sdrvb771
C                                             that output points have   sdrvb772
C                                             been or will be overtaken.sdrvb773
      IF (IWORK(ICNVRG) .EQ. 1) THEN                                    sdrvb774
        CONVRG = .TRUE.                                                 sdrvb775
      ELSE                                                              sdrvb776
        CONVRG = .FALSE.                                                sdrvb777
      END IF                                                            sdrvb778
      t = WORK(IT)                                                      sdrvb779
! *****
      write(unit = *, fmt = "(a10, a8, l2, a10, i2, a5, 1pe8.2)")
     1  " sdrvb779:", " Convrg =", Convrg, ", NState =", NState,
     2  ", t =", t
! *****
      H = WORK(IH)                                                      sdrvb780
      HSIGN = WORK(IHSIGN)                                              sdrvb781
      GO TO 180                                                         sdrvb769
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180                               sdrvb782
C                                                                       sdrvb783
C                                   IWORK(IJROOT) flags unreported      sdrvb784
C                                   roots, and is set to the value of   sdrvb785
C                                   NTASK when a root was last selected.sdrvb786
C                                   It is set to zero when all roots    sdrvb787
C                                   have been reported.  IWORK(INROOT)  sdrvb788
C                                   contains the index and WORK(ITOUT)  sdrvb789
C                                   contains the value of the root last sdrvb790
C                                   selected to be reported.            sdrvb791
C                                   IWORK(INRTLD) contains the value of sdrvb792
C                                   NROOT and IWORK(INDTRT) contains    sdrvb793
C                                   the value of ITROOT when the array  sdrvb794
C                                   of roots was last calculated.       sdrvb795
      IF(NROOT .NE. 0) THEN                                             sdrvb796
        JROOT = IWORK(IJROOT)                                           sdrvb797
        IF (JROOT .GT. 0) THEN                                          sdrvb798
C                                      TOUT has just been reported.     sdrvb799
C                                      If TROOT .LE. TOUT, report TROOT.sdrvb800
          IF (NSTATE .NE. 5) THEN                                       sdrvb801
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN                 sdrvb802
              TROOT = WORK(ITOUT)                                       sdrvb803
              CALL SDNTPB(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y) sdrvb804
              T = TROOT                                                 sdrvb805
              NSTATE = 5                                                sdrvb806
              GO TO 580                                                 sdrvb807
            END IF                                                      sdrvb808
C                                         A root has just been reported.sdrvb809
C                                         Select the next root.         sdrvb810
          ELSE                                                          sdrvb811
            TROOT = T                                                   sdrvb812
            IROOT = 0                                                   sdrvb813
            DO 50 I = 1,IWORK(INRTLD)                                   sdrvb814
              JTROOT = IWORK(INDTRT) + I - 1                            sdrvb815
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN             sdrvb816
C                                                                       sdrvb817
C                                              Check for multiple roots.sdrvb818
C                                                                       sdrvb819
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.                 sdrvb820
     8          I .GT. IWORK(INROOT)) THEN                              sdrvb821
                  IROOT = I                                             sdrvb822
                  TROOT = WORK(JTROOT)                                  sdrvb823
                  GO TO 60                                              sdrvb824
                END IF                                                  sdrvb825
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN     sdrvb826
                  IROOT = I                                             sdrvb827
                  TROOT = WORK(JTROOT)                                  sdrvb828
                END IF                                                  sdrvb829
              END IF                                                    sdrvb830
 50           CONTINUE                                                  sdrvb831
 60         IWORK(INROOT) = IROOT                                       sdrvb832
            WORK(ITOUT) = TROOT                                         sdrvb833
            IWORK(IJROOT) = NTASK                                       sdrvb834
            IF (NTASK .EQ. 1) THEN                                      sdrvb835
              IF (IROOT .EQ. 0) THEN                                    sdrvb836
                IWORK(IJROOT) = 0                                       sdrvb837
              ELSE                                                      sdrvb838
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN                   sdrvb839
                  CALL SDNTPB(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)sdrvb840
                  NSTATE = 5                                            sdrvb841
                  T = TROOT                                             sdrvb842
                  GO TO 580                                             sdrvb843
                END IF                                                  sdrvb844
              END IF                                                    sdrvb845
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN               sdrvb846
C                                                                       sdrvb847
C                                     If there are no more roots, or thesdrvb848
C                                     user has altered TOUT to be less  sdrvb849
C                                     than a root, set IJROOT to zero.  sdrvb850
C                                                                       sdrvb851
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN sdrvb852
                IWORK(IJROOT) = 0                                       sdrvb853
              ELSE                                                      sdrvb854
                CALL SDNTPB(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)sdrvb855
                NSTATE = 5                                              sdrvb856
                T = TROOT                                               sdrvb857
                GO TO 580                                               sdrvb858
              END IF                                                    sdrvb859
            END IF                                                      sdrvb860
          END IF                                                        sdrvb861
        END IF                                                          sdrvb862
      END IF                                                            sdrvb863
C                                                                       sdrvb864
      IF (NTASK .EQ. 1) THEN                                            sdrvb865
        NSTATE = 2                                                      sdrvb866
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN                               sdrvb867
          CALL SDNTPB (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)     sdrvb868
          T = TOUT                                                      sdrvb869
          GO TO 580                                                     sdrvb870
        END IF                                                          sdrvb871
      ELSE IF (NTASK .EQ. 2) THEN                                       sdrvb872
C                                                      Check if TOUT hassdrvb873
C                                                      been reset .LT. Tsdrvb874
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN                               sdrvb875
          WRITE(MSG, '(''SDRVB32WRN With NTASK='', I1, '' on input, '', sdrvb876
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',sdrvb877
     8    '' obtained by interpolation.'')') NTASK, T, TOUT             sdrvb878
C         CALL XERROR(MSG, 124, 2, 0)                                   sdrvb879
          CALL SDNTPB (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)     sdrvb880
          T = TOUT                                                      sdrvb881
          NSTATE = 2                                                    sdrvb882
          GO TO 580                                                     sdrvb883
        END IF                                                          sdrvb884
C                                   Determine if TOUT has been overtakensdrvb885
C                                                                       sdrvb886
        CMPR = MAX(ABS(T), ABS(TOUT)) + RELOUT*ABS(TOUT - T)            sdrvb887
!      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        IF (CMPR .EQ. MAX(ABS(T), ABS(TOUT))) THEN                      sdrvb888
          T = TOUT                                                      sdrvb889
          NSTATE = 2                                                    sdrvb890
          GO TO 560                                                     sdrvb891
        END IF                                                          sdrvb892
C                                             If there are no more rootssdrvb893
C                                             to report, report T.      sdrvb894
        IF (NSTATE .EQ. 5) THEN                                         sdrvb895
          NSTATE = 2                                                    sdrvb896
          GO TO 560                                                     sdrvb897
        END IF                                                          sdrvb898
        NSTATE = 2                                                      sdrvb899
C                                                       See if TOUT willsdrvb900
C                                                       be overtaken.   sdrvb901
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN                         sdrvb902
          H = TOUT - T                                                  sdrvb903
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND) sdrvb904
          WORK(IH) = H                                                  sdrvb905
          IF (H .EQ. 0.E0) GO TO 670                                    sdrvb906
          IWORK(IJTASK) = -1                                            sdrvb907
        END IF                                                          sdrvb908
      ELSE IF (NTASK .EQ. 3) THEN                                       sdrvb909
        NSTATE = 2                                                      sdrvb910
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN                               sdrvb911
          WRITE(MSG, '(''SDRVB32WRN With NTASK='', I1, '' on input, '', sdrvb912
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',sdrvb913
     8    '' obtained by interpolation.'')') NTASK, T, TOUT             sdrvb914
C         CALL XERROR(MSG, 124, 2, 0)                                   sdrvb915
          CALL SDNTPB (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)     sdrvb916
          T = TOUT                                                      sdrvb917
          GO TO 580                                                     sdrvb918
        END IF                                                          sdrvb919
        CMPR = MAX(ABS(T), ABS(TOUT)) + RELOUT*ABS(TOUT - T)            sdrvb920
        IF (CMPR .EQ. MAX(ABS(T), ABS(TOUT))) THEN                      sdrvb921
          T = TOUT                                                      sdrvb922
          GO TO 560                                                     sdrvb923
        END IF                                                          sdrvb924
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN                         sdrvb925
          H = TOUT - T                                                  sdrvb926
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND) sdrvb927
          WORK(IH) = H                                                  sdrvb928
          IF (H .EQ. 0.E0) GO TO 670                                    sdrvb929
          IWORK(IJTASK) = -1                                            sdrvb930
        END IF                                                          sdrvb931
      END IF                                                            sdrvb932
C                         Implement changes in MINT, MITER, and/or HMAX.sdrvb933
C                                                                       sdrvb934
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND. sdrvb935
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1      sdrvb936
      IF (HMAX .NE. WORK(IHMAX)) THEN                                   sdrvb937
        H = SIGN(MIN(ABS(H), HMAX), H)                                  sdrvb938
        IF (H .NE. WORK(IH)) THEN                                       sdrvb939
          IWORK(IJTASK) = -1                                            sdrvb940
          WORK(IH) = H                                                  sdrvb941
        END IF                                                          sdrvb942
        WORK(IHMAX) = HMAX                                              sdrvb943
      END IF                                                            sdrvb944
C                                                                       sdrvb945
 180  NStepL = IWORK(iNStep)                                            sdrvb946
      do 190 i = 1, N                                                    sdrvb947
        JYH = IYH + I - 1                                               sdrvb948
 190    Y(I) = WORK(JYH)                                                sdrvb949
! *****
!      write(unit = 9, fmt = *) " NRoot = ", NRoot, " N = ", N, 
!     1  ", i = ", i, ", Y = ", Y(i)
! *****
      if (NRoot .ne. 0) then                                            sdrvb950
        DO 200 I = 1, NROOT                                             sdrvb951
          JGNOW = IGNOW + I - 1                                         sdrvb952
 200      WORK(JGNOW) = G (N, T, Y, I)                                  sdrvb953
      else
      end if                                                            sdrvb954
      IF (IERROR .EQ. 1) THEN                                           sdrvb955
        DO 230 I = 1,N                                                  sdrvb956
          JYWT = I + IYWT - 1                                           sdrvb957
 230      WORK(JYWT) = 1.E0                                             sdrvb958
        GO TO 410                                                       sdrvb959
      ELSE IF (IERROR .EQ. 5) THEN                                      sdrvb960
        DO 250 I = 1,N                                                  sdrvb961
          JYWT = I + IYWT - 1                                           sdrvb962
 250      WORK(JYWT) = EWT(I)                                           sdrvb963
        GO TO 410                                                       sdrvb964
      end if                                                            sdrvb965
C                                       Reset YWT array.  Looping point.sdrvb966
 260  IF (IERROR .EQ. 2) THEN                                           sdrvb967
        dO 280 i = 1, N                                                  sdrvb968
! *****
      write(unit = 9, fmt = *) " sdrvb: 968"
      write(unit = 9, fmt = *) " i = ", i, ", Y = ", Y(i)
      write(unit = 9, fmt = *) " IJTASK = ", IJTASK, 
     1  " IWork(IJTASK) = ", IWork(IJTASK)
!         pause
      write(unit = 9, fmt = *) " sdrvb969:  Y(i) = ", Y(i)
! *****
          IF (Y(I) .EQ. 0.E0) GO TO 290                                 sdrvb969
          JYWT = I + IYWT - 1                                           sdrvb970
 280      WORK(JYWT) = ABS(Y(I))                                        sdrvb971
        GO TO 410                                                       sdrvb972
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (N, T, Y, WORK(ISAVE2))                                sdrvb974
          IWORK(INFE) = IWORK(INFE) + 1                                 sdrvb975
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN                      sdrvb976
            IFLAG = 0                                                   sdrvb977
C           CALL USERS(Y,WORK(IYH),WORK(IYWT),WORK(ISAVE1),WORK(ISAVE2),sdrvb978
C    8                 T, H, WORK(IEL), IMPL, N, NDECOM, IFLAG)         sdrvb979
          ELSE IF (IMPL .EQ. 1) THEN                                    sdrvb980
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN                    sdrvb981
              CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)       sdrvb982
              CALL SGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)     sdrvb983
              IF (INFO .NE. 0) GO TO 690                                sdrvb984
              CALL SGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)sdrvb985
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN               sdrvb986
              JAML = IA + ML                                            sdrvb987
              CALL FA (N, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)     sdrvb988
              CALL SGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)   sdrvb989
              IF (INFO .NE. 0) GO TO 690                                sdrvb990
              CALL SGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),   sdrvb991
     8                    WORK(ISAVE2), 0)                              sdrvb992
            END IF                                                      sdrvb993
          ELSE IF (IMPL .EQ. 2) THEN                                    sdrvb994
            CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)         sdrvb995
            DO 340 I = 1,NDECOM                                         sdrvb996
              JA = I + IA - 1                                           sdrvb997
              JSAVE2 = I + ISAVE2 - 1                                   sdrvb998
              IF(WORK(JA) .EQ. 0.E0) GO TO 690                          sdrvb999
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)                      sdrv1000
          END IF                                                        sdrv1001
        END IF                                                          sdrv1002
        DO 360 J = I,N                                                  sdrv1003
          JYWT = J + IYWT - 1                                           sdrv1004
          IF (Y(J) .NE. 0.E0) THEN                                      sdrv1005
            WORK(JYWT) = ABS(Y(J))                                      sdrv1006
          ELSE                                                          sdrv1007
            IF (IWORK(IJTASK) .EQ. 0) THEN                              sdrv1008
              JSAVE2 = J + ISAVE2 - 1                                   sdrv1009
              WORK(JYWT) = ABS(H*WORK(JSAVE2))                          sdrv1010
            ELSE                                                        sdrv1011
              JHYP = J + IYH + N - 1                                    sdrv1012
              WORK(JYWT) = ABS(WORK(JHYP))                              sdrv1013
            END IF                                                      sdrv1014
          END IF                                                        sdrv1015
          IF (WORK(JYWT) .EQ. 0.E0) WORK(JYWT) = UROUND                 sdrv1016
 360      CONTINUE                                                      sdrv1017
      ELSE IF (IERROR .EQ. 3) THEN                                      sdrv1018
        DO 380 I = 1, N                                                  sdrv1019
          JYWT = I + IYWT - 1                                           sdrv1020
 380      WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))                           sdrv1021
      else if (iError .eq. 4) then                                      sdrv1022
        do 400 i = 1, N                                                 sdrv1023
          JYWT = i + IYWT - 1                                           sdrv1024
 400      Work(JYWT) = max(EWT(i), abs(Y(i)))                           sdrv1025
      end if                                                            sdrv1026
C                                                                       sdrv1027
 410  do 420 i = 1, N                                                    sdrv1028
        JYWT = i + IYWT - 1                                             sdrv1029
        JSAVE2 = i + ISAVE2 - 1                                         sdrv1030
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)                                  sdrv1031
      SUM = SNRM2(N, WORK(ISAVE2), 1)/SQRT(REAL(N))                     sdrv1032
      CMPR = SUM + EPS                                                  sdrv1033
      IF (CMPR .EQ. SUM) THEN                                           sdrv1034
        EPS = CMPR*UROUND*(1.0E+00 + 10.0E+00*UROUND)                   sdrv1035
        WRITE(MSG, '(''SDRVB34REC At T,'', E16.8, '', the requested '', sdrv1036
     8  ''accuracy, EPS, was not obtainable with the machine '',        sdrv1037
     8  ''precision.  EPS has been increased to'')') T                  sdrv1038
        WRITE(MSG(137:), '(E16.8)') EPS                                 sdrv1039
C       CALL XERROR(MSG, 152, 4, 1)                                     sdrv1040
        NSTATE = 4                                                      sdrv1041
        GO TO 560                                                       sdrv1042
      END IF                                                            sdrv1043
!        write(unit = *, fmt = *) " SDRVB3: 1043:  T=", T, " H=", H, 
!     1    " NTask=", NTask
!        pause
      CMPR = T + H                                                      sdrv1044
      IF (CMPR .NE. T) THEN                                             sdrv1045
        IPRNT = 0                                                       sdrv1046
      ELSE IF (IPRNT .EQ. 0) THEN                                       sdrv1047
        WRITE(MSG, '(''SDRVB35WRN At T,'', E16.8, '', the step size,'', sdrv1048
     8  E16.8, '', is smaller than the roundoff level of T.  '')') T, H sdrv1049
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',   sdrv1050
     8  ''change in the right hand side of the differential '',         sdrv1051
     8  ''equations.'')')                                               sdrv1052
C       CALL XERROR(MSG, 205, 5, 0)                                     sdrv1053
        IPRNT = 1                                                       sdrv1054
      END IF                                                            sdrv1055
      IF (NTASK.NE.2) THEN                                              sdrv1056
        if((IWORK(INStep) - NSTEPL) .GT. MXSTEP) THEN                   sdrv1057
          WRITE(MSG, '(''SDRVB33WRN At T,'', E16.8, '', '', I8,         sdrv1058
     8    '' steps have been taken without reaching TOUT,'', E16.8)')   sdrv1059
     8    T, MXSTEP, TOUT                                               sdrv1060
C         CALL XERROR(MSG, 103, 3, 0)                                   sdrv1061
          NSTATE = 3                                                    sdrv1062
          GO TO 560                                                     sdrv1063
        END IF                                                          sdrv1064
      END IF                                                            sdrv1065
C                                                                       sdrv1066
C     CALL SDSTPB (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,      sdrv1067
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND,             sdrv1068
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,        sdrv1069
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, ,        sdrv1070
C    8            DFDY, EL, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,          sdrv1071
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV) sdrv1072
C                                                                       sdrv1073
      CALL SDSTPB (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,       sdrv1074
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, N,   sdrv1075
     8            NDECOM, WORK(IYWT), UROUND, WORK(IAVGH), WORK(IAVGRD),sdrv1076
     8            WORK(IH), WORK(IHUSED), IWORK(IJTASK), IWORK(IMNTLD), sdrv1077
     8            IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),              sdrv1078
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH), sdrv1079
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IHOLD),sdrv1080
     8            IWORK(INDPVT), JSTATE, IWORK(INQ), IWORK(INWAIT),     sdrv1081
     8            WORK(IRC), WORK(IRMAX), WORK(ISAVE1), WORK(ISAVE2),   sdrv1082
     8            WORK(ITQ), WORK(ITREND), MINT, IWORK(IMTRSV),         sdrv1083
     8            IWORK(IMXRDS))                                        sdrv1084
      T = WORK(IT)                                                      sdrv1085
      H = WORK(IH)                                                      sdrv1086
      GO TO (470, 670, 680, 690), JSTATE                                sdrv1087
 470  IWORK(IJTASK) = 1                                                 sdrv1088
C                                 Determine if a root has been overtakensdrv1089
      CMPR = T + H                                                      sdrv1090
      IF (CMPR .NE. T .AND. NROOT .NE. 0) THEN                          sdrv1091
        IROOT = 0                                                       sdrv1092
        DO 500 I = 1,NROOT                                              sdrv1093
          JTROOT = ITROOT + I - 1                                       sdrv1094
          JGNOW = IGNOW + I - 1                                         sdrv1095
          GLAST = WORK(JGNOW)                                           sdrv1096
          WORK(JGNOW) = G (N, T, Y, I)                                  sdrv1097
          IF (GLAST*WORK(JGNOW) .GT. 0.E0) THEN                         sdrv1098
            WORK(JTROOT) = T + H                                        sdrv1099
          ELSE                                                          sdrv1100
            IF (WORK(JGNOW) .EQ. 0.E0) THEN                             sdrv1101
              WORK(JTROOT) = T                                          sdrv1102
              IROOT = I                                                 sdrv1103
            ELSE                                                        sdrv1104
              IF (GLAST .EQ. 0.E0) THEN                                 sdrv1105
                WORK(JTROOT) = T + H                                    sdrv1106
              ELSE                                                      sdrv1107
                TLAST = T - WORK(IHUSED)                                sdrv1108
                IROOT = I                                               sdrv1109
                TROOT = T                                               sdrv1110
                CALL SDZROB (AE, G, H, N, IWORK(INQ), IROOT, RE, T,     sdrv1111
     8                      WORK(IYH), UROUND,  TROOT, TLAST,           sdrv1112
     8                      WORK(JGNOW), GLAST,  WORK(ISAVE1))          sdrv1113
                WORK(JTROOT) = TROOT                                    sdrv1114
              END IF                                                    sdrv1115
            END IF                                                      sdrv1116
          END IF                                                        sdrv1117
 500      CONTINUE                                                      sdrv1118
        IF (IROOT .EQ. 0) THEN                                          sdrv1119
          IWORK(IJROOT) = 0                                             sdrv1120
C                                                  Select the first rootsdrv1121
        ELSE                                                            sdrv1122
          IWORK(IJROOT) = NTASK                                         sdrv1123
          IWORK(INRTLD) = NROOT                                         sdrv1124
          IWORK(INDTRT) = ITROOT                                        sdrv1125
          TROOT = T + H                                                 sdrv1126
          DO 510 I = 1,NROOT                                            sdrv1127
            JTROOT = ITROOT + I - 1                                     sdrv1128
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN               sdrv1129
              TROOT = WORK(JTROOT)                                      sdrv1130
              IROOT = I                                                 sdrv1131
            END IF                                                      sdrv1132
 510        CONTINUE                                                    sdrv1133
          IWORK(INROOT) = IROOT                                         sdrv1134
          WORK(ITOUT) = TROOT                                           sdrv1135
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN                         sdrv1136
            CALL SDNTPB (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)  sdrv1137
            NSTATE = 5                                                  sdrv1138
            T = TROOT                                                   sdrv1139
            GO TO 580                                                   sdrv1140
          END IF                                                        sdrv1141
        END IF                                                          sdrv1142
      END IF                                                            sdrv1143
C                               Test for NTASK condition to be satisfiedsdrv1144
      NSTATE = 2                                                        sdrv1145
      IF (NTASK .EQ. 1) THEN                                            sdrv1146
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260                          sdrv1147
        CALL SDNTPB (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)       sdrv1148
        T = TOUT                                                        sdrv1149
        GO TO 580                                                       sdrv1150
C                               TOUT is assumed to have been attained   sdrv1151
C                               exactly if T is within twenty roundoff  sdrv1152
C                               units of TOUT, relative to max(TOUT, T).sdrv1153
      ELSE IF (NTASK .EQ. 2) THEN                                       sdrv1154
        CMPR = MAX(ABS(T), ABS(TOUT)) + RELOUT*ABS(TOUT - T)            sdrv1155
        IF (CMPR .EQ. MAX(ABS(T), ABS(TOUT))) THEN                      sdrv1156
          T = TOUT                                                      sdrv1157
        ELSE                                                            sdrv1158
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN                       sdrv1159
            H = TOUT - T                                                sdrv1160
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND) sdrv1161
            WORK(IH) = H                                                sdrv1162
            IF (H .EQ. 0.E0) GO TO 670                                  sdrv1163
            IWORK(IJTASK) = -1                                          sdrv1164
          END IF                                                        sdrv1165
        END IF                                                          sdrv1166
      ELSE IF (NTASK .EQ. 3) THEN                                       sdrv1167
        CMPR = MAX(ABS(T), ABS(TOUT)) + RELOUT*ABS(TOUT - T)            sdrv1168
        IF (CMPR .EQ. MAX(ABS(T), ABS(TOUT))) THEN                      sdrv1169
          T = TOUT                                                      sdrv1170
        ELSE                                                            sdrv1171
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN                       sdrv1172
            H = TOUT - T                                                sdrv1173
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND) sdrv1174
            WORK(IH) = H                                                sdrv1175
            IF (H .EQ. 0.E0) GO TO 670                                  sdrv1176
            IWORK(IJTASK) = -1                                          sdrv1177
          END IF                                                        sdrv1178
          GO TO 260                                                     sdrv1179
        END IF                                                          sdrv1180
      END IF                                                            sdrv1181
C                                      All returns are made through thissdrv1182
C                                      section.  IMXERR is determined.  sdrv1183
 560  DO 570 I = 1,N                                                    sdrv1184
        JYH = I + IYH - 1                                               sdrv1185
 570    Y(I) = WORK(JYH)                                                sdrv1186
! *****
  !   write(unit = 9, fmt = *) " Line 1671:  N =", N, " IYH =", IYH, 
  !  1  " JYH =", JYH
  !   write(unit = 9, fmt = "(i4, 1pe10.3)") (i, Y(i), i = 1, 157)
 580  write(unit = *, fmt = *) "sdrv1187", " Convrg =", Convrg
      IF (CONVRG) THEN                                                  sdrv1187
        IWORK(ICNVRG) = 1                                               sdrv1188
      ELSE                                                              sdrv1189
        IWORK(ICNVRG) = 0                                               sdrv1190
      END IF                                                            sdrv1191
      IF (IWORK(IJTASK) .EQ. 0) RETURN                                  sdrv1192
      BIG = 0.0E+00                                                     sdrv1193
      IMXERR = 1                                                        sdrv1194
      IWORK(INDMXR) = IMXERR                                            sdrv1195
      DO  590 I = 1,N                                                   sdrv1196
C                                            SIZE = ABS(ERROR(I)/YWT(I))sdrv1197
        JYWT = I + IYWT - 1                                             sdrv1198
        JERROR = I + ISAVE1 - 1                                         sdrv1199
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))                             sdrv1200
        IF (BIG .LT. SIZE) THEN                                         sdrv1201
          BIG = SIZE                                                    sdrv1202
          IMXERR = I                                                    sdrv1203
          IWORK(INDMXR) = IMXERR                                        sdrv1204
        END IF                                                          sdrv1205
 590    CONTINUE                                                        sdrv1206
      RETURN                                                            sdrv1207
! **********
! END REVISION
! **********
C                                        Fatal errors are processed heresdrv1208
C                                                                       sdrv1209
 670  WRITE(MSG, '(''SDRVB311FE At T,'', E16.8, '', the attempted '',   sdrv1210
     8  ''step size has gone to zero.  Often this occurs if the '',     sdrv1211
     8  ''problem setup is incorrect.'')') T                            sdrv1212
C     CALL XERROR(MSG, 129, 11, 2)                                      sdrv1213
      RETURN                                                            sdrv1214
C                                                                       sdrv1215
 680  WRITE(MSG, '(''SDRVB312FE At T,'', E16.8, '', the step size has'',sdrv1216
     8  '' been reduced about 50 times without advancing the '')') T    sdrv1217
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',      sdrv1218
     8  ''problem setup is incorrect.'')')                              sdrv1219
C     CALL XERROR(MSG, 165, 12, 2)                                      sdrv1220
      RETURN                                                            sdrv1221
C                                                                       sdrv1222
 690  WRITE(MSG, '(''SDRVB313FE At T,'', E16.8, '', while solving'',    sdrv1223
     8  '' A*YDOT = F, A is singular.'')') T                            sdrv1224
C     CALL XERROR(MSG, 74, 13, 2)                                       sdrv1225
! *****
      write(unit = 9, fmt = *) " After return from sdrvb3"
      write(unit = 9, fmt = "(i4, 1pe10.3)") (i, Y(i), 
     1  i = 1, 2)
!      pause
! *****
      return                                                            sdrv1226
      end                                                               sdrv1227
!
! FChem & Rates go here
!
!

c
c
c
c These are now compiled separately and linked into the executable
c
c
c
c


      SUBROUTINE SDZROB (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)  sdzrob 3
C***BEGIN PROLOGUE  SDZROB                                              sdzrob 4
C***REFER TO  B3                                                        sdzrob 5
C     This is a special purpose version of ZEROIN, modified for use withsdzrob 6
C     the B1 package.                                                   sdzrob 7
C                                                                       sdzrob 8
C     Sandia Mathematical Program Library                               sdzrob 9
C     Mathematical Computing Services Division 5422                     sdzrob10
C     Sandia Laboratories                                               sdzrob11
C     P. O. Box 5800                                                    sdzrob12
C     Albuquerque, New Mexico  87115                                    sdzrob13
C     Control Data 6600 Version 4.5, 1 November 1971                    sdzrob14
C                                                                       sdzrob15
C     ABSTRACT                                                          sdzrob16
C        ZEROIN searches for a zero of a function F(N,T,Y,IROOT) betweensdzrob17
C        the given values B and C until the width of the interval       sdzrob18
C        (B,C) has collapsed to within a tolerance specified by         sdzrob19
C        the stopping criterion, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).       sdzrob20
C                                                                       sdzrob21
C     Description of parameters                                         sdzrob22
C        F     - Name of the external function, which returns a         sdzrob23
C                real result.  This name must be in an                  sdzrob24
C                EXTERNAL statement in the calling program.             sdzrob25
C        B     - One end of the interval (B,C).  The value returned for sdzrob26
C                B usually is the better approximation to a zero of F.  sdzrob27
C        C     - The other end of the interval (B,C)                    sdzrob28
C        RE    - Relative error used for RW in the stopping criterion.  sdzrob29
C                If the requested RE is less than machine precision,    sdzrob30
C                then RW is set to approximately machine precision.     sdzrob31
C        AE    - Absolute error used in the stopping criterion.  If the sdzrob32
C                given interval (B,C) contains the origin, then a       sdzrob33
C                nonzero value should be chosen for AE.                 sdzrob34
C                                                                       sdzrob35
C     REFERENCES                                                        sdzrob36
C       1.  L F Shampine and H A Watts, ZEROIN, A Root-Solving Routine, sdzrob37
C           SC-TM-70-631, Sept 1970.                                    sdzrob38
C       2.  T J Dekker, Finding a Zero by Means of Successive Linear    sdzrob39
C           Interpolation, "Constructive Aspects of the Fundamental     sdzrob40
C           Theorem of Algebra", edited by B Dejon and P Henrici, 1969. sdzrob41
C***ROUTINES CALLED  SDNTPB                                             sdzrob42
C***DATE WRITTEN   790601   (YYMMDD)                                    sdzrob43
C***REVISION DATE  830920   (YYMMDD)                                    sdzrob44
C***CATEGORY NO.  I1A2,I1A1B                                            sdzrob45
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                sdzrob46
C             INITIAL VALUE PROBLEMS,GEAR'S OD                      sdzrob47
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdzrob48
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdzrob49
C***END PROLOGUE  SDZROB                                                sdzrob50
!      REAL A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,             sdzrob51
      REAL :: A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,             sdzrob51
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)               sdzrob52
C***FIRST EXECUTABLE STATEMENT  SDZROB                                  sdzrob53
      ER = 4.E0*UROUND                                                  sdzrob54
      RW=MAX(RE,ER)                                                     sdzrob55
      IC=0                                                              sdzrob56
      ACBS=ABS(B-C)                                                     sdzrob57
      A=C                                                               sdzrob58
      FA = FC                                                           sdzrob59
      KOUNT = 0                                                         sdzrob60
    1 IF (ABS(FC) .LT. ABS(FB)) THEN                                    sdzrob61
C                                                    Perform interchangesdzrob62
        A=B                                                             sdzrob63
        FA=FB                                                           sdzrob64
        B=C                                                             sdzrob65
        FB=FC                                                           sdzrob66
        C=A                                                             sdzrob67
        FC=FA                                                           sdzrob68
      END IF                                                            sdzrob69
      CMB=0.5E0*(C-B)                                                   sdzrob70
      ACMB=ABS(CMB)                                                     sdzrob71
      TOL=RW*ABS(B)+AE                                                  sdzrob72
C                                                Test stopping criterionsdzrob73
      IF (ACMB .LE. TOL) RETURN                                         sdzrob74
      IF (KOUNT .GT. 50) RETURN                                         sdzrob75
C                               Calculate new iterate implicitly as     sdzrob76
C                               B+P/Q, where we arrange P .GE. 0.       sdzrob77
C                         The implicit form is used to prevent overflow.sdzrob78
      P=(B-A)*FB                                                        sdzrob79
      Q=FA-FB                                                           sdzrob80
      IF (P .LT. 0.E0) THEN                                             sdzrob81
        P=-P                                                            sdzrob82
        Q=-Q                                                            sdzrob83
      END IF                                                            sdzrob84
C                         Update A and check for satisfactory reduction sdzrob85
C                         in the size of our bounding interval.         sdzrob86
      A=B                                                               sdzrob87
      FA=FB                                                             sdzrob88
      IC=IC+1                                                           sdzrob89
      IF (IC .GE. 4) THEN                                               sdzrob90
        IF (8.E0*ACMB .GE. ACBS) GO TO 6                                sdzrob91
        IC=0                                                            sdzrob92
      END IF                                                            sdzrob93
      ACBS=ACMB                                                         sdzrob94
C                                            Test for too small a changesdzrob95
      IF (P .LE. ABS(Q)*TOL) THEN                                       sdzrob96
C                                                 Increment by tolerancesdzrob97
        B=B+SIGN(TOL,CMB)                                               sdzrob98
        GO TO 7                                                         sdzrob99
      END IF                                                            sdzro100
C                                Root ought to be between B and (C+B)/2.sdzro101
      IF (P .LT. CMB*Q) THEN                                            sdzro102
C                                                            Interpolatesdzro103
        B=B+P/Q                                                         sdzro104
        GO TO 7                                                         sdzro105
      END IF                                                            sdzro106
    6 B=0.5E0*(C+B)                                                     sdzro107
C                                                                 Bisectsdzro108
C                                                                       sdzro109
C                           Have completed computation for new iterate Bsdzro110
    7 CALL SDNTPB (H, 0, N, NQ, T, B, YH,  Y)                           sdzro111
! *****
!      write(= *, fmt = *) " SDNTPB 1"
!      pause
! *****
      FB = F(N, B, Y, IROOT)                                            sdzro112
      IF (FB .EQ. 0.E0) RETURN                                          sdzro113
      KOUNT = KOUNT + 1                                                 sdzro114
C                                                                       sdzro115
C             Decide whether next step is interpolation or extrapolationsdzro116
C                                                                       sdzro117
      IF (SIGN(1.0E0,FB) .EQ. SIGN(1.0E0,FC)) THEN                      sdzro118
        C=A                                                             sdzro119
        FC=FA                                                           sdzro120
      END IF                                                            sdzro121
      GO TO 1                                                           sdzro122
      END                                                               sdzro123
!
      SUBROUTINE SDSTPB (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,  sdstpb 3
     8   MITER,ML,MU,N,NDE,YWT,UROUND,AVGH,AVGORD,H,HUSED,JTASK,MNTOLD, sdstpb 4
     8   MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,HOLD,IPVT, sdstpb 5
     8   JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,MTRSV,     sdstpb 6
     8   MXRDSV)                                                        sdstpb 7
C***BEGIN PROLOGUE  SDSTPB                                              sdstpb 8
C***REFER TO  B3                                                        sdstpb 9
C  SDSTPB performs one step of the integration of an initial value      sdstpb10
C  problem for a system of ordinary differential equations.             sdstpb11
C  Communication with SDSTPB is done with the following variables:      sdstpb12
C                                                                       sdstpb13
C    YH      An N by MAXORD+1 array containing the dependent variables  sdstpb14
C              and their scaled derivatives.  MAXORD, the maximum order sdstpb15
C              used, is currently 12 for the Adams methods and 5 for thesdstpb16
C              Gear methods.  YH(I,J+1) contains the J-th derivative of sdstpb17
C              Y(I), scaled by H**J/factorial(J).  Only Y(I),           sdstpb18
C              1 .LE. I .LE. N, need be set by the calling program on   sdstpb19
C              the first entry.  The YH array should not be altered by  sdstpb20
C              the calling program.  When referencing YH as a           sdstpb21
C              2-dimensional array, use a column length of N, as this issdstpb22
C              the value used in SDSTPB.                                sdstpb23
C    DFDY    A block of locations used for partial derivatives if MITER sdstpb24
C              is not 0.  If MITER is 1 or 2 its length must be at leastsdstpb25
C              N*N.  If MITER is 4 or 5 its length must be at least     sdstpb26
C              (2*ML+MU+1)*N.                                           sdstpb27
C    YWT     An array of N locations used in convergence and error testssdstpb28
C    SAVE1                                                              sdstpb29
C    SAVE2   Arrays of length N used for temporary storage.             sdstpb30
C    IPVT    An integer array of length N used by the linear system     sdstpb31
C              solvers for the storage of row interchange information.  sdstpb32
C    A       A block of locations used to store the matrix A, when usingsdstpb33
C              the implicit method.  If IMPL is 1, A is a MATDIM by N   sdstpb34
C              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4sdstpb35
C              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N. sdstpb36
C    JTASK   An integer used on input.                                  sdstpb37
C              It has the following values and meanings:                sdstpb38
C                .EQ.0  Perform the first step.  This value enables     sdstpb39
C                         the subroutine to initialize itself.          sdstpb40
C                .GT.0  Take a new step continuing from the last.       sdstpb41
C                         Assumes the last step was successful and      sdstpb42
C                         user has not changed any parameters.          sdstpb43
C                .LT.0  Take a new step with a new value of H and/or    sdstpb44
C                         MINT and/or MITER.                            sdstpb45
C    JSTATE  A completion code with the following meanings:             sdstpb46
C                1  The step was successful.                            sdstpb47
C                2  A solution could not be obtained with H.NE.0.       sdstpb48
C                3  A solution was not obtained in MXTRY attempts.      sdstpb49
C                4  For IMPL.NE.0, the matrix A is singular.            sdstpb50
C              On a return with JSTATE.GT.1, the values of T and        sdstpb51
C              the YH array are as of the beginning of the last         sdstpb52
C              step, and H is the last step size attempted.             sdstpb53
C***ROUTINES CALLED  SDNTLB,SDPSTB,SDCORB,SDPSCB,SDSCLB,SNRM2           sdstpb54
C***DATE WRITTEN   790601   (YYMMDD)                                    sdstpb55
C***REVISION DATE  830920   (YYMMDD)                                    sdstpb56
C***CATEGORY NO.  I1A2,I1A1B                                            sdstpb57
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                sdstpb58
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD                      sdstpb59
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdstpb60
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdstpb61
C***END PROLOGUE  SDSTPB                                                sdstpb62
      EXTERNAL F, JACOBN, FA                                            sdstpb63
!      REAL A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,              sdstpb64
      REAL :: A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,           sdstpb64
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,    sdstpb65
     8     ERDN, ERUP, ETEST, H, HMAX, HN, HOLD, HS, HUSED, NUMER, RC,  sdstpb66
     8     RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM, SAVE1(*),   sdstpb67
     8     SAVE2(*), SNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, UROUND,   sdstpb68
     8     Y(*), YH(N,*), YWT(*), Y0NRM                                 sdstpb69
!      INTEGER IPVT(*)                                                   sdstpb70
      INTEGER :: IPVT(*)                                                sdstpb70
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH                       sdstpb71
      DATA BIAS1 /1.3E0/,  BIAS2 /1.2E0/, BIAS3 /1.4E0/, MXFAIL /3/,    sdstpb72
     8     MXITER /3/,     MXTRY /50/,    RCTEST /.3E0/, RMFAIL /2.E0/, sdstpb73
     8     RMNORM /10.E0/, TRSHLD /1.E0/                                sdstpb74
C***FIRST EXECUTABLE STATEMENT  SDSTPB                                  sdstpb75
      BND = 0.E0                                                        sdstpb76
      SWITCH = .FALSE.                                                  sdstpb77
      NTRY = 0                                                          sdstpb78
      TOLD = T                                                          sdstpb79
      NFAIL = 0                                                         sdstpb80
! *****
   !   write(unit = 8, fmt = *) " JTASK =", JTASK, "H =", H, 
   !  1  "IER =", IER, "NTRY =", NTRY, "MXTRY =", MXTRY
! *****
      IF (JTASK.LE.0) THEN                                              sdstpb81
        CALL SDNTLB (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,       sdstpb82
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,      sdstpb83
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,            sdstpb84
     8              A, Convrg, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,      sdstpb85
     8              SAVE2, TQ, TREND, ISWFLG)                           sdstpb86
        IF (H.EQ.0.E0) GO TO 400                                        sdstpb87
        IF (IER) GO TO 420                                              sdstpb88
      END IF                                                            sdstpb89
 100  NTRY = NTRY + 1                                                   sdstpb90
      IF (NTRY.GT.MXTRY) GO TO 410                                      sdstpb91
      T = T + H                                                         sdstpb92
      CALL SDPSCB (1, N, NQ,  YH)                                       sdstpb93
      EVALJC = ((ABS(RC - 1.E0).GT.RCTEST).AND.(MITER.NE.0))            sdstpb94
      EVALFA = .NOT.EVALJC                                              sdstpb95
C                                                                       sdstpb96
 110  ITER = 0                                                          sdstpb97
      DO 115 I = 1,N                                                    sdstpb98
 115    Y(I) = YH(I,1)                                                  sdstpb99
      CALL F (N, T, Y, SAVE2)                                           sdstp100
      NFE = NFE + 1                                                     sdstp101
      IF (EVALJC) THEN                                                  sdstp102
        CALL SDPSTB (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,     sdstp103
     8              MU, N, NDE, NQ, SAVE2, T, Y, YH, YWT, UROUND,       sdstp104
     8              NFE, NJE,  A, DFDY, IER, IPVT, SAVE1, ISWFLG, BND)  sdstp105
        IF (IER) GO TO 160                                              sdstp106
        CONVRG = .FALSE.                                                sdstp107
        RC = 1.E0                                                       sdstp108
      END IF                                                            sdstp109
      DO 125 I = 1,N                                                    sdstp110
 125    SAVE1(I) = 0.E0                                                 sdstp111
C                      Up to MXITER corrector iterations are taken.     sdstp112
C                      Convergence is tested by requiring the r.m.s.    sdstp113
C                      norm of changes to be less than EPS.  The sum of sdstp114
C                      the corrections is accumulated in the vector     sdstp115
C                      SAVE1(I).  It is approximately equal to the L-th sdstp116
C                      derivative of Y multiplied by                    sdstp117
C                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus      sdstp118
C                      proportional to the actual errors to the lowest  sdstp119
C                      power of H present (H**L).  The YH array is not  sdstp120
C                      altered in the correction loop.  The norm of the sdstp121
C                      iterate difference is stored in D.  If           sdstp122
C                      ITER.GT.0, an estimate of the convergence rate   sdstp123
C                      constant is stored in TREND, and this is used in sdstp124
C                      the convergence test.                            sdstp125
C                                                                       sdstp126
 130  CALL SDCORB (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,      sdstp127
     8            MU, N, NDE, NQ, T, Y, YH, YWT,  EVALFA, SAVE1,        sdstp128
     8            SAVE2,  A, D)                                         sdstp129
      IF (ISWFLG.EQ.3 .AND. MINT.EQ.1) THEN                             sdstp130
        IF (ITER.EQ.0) THEN                                             sdstp131
          NUMER = SNRM2(N, SAVE1, 1)                                    sdstp132
          DO 132 I = 1,N                                                sdstp133
 132        DFDY(1,I) = SAVE1(I)                                        sdstp134
          Y0NRM = SNRM2(N, YH, 1)                                       sdstp135
        ELSE                                                            sdstp136
          DENOM = NUMER                                                 sdstp137
          DO 134 I = 1,N                                                sdstp138
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)                            sdstp139
          NUMER = SNRM2(N, DFDY, MATDIM)                                sdstp140
          IF (EL(1,NQ)*NUMER.LE.100.E0*UROUND*Y0NRM) THEN               sdstp141
            IF (RMAX.EQ.RMFAIL) THEN                                    sdstp142
              SWITCH = .TRUE.                                           sdstp143
              GO TO 170                                                 sdstp144
            END IF                                                      sdstp145
          END IF                                                        sdstp146
          DO 136 I = 1,N                                                sdstp147
 136        DFDY(1,I) = SAVE1(I)                                        sdstp148
          IF (DENOM.NE.0.E0)                                            sdstp149
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))                 sdstp150
        END IF                                                          sdstp151
      END IF                                                            sdstp152
      IF (ITER.GT.0) TREND = MAX(.9E0*TREND, D/D1)                      sdstp153
      D1 = D                                                            sdstp154
      CTEST = MIN(2.E0*TREND, 1.E0)*D                                   sdstp155
      IF (CTEST.LE.EPS) GO TO 170                                       sdstp156
      ITER = ITER + 1                                                   sdstp157
      IF (ITER.LT.MXITER) THEN                                          sdstp158
        DO 140 I = 1,N                                                  sdstp159
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)                            sdstp160
        CALL F (N, T, Y, SAVE2)                                         sdstp161
        NFE = NFE + 1                                                   sdstp162
        GO TO 130                                                       sdstp163
      END IF                                                            sdstp164
C                     The corrector iteration failed to converge in     sdstp165
C                     MXITER tries.  If partials are involved but are   sdstp166
C                     not up to date, they are reevaluated for the next sdstp167
C                     try.  Otherwise the YH array is retracted to its  sdstp168
C                     values before prediction, and H is reduced, if    sdstp169
C                     possible.  If not, a no-convergence exit is taken.sdstp170
      IF (CONVRG) THEN                                                  sdstp171
        EVALJC = .TRUE.                                                 sdstp172
        EVALFA = .FALSE.                                                sdstp173
        GO TO 110                                                       sdstp174
      END IF                                                            sdstp175
 160  T = TOLD                                                          sdstp176
      CALL SDPSCB (-1, N, NQ,  YH)                                      sdstp177
      NWAIT = NQ + 2                                                    sdstp178
      IF (JTASK.NE.0 .AND. JTASK.NE.2) RMAX = RMFAIL                    sdstp179
      IF (ITER.EQ.0) THEN                                               sdstp180
        RH = .3E0                                                       sdstp181
      ELSE                                                              sdstp182
        RH = .9E0*(EPS/CTEST)**(.2E0)                                   sdstp183
      END IF                                                            sdstp184
      IF (RH*H.EQ.0.E0) GO TO 400                                       sdstp185
      CALL SDSCLB (HMAX, N, NQ, RMAX,  H, RC, RH, YH)                   sdstp186
      GO TO 100                                                         sdstp187
C                          The corrector has converged.  CONVRG is set  sdstp188
C                          to .TRUE. if partial derivatives were used,  sdstp189
C                          to indicate that they may need updating on   sdstp190
C                          subsequent steps.  The error test is made.   sdstp191
 170  CONVRG = (MITER.NE.0)                                             sdstp192
      DO 180 I = 1,NDE                                                  sdstp193
!      write(unit = 8, fmt = *) "SSSSSSSSSSSSSSSSSSSSSSS1"
 180    SAVE2(I) = SAVE1(I)/YWT(I)                                      sdstp194
      ETEST = SNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(REAL(NDE)))           sdstp195
C                                                                       sdstp196
C                           The error test failed.  NFAIL keeps track ofsdstp197
C                           multiple failures.  Restore T and the YH    sdstp198
C                           array to their previous values, and prepare sdstp199
C                           to try the step again.  Compute the optimum sdstp200
C                           step size for this or one lower order.      sdstp201
      IF (ETEST.GT.EPS) THEN                                            sdstp202
        T = TOLD                                                        sdstp203
        CALL SDPSCB (-1, N, NQ,  YH)                                    sdstp204
        NFAIL = NFAIL + 1                                               sdstp205
        IF (NFAIL.LT.MXFAIL) THEN                                       sdstp206
          IF (JTASK.NE.0 .AND. JTASK.NE.2) RMAX = RMFAIL                sdstp207
          RH2 = 1.E0/(BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))             sdstp208
          IF (NQ.GT.1) THEN                                             sdstp209
            DO 190 I = 1,NDE                                            sdstp210
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)                              sdstp211
            ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))      sdstp212
            RH1 = 1.E0/MAX(1.E0, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))     sdstp213
            IF (RH2.LT.RH1) THEN                                        sdstp214
              NQ = NQ - 1                                               sdstp215
              NWAIT = NQ + 2                                            sdstp216
              RC = RC*EL(1,NQ)/EL(1,NQ+1)                               sdstp217
              RH = RH1                                                  sdstp218
              IF (RH*H.EQ.0.E0) GO TO 400                               sdstp219
              CALL SDSCLB (HMAX, N, NQ, RMAX,  H, RC, RH, YH)           sdstp220
              GO TO 100                                                 sdstp221
            END IF                                                      sdstp222
          END IF                                                        sdstp223
          NWAIT = NQ + 2                                                sdstp224
          RH = RH2                                                      sdstp225
! *****
   !   write(unit = *, fmt = *) " At sdstp226:  RH*H = ", RH*H
!      pause
! *****          
          IF (RH*H .EQ. 0.E0) GO TO 400                                 sdstp226
          CALL SDSCLB (HMAX, N, NQ, RMAX,  H, RC, RH, YH)               sdstp227
          GO TO 100                                                     sdstp228
        END IF                                                          sdstp229
C                Control reaches this section if the error test has     sdstp230
C                failed MXFAIL or more times.  It is assumed that the   sdstp231
C                derivatives that have accumulated in the YH array have sdstp232
C                errors of the wrong order.  Hence the first derivative sdstp233
C                is recomputed, the order is set to 1, and the step is  sdstp234
C                retried.                                               sdstp235
        NFAIL = 0                                                       sdstp236
        JTASK = 2                                                       sdstp237
        DO 215 I = 1, N                                                 sdstp238
 215      Y(I) = YH(I, 1)
        CALL SDNTLB (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,       sdstp240
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,      sdstp241
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,            sdstp242
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,      sdstp243
     8              SAVE2, TQ, TREND, ISWFLG)                           sdstp244
        IF (H.EQ.0.E0) GO TO 400                                        sdstp245
        IF (IER) GO TO 420                                              sdstp246
        GO TO 100                                                       sdstp247
      END IF                                                            sdstp248
C                          After a successful step, update the YH array.sdstp249
      NStep = NStep + 1                                                 sdstp250
      HUSED = H                                                         sdstp251
      NQUSED = NQ                                                       sdstp252
      AVGH = (REAL(NStep - 1)*AVGH + H)/REAL(NSTEP)                       sdstp253
      AVGORD = (REAL(NStep - 1)*AVGORD + REAL(NQ))/REAL(NStep)            sdstp254
      DO 230 J = 1,NQ+1                                                 sdstp255
        DO 230 I = 1,N                                                  sdstp256
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)                         sdstp257
      DO 235 I = 1,N                                                    sdstp258
 235    Y(I) = YH(I,1)                                                  sdstp259
C                                          If ISWFLG is 3, consider     sdstp260
C                                          changing integration methods.sdstp261
C                                                                       sdstp262
      IF (ISWFLG.EQ.3) THEN                                             sdstp263
        IF (BND.NE.0.E0) THEN                                           sdstp264
          IF (MINT.EQ.1 .AND. NQ.LE.5) THEN                             sdstp265
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))     sdstp266
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))                      sdstp267
            HS = ABS(H)/MAX(UROUND,                                     sdstp268
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.E0/REAL(NQ+1)))                sdstp269
            IF (HS.GT.1.2E0*HN) THEN                                    sdstp270
              MINT = 2                                                  sdstp271
              MNTOLD = MINT                                             sdstp272
              MITER = MTRSV                                             sdstp273
              MTROLD = MITER                                            sdstp274
              MAXORD = MIN(MXRDSV, 5)                                   sdstp275
              RC = 0.E0                                                 sdstp276
              RMAX = RMNORM                                             sdstp277
              TREND = 1.E0                                              sdstp278
              CALL SDCSTB (MAXORD, MINT, ISWFLG, EL, TQ)                sdstp279
              NWAIT = NQ + 2                                            sdstp280
            END IF                                                      sdstp281
          ELSE IF (MINT.EQ.2) THEN                                      sdstp282
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))     sdstp283
            HN = ABS(H)/MAX(UROUND,                                     sdstp284
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.E0/REAL(NQ+1)))                  sdstp285
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))                      sdstp286
            IF (HN.GE.HS) THEN                                          sdstp287
              MINT = 1                                                  sdstp288
              MNTOLD = MINT                                             sdstp289
              MITER = 0                                                 sdstp290
              MTROLD = MITER                                            sdstp291
              MAXORD = MIN(MXRDSV, 12)                                  sdstp292
              RMAX = RMNORM                                             sdstp293
              TREND = 1.E0                                              sdstp294
              CONVRG = .FALSE.                                          sdstp295
              CALL SDCSTB (MAXORD, MINT, ISWFLG, EL, TQ)                sdstp296
              NWAIT = NQ + 2                                            sdstp297
            END IF                                                      sdstp298
          END IF                                                        sdstp299
        END IF                                                          sdstp300
      END IF                                                            sdstp301
      IF (SWITCH) THEN                                                  sdstp302
        MINT = 2                                                        sdstp303
        MNTOLD = MINT                                                   sdstp304
        MITER = MTRSV                                                   sdstp305
        MTROLD = MITER                                                  sdstp306
        MAXORD = MIN(MXRDSV, 5)                                         sdstp307
        NQ = MIN(NQ, MAXORD)                                            sdstp308
        RC = 0.E0                                                       sdstp309
        RMAX = RMNORM                                                   sdstp310
        TREND = 1.E0                                                    sdstp311
        CALL SDCSTB (MAXORD, MINT, ISWFLG, EL, TQ)                      sdstp312
        NWAIT = NQ + 2                                                  sdstp313
      END IF                                                            sdstp314
C                           Consider changing H if NWAIT = 1.  Otherwisesdstp315
C                           decrease NWAIT by 1.  If NWAIT is then 1 andsdstp316
C                           NQ.LT.MAXORD, then SAVE1 is saved for use insdstp317
C                           a possible order increase on the next step. sdstp318
C                                                                       sdstp319
      IF (JTASK.EQ.0 .OR. JTASK.EQ.2) THEN                              sdstp320
        RH = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))     sdstp321
        IF (RH.GT.TRSHLD) CALL SDSCLB (HMAX, N, NQ, RMAX, H, RC, RH, YH)sdstp322
      ELSE IF (NWAIT.GT.1) THEN                                         sdstp323
        NWAIT = NWAIT - 1                                               sdstp324
        IF (NWAIT.EQ.1 .AND. NQ.LT.MAXORD) THEN                         sdstp325
          DO 250 I = 1,NDE                                              sdstp326
 250        YH(I,MAXORD+1) = SAVE1(I)                                   sdstp327
        END IF                                                          sdstp328
C             If a change in H is considered, an increase or decrease insdstp329
C             order by one is considered also.  A change in H is made   sdstp330
C             only if it is by a factor of at least TRSHLD.  Factors    sdstp331
C             RH1, RH2, and RH3 are computed, by which H could be       sdstp332
C             multiplied at order NQ - 1, order NQ, or order NQ + 1,    sdstp333
C             respectively.  The largest of these is determined and the sdstp334
C             new order chosen accordingly.  If the order is to be      sdstp335
C             increased, we compute one additional scaled derivative.   sdstp336
C             If there is a change of order, reset NQ and the           sdstp337
C             coefficients.  In any case H is reset according to RH and sdstp338
C             the YH array is rescaled.                                 sdstp339
      ELSE                                                              sdstp340
        IF (NQ.EQ.1) THEN                                               sdstp341
          RH1 = 0.E0                                                    sdstp342
        ELSE                                                            sdstp343
          DO 270 I = 1,NDE                                              sdstp344
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)                                sdstp345
          ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))        sdstp346
          RH1 = 1.E0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))     sdstp347
        END IF                                                          sdstp348
        RH2 = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))    sdstp349
        IF (NQ.EQ.MAXORD) THEN                                          sdstp350
          RH3 = 0.E0                                                    sdstp351
        ELSE                                                            sdstp352
          DO 290 I = 1,NDE                                              sdstp353
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)               sdstp354
          ERUP = SNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(REAL(NDE)))        sdstp355
          RH3 = 1.E0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.E0/REAL(NQ+2)))   sdstp356
        END IF                                                          sdstp357
        IF (RH1.GT.RH2 .AND. RH1.GE.RH3) THEN                           sdstp358
          RH = RH1                                                      sdstp359
          IF (RH.LE.TRSHLD) GO TO 380                                   sdstp360
          NQ = NQ - 1                                                   sdstp361
          RC = RC*EL(1,NQ)/EL(1,NQ+1)                                   sdstp362
        ELSE IF (RH2.GE.RH1 .AND. RH2.GE.RH3) THEN                      sdstp363
          RH = RH2                                                      sdstp364
          IF (RH.LE.TRSHLD) GO TO 380                                   sdstp365
        ELSE                                                            sdstp366
          RH = RH3                                                      sdstp367
          IF (RH.LE.TRSHLD) GO TO 380                                   sdstp368
          DO 360 I = 1,N                                                sdstp369
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/REAL(NQ+1)                sdstp370
          NQ = NQ + 1                                                   sdstp371
          RC = RC*EL(1,NQ)/EL(1,NQ-1)                                   sdstp372
        END IF                                                          sdstp373
        IF (ISWFLG.EQ.3 .AND. MINT.EQ.1) THEN                           sdstp374
          IF (BND.NE.0.E0) RH = MIN(RH, 1.E0/(2.E0*EL(1,NQ)*BND*ABS(H)))sdstp375
        END IF                                                          sdstp376
        CALL SDSCLB (HMAX, N, NQ, RMAX,  H, RC, RH, YH)                 sdstp377
        RMAX = RMNORM                                                   sdstp378
 380    NWAIT = NQ + 2                                                  sdstp379
      END IF                                                            sdstp380
C               All returns are made through this section.  H is saved  sdstp381
C               in HOLD to allow the caller to change H on the next stepsdstp382
      JSTATE = 1                                                        sdstp383
      HOLD = H                                                          sdstp384
      RETURN                                                            sdstp385
C                                                                       sdstp386
 400  JSTATE = 2                                                        sdstp387
      HOLD = H                                                          sdstp388
      DO 405 I = 1,N                                                    sdstp389
 405    Y(I) = YH(I,1)                                                  sdstp390
      RETURN                                                            sdstp391
C                                                                       sdstp392
 410  JSTATE = 3                                                        sdstp393
      HOLD = H                                                          sdstp394
      RETURN                                                            sdstp395
C                                                                       sdstp396
 420  JSTATE = 4                                                        sdstp397
      HOLD = H                                                          sdstp398
      RETURN                                                            sdstp399
      END                                                               sdstp400
!
      REAL FUNCTION SNRM2(N,SX,INCX)                                    SNRM2  4
C***BEGIN PROLOGUE  SNRM2                                               SNRM2  5
C***DATE WRITTEN   791001   (YYMMDD)                                    SNRM2  6
C***REVISION DATE  820801   (YYMMDD)                                    SNRM2  7
C***CATEGORY NO.  D1A3B                                                 SNRM2  8
C***KEYWORDS  BLAS,EUCLIDEAN,L2,LENGTH,LINEAR ALGEBRA,NORM,VECTOR       SNRM2  9
C***AUTHOR  LAWSON, C. L., (JPL)                                        SNRM2 10
C           HANSON, R. J., (SNLA)                                       SNRM2 11
C           KINCAID, D. R., (U. OF TEXAS)                               SNRM2 12
C           KROGH, F. T., (JPL)                                         SNRM2 13
C***PURPOSE  EUCLIDEAN LENGTH (L2 NORM) OF S.P. VECTOR                  SNRM2 14
C***DESCRIPTION                                                         SNRM2 15
C                                                                       SNRM2 16
C                B L A S  SUBPROGRAM                                    SNRM2 17
C    DESCRIPTION OF PARAMETERS                                          SNRM2 18
C                                                                       SNRM2 19
C     --INPUT--                                                         SNRM2 20
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       SNRM2 21
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     SNRM2 22
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      SNRM2 23
C                                                                       SNRM2 24
C     --OUTPUT--                                                        SNRM2 25
C    SNRM2  SINGLE PRECISION RESULT (ZERO IF N .LE. 0)                  SNRM2 26
C                                                                       SNRM2 27
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE        SNRM2 28
C     INCREMENT INCX .                                                  SNRM2 29
C     IF N .LE. 0, RETURN WITH RESULT = 0.                              SNRM2 30
C     IF N .GE. 1, THEN INCX MUST BE .GE. 1                             SNRM2 31
C                                                                       SNRM2 32
C           C. L. LAWSON, 1978 JAN 08                                   SNRM2 33
C                                                                       SNRM2 34
C     FOUR PHASE OD     USING TWO BUILT-IN CONSTANTS THAT ARE       SNRM2 35
C     HOPEFULLY APPLICABLE TO ALL MACHINES.                             SNRM2 36
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.     SNRM2 37
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.     SNRM2 38
C     WHERE                                                             SNRM2 39
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.                 SNRM2 40
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)               SNRM2 41
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)               SNRM2 42
C                                                                       SNRM2 43
C     BRIEF OUTLINE OF ALGORITHM..                                      SNRM2 44
C                                                                       SNRM2 45
C     PHASE 1 SCANS ZERO COMPONENTS.                                    SNRM2 46
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO        SNRM2 47
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO                    SNRM2 48
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M                  SNRM2 49
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.                 SNRM2 50
C                                                                       SNRM2 51
C     VALUES FOR CUTLO AND CUTHI..                                      SNRM2 52
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER    SNRM2 53
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..                     SNRM2 54
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARESNRM2 55
C                   UNIVAC AND DEC AT 2**(-103)                         SNRM2 56
C                   THUS CUTLO = 2**(-51) = 4.44089E-16                 SNRM2 57
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.          SNRM2 58
C                   THUS CUTHI = 2**(63.5) = 1.30438E19                 SNRM2 59
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.             SNRM2 60
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11               SNRM2 61
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19                    SNRM2 62
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /                        SNRM2 63
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /                        SNRM2 64
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     SNRM2 65
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, SNRM2 66
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       SNRM2 67
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 SNRM2 68
C***ROUTINES CALLED  (NONE)                                             SNRM2 69
C***END PROLOGUE  SNRM2                                                 SNRM2 70
      INTEGER          NEXT                                             SNRM2 71
!      REAL   SX(*),  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE         SNRM2 72
      REAL :: SX(*),  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE         SNRM2 72
      DATA   ZERO, ONE /0.0E0, 1.0E0/                                   SNRM2 73
C                                                                       SNRM2 74
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /                        SNRM2 75
C***FIRST EXECUTABLE STATEMENT  SNRM2                                   SNRM2 76
      IF(N .GT. 0) GO TO 10                                             SNRM2 77
         SNRM2  = ZERO                                                  SNRM2 78
         GO TO 300                                                      SNRM2 79
C                                                                       SNRM2 80
   10 ASSIGN 30 TO NEXT                                                 SNRM2 81
      SUM = ZERO                                                        SNRM2 82
      NN = N * INCX                                                     SNRM2 83
C                                                 BEGIN MAIN LOOP       SNRM2 84
      I = 1                                                             SNRM2 85
   20    GO TO NEXT,(30, 50, 70, 110)                                   SNRM2 86
   30 IF( ABS(SX(I)) .GT. CUTLO) GO TO 85                               SNRM2 87
      ASSIGN 50 TO NEXT                                                 SNRM2 88
      XMAX = ZERO                                                       SNRM2 89
C                                                                       SNRM2 90
C                        PHASE 1.  SUM IS ZERO                          SNRM2 91
C                                                                       SNRM2 92
   50 IF( SX(I) .EQ. ZERO) GO TO 200                                    SNRM2 93
      IF( ABS(SX(I)) .GT. CUTLO) GO TO 85                               SNRM2 94
C                                                                       SNRM2 95
C                                PREPARE FOR PHASE 2.                   SNRM2 96
      ASSIGN 70 TO NEXT                                                 SNRM2 97
      GO TO 105                                                         SNRM2 98
C                                                                       SNRM2 99
C                                PREPARE FOR PHASE 4.                   SNRM2100
C                                                                       SNRM2101
  100 I = J                                                             SNRM2102
      ASSIGN 110 TO NEXT                                                SNRM2103
      SUM = (SUM / SX(I)) / SX(I)                                       SNRM2104
  105 XMAX = ABS(SX(I))                                                 SNRM2105
      GO TO 115                                                         SNRM2106
C                                                                       SNRM2107
C                   PHASE 2.  SUM IS SMALL.                             SNRM2108
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.     SNRM2109
C                                                                       SNRM2110
   70 IF( ABS(SX(I)) .GT. CUTLO ) GO TO 75                              SNRM2111
C                                                                       SNRM2112
C                     COMMON CODE FOR PHASES 2 AND 4.                   SNRM2113
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.SNRM2114
C                                                                       SNRM2115
  110 IF( ABS(SX(I)) .LE. XMAX ) GO TO 115                              SNRM2116
         SUM = ONE + SUM * (XMAX / SX(I))**2                            SNRM2117
         XMAX = ABS(SX(I))                                              SNRM2118
         GO TO 200                                                      SNRM2119
C                                                                       SNRM2120
  115 SUM = SUM + (SX(I)/XMAX)**2                                       SNRM2121
      GO TO 200                                                         SNRM2122
C                                                                       SNRM2123
C                                                                       SNRM2124
C                  PREPARE FOR PHASE 3.                                 SNRM2125
C                                                                       SNRM2126
   75 SUM = (SUM * XMAX) * XMAX                                         SNRM2127
C                                                                       SNRM2128
C                                                                       SNRM2129
C     FOR REAL OR D.P. SET HITEST = CUTHI/N                             SNRM2130
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)                         SNRM2131
C                                                                       SNRM2132
   85 HITEST = CUTHI/FLOAT( N )                                         SNRM2133
C                                                                       SNRM2134
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.            SNRM2135
C                                                                       SNRM2136
      DO 95 J =I,NN,INCX                                                SNRM2137
      IF(ABS(SX(J)) .GE. HITEST) GO TO 100                              SNRM2138
   95    SUM = SUM + SX(J)**2                                           SNRM2139
      SNRM2 = SQRT( SUM )                                               SNRM2140
      GO TO 300                                                         SNRM2141
C                                                                       SNRM2142
  200 CONTINUE                                                          SNRM2143
      I = I + INCX                                                      SNRM2144
      IF ( I .LE. NN ) GO TO 20                                         SNRM2145
C                                                                       SNRM2146
C              END OF MAIN LOOP.                                        SNRM2147
C                                                                       SNRM2148
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.              SNRM2149
C                                                                       SNRM2150
      SNRM2 = XMAX * SQRT(SUM)                                          SNRM2151
  300 CONTINUE                                                          SNRM2152
      RETURN                                                            SNRM2153
      END                                                               SNRM2154
!
      SUBROUTINE SGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)                      SGBSL  3
C***BEGIN PROLOGUE  SGBSL                                               SGBSL  4
C***DATE WRITTEN   780814   (YYMMDD)                                    SGBSL  5
C***REVISION DATE  820801   (YYMMDD)                                    SGBSL  6
C***CATEGORY NO.  D2A2                                                  SGBSL  7
C***KEYWORDS  BANDED,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE                SGBSL  8
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                            SGBSL  9
C***PURPOSE  SOLVES THE REAL BAND SYSTEM A*X=B OR TRANS(A)*X=B          SGBSL 10
C            USING THE FACTORS COMPUTED BY SGBCO OR SGBFA.              SGBSL 11
C***DESCRIPTION                                                         SGBSL 12
C                                                                       SGBSL 13
C     SGBSL SOLVES THE REAL BAND SYSTEM                                 SGBSL 14
C     A * X = B  OR  TRANS(A) * X = B                                   SGBSL 15
C     USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.                     SGBSL 16
C                                                                       SGBSL 17
C     ON ENTRY                                                          SGBSL 18
C                                                                       SGBSL 19
C        ABD     REAL(LDA, N)                                           SGBSL 20
C                THE OUTPUT FROM SBGCO OR SGBFA.                        SGBSL 21
C                                                                       SGBSL 22
C        LDA     INTEGER                                                SGBSL 23
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              SGBSL 24
C                                                                       SGBSL 25
C        N       INTEGER                                                SGBSL 26
C                THE ORDER OF THE ORIGINAL MATRIX.                      SGBSL 27
C                                                                       SGBSL 28
C        ML      INTEGER                                                SGBSL 29
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           SGBSL 30
C                                                                       SGBSL 31
C        MU      INTEGER                                                SGBSL 32
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           SGBSL 33
C                                                                       SGBSL 34
C        IPVT    INTEGER(N)                                             SGBSL 35
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.                  SGBSL 36
C                                                                       SGBSL 37
C        B       REAL(N)                                                SGBSL 38
C                THE RIGHT HAND SIDE VECTOR.                            SGBSL 39
C                                                                       SGBSL 40
C        JOB     INTEGER                                                SGBSL 41
C                = 0         TO SOLVE  A*X = B ,                        SGBSL 42
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE           SGBSL 43
C                            TRANS(A)  IS THE TRANSPOSE.                SGBSL 44
C                                                                       SGBSL 45
C     ON RETURN                                                         SGBSL 46
C                                                                       SGBSL 47
C        B       THE SOLUTION VECTOR  X .                               SGBSL 48
C                                                                       SGBSL 49
C     ERROR CONDITION                                                   SGBSL 50
C                                                                       SGBSL 51
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   SGBSL 52
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,SGBSL 53
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       SGBSL 54
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     SGBSL 55
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0           SGBSL 56
C        OR SGBFA HAS SET INFO .EQ. 0 .                                 SGBSL 57
C                                                                       SGBSL 58
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 SGBSL 59
C     WITH  P  COLUMNS                                                  SGBSL 60
C           CALL SBGCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    SGBSL 61
C           IF (RCOND IS TOO SMALL) GO TO ...                           SGBSL 62
C           DO 10 J = 1, P                                              SGBSL 63
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)                SGBSL 64
C        10 CONTINUE                                                    SGBSL 65
C                                                                       SGBSL 66
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           SGBSL 67
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      SGBSL 68
C                                                                       SGBSL 69
C     SUBROUTINES AND FUNCTIONS                                         SGBSL 70
C                                                                       SGBSL 71
C     BLAS SAXPY,SDOT                                                   SGBSL 72
C     FORTRAN MIN0                                                      SGBSL 73
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    SGBSL 74
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   SGBSL 75
C***ROUTINES CALLED  SAXPY,SDOT                                         SGBSL 76
C***END PROLOGUE  SGBSL                                                 SGBSL 77
      INTEGER LDA,N,ML,MU,IPVT(*),JOB                                   SGBSL 78
!      REAL ABD(LDA,*),B(*)                                              SGBSL 79
      REAL :: ABD(LDA,*),B(*)                                              SGBSL 79
C                                                                       SGBSL 80
!      REAL SDOT,T                                                       SGBSL 81
      REAL :: SDOT,T                                                       SGBSL 81
!      INTEGER ki,KB,L,LA,LB,LM,M,NM1                                     SGBSL 82
      INTEGER :: ki,KB,L,LA,LB,LM,M,NM1                                     SGBSL 82
C***FIRST EXECUTABLE STATEMENT  SGBSL                                   SGBSL 83
      M = MU + ML + 1                                                   SGBSL 84
      NM1 = N - 1                                                       SGBSL 85
      IF (JOB .NE. 0) GO TO 50                                          SGBSL 86
C                                                                       SGBSL 87
C        JOB = 0 , SOLVE  A * X = B                                     SGBSL 88
C        FIRST SOLVE L*Y = B                                            SGBSL 89
C                                                                       SGBSL 90
         IF (ML .EQ. 0) GO TO 30                                        SGBSL 91
         IF (NM1 .LT. 1) GO TO 30                                       SGBSL 92
            DO 20 ki = 1, NM1                                            SGBSL 93
               LM = MIN0(ML,N-ki)                                        SGBSL 94
               L = IPVT(ki)                                              SGBSL 95
               T = B(L)                                                 SGBSL 96
               IF (L .EQ. ki) GO TO 10                                   SGBSL 97
                  B(L) = B(ki)                                           SGBSL 98
                  B(ki) = T                                              SGBSL 99
   10          CONTINUE                                                 SGBSL100
               CALL SAXPY(LM,T,ABD(M+1,ki),1,B(ki+1),1)                   SGBSL101
   20       CONTINUE                                                    SGBSL102
   30    CONTINUE                                                       SGBSL103
C                                                                       SGBSL104
C        NOW SOLVE  U*X = Y                                             SGBSL105
C                                                                       SGBSL106
         DO 40 KB = 1, N                                                SGBSL107
            ki = N + 1 - KB                                              SGBSL108
            B(ki) = B(ki)/ABD(M,ki)                                        SGBSL109
            LM = MIN0(ki,M) - 1                                          SGBSL110
            LA = M - LM                                                 SGBSL111
            LB = ki - LM                                                 SGBSL112
            T = -B(ki)                                                   SGBSL113
            CALL SAXPY(LM,T,ABD(LA,ki),1,B(LB),1)                        SGBSL114
   40    CONTINUE                                                       SGBSL115
      GO TO 100                                                         SGBSL116
   50 CONTINUE                                                          SGBSL117
C                                                                       SGBSL118
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         SGBSL119
C        FIRST SOLVE  TRANS(U)*Y = B                                    SGBSL120
C                                                                       SGBSL121
         DO 60 ki = 1, N                                                 SGBSL122
            LM = MIN0(ki,M) - 1                                          SGBSL123
            LA = M - LM                                                 SGBSL124
            LB = ki - LM                                                 SGBSL125
            T = SDOT(LM,ABD(LA,ki),1,B(LB),1)                            SGBSL126
            B(ki) = (B(ki) - T)/ABD(M,ki)                                  SGBSL127
   60    CONTINUE                                                       SGBSL128
C                                                                       SGBSL129
C        NOW SOLVE TRANS(L)*X = Y                                       SGBSL130
C                                                                       SGBSL131
         IF (ML .EQ. 0) GO TO 90                                        SGBSL132
         IF (NM1 .LT. 1) GO TO 90                                       SGBSL133
            DO 80 KB = 1, NM1                                           SGBSL134
               ki = N - KB                                               SGBSL135
               LM = MIN0(ML,N-ki)                                        SGBSL136
               B(ki) = B(ki) + SDOT(LM,ABD(M+1,ki),1,B(ki+1),1)             SGBSL137
               L = IPVT(ki)                                              SGBSL138
               IF (L .EQ. ki) GO TO 70                                   SGBSL139
                  T = B(L)                                              SGBSL140
                  B(L) = B(ki)                                           SGBSL141
                  B(ki) = T                                              SGBSL142
   70          CONTINUE                                                 SGBSL143
   80       CONTINUE                                                    SGBSL144
   90    CONTINUE                                                       SGBSL145
  100 CONTINUE                                                          SGBSL146
      RETURN                                                            SGBSL147
      END                                                               SGBSL148
!
      SUBROUTINE SGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                       SGBFA  3
C***BEGIN PROLOGUE  SGBFA                                               SGBFA  4
C***DATE WRITTEN   780814   (YYMMDD)                                    SGBFA  5
C***REVISION DATE  820801   (YYMMDD)                                    SGBFA  6
C***CATEGORY NO.  D2A2                                                  SGBFA  7
C***KEYWORDS  BANDED,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX               SGBFA  8
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                            SGBFA  9
C***PURPOSE  FACTORS A REAL BAND MATRIX BY ELIMINATION.                 SGBFA 10
C***DESCRIPTION                                                         SGBFA 11
C                                                                       SGBFA 12
C     SGBFA FACTORS A REAL BAND MATRIX BY ELIMINATION.                  SGBFA 13
C                                                                       SGBFA 14
C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED            SGBFA 15
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          SGBFA 16
C                                                                       SGBFA 17
C     ON ENTRY                                                          SGBFA 18
C                                                                       SGBFA 19
C        ABD     REAL(LDA, N)                                           SGBFA 20
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      SGBFA 21
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   SGBFA 22
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         SGBFA 23
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       SGBFA 24
C                SEE THE COMMENTS BELOW FOR DETAILS.                    SGBFA 25
C                                                                       SGBFA 26
C        LDA     INTEGER                                                SGBFA 27
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              SGBFA 28
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       SGBFA 29
C                                                                       SGBFA 30
C        N       INTEGER                                                SGBFA 31
C                THE ORDER OF THE ORIGINAL MATRIX.                      SGBFA 32
C                                                                       SGBFA 33
C        ML      INTEGER                                                SGBFA 34
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           SGBFA 35
C                0 .LE. ML .LT. N .                                     SGBFA 36
C                                                                       SGBFA 37
C        MU      INTEGER                                                SGBFA 38
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           SGBFA 39
C                0 .LE. MU .LT. N .                                     SGBFA 40
C                MORE EFFICIENT IF  ML .LE. MU .                        SGBFA 41
C     ON RETURN                                                         SGBFA 42
C                                                                       SGBFA 43
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         SGBFA 44
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          SGBFA 45
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE      SGBFA 46
C                L  IS A PRODUCT OF PERMUTATION AND  LOWER              SGBFA 47
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       SGBFA 48
C                                                                       SGBFA 49
C        IPVT    INTEGER(N)                                             SGBFA 50
C                AN INTEGER VECTOR OF PIVOT INDICES.                    SGBFA 51
C                                                                       SGBFA 52
C        INFO    INTEGER                                                SGBFA 53
C                = 0  NORMAL VALUE.                                     SGBFA 54
C                = ki  IF  U(ki,ki) .EQ. 0.0 .  THIS IS NOT AN ERROR     SGBFA 55
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        SGBFA 56
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF        SGBFA 57
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE      SGBFA 58
C                     INDICATION OF SINGULARITY.                        SGBFA 59
C                                                                       SGBFA 60
C     BAND STORAGE                                                      SGBFA 61
C                                                                       SGBFA 62
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      SGBFA 63
C           WILL SET UP THE INPUT.                                      SGBFA 64
C                                                                       SGBFA 65
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                SGBFA 66
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                SGBFA 67
C                   M = ML + MU + 1                                     SGBFA 68
C                   DO 20 J = 1, N                                      SGBFA 69
C                      I1 = MAX0(1, J-MU)                               SGBFA 70
C                      I2 = MIN0(N, J+ML)                               SGBFA 71
C                      DO 10 I = I1, I2                                 SGBFA 72
C                         ki = I - J + M                                SGBFA 73
C                         ABD(ki,J) = A(I,J)                             SGBFA 74
C                10    CONTINUE                                         SGBFA 75
C                20 CONTINUE                                            SGBFA 76
C                                                                       SGBFA 77
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         SGBFA 78
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      SGBFA 79
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            SGBFA 80
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    SGBFA 81
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            SGBFA 82
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          SGBFA 83
C                                                                       SGBFA 84
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           SGBFA 85
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      SGBFA 86
C                                                                       SGBFA 87
C     SUBROUTINES AND FUNCTIONS                                         SGBFA 88
C                                                                       SGBFA 89
C     BLAS SAXPY,SSCAL,ISAMAX                                           SGBFA 90
C     FORTRAN MAX0,MIN0                                                 SGBFA 91
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    SGBFA 92
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   SGBFA 93
C***ROUTINES CALLED  ISAMAX,SAXPY,SSCAL                                 SGBFA 94
C***END PROLOGUE  SGBFA                                                 SGBFA 95
      INTEGER LDA,N,ML,MU,IPVT(*),INFO                                  SGBFA 96
!      REAL ABD(LDA,*)                                                   SGBFA 97
      REAL :: ABD(LDA,*)                                                   SGBFA 97
C                                                                       SGBFA 98
!      REAL T                                                            SGBFA 99
      REAL :: T                                                            SGBFA 99
!      INTEGER I,ISAMAX,I0,J,JU,JZ,J0,J1,ki,KP1,L,LM,M,MM,NM1             SGBFA100
      INTEGER :: I,ISAMAX,I0,J,JU,JZ,J0,J1,ki,KP1,L,LM,M,MM,NM1             SGBFA100
C                                                                       SGBFA101
C***FIRST EXECUTABLE STATEMENT  SGBFA                                   SGBFA102
      M = ML + MU + 1                                                   SGBFA103
      INFO = 0                                                          SGBFA104
C                                                                       SGBFA105
C     ZERO INITIAL FILL-IN COLUMNS                                      SGBFA106
C                                                                       SGBFA107
      J0 = MU + 2                                                       SGBFA108
      J1 = MIN0(N,M) - 1                                                SGBFA109
      IF (J1 .LT. J0) GO TO 30                                          SGBFA110
      DO 20 JZ = J0, J1                                                 SGBFA111
         I0 = M + 1 - JZ                                                SGBFA112
         DO 10 I = I0, ML                                               SGBFA113
            ABD(I,JZ) = 0.0E0                                           SGBFA114
   10    CONTINUE                                                       SGBFA115
   20 CONTINUE                                                          SGBFA116
   30 CONTINUE                                                          SGBFA117
      JZ = J1                                                           SGBFA118
      JU = 0                                                            SGBFA119
C                                                                       SGBFA120
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        SGBFA121
C                                                                       SGBFA122
      NM1 = N - 1                                                       SGBFA123
      IF (NM1 .LT. 1) GO TO 130                                         SGBFA124
      DO 120 ki = 1, NM1                                                 SGBFA125
         KP1 = ki + 1                                                    SGBFA126
C                                                                       SGBFA127
C        ZERO NEXT FILL-IN COLUMN                                       SGBFA128
C                                                                       SGBFA129
         JZ = JZ + 1                                                    SGBFA130
         IF (JZ .GT. N) GO TO 50                                        SGBFA131
         IF (ML .LT. 1) GO TO 50                                        SGBFA132
            DO 40 I = 1, ML                                             SGBFA133
               ABD(I,JZ) = 0.0E0                                        SGBFA134
   40       CONTINUE                                                    SGBFA135
   50    CONTINUE                                                       SGBFA136
C                                                                       SGBFA137
C        FIND L = PIVOT INDEX                                           SGBFA138
C                                                                       SGBFA139
         LM = MIN0(ML,N-ki)                                              SGBFA140
         L = ISAMAX(LM+1,ABD(M,ki),1) + M - 1                            SGBFA141
         IPVT(ki) = L + ki - M                                            SGBFA142
C                                                                       SGBFA143
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          SGBFA144
C                                                                       SGBFA145
         IF (ABD(L,ki) .EQ. 0.0E0) GO TO 100                             SGBFA146
C                                                                       SGBFA147
C           INTERCHANGE IF NECESSARY                                    SGBFA148
C                                                                       SGBFA149
            IF (L .EQ. M) GO TO 60                                      SGBFA150
               T = ABD(L,ki)                                             SGBFA151
               ABD(L,ki) = ABD(M,ki)                                      SGBFA152
               ABD(M,ki) = T                                             SGBFA153
   60       CONTINUE                                                    SGBFA154
C                                                                       SGBFA155
C           COMPUTE MULTIPLIERS                                         SGBFA156
C                                                                       SGBFA157
            T = -1.0E0/ABD(M,ki)                                         SGBFA158
            CALL SSCAL(LM,T,ABD(M+1,ki),1)                               SGBFA159
C                                                                       SGBFA160
C           ROW ELIMINATION WITH COLUMN INDEXING                        SGBFA161
C                                                                       SGBFA162
            JU = MIN0(MAX0(JU,MU+IPVT(ki)),N)                            SGBFA163
            MM = M                                                      SGBFA164
            IF (JU .LT. KP1) GO TO 90                                   SGBFA165
            DO 80 J = KP1, JU                                           SGBFA166
               L = L - 1                                                SGBFA167
               MM = MM - 1                                              SGBFA168
               T = ABD(L,J)                                             SGBFA169
               IF (L .EQ. MM) GO TO 70                                  SGBFA170
                  ABD(L,J) = ABD(MM,J)                                  SGBFA171
                  ABD(MM,J) = T                                         SGBFA172
   70          CONTINUE                                                 SGBFA173
               CALL SAXPY(LM,T,ABD(M+1,ki),1,ABD(MM+1,J),1)              SGBFA174
   80       CONTINUE                                                    SGBFA175
   90       CONTINUE                                                    SGBFA176
         GO TO 110                                                      SGBFA177
  100    CONTINUE                                                       SGBFA178
            INFO = ki                                                    SGBFA179
  110    CONTINUE                                                       SGBFA180
  120 CONTINUE                                                          SGBFA181
  130 CONTINUE                                                          SGBFA182
      IPVT(N) = N                                                       SGBFA183
      IF (ABD(M,N) .EQ. 0.0E0) INFO = N                                 SGBFA184
      RETURN                                                            SGBFA185
      END                                                               SGBFA186
!
      SUBROUTINE SGESL(A,LDA,N,IPVT,B,JOB)                              SGESL  4
C***BEGIN PROLOGUE  SGESL                                               SGESL  5
C***DATE WRITTEN   780814   (YYMMDD)                                    SGESL  6
C***REVISION DATE  820801   (YYMMDD)                                    SGESL  7
C***CATEGORY NO.  D2A1                                                  SGESL  8
C***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE                       SGESL  9
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                            SGESL 10
C***PURPOSE  SOLVES THE REAL SYSTEM A*X=B OR TRANS(A)*X=B               SGESL 11
C            USING THE FACTORS OF SGECO OR SGEFA                        SGESL 12
C***DESCRIPTION                                                         SGESL 13
C                                                                       SGESL 14
C     SGESL SOLVES THE REAL SYSTEM                                      SGESL 15
C     A * X = B  OR  TRANS(A) * X = B                                   SGESL 16
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.                     SGESL 17
C                                                                       SGESL 18
C     ON ENTRY                                                          SGESL 19
C                                                                       SGESL 20
C        A       REAL(LDA, N)                                           SGESL 21
C                THE OUTPUT FROM SGECO OR SGEFA.                        SGESL 22
C                                                                       SGESL 23
C        LDA     INTEGER                                                SGESL 24
C                THE LEADING DIMENSION OF THE ARRAY  A .                SGESL 25
C                                                                       SGESL 26
C        N       INTEGER                                                SGESL 27
C                THE ORDER OF THE MATRIX  A .                           SGESL 28
C                                                                       SGESL 29
C        IPVT    INTEGER(N)                                             SGESL 30
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.                  SGESL 31
C                                                                       SGESL 32
C        B       REAL(N)                                                SGESL 33
C                THE RIGHT HAND SIDE VECTOR.                            SGESL 34
C                                                                       SGESL 35
C        JOB     INTEGER                                                SGESL 36
C                = 0         TO SOLVE  A*X = B ,                        SGESL 37
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            SGESL 38
C                            TRANS(A)  IS THE TRANSPOSE.                SGESL 39
C                                                                       SGESL 40
C     ON RETURN                                                         SGESL 41
C                                                                       SGESL 42
C        B       THE SOLUTION VECTOR  X .                               SGESL 43
C                                                                       SGESL 44
C     ERROR CONDITION                                                   SGESL 45
C                                                                       SGESL 46
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   SGESL 47
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDES SINGULARITY,SGESL 48
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       SGESL 49
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     SGESL 50
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0           SGESL 51
C        OR SGEFA HAS SET INFO .EQ. 0 .                                 SGESL 52
C                                                                       SGESL 53
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 SGESL 54
C     WITH  P  COLUMNS                                                  SGESL 55
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)                            SGESL 56
C           IF (RCOND IS TOO SMALL) GO TO ...                           SGESL 57
C           DO 10 J = 1, P                                              SGESL 58
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)                        SGESL 59
C        10 CONTINUE                                                    SGESL 60
C                                                                       SGESL 61
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           SGESL 62
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      SGESL 63
C                                                                       SGESL 64
C     SUBROUTINES AND FUNCTIONS                                         SGESL 65
C                                                                       SGESL 66
C     BLAS SAXPY,SDOT                                                   SGESL 67
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    SGESL 68
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   SGESL 69
C***ROUTINES CALLED  SAXPY,SDOT                                         SGESL 70
C***END PROLOGUE  SGESL                                                 SGESL 71
      INTEGER LDA,N,IPVT(*),JOB                                         SGESL 72
!      REAL A(LDA,*),B(*)                                                SGESL 73
      REAL :: A(LDA,*),B(*)                                                SGESL 73
C                                                                       SGESL 74
!      REAL SDOT,T                                                       SGESL 75
      REAL :: SDOT,T                                                       SGESL 75
!      INTEGER ki,KB,L,NM1                                                SGESL 76
      INTEGER :: ki,KB,L,NM1                                                SGESL 76
C***FIRST EXECUTABLE STATEMENT  SGESL                                   SGESL 77
      NM1 = N - 1                                                       SGESL 78
      IF (JOB .NE. 0) GO TO 50                                          SGESL 79
C                                                                       SGESL 80
C        JOB = 0 , SOLVE  A * X = B                                     SGESL 81
C        FIRST SOLVE  L*Y = B                                           SGESL 82
C                                                                       SGESL 83
         IF (NM1 .LT. 1) GO TO 30                                       SGESL 84
         DO 20 ki = 1, NM1                                               SGESL 85
            L = IPVT(ki)                                                 SGESL 86
            T = B(L)                                                    SGESL 87
            IF (L .EQ. ki) GO TO 10                                      SGESL 88
               B(L) = B(ki)                                              SGESL 89
               B(ki) = T                                                 SGESL 90
   10       CONTINUE                                                    SGESL 91
            CALL SAXPY(N-ki,T,A(ki+1,ki),1,B(ki+1),1)                       SGESL 92
   20    CONTINUE                                                       SGESL 93
   30    CONTINUE                                                       SGESL 94
C                                                                       SGESL 95
C        NOW SOLVE  U*X = Y                                             SGESL 96
C                                                                       SGESL 97
         DO 40 KB = 1, N                                                SGESL 98
            ki = N + 1 - KB                                              SGESL 99
            B(ki) = B(ki)/A(ki,ki)                                          SGESL100
            T = -B(ki)                                                   SGESL101
            CALL SAXPY(ki-1,T,A(1,ki),1,B(1),1)                           SGESL102
   40    CONTINUE                                                       SGESL103
      GO TO 100                                                         SGESL104
   50 CONTINUE                                                          SGESL105
C                                                                       SGESL106
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         SGESL107
C        FIRST SOLVE  TRANS(U)*Y = B                                    SGESL108
C                                                                       SGESL109
         DO 60 ki = 1, N                                                 SGESL110
            T = SDOT(ki-1,A(1,ki),1,B(1),1)                               SGESL111
            B(ki) = (B(ki) - T)/A(ki,ki)                                    SGESL112
   60    CONTINUE                                                       SGESL113
C                                                                       SGESL114
C        NOW SOLVE TRANS(L)*X = Y                                       SGESL115
C                                                                       SGESL116
         IF (NM1 .LT. 1) GO TO 90                                       SGESL117
         DO 80 KB = 1, NM1                                              SGESL118
            ki = N - KB                                                  SGESL119
            B(ki) = B(ki) + SDOT(N-ki,A(ki+1,ki),1,B(ki+1),1)                 SGESL120
            L = IPVT(ki)                                                 SGESL121
            IF (L .EQ. ki) GO TO 70                                      SGESL122
               T = B(L)                                                 SGESL123
               B(L) = B(ki)                                              SGESL124
               B(ki) = T                                                 SGESL125
   70       CONTINUE                                                    SGESL126
   80    CONTINUE                                                       SGESL127
   90    CONTINUE                                                       SGESL128
  100 CONTINUE                                                          SGESL129
      RETURN                                                            SGESL130
      END                                                               SGESL131
!
      SUBROUTINE SGEFA(A,LDA,N,IPVT,INFO)                               SGEFA  4
C***BEGIN PROLOGUE  SGEFA                                               SGEFA  5
C***DATE WRITTEN   780814   (YYMMDD)                                    SGEFA  6
C***REVISION DATE  820801   (YYMMDD)                                    SGEFA  7
C***CATEGORY NO.  D2A1                                                  SGEFA  8
C***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX                      SGEFA  9
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                            SGEFA 10
C***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.             SGEFA 11
C***DESCRIPTION                                                         SGEFA 12
C                                                                       SGEFA 13
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.              SGEFA 14
C                                                                       SGEFA 15
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED            SGEFA 16
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          SGEFA 17
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .                   SGEFA 18
C                                                                       SGEFA 19
C     ON ENTRY                                                          SGEFA 20
C                                                                       SGEFA 21
C        A       REAL(LDA, N)                                           SGEFA 22
C                THE MATRIX TO BE FACTORED.                             SGEFA 23
C                                                                       SGEFA 24
C        LDA     INTEGER                                                SGEFA 25
C                THE LEADING DIMENSION OF THE ARRAY  A .                SGEFA 26
C                                                                       SGEFA 27
C        N       INTEGER                                                SGEFA 28
C                THE ORDER OF THE MATRIX  A .                           SGEFA 29
C                                                                       SGEFA 30
C     ON RETURN                                                         SGEFA 31
C                                                                       SGEFA 32
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         SGEFA 33
C                WHICH WERE USED TO OBTAIN IT.                          SGEFA 34
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE      SGEFA 35
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          SGEFA 36
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       SGEFA 37
C                                                                       SGEFA 38
C        IPVT    INTEGER(N)                                             SGEFA 39
C                AN INTEGER VECTOR OF PIVOT INDICES.                    SGEFA 40
C                                                                       SGEFA 41
C        INFO    INTEGER                                                SGEFA 42
C                = 0  NORMAL VALUE.                                     SGEFA 43
C                = ki  IF  U(ki,ki) .EQ. 0.0 .  THIS IS NOT AN ERROR       SGEFA 44
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        SGEFA 45
C                     INDE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO  SGEFA 46
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE   SGEFA 47
C                     INDION OF SINGULARITY.                        SGEFA 48
C                                                                       SGEFA 49
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           SGEFA 50
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      SGEFA 51
C                                                                       SGEFA 52
C     SUBROUTINES AND FUNCTIONS                                         SGEFA 53
C                                                                       SGEFA 54
C     BLAS SAXPY,SSCAL,ISAMAX                                           SGEFA 55
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    SGEFA 56
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   SGEFA 57
C***ROUTINES CALLED  ISAMAX,SAXPY,SSCAL                                 SGEFA 58
C***END PROLOGUE  SGEFA                                                 SGEFA 59
      INTEGER LDA,N,IPVT(*),INFO                                        SGEFA 60
!      REAL A(LDA,*)                                                     SGEFA 61
      REAL :: A(LDA,*)                                                     SGEFA 61
C                                                                       SGEFA 62
!      REAL T                                                            SGEFA 63
      REAL :: T                                                            SGEFA 63
!      INTEGER ISAMAX,J,ki,KP1,L,NM1                                      SGEFA 64
      INTEGER :: ISAMAX,J,ki,KP1,L,NM1                                      SGEFA 64
C                                                                       SGEFA 65
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        SGEFA 66
C                                                                       SGEFA 67
C***FIRST EXECUTABLE STATEMENT  SGEFA                                   SGEFA 68
      INFO = 0                                                          SGEFA 69
      NM1 = N - 1                                                       SGEFA 70
      IF (NM1 .LT. 1) GO TO 70                                          SGEFA 71
      DO 60 ki = 1, NM1                                                  SGEFA 72
         KP1 = ki + 1                                                    SGEFA 73
C                                                                       SGEFA 74
C        FIND L = PIVOT INDEX                                           SGEFA 75
C                                                                       SGEFA 76
         L = ISAMAX(N-ki+1,A(ki,ki),1) + ki - 1                             SGEFA 77
         IPVT(ki) = L                                                    SGEFA 78
C                                                                       SGEFA 79
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          SGEFA 80
C                                                                       SGEFA 81
         IF (A(L,ki) .EQ. 0.0E0) GO TO 40                                SGEFA 82
C                                                                       SGEFA 83
C           INTERCHANGE IF NECESSARY                                    SGEFA 84
C                                                                       SGEFA 85
            IF (L .EQ. ki) GO TO 10                                      SGEFA 86
               T = A(L,ki)                                               SGEFA 87
               A(L,ki) = A(ki,ki)                                          SGEFA 88
               A(ki,ki) = T                                               SGEFA 89
   10       CONTINUE                                                    SGEFA 90
C                                                                       SGEFA 91
C           COMPUTE MULTIPLIERS                                         SGEFA 92
C                                                                       SGEFA 93
            T = -1.0E0/A(ki,ki)                                           SGEFA 94
            CALL SSCAL(N-ki,T,A(ki+1,ki),1)                                SGEFA 95
C                                                                       SGEFA 96
C           ROW ELIMINATION WITH COLUMN INDEXING                        SGEFA 97
C                                                                       SGEFA 98
            DO 30 J = KP1, N                                            SGEFA 99
               T = A(L,J)                                               SGEFA100
               IF (L .EQ. ki) GO TO 20                                   SGEFA101
                  A(L,J) = A(ki,J)                                       SGEFA102
                  A(ki,J) = T                                            SGEFA103
   20          CONTINUE                                                 SGEFA104
               CALL SAXPY(N-ki,T,A(ki+1,ki),1,A(ki+1,J),1)                  SGEFA105
   30       CONTINUE                                                    SGEFA106
         GO TO 50                                                       SGEFA107
   40    CONTINUE                                                       SGEFA108
            INFO = ki                                                    SGEFA109
   50    CONTINUE                                                       SGEFA110
   60 CONTINUE                                                          SGEFA111
   70 CONTINUE                                                          SGEFA112
      IPVT(N) = N                                                       SGEFA113
      IF (A(N,N) .EQ. 0.0E0) INFO = N                                   SGEFA114
      RETURN                                                            SGEFA115
      END                                                               SGEFA116
!
      SUBROUTINE SDNTPB (H,ki,N,NQ,T,TOUT,YH,Y)                          sdntpb 3
C***BEGIN PROLOGUE  SDNTPB                                              sdntpb 4
C***REFER TO  B3                                                        sdntpb 5
C   Subroutine SDNTPB interpolates the ki-th derivative of Y at TOUT,    sdntpb 6
C   using the data in the YH array.  If ki has a value greater than NQ,  sdntpb 7
C   the NQ-th derivative is calculated.                                 sdntpb 8
C***ROUTINES CALLED  (NONE)                                             sdntpb 9
C***DATE WRITTEN   790601   (YYMMDD)                                    sdntpb10
C***REVISION DATE  830920   (YYMMDD)                                    sdntpb11
C***CATEGORY NO.  I1A2,I1A1B                                            sdntpb12
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                sdntpb13
C             INITIAL VALUE PROBLEMS,GEAR'S OD                      sdntpb14
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdntpb15
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdntpb16
C***END PROLOGUE  SDNTPB                                                sdntpb17
!      REAL FACTOR, H, R, T, TOUT, Y(*), YH(N,*)                         sdntpb18
      REAL :: FACTOR, H, R, T, TOUT, Y(*), YH(N,*)                         sdntpb18
C***FIRST EXECUTABLE STATEMENT  SDNTPB                                  sdntpb19
      KUSED = MIN(ki, NQ)                                                sdntpb20
      IF (KUSED.EQ.0) THEN                                              sdntpb21
        DO 10 I = 1,N                                                   sdntpb22
 10       Y(I) = YH(I,NQ+1)                                             sdntpb23
        R = ((TOUT - T)/H)                                              sdntpb24
        DO 20 JJ = 1,NQ                                                 sdntpb25
          J = NQ + 1 - JJ                                               sdntpb26
          DO 20 I = 1,N                                                 sdntpb27
 20         Y(I) = YH(I,J) + R*Y(I)                                     sdntpb28
      ELSE                                                              sdntpb29
        FACTOR = 1.E0                                                   sdntpb30
        DO 40 KK = 1,KUSED                                              sdntpb31
 40       FACTOR = FACTOR*REAL(NQ+1-KK)                                 sdntpb32
        DO 50 I = 1,N                                                   sdntpb33
 50       Y(I) = FACTOR*YH(I,NQ+1)                                      sdntpb34
        IF (KUSED.NE.NQ) THEN                                           sdntpb35
          R = ((TOUT - T)/H)                                            sdntpb36
          DO 80 JJ = KUSED+1,NQ                                         sdntpb37
            J = ki + 1 + NQ - JJ                                         sdntpb38
            FACTOR = 1.E0                                               sdntpb39
            DO 60 KK = 1,KUSED                                          sdntpb40
 60           FACTOR = FACTOR*REAL(J-KK)                                sdntpb41
            DO 70 I = 1,N                                               sdntpb42
 70           Y(I) = FACTOR*YH(I,J) + R*Y(I)                            sdntpb43
 80         CONTINUE                                                    sdntpb44
        END IF                                                          sdntpb45
        do 100 i = 1, N                                                  sdntpb46
 100      Y(i) = Y(i)*H**(-KUSED)                                       sdntpb47
      end if                                                            sdntpb48
      end                                                               sdntpb49
!
      SUBROUTINE SDCSTB (MAXORD,MINT,ISWFLG,EL,TQ)                      sdcstb 3
C***BEGIN PROLOGUE  SDCSTB                                              sdcstb 4
C***REFER TO  B3                                                        sdcstb 5
C  SDCSTB is called by SDNTLB and sets coefficients used by the core    sdcstb 6
C  integrator SDSTPB.  The array EL determines the basic od.        sdcstb 7
C  The array TQ is involved in adjusting the step size in relation      sdcstb 8
C  to truncation error.  EL and TQ depend upon MINT, and are calculated sdcstb 9
C  for orders 1 to MAXORD(.LE.12).  For each order NQ, the coefficients sdcstb10
C  EL are calculated from the generating polynomial:                    sdcstb11
C    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.            sdcstb12
C  For the implicit Adams ods, L(T) is given by                     sdcstb13
C    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/ki,   L(-1) = 0,                 sdcstb14
C    where      ki = factorial(NQ-1).                                    sdcstb15
C  For the Gear ods,                                                sdcstb16
C    L(T) = (1+T)*(2+T)* ... *(NQ+T)/ki,                                 sdcstb17
C    where      ki = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).               sdcstb18
C  For each order NQ, there are three components of TQ.                 sdcstb19
C***ROUTINES CALLED  (NONE)                                             sdcstb20
C***DATE WRITTEN   790601   (YYMMDD)                                    sdcstb21
C***REVISION DATE  830920   (YYMMDD)                                    sdcstb22
C***CATEGORY NO.  I1A2,I1A1B                                            sdcstb23
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                sdcstb24
C             INITIAL VALUE PROBLEMS,GEAR'S OD                      sdcstb25
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdcstb26
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdcstb27
C***END PROLOGUE  SDCSTB                                                sdcstb28
!      REAL EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)              sdcstb29
      REAL :: EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)              sdcstb29
C***FIRST EXECUTABLE STATEMENT  SDCSTB                                  sdcstb30
      FACTRL(1) = 1.E0                                                  sdcstb31
      IF (MAXORD.GE.2) THEN                                             sdcstb32
        DO 10 I = 2,MAXORD                                              sdcstb33
 10       FACTRL(I) = REAL(I)*FACTRL(I-1)                               sdcstb34
      END IF                                                            sdcstb35
C                                             COMPUTE ADAMS COEFFICIENTSsdcstb36
      IF (MINT.EQ.1) THEN                                               sdcstb37
        GAMMA(1) = 1.E0                                                 sdcstb38
        DO 40 I = 1,MAXORD+1                                            sdcstb39
          SUM = 0.E0                                                    sdcstb40
          DO 30 J = 1,I                                                 sdcstb41
 30         SUM = SUM - GAMMA(J)/REAL(I-J+2)                            sdcstb42
 40       GAMMA(I+1) = SUM                                              sdcstb43
        EL(1,1) = 1.E0                                                  sdcstb44
        EL(2,1) = 1.E0                                                  sdcstb45
        EL(2,2) = 1.E0                                                  sdcstb46
        EL(3,2) = 1.E0                                                  sdcstb47
        IF (MAXORD.GE.3) THEN                                           sdcstb48
          DO 60 J = 3,MAXORD                                            sdcstb49
            EL(2,J) = FACTRL(J-1)                                       sdcstb50
            DO 50 I = 3,J                                               sdcstb51
 50           EL(I,J) = REAL(J-1)*EL(I,J-1) + EL(I-1,J-1)               sdcstb52
 60         EL(J+1,J) = 1.E0                                            sdcstb53
        END IF                                                          sdcstb54
        IF (MAXORD.GE.2) THEN                                           sdcstb55
          DO 80 J = 2,MAXORD                                            sdcstb56
            EL(1,J) = EL(1,J-1) + GAMMA(J)                              sdcstb57
            EL(2,J) = 1.E0                                              sdcstb58
            DO 80 I = 3,J+1                                             sdcstb59
 80           EL(I,J) = EL(I,J)/(REAL(I-1)*FACTRL(J-1))                 sdcstb60
        END IF                                                          sdcstb61
        DO 100 J = 1,MAXORD                                             sdcstb62
          TQ(1,J) = -1.E0/(FACTRL(J)*GAMMA(J))                          sdcstb63
          TQ(2,J) = -1.E0/GAMMA(J+1)                                    sdcstb64
 100      TQ(3,J) = -1.E0/GAMMA(J+2)                                    sdcstb65
C                                              COMPUTE GEAR COEFFICIENTSsdcstb66
      ELSE IF (MINT.EQ.2) THEN                                          sdcstb67
        EL(1,1) = 1.E0                                                  sdcstb68
        EL(2,1) = 1.E0                                                  sdcstb69
        IF (MAXORD.GE.2) THEN                                           sdcstb70
          DO 130 J = 2,MAXORD                                           sdcstb71
            EL(1,J) = FACTRL(J)                                         sdcstb72
            DO 120 I = 2,J                                              sdcstb73
 120          EL(I,J) = REAL(J)*EL(I,J-1) + EL(I-1,J-1)                 sdcstb74
 130        EL(J+1,J) = 1.E0                                            sdcstb75
          SUM = 1.E0                                                    sdcstb76
          DO 150 J = 2,MAXORD                                           sdcstb77
            SUM = SUM + 1.E0/REAL(J)                                    sdcstb78
            DO 150 I = 1,J+1                                            sdcstb79
 150          EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)                         sdcstb80
        END IF                                                          sdcstb81
        DO 170 J = 1,MAXORD                                             sdcstb82
          IF (J.GT.1) TQ(1,J) = 1.E0/FACTRL(J-1)                        sdcstb83
          TQ(2,J) = REAL(J+1)/EL(1,J)                                   sdcstb84
 170      TQ(3,J) = REAL(J+2)/EL(1,J)                                   sdcstb85
      END IF                                                            sdcstb86
C                          Compute constants used in the stiffness test.sdcstb87
C                          These are the ratio of TQ(2,NQ) for the Gear sdcstb88
C                          ods to those for the Adams ods.      sdcstb89
      IF (ISWFLG.EQ.3) THEN                                             sdcstb90
        MXRD = MIN(MAXORD, 5)                                           sdcstb91
        IF (MINT.EQ.2) THEN                                             sdcstb92
          GAMMA(1) = 1.E0                                               sdcstb93
          DO 190 I = 1,MXRD                                             sdcstb94
            SUM = 0.E0                                                  sdcstb95
            DO 180 J = 1,I                                              sdcstb96
 180          SUM = SUM - GAMMA(J)/REAL(I-J+2)                          sdcstb97
 190        GAMMA(I+1) = SUM                                            sdcstb98
        END IF                                                          sdcstb99
        IF (MXRD.GE.2) THEN                                             sdcst100
          SUM = 1.E0                                                    sdcst101
          DO 200 I = 2,MXRD                                             sdcst102
            SUM = SUM + 1.E0/REAL(I)                                    sdcst103
 200        EL(1+I,1) = -REAL(I+1)*SUM*GAMMA(I+1)                       sdcst104
        END IF                                                          sdcst105
      END IF                                                            sdcst106
      END                                                               sdcst107
!
      SUBROUTINE SDSCLB (HMAX,N,NQ,RMAX,H,RC,RH,YH)                     
C***BEGIN PROLOGUE  SDSCLB                                              
C***REFER TO  B3                                                    
C   THIS SUBROUTINE RESCALES THE YH ARRAY WHENEVER THE STEP SIZE        
C   IS CHANGED.                                                         
C***ROUTINES CALLED  (NONE)                                             
C***DATE WRITTEN   790601   (YYMMDD)                                    
C***REVISION DATE  830920   (YYMMDD)                                    
C***CATEGORY NO.  I1A2,I1A1B                                            
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD                      
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           
C***END PROLOGUE  SDSCLB                                                
!      REAL H, HMAX, RC, RH, RMAX, R1, YH(N,*)                           
      REAL :: H, HMAX, RC, RH, RMAX, R1, YH(N,*)                           
C***FIRST EXECUTABLE STATEMENT  SDSCLB                                  
      RH = MIN(ABS(H)*RH, HMAX, ABS(H)*RMAX)/ABS(H)                     
      R1 = 1.E0                                                         
      DO 10 J = 1,NQ                                                    
        R1 = R1*RH                                                      
        DO 10 I = 1,N                                                   
 10       YH(I,J+1) = YH(I,J+1)*R1                                      
      H = H*RH                                                          
      RC = RC*RH                                                        
      END                                                               
!
      SUBROUTINE SDCORB (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,   
     8   NDE,NQ,T,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D)                      
C***BEGIN PROLOGUE  SDCORB                                              
C***REFER TO  B3                                                    
C  SUBROUTINE SDCORB IS CALLED TO COMPUTE CORRECTIONS TO THE Y ARRAY.   
C  IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM THE      
C  RESULT OF THE LAST CALL TO F.                                        
C  IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR AND     
C  SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT HAND SIDE AND DFDY AS     
C  COEFFICIENT MATRIX, USING THE LU DECOMPOSITION IF MITER IS 1, 2, 4,  
C  OR 5.                                                                
C***ROUTINES CALLED  SGESL,SGBSL,SNRM2                                  
C***DATE WRITTEN   790601   (YYMMDD)                                    
C***REVISION DATE  830920   (YYMMDD)                                    
C***CATEGORY NO.  I1A2,I1A1B                                            
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD                      
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           
C***END PROLOGUE  SDCORB                                                
!      REAL A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,                
      REAL :: A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,                
     8     SAVE1(*), SAVE2(*), SNRM2, T, Y(*), YH(N,*), YWT(*)          
!      INTEGER IPVT(*)                                                   
      INTEGER :: IPVT(*)                                                   
!      LOGICAL EVALFA                                                    
      LOGICAL :: EVALFA                                                    
C***FIRST EXECUTABLE STATEMENT  SDCORB                                  
      IF (MITER.EQ.0) THEN                                              
        DO 100 I = 1, N                                                  
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)           
        D = SNRM2(N, SAVE1, 1)/SQRT(REAL(N))                            
        DO 105 I = 1,N                                                  
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)                               
      ELSE IF (MITER.EQ.1 .OR. MITER.EQ.2) THEN                         
        IF (IMPL.EQ.0) THEN                                             
! *****
      write(unit = 8, fmt = "(a23, 1pe10.3, a4, i4)") " sdcorb near star
     1t: H =", H, " N =", N
      write(unit = 8, fmt = *) "  i    Save2       YH        Save1"
! *****
          DO 130 I = 1,N                                                
! *****
   !   WRITE(unit = 8, FMT = "(i4, 3e11.3)") i, Save2(i), YH(i,2), 
   !  1  Save1(i)
! *****
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)                  
! *****
   !   WRITE(unit = 8, FMT = "(i4, a8, 1pe11.3)") i, " Save2 =", Save2(i)
! *****
        ELSE IF (IMPL.EQ.1) THEN                                        
          IF (EVALFA) THEN                                              
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                   
          ELSE                                                          
            EVALFA = .TRUE.                                             
          END IF                                                        
          dO 150 i = 1, N                                                
 150        SAVE2(I) = H*SAVE2(I)                                       
          DO 160 J = 1,N                                                
            DO 160 I = 1,N                                              
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))         
        ELSE IF (IMPL.EQ.2) THEN                                        
          IF (EVALFA) THEN                                              
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                   
          ELSE                                                          
            EVALFA = .TRUE.                                             
          END IF                                                        
          DO 180 I = 1,N                                                
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))         
        END IF                                                          
        CALL SGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)                    
! *****
   !   write(unit = 8, fmt = *) " sdcorb after '180':  N =", N
   !   write(unit = 8, fmt = *) "  i   Save1     Save2      YWT"
! *****
        DO 200 I = 1,N                                                  
! *****
   !   WRITE(unit = 8, FMT = "(i4, 3e10.3)") i, Save1(i), Save2(i), 
   !  1  YWT(i)
! *****
          SAVE1(I) = SAVE1(I) + SAVE2(I)                                
 200      SAVE2(I) = SAVE2(I)/YWT(I)                                    
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))                            
      ELSE IF (MITER.EQ.4 .OR. MITER.EQ.5) THEN                         
        IF (IMPL.EQ.0) THEN                                             
          DO 230 I = 1,N                                                
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)                  
        ELSE IF (IMPL.EQ.1) THEN                                        
          IF (EVALFA) THEN                                              
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)           
          ELSE                                                          
            EVALFA = .TRUE.                                             
          END IF                                                        
          DO 250 I = 1,N                                                
 250        SAVE2(I) = H*SAVE2(I)                                       
          MW = ML + 1 + MU                                              
          DO 260 J = 1,N                                                
            I1 = MAX(ML+1, MW+1-J)                                      
            I2 = MIN(MW+N-J, MW+ML)                                     
            DO 260 I = I1,I2                                            
              I3 = I + J - MW                                           
 260          SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))       
        ELSE IF (IMPL.EQ.2) THEN                                        
          IF (EVALFA) THEN                                              
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                   
          ELSE                                                          
            EVALFA = .TRUE.                                             
          END IF                                                        
          DO 280 I = 1,N                                                
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))         
        END IF                                                          
        CALL SGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)            
        DO 300 I = 1,N                                                  
          SAVE1(I) = SAVE1(I) + SAVE2(I)                                
 300      SAVE2(I) = SAVE2(I)/YWT(I)                                    
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))                            
      ELSE IF (MITER.EQ.3) THEN                                         
        IFLAG = 2                                                       
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)                                      
        DO 320 I = 1,N                                                  
          SAVE1(I) = SAVE1(I) + SAVE2(I)                                
 320      SAVE2(I) = SAVE2(I)/YWT(I)                                    
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))                            
      END IF                                                            
      END                                                               
!
      SUBROUTINE SDPSTB (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
     8   NQ,SAVE2,T,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,IER,IPVT,SAVE1,      
     8   ISWFLG,BND)                                                    
C***BEGIN PROLOGUE  SDPSTB                                              
C***REFER TO  B3                                                    
C  SUBROUTINE SDPSTB IS CALLED TO REEVALUATE THE PARTIALS.              
C  IF MITER IS 1, 2, 4, OR 5, THE MATRIX                                
C  P = I - L(0)*H*JACOBIAN IS STORED IN DFDY AND SUBJECTED TO LU        
C  DECOMPOSITION, WITH THE RESULTS ALSO STORED IN DFDY.                 
C***ROUTINES CALLED  SGEFA,SGBFA,SNRM2                                  
C***DATE WRITTEN   790601   (YYMMDD)                                    
C***REVISION DATE  830920   (YYMMDD)                                    
C***CATEGORY NO.  I1A2,I1A1B                                            
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD                      
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           
C***END PROLOGUE  SDPSTB                                                
!      REAL A(MATDIM,*), BND, DFDY(MATDIM,*), DFDYMX, DY,                
      REAL :: A(MATDIM,*), BND, DFDY(MATDIM,*), DFDYMX, DY,                
     8     EL(13,12), EPSJ, EPSJY, FACTOR, H, RDELP, R0, SAVE1(*),      
     8     SAVE2(*), SNRM2, T, UROUND, Y(*), YH(N,*), YJ, YWT(*)        
!      INTEGER IPVT(*)                                                   
      INTEGER :: IPVT(*)                                                   
!      LOGICAL IER                                                       
      LOGICAL :: IER                                                       
      DATA RDELP /1000.E0/                                              
C***FIRST EXECUTABLE STATEMENT  SDPSTB                                  
      NJE = NJE + 1                                                     
      IER = .FALSE.                                                     
      IF (MITER.EQ.1 .OR. MITER.EQ.2) THEN                              
        IF (MITER.EQ.1) THEN                                            
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)                   
          IF (ISWFLG.EQ.3) BND = SNRM2(N*N, DFDY, 1)                    
          FACTOR = -EL(1,NQ)*H                                          
          DO 110 J = 1,N                                                
            DO 110 I = 1,N                                              
 110          DFDY(I,J) = FACTOR*DFDY(I,J)                              
        ELSE IF (MITER.EQ.2) THEN                                       
          R0 = ABS(H)*RDELP*UROUND*SNRM2(N, SAVE2, 1)                   
          EPSJ = SQRT(UROUND)                                           
          DO 170 J = 1,N                                                
            EPSJY = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))                    
            IF (EPSJY.GT.R0) THEN                                       
              DY = EPSJY                                                
            ELSE IF (R0.NE.0.E0) THEN                                   
              DY = R0                                                   
            ELSE                                                        
              DY = UROUND                                               
            END IF                                                      
            YJ = Y(J)                                                   
            Y(J) = Y(J) + DY                                            
            CALL F (N, T, Y, SAVE1)                                     
            Y(J) = YJ                                                   
            FACTOR = -EL(1,NQ)*H/DY                                     
            DO 140 I = 1,N                                              
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR                  
 170        CONTINUE                                                    
          IF (ISWFLG.EQ.3) BND = SNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)      
          NFE = NFE + N                                                 
        END IF                                                          
        IF (IMPL.EQ.0) THEN                                             
          DO 190 I = 1,N                                                
 190        DFDY(I,I) = DFDY(I,I) + 1.E0                                
        ELSE IF (IMPL.EQ.1) THEN                                        
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                     
          DO 210 J = 1,N                                                
            DO 210 I = 1,N                                              
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)                            
        ELSE IF (IMPL.EQ.2) THEN                                        
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                     
          DO 230 I = 1,NDE                                              
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)                              
        END IF                                                          
        CALL SGEFA (DFDY, MATDIM, N, IPVT, INFO)                        
        IF (INFO.NE.0) IER = .TRUE.                                     
      ELSE IF (MITER.EQ.4 .OR. MITER.EQ.5) THEN                         
        IF (MITER.EQ.4) THEN                                            
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)           
          FACTOR = -EL(1,NQ)*H                                          
          MW = ML + MU + 1                                              
          DO 260 J = 1,N                                                
            I1 = MAX(ML+1, MW+1-J)                                      
            I2 = MIN(MW+N-J, MW+ML)                                     
            DO 260 I = I1,I2                                            
 260          DFDY(I,J) = FACTOR*DFDY(I,J)                              
        ELSE IF (MITER.EQ.5) THEN                                       
          R0 = ABS(H)*RDELP*UROUND*SNRM2(N, SAVE2, 1)                   
          EPSJ = SQRT(UROUND)                                           
          MW = ML + MU + 1                                              
          J2 = MIN(MW, N)                                               
          DO 340 J = 1,J2                                               
            DO 290 ki = J,N,MW                                           
              EPSJY = EPSJ*MAX(ABS(YWT(ki)), ABS(Y(ki)))                  
              IF (EPSJY.GT.R0) THEN                                     
                DY = EPSJY                                              
              ELSE IF (R0.NE.0.E0) THEN                                 
                DY = R0                                                 
              ELSE                                                      
                DY = UROUND                                             
              END IF                                                    
              DFDY(MW,ki) = Y(ki)                                         
 290          Y(ki) = Y(ki) + DY                                          
            CALL F (N, T, Y, SAVE1)                                     
            DO 330 ki = J,N,MW                                           
              Y(ki) = DFDY(MW,ki)                                         
              EPSJY = EPSJ*MAX(ABS(YWT(ki)), ABS(Y(ki)))                  
              IF (EPSJY.GT.R0) THEN                                     
                DY = EPSJY                                              
              ELSE IF (R0.NE.0.E0) THEN                                 
                DY = R0                                                 
              ELSE                                                      
                DY = UROUND                                             
              END IF                                                    
              FACTOR = -EL(1,NQ)*H/DY                                   
              I1 = MAX(ML+1, MW+1-ki)                                    
              I2 = MIN(MW+N-ki, MW+ML)                                   
              DO 300 I = I1,I2                                          
                I3 = ki + I - MW                                         
 300            DFDY(I,ki) = FACTOR*(SAVE1(I3) - SAVE2(I3))              
 330          CONTINUE                                                  
 340        CONTINUE                                                    
          NFE = NFE + J2                                                
        END IF                                                          
        IF (ISWFLG.EQ.3) THEN                                           
          DFDYMX = 0.E0                                                 
          DO 345 J = 1,N                                                
            I1 = MAX(ML+1, MW+1-J)                                      
            I2 = MIN(MW+N-J, MW+ML)                                     
            DO 345 I = I1,I2                                            
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))                      
          BND = 0.E0                                                    
          IF (DFDYMX.NE.0.E0) THEN                                      
            DO 350 J = 1,N                                              
              I1 = MAX(ML+1, MW+1-J)                                    
              I2 = MIN(MW+N-J, MW+ML)                                   
              DO 350 I = I1,I2                                          
 350            BND = BND + (DFDY(I,J)/DFDYMX)**2                       
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)                        
          END IF                                                        
        END IF                                                          
        IF (IMPL.EQ.0) THEN                                             
          DO 360 J = 1,N                                                
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0                              
        ELSE IF (IMPL.EQ.1) THEN                                        
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)             
          DO 380 J = 1,N                                                
            I1 = MAX(ML+1, MW+1-J)                                      
            I2 = MIN(MW+N-J, MW+ML)                                     
            DO 380 I = I1,I2                                            
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)                            
        ELSE IF (IMPL.EQ.2) THEN                                        
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                     
          DO 400 J = 1,NDE                                              
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)                           
        END IF                                                          
        CALL SGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)                
        IF (INFO.NE.0) IER = .TRUE.                                     
      ELSE IF (MITER.EQ.3) THEN                                         
        IFLAG = 1                                                       
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)                                      
      END IF                                                            
      END                                                               
!
      SUBROUTINE SDPSCB (KSGN,N,NQ,YH)                                  
C***BEGIN PROLOGUE  SDPSCB                                              
C***REFER TO  B3                                                    
C     THIS SUBROUTINE COMPUTES THE PREDICTED YH VALUES BY EFFECTIVELY   
C     MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX WHEN KSGN  
C     IS +1, AND PERFORMS THE INVERSE FUNCTION WHEN KSGN IS -1.         
C***ROUTINES CALLED  (NONE)                                             
C***DATE WRITTEN   790601   (YYMMDD)                                    
C***REVISION DATE  830920   (YYMMDD)                                    
C***CATEGORY NO.  I1A2,I1A1B                                            
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,                
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD                      
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           
C***END PROLOGUE  SDPSCB                                                
!      REAL YH(N,*)                                                      
      REAL :: YH(N,*)                                                      
C***FIRST EXECUTABLE STATEMENT  SDPSCB                                  
      IF (KSGN.GT.0) THEN                                               
        DO 10 J1 = 1,NQ                                                 
          DO 10 J2 = J1,NQ                                              
            J = NQ - J2 + J1                                            
            DO 10 I = 1,N                                               
 10           YH(I,J) = YH(I,J) + YH(I,J+1)                             
      ELSE                                                              
        DO 30 J1 = 1,NQ                                                 
          DO 30 J2 = J1,NQ                                              
            J = NQ - J2 + J1                                            
            DO 30 I = 1,N                                               
 30           YH(I,J) = YH(I,J) - YH(I,J+1)                             
      END IF                                                            
      END                                                               
!
      SUBROUTINE SDNTLB (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,   sdntlb 3
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,Y,YWT,H,MNTOLD,MTROLD,NFE,RC,YH,sdntlb 4
     8   A,CONVRG,EL,IER,IPVT,NQ,NWAIT,RH,RMAX,SAVE2,TQ,TREND,ISWFLG)   sdntlb 5
C***BEGIN PROLOGUE  SDNTLB                                              sdntlb 6
C***REFER TO  B3                                                        sdntlb 7
C  Subroutine SDNTLB is called to set parameters on the first call      sdntlb 8
C  to SDSTPB, on an internal restart, or when the user has altered      sdntlb 9
C  MINT, MITER, and/or H.                                               sdntlb10
C  On the first call, the order is set to 1 and the initial derivatives sdntlb11
C  are calculated.  RMAX is the maximum ratio by which H can be         sdntlb12
C  increased in one step.  It is initially RMINIT to compensate         sdntlb13
C  for the small initial H, but then is normally equal to RMNORM.       sdntlb14
C  If a failure occurs (in corrector convergence or error test), RMAX   sdntlb15
C  is set at RMFAIL for the next increase.                              sdntlb16
C  If the caller has changed MINT, or if JTASK = 0, SDCSTB is called    sdntlb17
C  to set the coefficients of the method.  If the caller has changed H, sdntlb18
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is        sdntlb19
C  reset to NQ + 2 to prevent further increases in H for that many      sdntlb20
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of  sdntlb21
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is      sdntlb22
C  set to 0 to force the partials to be updated, if partials are used.  sdntlb23
C***ROUTINES CALLED  SDCSTB,SDSCLB,SGEFA,SGESL,SGBFA,SGBSL,SNRM2        sdntlb24
C***DATE WRITTEN   790601   (YYMMDD)                                    sdntlb25
C***REVISION DATE  841119   (YYMMDD)                                    sdntlb26
C***CATEGORY NO.  I1A2,I1A1B                                            sdntlb27
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,               sdntlb28
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY           sdntlb29
C***END PROLOGUE  SDNTLB                                                sdntlb30
!      REAL A(MATDIM,*), EL(13,12), EPS, H, HMAX, HNEW, HOLD,            sdntlb31
      REAL :: A(MATDIM,*), EL(13,12), EPS, H, HMAX, HNEW, HOLD,            sdntlb31
     8     OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX, SMIN, sdntlb32
     8     SNRM2, SUM, SUM0, T, TQ(3,12), TREND, Y(*), YH(N,*), YWT(*)  sdntlb33
!      INTEGER IPVT(*)                                                   sdntlb34
      INTEGER :: IPVT(*)                                                   sdntlb34
!      LOGICAL CONVRG, IER                                               sdntlb35
      LOGICAL :: CONVRG, IER                                               sdntlb35
      DATA RMINIT /10000.E0/                                            sdntlb36
! *****
      NSpeci = 329           ! 329 = NSpeci (see write statement below)
! *****
C***FIRST EXECUTABLE STATEMENT  SDNTLB                                  sdntlb37
      IER = .FALSE.                                                     sdntlb38
! *****
   !   write(unit = 8, fmt = *) "At sdntlb 39:  JTASK =", JTASK, " call s
   !  1dcstb"
! *****
! *****
   !   write(unit = 8, fmt = *) " JTASK = ", JTASK
   !   write(unit = *, fmt = *) " JTASK = ", JTASK
!      pause
! *****
      IF (JTASK .GE. 0) THEN                                            sdntlb39
        IF (JTASK .EQ. 0) CALL SDCSTB (MAXORD, MINT, ISWFLG,  EL, TQ)   sdntlb40
        RC = 0.E0                                                       sdntlb41
        CONVRG = .FALSE.                                                sdntlb42
        TREND = 1.E0                                                    sdntlb43
        RMAX = RMINIT                                                   sdntlb44
        NQ = 1                                                          sdntlb45
        NWAIT = 3                                                       sdntlb46
! *****
   !   write(unit = *, fmt = *) " sdntlb46:  call F"
!      pause
   !   write(unit = 8, fmt = *) " sdntlb 46: befor call F:  RC =", RC, 
   !  1  " Convrg =", Convrg, " Trend =", Trend, " Rmax =", Rmax, 
   !  2  " NQ =", NQ, " NWait =", NWait
! *****
! *****
   !   write(unit = 8, fmt = *) " N =", N, " T =", T
! *****
      write(unit = 8, fmt = *) " sdntlb 46: before call F"
      write(unit = 8, fmt = *) "       Y      Save2 = YP"
 !     write(unit = *, fmt = *) " NSpeci = ", NSpeci
      do i = 1, NSpeci          ! XXXXXXXXXXXXXXXX Should be NSpeci? XXXXXXXXXX
        write(unit = 8, fmt = "(i4, 1x, 1pe10.3, e10.3)") i, Y(i), 
     1    Save2(i)
      end do
! *****
   !   write(unit = *, fmt = *)  " N = ", N, " T = ", T!, " Y = ", Y, 
!     1  " Save2 = ", Save2
        CALL F (N, T, Y, SAVE2)                                         sdntlb47
   !   write(unit = *, fmt = *)  " After sdntlb47:  Call F" 
! *****
      write(unit = 8, fmt = *) " sdntlb 47: after call F"
! *****
        NFE = NFE + 1                                                   sdntlb48
! *****
   !   write(unit = *, fmt = *) " sdntlb48: IMPL =", IMPL, " MITER =", 
   !  1   MITER, " NFE =", NFE
! *****
        IF (IMPL .NE. 0) THEN                                           sdntlb49
          IF (MITER .EQ. 3) THEN                                        sdntlb50
            IFLAG = 0                                                   sdntlb51
C           CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,    sdntlb52
C    8                  NDE, IFLAG)                                     sdntlb53
          ELSE IF (IMPL .EQ. 1) THEN                                    sdntlb54
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN                    sdntlb55
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                 sdntlb56
              CALL SGEFA (A, MATDIM, N, IPVT, INFO)                     sdntlb57
              IF (INFO .NE. 0) THEN                                     sdntlb58
                IER = .TRUE.                                            sdntlb59
                RETURN                                                  sdntlb60
              END IF                                                    sdntlb61
              CALL SGESL (A, MATDIM, N, IPVT, SAVE2, 0)                 sdntlb62
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN               sdntlb63
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)         sdntlb64
              CALL SGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)             sdntlb65
              IF (INFO .NE. 0) THEN                                     sdntlb66
                IER = .TRUE.                                            sdntlb67
                RETURN                                                  sdntlb68
              END IF                                                    sdntlb69
              CALL SGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)         sdntlb70
            END IF                                                      sdntlb71
          ELSE IF (IMPL .EQ. 2) THEN                                    sdntlb72
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)                   sdntlb73
            DO 150 I = 1,NDE                                            sdntlb74
              IF(A(I,1) .EQ. 0.E0) THEN                                 sdntlb75
                IER = .TRUE.                                            sdntlb76
                RETURN                                                  sdntlb77
              ELSE                                                      sdntlb78
                SAVE2(I) = SAVE2(I)/A(I,1)                              sdntlb79
              END IF                                                    sdntlb80
 150          CONTINUE                                                  sdntlb81
            DO 155 I = NDE+1,N                                          sdntlb82
 155          A(I,1) = 0.E0                                             sdntlb83
          END IF                                                        sdntlb84
        END IF                                                          sdntlb85
! *****
   !   write(unit = 8, fmt = *) " sdntlb 85: Here IMPL = 0.  NDE =", NDE
! *****
        DO 170 I = 1,NDE                                                sdntlb86
 170      SAVE1(I) = SAVE2(I)/YWT(I)                                    sdntlb87
        SUM = SNRM2(NDE, SAVE1, 1)                                      sdntlb88
        SUM0 = 1.E0/MAX(1.E0, ABS(T))                                   sdntlb89
        SMAX = MAX(SUM0, SUM)                                           sdntlb90
        SMIN = MIN(SUM0, SUM)                                           sdntlb91
        SUM = SMAX*SQRT(1.E0 + (SMIN/SMAX)**2)/SQRT(REAL(NDE))          sdntlb92
        H = SIGN(MIN(2.E0*EPS/SUM, ABS(H)), H)                          sdntlb93
! *****
   !   write(unit = 8, fmt = "(a16, 1pe10.3, a4, i4)") " sdntlb 94: H =",
   !  1  H, " N =", N
! *****
        DO 180 I = 1,N                                                  sdntlb94
 180      YH(I,2) = H*SAVE2(I)                                          sdntlb95
! *****
   !   do i = 1, N
   !   write(unit = 8, fmt = "(a4, i4, a5, 1pe10.3, a8, e10.3)") 
   !  1  " i =", i, " YH =", YH(i,2), " Save2 =", Save2(i)
   !   end do
   !   write(unit = 8, fmt = *)
! *****
      ELSE                                                              sdntlb96
        IF (MITER .NE. MTROLD) THEN                                     sdntlb97
          MTROLD = MITER                                                sdntlb98
          RC = 0.E0                                                     sdntlb99
          CONVRG = .FALSE.                                              sdntl100
        END IF                                                          sdntl101
        IF (MINT .NE. MNTOLD) THEN                                      sdntl102
          MNTOLD = MINT                                                 sdntl103
          OLDL0 = EL(1,NQ)                                              sdntl104
          CALL SDCSTB (MAXORD, MINT, ISWFLG,  EL, TQ)                   sdntl105
          RC = RC*EL(1,NQ)/OLDL0                                        sdntl106
          NWAIT = NQ + 2                                                sdntl107
        END IF                                                          sdntl108
        IF (H .NE. HOLD) THEN                                           sdntl109
          NWAIT = NQ + 2                                                sdntl110
          HNEW = H                                                      sdntl111
          RH = H/HOLD                                                   sdntl112
          H = HOLD                                                      sdntl113
          CALL SDSCLB (HMAX, N, NQ, RMAX,  H, RC, RH, YH)               sdntl114
          H = SIGN(MIN(ABS(H), ABS(HNEW)), H)                           sdntl115
        END IF                                                          sdntl116
      END IF                                                            sdntl117
      END                                                               sdntl118
!
      function sdot(n,x,n1,y,n2)
      dimension x(*),y(*)
      if(n1.ne.1.or.n2.ne.1) stop 25
      sdot = 0.0
      do 10 i=1,n
      sdot = sdot + x(i)*y(i)
   10 continue
      return
      end
!
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)                            SAXPY  4
C***BEGIN PROLOGUE  SAXPY                                               SAXPY  5
C***DATE WRITTEN   791001   (YYMMDD)                                    SAXPY  6
C***REVISION DATE  820801   (YYMMDD)                                    SAXPY  7
C***CATEGORY NO.  D1A7                                                  SAXPY  8
C***KEYWORDS  BLAS,LINEAR ALGEBRA,TRIAD,VECTOR                          SAXPY  9
C***AUTHOR  LAWSON, C. L., (JPL)                                        SAXPY 10
C           HANSON, R. J., (SNLA)                                       SAXPY 11
C           KINCAID, D. R., (U. OF TEXAS)                               SAXPY 12
C           KROGH, F. T., (JPL)                                         SAXPY 13
C***PURPOSE  S.P. COMPUTATION Y = A*X + Y                               SAXPY 14
C***DESCRIPTION                                                         SAXPY 15
C                                                                       SAXPY 16
C                B L A S  SUBPROGRAM                                    SAXPY 17
C    DESCRIPTION OF PARAMETERS                                          SAXPY 18
C                                                                       SAXPY 19
C     --INPUT--                                                         SAXPY 20
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       SAXPY 21
C       SA  SINGLE PRECISION SCALAR MULTIPLIER                          SAXPY 22
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     SAXPY 23
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      SAXPY 24
C       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS                     SAXPY 25
C     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY                      SAXPY 26
C                                                                       SAXPY 27
C     --OUTPUT--                                                        SAXPY 28
C       SY  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)             SAXPY 29
C                                                                       SAXPY 30
C     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.    SAXPY 31
C     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +  SAXPY 32
C       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N SAXPY 33
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.                  SAXPY 34
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     SAXPY 35
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, SAXPY 36
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       SAXPY 37
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 SAXPY 38
C***ROUTINES CALLED  (NONE)                                             SAXPY 39
C***END PROLOGUE  SAXPY                                                 SAXPY 40
C                                                                       SAXPY 41
!      REAL SX(*),SY(*),SA                                               SAXPY 42
      REAL :: SX(*),SY(*),SA                                               SAXPY 42
C***FIRST EXECUTABLE STATEMENT  SAXPY                                   SAXPY 43
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN                                   SAXPY 44
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60                               SAXPY 45
    5 CONTINUE                                                          SAXPY 46
C                                                                       SAXPY 47
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.                   SAXPY 48
C                                                                       SAXPY 49
      IX = 1                                                            SAXPY 50
      IY = 1                                                            SAXPY 51
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 SAXPY 52
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 SAXPY 53
      DO 10 I = 1,N                                                     SAXPY 54
        SY(IY) = SY(IY) + SA*SX(IX)                                     SAXPY 55
        IX = IX + INCX                                                  SAXPY 56
        IY = IY + INCY                                                  SAXPY 57
   10 CONTINUE                                                          SAXPY 58
      RETURN                                                            SAXPY 59
C                                                                       SAXPY 60
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            SAXPY 61
C                                                                       SAXPY 62
C                                                                       SAXPY 63
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.   SAXPY 64
C                                                                       SAXPY 65
   20 M = MOD(N,4)                                                      SAXPY 66
      IF( M .EQ. 0 ) GO TO 40                                           SAXPY 67
      DO 30 I = 1,M                                                     SAXPY 68
        SY(I) = SY(I) + SA*SX(I)                                        SAXPY 69
   30 CONTINUE                                                          SAXPY 70
      IF( N .LT. 4 ) RETURN                                             SAXPY 71
   40 MP1 = M + 1                                                       SAXPY 72
      DO 50 I = MP1,N,4                                                 SAXPY 73
        SY(I) = SY(I) + SA*SX(I)                                        SAXPY 74
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)                            SAXPY 75
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)                            SAXPY 76
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)                            SAXPY 77
   50 CONTINUE                                                          SAXPY 78
      RETURN                                                            SAXPY 79
C                                                                       SAXPY 80
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.                  SAXPY 81
C                                                                       SAXPY 82
   60 CONTINUE                                                          SAXPY 83
      NS = N*INCX                                                       SAXPY 84
          DO 70 I=1,NS,INCX                                             SAXPY 85
          SY(I) = SA*SX(I) + SY(I)                                      SAXPY 86
   70     CONTINUE                                                      SAXPY 87
      RETURN                                                            SAXPY 88
      END                                                               SAXPY 89
!
      SUBROUTINE SSCAL(N,SA,SX,INCX)                                    SSCAL  4
C***BEGIN PROLOGUE  SSCAL                                               SSCAL  5
C***DATE WRITTEN   791001   (YYMMDD)                                    SSCAL  6
C***REVISION DATE  820801   (YYMMDD)                                    SSCAL  7
C***CATEGORY NO.  D1A6                                                  SSCAL  8
C***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR                          SSCAL  9
C***AUTHOR  LAWSON, C. L., (JPL)                                        SSCAL 10
C           HANSON, R. J., (SNLA)                                       SSCAL 11
C           KINCAID, D. R., (U. OF TEXAS)                               SSCAL 12
C           KROGH, F. T., (JPL)                                         SSCAL 13
C***PURPOSE  S.P. VECTOR SCALE X = A*X                                  SSCAL 14
C***DESCRIPTION                                                         SSCAL 15
C                                                                       SSCAL 16
C                B L A S  SUBPROGRAM                                    SSCAL 17
C    DESCRIPTION OF PARAMETERS                                          SSCAL 18
C                                                                       SSCAL 19
C     --INPUT--                                                         SSCAL 20
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       SSCAL 21
C       SA  SINGLE PRECISION SCALE FACTOR                               SSCAL 22
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     SSCAL 23
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      SSCAL 24
C                                                                       SSCAL 25
C     --OUTPUT--                                                        SSCAL 26
C       SX  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)             SSCAL 27
C                                                                       SSCAL 28
C     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.            SSCAL 29
C     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)    SSCAL 30
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     SSCAL 31
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, SSCAL 32
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       SSCAL 33
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 SSCAL 34
C***ROUTINES CALLED  (NONE)                                             SSCAL 35
C***END PROLOGUE  SSCAL                                                 SSCAL 36
C                                                                       SSCAL 37
!      REAL SA,SX(*)                                                     SSCAL 38
      REAL :: SA,SX(*)                                                     SSCAL 38
C***FIRST EXECUTABLE STATEMENT  SSCAL                                   SSCAL 39
      IF(N.LE.0)RETURN                                                  SSCAL 40
      IF(INCX.EQ.1)GOTO 20                                              SSCAL 41
C                                                                       SSCAL 42
C        CODE FOR INCREMENTS NOT EQUAL TO 1.                            SSCAL 43
C                                                                       SSCAL 44
      NS = N*INCX                                                       SSCAL 45
          DO 10 I = 1,NS,INCX                                           SSCAL 46
          SX(I) = SA*SX(I)                                              SSCAL 47
   10     CONTINUE                                                      SSCAL 48
      RETURN                                                            SSCAL 49
C                                                                       SSCAL 50
C        CODE FOR INCREMENTS EQUAL TO 1.                                SSCAL 51
C                                                                       SSCAL 52
C                                                                       SSCAL 53
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.   SSCAL 54
C                                                                       SSCAL 55
   20 M = MOD(N,5)                                                      SSCAL 56
      IF( M .EQ. 0 ) GO TO 40                                           SSCAL 57
      DO 30 I = 1,M                                                     SSCAL 58
        SX(I) = SA*SX(I)                                                SSCAL 59
   30 CONTINUE                                                          SSCAL 60
      IF( N .LT. 5 ) RETURN                                             SSCAL 61
   40 MP1 = M + 1                                                       SSCAL 62
      DO 50 I = MP1,N,5                                                 SSCAL 63
        SX(I) = SA*SX(I)                                                SSCAL 64
        SX(I + 1) = SA*SX(I + 1)                                        SSCAL 65
        SX(I + 2) = SA*SX(I + 2)                                        SSCAL 66
        SX(I + 3) = SA*SX(I + 3)                                        SSCAL 67
        SX(I + 4) = SA*SX(I + 4)                                        SSCAL 68
   50 CONTINUE                                                          SSCAL 69
      RETURN                                                            SSCAL 70
      END                                                               SSCAL 71
!
      INTEGER FUNCTION ISAMAX(N,SX,INCX)                                ISAMAX 4
C***BEGIN PROLOGUE  ISAMAX                                              ISAMAX 5
C***DATE WRITTEN   791001   (YYMMDD)                                    ISAMAX 6
C***REVISION DATE  820801   (YYMMDD)                                    ISAMAX 7
C***CATEGORY NO.  D1A2                                                  ISAMAX 8
C***KEYWORDS  BLAS,LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR              ISAMAX 9
C***AUTHOR  LAWSON, C. L., (JPL)                                        ISAMAX10
C           HANSON, R. J., (SNLA)                                       ISAMAX11
C           KINCAID, D. R., (U. OF TEXAS)                               ISAMAX12
C           KROGH, F. T., (JPL)                                         ISAMAX13
C***PURPOSE  FIND LARGEST COMPONENT OF S.P. VECTOR                      ISAMAX14
C***DESCRIPTION                                                         ISAMAX15
C                                                                       ISAMAX16
C                B L A S  SUBPROGRAM                                    ISAMAX17
C    DESCRIPTION OF PARAMETERS                                          ISAMAX18
C                                                                       ISAMAX19
C     --INPUT--                                                         ISAMAX20
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)                       ISAMAX21
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS                     ISAMAX22
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX                      ISAMAX23
C                                                                       ISAMAX24
C     --OUTPUT--                                                        ISAMAX25
C   ISAMAX  SMALLEST INDEX (ZERO IF N .LE. 0)                           ISAMAX26
C                                                                       ISAMAX27
C     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.  ISAMAX28
C     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX) ISAMAX29
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,     ISAMAX30
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*, ISAMAX31
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL       ISAMAX32
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323 ISAMAX33
C***ROUTINES CALLED  (NONE)                                             ISAMAX34
C***END PROLOGUE  ISAMAX                                                ISAMAX35
C                                                                       ISAMAX36
!      REAL SX(*),SMAX,XMAG                                              ISAMAX37
      REAL :: SX(*),SMAX,XMAG                                              ISAMAX37
C***FIRST EXECUTABLE STATEMENT  ISAMAX                                  ISAMAX38
      ISAMAX = 0                                                        ISAMAX39
      IF(N.LE.0) RETURN                                                 ISAMAX40
      ISAMAX = 1                                                        ISAMAX41
      IF(N.LE.1)RETURN                                                  ISAMAX42
      IF(INCX.EQ.1)GOTO 20                                              ISAMAX43
C                                                                       ISAMAX44
C        CODE FOR INCREMENTS NOT EQUAL TO 1.                            ISAMAX45
C                                                                       ISAMAX46
      SMAX = ABS(SX(1))                                                 ISAMAX47
      NS = N*INCX                                                       ISAMAX48
      II = 1                                                            ISAMAX49
          DO 10 I=1,NS,INCX                                             ISAMAX50
          XMAG = ABS(SX(I))                                             ISAMAX51
          IF(XMAG.LE.SMAX) GO TO 5                                      ISAMAX52
          ISAMAX = II                                                   ISAMAX53
          SMAX = XMAG                                                   ISAMAX54
    5     II = II + 1                                                   ISAMAX55
   10     CONTINUE                                                      ISAMAX56
      RETURN                                                            ISAMAX57
C                                                                       ISAMAX58
C        CODE FOR INCREMENTS EQUAL TO 1.                                ISAMAX59
C                                                                       ISAMAX60
   20 SMAX = ABS(SX(1))                                                 ISAMAX61
      DO 30 I = 2,N                                                     ISAMAX62
         XMAG = ABS(SX(I))                                              ISAMAX63
         IF(XMAG.LE.SMAX) GO TO 30                                      ISAMAX64
         ISAMAX = I                                                     ISAMAX65
         SMAX = XMAG                                                    ISAMAX66
   30 CONTINUE                                                          ISAMAX67
      RETURN                                                            ISAMAX68
      END                                                               ISAMAX69
!
      SUBROUTINE USERS(R1,R2,R3,R4,R5,R6,R7,R8,I1,I2,I3,I4)
      write(8, 1)
    1 FORMAT(' ','Subroutine Users called')
      RETURN
      END
!

