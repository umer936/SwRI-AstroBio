      program EqWrt
!
! This program needs a blank file 'SpecReac' in Debug at the start of 
! any new project.
!
! This program generates FORTRAN source code for subroutines 'FChem', 
! 'Jacobn', and 'Rates', and writes them to FORTRAN file 'EWS.f' (EWS 
! = Equation Write Subroutines) in the projects for the chemistry of 
! the Solar Nebula, Comet Comae, or Titan Atmosphere, depending on the 
! project name at the beginning of the 'Species' and 'Reactn' input 
! files.  Project names should be 6 characters long.  
!
! Quantities in common blocks /LC1/ and /LC2/ are written to the 
! file 'DatFil' for use in subroutines 'RateCn' and 'Photo' in some 
! user programs, such as 'ComChem'.  'Datfil' is not needed in Solar 
! Nebula program.
!
! Input consists of the number of species participating in the 
! problem, 'NSpeci', and the number of ACTIVE chemical species = the 
! number of differential equations, 'NLIWDE', followed by the set of all 
! symbols, together with their stoichiometric coefficients in the file 
! 'Species'.  Next is the set of reactions, in the file 'Reactn'.  If 
! the number of reactions and/or the number of species is increased, 
! the dimensions (arrays) of files 'iWork(LiW)' and 'CWork(LCW)' should 
! be increased (see: 'integer, parameter' below).  Also 'iCRTbl(LTabl)' 
! may have to be increased (see 'integer, parameter').  The size of 
! 'iCRTbl' is more difficult to predict.  Increase 'LTabl', then run 
! the program and chck the 'Length of the Cross Reference Table' near 
! the end of file 'EqnOut' (or in 'DatFil' (but it is not well labled 
! there; it is the 5th value on first line).  Then reduce 'LTabl' if 
! necessary.
!
! Make sure there are no blank lines at the end of files 'Species' 
! and 'Reactn'!
!
! NSpeci = Number of chemical species, including catalysts (M). Hot 
!          hydrogens [H(h), H2(h)] should be included if they are in a 
!          chemical reaction, however, if escape of hot hydrogen is 
!          calculated only from the tail part of the speed distribution 
!          function then they should be excluded [Huebner and Keady, 
!          First-Flight Escape from Spheres with R**(-2) Density 
!          Distribution.  Astron. Astrophys. 135, 177-180 (1984)].
! NDE    = Number of differential equations to be solved: 
!          excludes hot hydrogens and catalysts (M).
! NReac  = Number of chemical Reactions to be solved.
! NSym   = Number of reactants and products in a reaction (typically 
!          NSym = 7, 3 reactants + 4 products).
           
!
! Expanded, modified, and up-dated:  June 2007 and November 2014 
! by W. F. Huebner.
!
      integer, parameter :: MxReac = 8000, NSym = 7, MaxSpec = 600, 
     1                      MxStoich = 18, LenTbl = 19838, LIW = 160000,
     2                      LCW = 120000 
! Set LIW = 20*NReac, LCW = 15*NReac
! Set LIW = 20*MxReac, LCW = 15*MxReac
      dimension WtM(MaxSpec), gma(MaxSpec), EX(MxReac), Numb(MxReac)
      dimension NAt(MaxSpec), NDegF(MaxSpec), T3(MaxSpec), P3(MaxSpec) 
      dimension TLo(MxReac), THi(MxReac), iEqns(7)
!      dimension BF(MxReac), CF(MxReac), iEqns(7)
      character SYMBP*8, Eqns*8, iRMLim*9, Quelle*8
      character Cat*4
      character*18 Comment
      character CWork(LCW)*1
      integer   jStoich(MxStoich, MaxSpec), iSum(MxStoich), iWork(LIW)
      data LOut /7/, NStoich /18/
      common /LC1/ A(MxReac), B(MxReac), C(MxReac), RMLim(MxReac),
     8       Cat(MxReac), iStoich(MxStoich, MaxSpec), iCRTbl(LenTbl)
      common /LC2/ SymbP(MaxSpec), Eqns(NSym, MxReac), Quelle(MxReac),
     8       iRMLim(MxReac)
!      Common /FC/ MaxSpec, MxReac
!
      character(len = 8) :: Project, ProjR
      character(len = 92) :: ColHeads
      character(len = 7) :: FName
  !    character(len = 63) :: Fmt1      ! \Documents and Settings[ ... \Codes\
  !    character(len = 23) :: Fmt2      ! Enceladus\Torus\.....or equivalent 
  !    character(len = 34) :: Fmt3      ! Enceladus\Torus\bin\Debug\.....
  !    character(len = 82) :: FmtEWS    ! EWS = Equation Write Subroutines:  
!                                          FChem, Jacobn, Rates
  !    character(len = 96) :: FMTSRT    ! SRT = Species, Reactn, iCRTbl file 
!                                          (input file for equation solver 
!                                          program)
      character(len = 102) :: RColHead1
      character*6 Prog, Speci, Symbol, SpecIn, x
      character(len = 46) :: Dummy
      character(len = 92) :: RColHead2
!
      open(unit =  1, file = "SpeciesSN1", status = "old")
      open(unit =  2, fiLE = "ReactnSN1",  status = "old")
! 8888888888888888
! unit = 6 IS THE SCREEN!
! 8888888888888888
c
c This does not change the default unit (*) in some versions of
c the runtime system. A better solution for the existing code is
c command line redirection.
c      open(unit = 6, file = "EqnOut",  status = "replace")
c
      open(unit = 60, file = "DatFil",  status = "replace") ! Not needed for SolNeb
      write(unit = *, fmt = *) " LOut = ", LOut
! ***
  !    open(file = "FmtEWS", unit = LOut, status = "replace") ! See lower down
c      open(file = "EWS", unit = LOut, status = "replace") ! See lower down
!      open(unit = 3, file = "FmtSRT", status = "replace")
      open(unit = 3, file = "SRT", status = "replace")
! ***
!
      read(unit = 1, fmt = "(a8)") Project                 ! 'Species'
      write(unit = *, fmt = *) "Project = ", Project
      read(unit = 1, fmt = "(a92)") ColHeads               ! 'Species'
      write(unit = *, fmt = *) ColHeads
      read(unit = 2, fmt = "(a8)") ProjR                   ! 'Reactn'
      write(unit = *, fmt = *) "ProjR = ", ProjR
      if(Project .ne. ProjR) then
        write(unit = *, fmt = *) "'Reactn' file = ", Project, "is not ",
     1    "consistent with 'Species' file = ", ProjR
      stop
      else
      end if
!
 !     fmt1 = "C:\Documents and Settings\Walter F. Huebner\My Documents\C
 !    1odes\"
      if(Project == "SolNeb01" .or. Project == "SolNeb02") then     ! Solar Nebula
        FName = "EWS.f"                ! EWS = EqWrt generated subroutines 
  !      fmt2 = "SolNeb\"
  !      fmt2(8:12) = FName
  !      fmt3 = "SolNeb\bin\Debug\SpecReac"
      else
        if(Project == "Comets") then     ! Comet Coma
          FName = "EWS.f"                ! EWS = EqWrt generated subroutines 
  !        fmt2 = 
  !   1      "Comets\"
  !        fmt2(8:12) = FName
  !        fmt3 = "Comets\bin\Debug\SpecReac"
        else
          if(Project == "EncTor") then     ! Enceladus Torus
            FName = "EWS.f"                ! EWS = EqWrt generated subroutines 
  !          fmt2 = "EncTor\"
  !          fmt2(8:12) = FName
  !          fmt3 = "EncTor\bin\Debug\SpecReac" 
          else
            write(unit = *, fmt = *) "Need project names at beginning ", 
     1        "of files 'Reactn' and 'Species'."
            stop
          end if
        end if
      end if
!
  !    write(unit = *, fmt = *) "fmt1 = ", fmt1
  !    write(unit = *, fmt = *) "fmt2 = ", fmt2
  !    write(unit = *, fmt = *) "fmt3 = ", fmt3
!
  !    FmtEWS(1:63) = fmt1
  !    FmtEWS(64:84) = fmt2
c
c      write(unit = *, fmt = *) "FmtEWS = ", FmtEWS
c
!      FmtSRT(62:96) = fmt3
  !    FmtSRT(62:86) = fmt3
!      write(unit = *, fmt = *) "     * * * * * * * * * * "
!      write(unit = *, fmt = *) "Subroutines FChem, Jacobn, and Rates ",
!     1  "have been written to file: ", FmtEWS
!      write(unit = *, fmt = *) "     * * * * * * * * * * "
!      write(unit = *, fmt = *) "A file containing 'Species', ", 
!     1  "'Reactn', and 'iCRTbl' (cross reference table), has been ", 
!     2  "written to file: ", FmtSRT
!      write(unit = *, fmt = *) "     * * * * * * * * * * "
  !    write(unit = *, fmt = *) "FmtEWS = ", FmtEWS
! ***      
  !    open(unit = LOut, file = "FmtEWS", status = "replace")
      open(unit = LOut, file = FName, status = "replace")
! ***
!
! Read the set of all possible species symbols and the corresponding 
! stoichiometric coefficients.
!
      read(unit = 1, fmt = "(5i4, f5.2)") NReac, NSpeci, NDE, NHot, 
     1  iPrnt, rh                                                  ! 'Species'
!
      write(unit = *, fmt = "(1x, 5i4, f5.2)") NReac, NSpeci, NDE, 
     1  NHot, iPrnt, rh                                            ! 'Species'
!
      write(unit = 6, fmt = "(5x, 2i4)") NSpeci, NDE               ! 'EqnOut'
!
!      write(unit = *, fmt = "(5x, 2i4)") NSpeci, NDE              ! 'EqnOut'
!
      write(unit = 3, fmt = "(5i4, f5.2)") NReac, NSpeci, NDE, NHot, 
     1  iPrnt, rh
!
c      rewind(unit = 3)
c      nReac = 0
c      read(unit = 3, fmt = "(5i4, f5.2)") NReac, NSpeci, NDE, NHot, 
c     1  iPrnt, rh
c      write(unit = *, fmt = "(5i4, f5.2)") NReac, NSpeci, NDE, NHot, 
c     1  iPrnt, rh
c      pause
      i = 1
! Read Chemical species, charge, stochiometric coefficients, molecular mass, 
!     number of atoms in molecule, degrees of freedom at low temperature, 
!     triple point temperature (T3 [K]) and pressuren (P3 [kPa]).
   10 read(unit = 1, fmt = "(a8, 18i3, f9.4, i3, i2, f7.2, e10.3)",  ! Species
     1  end = 15)   !     CHECK
     2  SymbP(i), (jStoich(j, i), j = 1, NStoich), WtM(i), NAt(i), 
     3  NDegF(i), T3(i), P3(i)
!     gma(i) deleted, 'Species'
      write(unit = 3, fmt = "(a8, 18i3, f9.4, i3, i2, f7.2, e10.3)") 
     1  SymbP(i), (jStoich(j, i), j = 1, NStoich), WtM(i), NAt(i), 
     3  NDegF(i), T3(i), P3(i) 
c
c      write(unit = *, fmt = "(a6, 18i3, f7.3, f9.4)") 
c     1  SymbP(i), (JStoich(j, i), j = 1, NStoich), gma(i), WTM(i)
c      pause
c   
      iSymbP = iChar(SymbP(i)(1:1))
C
!      if(i .gt. 20 .and. i .lt. 30) then
!        write(unit = *, fmt = "(1x, a8, i4)") SymbP(i), iSymbP
!        pause
!      else
!      end if
C
      write(unit = 6, fmt = "(1x, i6, 18i3, f9.4)") iSymbP,
     1  (JStoich(j, i), j = 1, NStoich), WtM(i)            ! gma(i) deleted, 'EqnOut'
      i = i + 1
      go to 10
   15 NSPTOT = i - 1
C
      write(unit = *, fmt = *) "NSpTot = ", NSpTot
C
!
! Read the set of reactions
!
      read(unit = 2, fmt = "(a102)")  RColHead1            ! Reactn
      write(unit = *, fmt = "(1x, a102)")  RColHead1
!      pause
      read(unit = 2, fmt = "(a92)")  RColHead2             ! Reactn
      write(unit = *, fmt = "(1x, a92)")  RColHead2
!      pause
      iEr = 0
!
      i = 1
!   60 read(unit = 2, fmt = "(7a8, 1pe8.2, 0pf8.4, 0pf7.1, 0pf7.4, 
!     1  0pf7.4, 0pf8.2, i5, a1, i2)", end = 230)
  !    A        B      C        Lo     Hi  
  !    5.00E-10  0.000     0.9     10. 41000. AD  UMIST12
   60 read(unit = 2, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
     1  0pf7.0, a4, a8)", end = 230)                                 ! 'Reactn'
     2  (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i)
 !     write(unit = *, fmt = "(1x, 7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
 !    1  0pf7.0, a4, a8)")                                            ! 'Reactn'
 !    2  (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
 !    3  Cat(i), Quelle(i)
!     3  iCat(i), Quelle(i)
!   60 read(unit = 2, fmt = "(7a8, 1pe8.2, 0pf8.4, 0pf7.1, 0pf8.0, 
!     1  0pf7.0, a4, a8)", end = 230)
!     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  EX(i), NUMB(i), Quelle(i), ICat(i)                         ! 'Reactn'
! *****************
! Note different format for "C" coefficient (0pf7.0 vs. 0pf7.1) in the 
! following write statements: 
! *****************
  !    write(unit = *, fmt = *) i, Cat(i)
  !    pause
!      if(icat(i) == 28 .or. icat(i) == 41 .or. icat(i) == 45 .or. 
!     1  icat(i) == 46 .or. icat(i) == 47) then
!      if(cat(i) == 28 .or. cat(i) == 41 .or. cat(i) == 45 .or. 
!     1  cat(i) == 46 .or. cat(i) == 47) then
      if(cat(i) == "28" .or. cat(i) == "41" .or. cat(i) == "45" .or. 
     1  cat(i) == "46" .or. cat(i) == "47") then
!        write(unit = 3, fmt = "(7A6, 1PE8.2, 0PF8.4, 0pF7.0, 0PF7.4, 
!     1  0PF7.4, 0PF8.2, I5, A1, I2)")
        write(unit = 3, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
     1  0pf7.0, a4, a8)")
     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i)
      else
  !      write(unit = 3, fmt = "(7a8, 1PE8.2, 0PF8.4, 0pF7.1, 0PF8.0, 
  !   1  0PF7.0, a4, a8)")
!        write(unit = 3, fmt = "(7A6, 1PE8.2, 0PF8.4, 0pF7.1, 0PF8.0, 
!     1  0PF7.4, 0PF8.2, I5, A1, I2)")
  !   2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
  !   3  CAT(i), Quelle(i)
!     3  EX(i), NUMB(i), Quelle(i), ICAT(i)
        write(unit = 3, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
     1  0pf7.0, a4, a8)")
     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i)
      end if
      do 14 j = 1, NSym
        iEqns(j) = iChar(Eqns(j, i)(1:1))
   14 continue
      write(unit = 6, fmt = "(4x, 7(1x, i8))") (iEqNs(j), j = 1, NSym) !'EqnOut'
!
! Missing symbol (typo error) check.  Element and charge balance check.
!
c
C      write(unit = *, fmt = *) ' NStoich = ', NStoich, ' NSym = ', NSym, 
C     1  ' NSpTot = ', NSpTot
c
      do 80 j = 4, NSym 
        if(Eqns(j, i) .ne. "        " ) go to 110
   80 continue
      do 100 j = 1, 3
        if(Eqns(j, i) .ne. "        ") then
          do 90 l = 1, NSpTot
            if(Eqns(j, i) .eq. SymbP(l)) then
              go to 100
            else
            end if
   90     continue
          write(unit = 6, fmt = "(a37)") "Symbol not in the input symbol
     1 table:"
          write(unit = 6, fmt = "(i4, 7(1x, a8), i2)") i,              ! 'EqnOut'
     1      (Eqns(l, i), l = 1, NSym), 1
          iEr = 1
          write (unit = *, fmt = *) " iER = 1.1"
          stop
          go to 60
        end if
  100 continue
      go to 180
!
  110 do 120 k = 1, NStoich
  120   iSum(k) = 0
      do 160 j = 1, NSym
        if(Eqns(j, i) .ne. "        ") then
          do 130 l = 1, NSpTot
            if(Eqns(j, i) .eq. SymbP(l)) then
              go to 140
            else
            end if
  130     continue
          write(unit = 6, fmt = "(a37)") "Symbol not in the input symbol
     1 table:"
          write(unit = *, fmt = "(i4, 7(1x, a8), i2)") i,              ! 'EqnOut'
     1      (Eqns(l, i), l = 1, NSym), 2
c
          iEr = 1
          write (unit = *, fmt = *) " iER = 1.2"
          stop
          go to 60
  140     do 150 k = 1, NStoich
            if(j .le. 3) then
              ISum(k) = ISum(k) + JStoich(k, l)
            else
              ISUM(K) = ISUM(K) - JStoich(K,L)
            end if
  150     continue
        end if
  160 continue
      do 170 k = 1, NStoich
        if(ISum(k) .ne. 0) then
          write(unit = 6, fmt = "(a15)") "Atom imbalance:"
          write(unit = 6, fmt = "(i4, 7(1x, a8))") i, 
     1      (EqNs(l, i), l = 1, NSym)                               ! 'EqnOut'
!          write(6, '(/ '' Atom imbalance: '' / I8, 8(1X, A6))')
!     1      i, (EQNS(L, i), L = 1, NSym)                           ! 'EqnOut'

          IER = 1
          write (unit = *, fmt = *) " iER = 1.3"
          stop
          go to 60
        end if
  170 continue
!
  180 do 200 j = 1, NSym
        if(Eqns(j, i) .ne. ' ') then
          do 190 L = 1, NSPECI
            if(EQNS(j, i) .EQ. SYMBP(L)) go to 200
  190     continue
          GO TO 60
        end if
  200 continue
!
!      if(ICat(i) .eq. 27 .or. ICat(i) .eq. 28) then
!      if(Cat(i) .eq. 27 .or. Cat(i) .eq. 28) then
      if(Cat(i) .eq. "27" .or. Cat(i) .eq. "28") then
        if(RMLIM(i) .eq. 0.) then
          IRMLIM(i) = ' '
!     R1MACH(2) temporarily ignored
          RMLIM(i) = .5*R1MACH(2)
!
        else
          write(IRMLIM(i), '(''('', 1pe7.1, '')'')') RMLIM(i)
        end if
      else
        IRMLIM(i) = ' '
      end if
! Note different format for "C" coefficient in the following write statements: 
!      if(ICat(i) .ne. 5) then
!      if(Cat(i) .ne. 5) then
      if(Cat(i) .ne. "5") then
!     "(7a8, 1pe8.2, 0pf8.4, 0pf7.1, 0pf8.0, 
!     1  0pf7.0, a4, a8)", end = 230)
!     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  EX(i), NUMB(i), Quelle(i), ICat(i) 
        write(unit = 6, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
     1  0pf7.0, a4, a8)")
     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i)
!        write(unit = 6, fmt = "(7(1x, a8), 1pe8.2, 0pf8.4, 0pf7.1, 
!     1    0pf8.0, 0pf7.0, a4, a8)")
!     2    (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3    Cat(i), Quelle(i)
!      write(6, '(7a6, 1pe8.2, 0pf8.4, 0pf8.1, 0pf7.4, 0pf7.4, 0pf8.2, 
!     1  i5, a1, i2)')
!     2  (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  EX(i), Numb(i), Quelle(i), iCat(i)                         ! 'EqnOut'
      else
        write(unit = 6, fmt = "(7a8, 1pe8.2, 0pf7.3, 0pf9.2, 0pf7.0, 
     1  0pf7.0, a4, a8)")
     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), TLo(i), THi(i), 
     3  Cat(i), Quelle(i)
!      write(unit = 6, fmt = "(7(1x, a8), 1pe8.2, 0pf8.4, 0pf7.1, 
!     1  0pf8.0, 0pf7.0, a4, a8)")
!     2  (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  Cat(i), Quelle(i)                                          ! 'EqnOut'
!      write(6, '(7a6, 1pe8.2, 0pf8.4, 0pf7.1, 0pf7.4, 0pf7.4, 0pf8.2, 
!     1  i5, a1, i2)')
!     2  (Eqns(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  EX(i), NUMB(i), Quelle(i), iCat(i)                         ! 'EqnOut'
      end if
      i = i + 1
      go to 60
!
  230 if(iEr .ne. 0) then
        write(unit = *, fmt = *) " iEr .ne. 0"
c       pause
      else
      end if
      NReac = i - 1
! **********
!      write(unit = *, fmt = "(1x, 7a8, 1pe8.2, 0pf8.4, 0pf7.1, 0pf8.0, 
!     1  0pf7.0, a4, a8)")
!     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
!     3  Cat(i), Quelle(i)                                           ! 'Reactn'
! **********
!     3  EX(i), NUMB(i), Quelle(i), ICat(i)                         ! 'Reactn'
!      Stop
! *****
!      write(unit = *, fmt = *) "NReac = ", NReac, " 1 "
! *****
!      stop
! XXXXXXX
      write(unit = 6, fmt = "(a24, i5)") " Number of reactions is ", 
     1  NReac                                                      ! 'EqnOut'
c
      write(unit = 6, fmt = *) " NSym = ", NSym, "NReac = ", NReac, 
     1  "NDE = ", NDE, "iCRTbl = ", iCRTbl, "LenTbl = ", LenTbl
      write(unit = 6, fmt = *) " call Diagnos", " iEr = ", iEr
      write(unit = *, fmt = *) " call Diagnos", " iEr = ", iEr
c     pause
c
      call Diagnos(Eqns, SymbP, NSym, NReac, NDE, iCRTbl, LenTbl)
! *****
        write(unit = 6, fmt = *) " SymbP = ", (SymbP(m), m = 8, 12),
     1  "after Diagnos" 
! *****
      write(unit = 6, fmt = "(a40, i5)") " LENGTH OF THE CROSS REFERENCE
     1 TABLE IS ", LenTbl                                          ! 'EqnOut'
      write(6, *) " Diagnos: ", NSpeci, NDE, NStoich, NReac, LenTbl
      write(unit = *, fmt = *)
     1  " Diagnos: ", NSpeci, NDE, NStoich, NReac, LenTbl
     
  !    write(60, *) NSpeci, NDE, NStoich, NReac, LenTbl
  !    write(unit = *, fmt = *)
!  'FilDat = tape 60' is not needed for SolNeb program
! !     write(60, '(12A6)') (SYMBP(i), i = 1, NSPECI)
  !    write(60, '(12a8)') (SymbP(i), i = 1, NSpeci)
  !    write(unit = *, fmt = "(1x, 12a8)")
  !   1  (SymbP(i), i = 1, 10)
  !    pause
  !   
  !    write(60, *) ((iStoich(i, j), i = 1, NStoich), j = 1, NSpeci)
  !    write(unit = *, fmt = *)
  !   1  ((iStoich(i, j), i = 1, 5), j = 1, 5)
  !    pause
  !   
  !    write(60, *) (A(i), i = 1, NReac)
  !    write(unit = *, fmt = *)
  !   1  (A(i), i = 1, 10)
  !    pause
  !   
  !    write(60, *) (B(i), i = 1, NReac)
  !    write(unit = *, fmt = *)
  !   1  (B(i), i = 1, 10)
  !    pause
  !   
  !    write(60, *) (C(i), i = 1, NReac)
  !    write(unit = *, fmt = *)
  !   1  (C(i), i = 1, 10)
  !    pause
  !   
  !    write(60, *) (RMLim(i), i = 1, NReac)
  !    write(unit = *, fmt = *)
  !   1  (RMLim(i), i = 1, 10)
  !    pause
  !   
  !    write(60, '(12a8)') ((Eqns(i, j), i = 1, NSym), j = 1, NReac)
!      write(60, '(12A6)') ((EQNS(i, j), i = 1, NSym), j = 1, NREAC)
  !    write(unit = *, fmt = "(1x, 12a8)")
  !   1  ((EQNS(i, j), i = 1, NSym), j = 1, 10)
  !    pause
  !   
  !    write(60, '(8A9)') (IRMLIM(i), i = 1, NREAC)
!      write(60, *) (ICAT(i), i = 1, NREAC)
  !    write(unit = *, fmt = "(1x, 8A9)")
  !   1  (IRMLIM(i), i = 1, 10)
  !   
  !    write(60, *) (Cat(i), i = 1, NReac)
  !    write(unit = *, fmt = *)
  !   1  (Cat(i), i = 1, 10)
  !    pause
  !   
  !    write(60, '(4a8)') (Quelle(i), i = 1, NReac)
! !     write(60, '(4A22)') (Quelle(i), i = 1, NREAC)
  !    write(unit = *, fmt = "(1x, 4a8)")
  !   1  (Quelle(i), i = 1, 10)
  !    pause
  !   
  !    write(60, '(i5)') LenTbl
  !    write(unit = *, fmt = "(1x, i5)")
  !   1  LenTbl
  !    pause
  !   
  !    write(60, '(16i5)') (iCRTbl(i), i = 1, LenTbl)
  !    write(unit = *, fmt = "(1x, 16i5)")
  !   1  (iCRTbl(i), i = 1, 10)
  !    pause
  !    
      write(unit = 3, fmt = "(i5)") LenTbl
      write(6, '(i5)') LenTbl                                      ! 'EqnOut'
!      write(LOut, '(i5)') LenTbl
      write(unit = 3, fmt = "(16i5)") (iCRTbl(i), i = 1, LenTbl)
      write(6, '(16I5)') (iCRTbl(i), i = 1, LENTBL)                ! 'EqnOut'
!      write(LOUT, '(16I5)') (iCRTbl(i), i=1,LENTBL)
      NProd = NSym - 3
c
c     pause
c      write(unit = *, fmt = "(a6, 16i3, f7.3, f9.4)") 
c     1  SymbP(i), (JStoich(j, i), j = 1, NStoich), gma(i), WTM(i)
c      pause
c        write(unit = *, fmt = "(7A6, 1PE8.2, 0PF8.4, 0pF7.0, 0PF7.4, 
c     1  0PF7.4, 0PF8.2, I5, A1, I2)")
c     2  (EQNS(j, i), j = 1, NSym), A(i), B(i), C(i), BF(i), CF(i), 
c     3  EX(i), NUMB(i), Quelle(i), ICAT(i)
c      write(unit = *, fmt = *) " call EqnWrt"
c      pause
c
! *****
!        write(unit = 6, fmt = *) " SymbP = ", (SymbP(m), m = 8, 12),
!     1  "before EQnWrt" 
! *****
!      stop
      call EqnWrt(NDE, NSpeci, NReac, NProd, SymbP, Eqns,
     8            LOut, iWork, LiW, CWork, LCW)
c
c     write(unit = *, fmt = *) " call Rates"
c
! *****
        write(unit = 6, fmt = *) " SymbP = ", (SymbP(m), m = 8, 12),
     1  "before Rates" 
! *****
      call Rates(Eqns, SymbP, NSym, NReac, NSpeci, LOut, CWork, LCW)
!      pause
      close(unit = 3)
      close(unit = 6)                                              ! 'EqnOut'
      close(unit = LOut)                                           ! unit = 7
!      close(unit = 7)
!      close(unit = 60)                                             ! 'DatFil'
      write(unit = *, fmt = *) " Finished"
      stop
      end program
!
      subroutine Diagnos(Eqns, SymbP, NSym, NReac, NDE, iCRTbl, LenTbl)
!      integer iCRTbl(*)
      integer iCRTbl(LenTbl)
!      write(unit = *, fmt = *) "iCRTbl = ", iCRTbl
!      pause
!      character SymbP(NDE)*6, Eqns(NSym, NReac)*6
      character SymbP(NDE)*8, Eqns(NSym, NReac)*8
!
      write(unit = 6, fmt = *) " In Diagnos"
! *****
        write(unit = 6, fmt = *) " NDE = ", NDE
        do iSy = 1, NDE
        write(unit = 6, fmt = "(a9)") SymbP(iSy)
!     1  "in Diagnos" 
        end do
! *****
!      write(unit = *, fmt = *) " Eqns = ", Eqns, "SymbP = ", SymbP, 
!     1  "NSym = ", NSym, "NReac = ", NReac, "NDE = ", NDE, "iCRTl = ", 
!     2  iCRTbl, "LenTbl = ", LenTbl
!      pause
      iPoint = 1
      do 40 i = 1, NDE
        NEntry = 0
        do 30 k = 1, NReac
!        do 30 k = 1, 1445
          iCount = 0
!          Write(unit = 6, fmt = *) " Eqns = ", (Eqns(j, k), j = 1, 7) 
          do 10 j = 1, 3
            if(Eqns(j, k) .eq. SymbP(i)) then
              iCount = iCount - 1
! *****
         write(unit = 6, fmt = *) "k= ", k, "Eqns(j, k) =", Eqns(j, k), 
     1   "Symbp(i) = ", (SymbP(in), in = 8, 12), " 1"
! *****
            else
            end if
   10     continue
          do 20 j = 4, NSym
            if(Eqns(j, k) .eq. SymbP(i)) then 
              iCount = iCount + 1
!              write(unit = *, fmt = *) " i = ", i, SymbP(i), " iCount"
!              pause
            else
            end if
   20     continue
          if(iCount .ne. 0) then
            NEntry = NEntry + 1
            iCRTbl(iPoint + NEntry) = sign(k, iCount)
!            write(unit = 6, fmt = *)  " k = ", k, " iCount = ", iCount, 
!     1        "iPoint = ", iPoint, "NEntry = ", NEntry
          else
!            write(unit = 6, fmt = *) " iCount = ", iCount
          end if
!          stop
   30   continue
        iCRTbl(iPoint) = NEntry
         write(unit = 6, fmt = "(3x, a6, i5)") SymbP(i), iCRTbl(iPoint)
         write(unit = *, fmt = *) SymbP(i), iCRTbl(iPoint)
   40   iPoint = iPoint + NEntry + 1
!     LenTbl = iPoint - 1
c     pause
      return
      end
!
      subroutine Rates(Eqns, SymbP, NSym, NReac, NSpeci, Lout, CWork, 
     8                 LCW)
!      subroutine Rates(MxReac, Eqns, SymbP, NSym, NReac, NSpeci, Lout, 
!     8                 CWork, LCW)
!      character Eqns(NSym, NReac)*6, SymbP(NSpeci)*6, CWork(LCW)*1
      character Eqns(NSym, NReac)*8, SymbP(NSpeci)*8, CWork(LCW)*1
!
      write(unit = LOut, fmt = "(6x, a19)") "subroutine Rates(y)"
!      write(unit = LOut, fmt = "(6x, a27)") "subroutine Rates(MxReac, y)
!     1"
! *****
!     Put MxReac in the calling sequence or in a common statement ?
        write(unit = 6, fmt = *) " SymbP = ", (SymbP(m), m = 8, 12)
! *****
      write(unit = LOut, fmt = "(6x, a25)") "parameter (MxReac = 8000)"
c      write(unit = LOut, fmt = "(6x, a25)") "parameter (MaxEqs = 600)"
!      write(unit = LOut, fmt = "(6x, a9, i4)") "MxReac = ", MxReac
      write(unit = LOut, fmt = "(6x, a14)") "dimension y(*)"
      write(unit = LOut, fmt = "(6x, a24)") "common/stcom9/ k(MxReac)"
      write(unit = LOut, fmt = "(6x, a25)") "common/stcom10/ R(MxReac)"
c      write(unit = LOut, fmt = "(6x, a24)") "common/stcom9/ k(MaxEqs)"
c      write(unit = LOut, fmt = "(6x, a25)") "common/stcom10/ R(MaxEqs)"
!      write(unit = LOut, '(6x, ''subroutine Rates(y, r)'')')
      IP = 1
c
c     write(unit = *, fmt = *) " call EqsPlt 1"
c
      call Eqsplt(CWork, IP, NSpeci, LCW)
      CWORK(IP) = ')'
      CWORK(IP+1) = ','
      CWORK(IP+2) = ' '
      CWORK(IP+3) = 'R'
      CWORK(IP+4) = '('
      IP = IP + 5
c
c     write(unit = *, fmt = *) " call EqsPlt 2"
c
      CALL EQSPLT (CWORK, IP, NREAC, LCW)
      CWORK(IP) = ')'
      write(unit = LOut, fmt = "(6x, a6/ a1)") "real k", "!"
!      write(unit = LOut, fmt = "(6x, a24/ a1)") 
!     !  "double precision R, k, y", "!"
!      write(LOUT, '(6X, ''REAL K, Y('', 56A1)') (CWORK(i), i=1,IP)
      IP = 1
c
c     write(unit = *, fmt = *) " call EqsPlt 3"
c
      CALL EqSPLT(CWORK, IP, NREAC, LCW)
      CWORK(IP) = ')'
!      write(LOUT, '(6X, ''common /SDRCOM/ K('', 48A1)')
!     8  (CWORK(i), i=1,IP)
      do 60 j = 1, NREAC
        CWORK(1) = 'R'
        CWORK(2) = '('
        IC = 3
c
c     write(unit = *, fmt = *) " call EqsPlt 4"
c
        call EqSPLT(CWORK, IC, j, LCW)
        CWORK(IC) = ')'
        CWORK(IC+1) = '='
        CWORK(IC+2) = 'K'
        CWORK(IC+3) = '('
        IC = IC + 4
c
c     write(unit = *, fmt = *) " call EqsPlt 5"
c
       CALL EQSPLT (CWORK, IC, j, LCW)
        CWORK(IC) = ')'
        do 40 k = 1, 3
          if(Eqns(k, j) .ne. ' ') then
            do 20 l = 1, NSpeci
              if(Eqns(k, j) .eq. SymbP(l)) go to 30
   20         continue
            write(unit = 6, fmt = *) "j = ", j, (Eqns(n, j), n = 1, 3), 
     1        " SymbP = ", (SymbP(m), m = 8, 12)
            write(6, '(// '' for reaction '', I4, '', there is a symb'',
     8      ''ol not in the input symbol table: '', a8 //)')
     8      j, Eqns(K, j)
            stop
   30       CWORK(IC+1) = '*'
            CWORK(IC+2) = 'Y'
            CWORK(IC+3) = '('
            IC = IC + 4
c
c     write(unit = *, fmt = *) " call EqsPlt 6"
c
            CALL EQSPLT (CWORK, IC, L, LCW)
            CWORK(IC) = ')'
          end if
   40     continue
   60   write(LOUT, '(6X, 66A1)') (CWORK(i), i=1,IC)
      write(LOut, '(6x, ''return'' / 6x, ''end'')')
      return
      end
!
      subroutine EqnWrt (NDE, NSpeci, NReac, NProd, Speci, Reac,
     8                   LOut, iWork, LiW, CWork, LCW)
!     8                   LOut, iWork, LeniW, CWork, LenCW)
! *****
!      write(unit = *, fmt = *) " NReac = ", NReac, "Test"
! *****
C  Begin prologue EqnWrt
C
C  THIS PROGRAM IS TO BE USED WITH THE FORTRAN 77 COMPILER.
C  C. D. SUTHERLAND, LOS ALAMOS NATIONAL LABORATORY, MARCH 21, 1983
C
C  i.   ABSTRACT  ............................................
C
C  THIS PROGRAM GENERATES TWO FORTRAN ROUTINES USED IN INTEGRATING THE
C  RATE EQUATIONS OF CHEMICAL KINETICS, NUCLEAR REACTIONS, ETC.  THE
C  ROUTINES ARE subroutine F, WHICH EVALUATES THE RIGHT HAND SIDE OF THE
C  DIFFERENTIAL EQUATIONS INVOLVED, AND subroutine JACOBN, WHICH
C  EVALUATES THE JACOBIAN MATRIX OF THE RIGHT HAND SIDE.  THE
C  subroutineS F AND JACOBN HAVE A CALL SEQUENCE COMPATIBLE WITH THE
C  REQUIREMENTS OF THE O.D.E. SOLVERS SDRIV1,2,3 (LANL PROGRAM
C  NUMBER D208).  THESE subroutineS CAN BE MADE COMPATIBLE WITH OTHER
C  O.D.E. SOLVERS BY SUITABLE MODIFICATION.
C
C  II.  PARAMETERS  ....................................
C
C    THE PARAMETERS IN THE CALL SEQUENCE ARE:
C
C    NDE    = (INPUT) THE NUMBER OF ACTIVE SPECIES, AND THUS THE NUMBER
C             OF DIFFERENTIAL EQUATIONS.
C
C    NSPECI = (INPUT) THE TOTAL NUMBER OF SPECIES APPEARING IN THE
C             REACTIONS, NSPECI .GE. NDE.  THE USER CAN REFER TO SPECIES
C             WHICH PARTICIPATE IN REACTIONS, BUT DO NOT CHANGE WITH
C             TIME, BY SETTING NSPECI GREATER THAN NDE.  THIS IS A
C             CONVENIENT WAY OF REFERRING TO THE THIRD BODY IN THREE
C             BODY REACTIONS.
C
C    NREAC  = (INPUT) THE NUMBER OF REACTIONS TO BE PROCESSED.
C
C    NPROD  = (INPUT) THE MAXIMUM NUMBER OF PRODUCT SPECIES IN ANY
C             REACTION.
C
C    SPECI  = (INPUT) A character ARRAY OF LENGTH NSPECI CONTAINING THE
C             SYMBOLS OF THE SPECIES THAT APPEAR IN THE REACTIONS.  THE
C             SYMBOLS FOR THE ACTIVE SPECIES MUST APPEAR FIRST, FOLLOWED
C             BY THE SYMBOLS FOR THE INACTIVE SPECIES.  THIS LIST OF
C             SPECIES DETERMINES AN ORDERING.  THUS IN THE OUTPUT
C             subroutineS F AND JACOBN, Y(i) AND YP(i) REFER TO THE i-TH
C             ELEMENT IN THE SPECI ARRAY.
C
C    REAC   = (INPUT) A (3+NPROD) BY NREAC character ARRAY CONTAINING
C             THE SYMBOLS OF THE REACTANT AND PRODUCT SPECIES FOR EACH
C             REACTION.  THE FIRST THREE COMPONENTS OF EACH COLUMN ARE
C             RESERVED FOR THE SYMBOLS OF THE REACTANTS.  THE NEXT NPROD
C             COMPONENTS (BEGINNING AT THE FOURTH POSITION) ARE RESERVED
C             FOR THE SYMBOLS OF THE PRODUCT SPECIES.  ANY COMPONENT
C             WHICH DOES NOT CONTAIN A SPECIES SYMBOL MUST CONTAIN A
C             BLANK.
C
C    LOUT   = (INPUT) LOGICAL unit NUMBER TO WHICH THE RESULTING FORTRAN
C             SOURCE PROGRAMS WILL BE WRITTEN.
C
C    IWORK  = (INPUT)
C    LENIW  = (INPUT)  replaced by next line
!    liW    = (INPUT)
C             IWORK IS AN integer ARRAY OF LENGTH LENIW USED FOR
C             TEMPORARY STORAGE.  THE USER MUST ALLOCATE SPACE FOR THIS
C             ARRAY IN THE CALLING PROGRAM BY A STATEMENT SUCH AS:
C                   integer IWORK(...)
C             AND SET LENIW TO THE VALUE USED IN THAT STATEMENT.  THE
C             STORAGE REQUIREMENTS CANNOT BE PREDICTED PRECISELY, BUT A
C             VALUE OF LENIW=20*NREAC IS A REASONABLE ESTIMATE.
C
C    CWORK  = (INPUT)
C    LENCW  = (INPUT)
!    LcWW   = (Input)
C             CWORK IS A character ARRAY OF LENGTH LENCW USED FOR
C             TEMPORARY STORAGE.  THE USER MUST ALLOCATE SPACE FOR THIS
C             ARRAY IN THE CALLING PROGRAM BY A STATEMENT SUCH AS:
C                   character CWORK(...)*1
C             AND SET LENCW TO THE VALUE USED IN THAT STATEMENT.  THE
C             STORAGE REQUIREMENTS CANNOT BE PREDICTED PRECISELY, BUT A
C             VALUE OF LENCW=15*NREAC IS A REASONABLE ESTIMATE.
C
C  III.  USAGE  ................................................
C
C    A. THE USER DEFINES A PROBLEM BY SELECTING THE SET OF SPECIES THAT
C       PARTICIPATE IN THE REACTIONS, AND THE SET OF REACTIONS.  THE
C       SYMBOLS FOR THE SET OF SPECIES ARE ENTERED IN THE ARRAY SPECI
C       AND THE SYMBOLS FOR THE REACTANT AND PRODUCT SPECIES ARE ENTERED
C       IN THE ARRAY REAC, WITH EACH REACTION STORED IN A SEPARATE
C       COLUMN.  AS AN EXAMPLE, CONSIDER THE FOUR SPECIES SYMBOLS O, O2,
C       O3, AND M, AND THE SIX REACTIONS:
C                    O   + O3        = O2  + O2
C                    O2  + O2        = O   + O3
C                    O   + O2  + M   = O3  + M
C                    O3  + M         = O   + O2  + M
C                    O2              = O   + O
C                    O3              = O   + O2
C       THE SETUP FOR RUNNING EQNWRT IS AS FOLLOWS:
C             PROGRAM EXAMPL
C             integer IWORK(90)
C             character SPECI(4)*2, REAC(6,6)*2, CWORK(120)*1
C             data (SPECI(i), i=1,4) /'O ', 'O2', 'O3', 'M'/
C             data ((REAC(i,j), i=1,6), j=1,6)
C            8    /'O ', 'O3', '  ',  'O2', 'O2', '  ',
C            8     'O2', 'O2', '  ',  'O ', 'O3', '  ',
C            8     'O ', 'O2', 'M ',  'O3', 'M ', '  ',
C            8     'O3', 'M ', '  ',  'O ', 'O2', 'M ',
C            8     'O2', '  ', '  ',  'O ', 'O ', '  ',
C            8     'O3', '  ', '  ',  'O ', 'O2', '  '/
C             NDE = 3
C             NSPECI = 4
C             NREAC = 6
C             NPROD = 3
C             LOUT = 6
C             open(file = 'TAPE6', unit=LOUT, status = 'new')
C             LENCW = 120
C             LENIW = 90
C             CALL EQNWRT (NDE, NSPECI, NREAC, NPROD, SPECI, REAC,
C            8             LOUT, IWORK, LENIW, CWORK, LENCW)
C             end
C       IN ACTUAL PRACTICE, LOADING OF THE ARRAYS SPECI AND REAC WOULD
C       BE DONE BY READING FROM A FILE INSTEAD OF BY data STATEMENTS.
C
C    B. FOR THE ABOVE EXAMPLE THE subroutineS F AND JACOBN, PRODUCED AS
C       OUTPUT ON THE FILE LOUT, ARE:
C
C     subroutine F (NDE, T, Y, YP)
C     REAL K, T, Y(NDE), YP(NDE), Q(3)
C     common /SDRCOM/ K(6)
C     data (Q(i), i=1,3)/3*0./
C
C     RIGHT HAND SIDE FOR SPECIES O
C     YP(1)=Y(2)*(+K(2)*Y(2)-K(3)*Y(1)*Y(4)+K(5)+K(5))+Y(3)*(-K(1)*Y(1)+
C    8 K(4)*Y(4)+K(6))+Q(1)
C     RIGHT HAND SIDE FOR SPECIES O2
C     YP(2)=Y(2)*(Y(2)*(-K(2)-K(2))-K(3)*Y(1)*Y(4)-K(5))+Y(3)*(Y(1)*(+K(
C    8 1)+K(1))+K(4)*Y(4)+K(6))+Q(2)
C     RIGHT HAND SIDE FOR SPECIES O3
C     YP(3)=Y(3)*(-K(1)*Y(1)-K(4)*Y(4)-K(6))+Y(2)*(+K(2)*Y(2)+K(3)*Y(1)*
C    8 Y(4))+Q(3)
C     return
C     end
C     subroutine JACOBN (NDE, T, Y, P, MATDIM, ML, MU)
C     REAL K, P(MATDIM,NDE), T, Y(NDE)
C     common /SDRCOM/ K(6)
C
C     DO 10 j = 1,3
C       DO 10 i = 1,3
C10       P(i,j) = 0.
C   P(i,j) CONTAINS THE PARTIAL DERIVATIVE OF YP(i) WITH RESPECT TO Y(j)
C
C     P(1,1)=-K(1)*Y(3)-K(3)*Y(2)*Y(4)
C     P(1,2)=Y(2)*(+K(2)+K(2))-K(3)*Y(1)*Y(4)+K(5)+K(5)
C     P(1,3)=-K(1)*Y(1)+K(4)*Y(4)+K(6)
C
C     P(2,1)=Y(3)*(+K(1)+K(1))-K(3)*Y(2)*Y(4)
C     P(2,2)=Y(2)*(-K(2)-K(2)-K(2)-K(2))-K(3)*Y(1)*Y(4)-K(5)
C     P(2,3)=Y(1)*(+K(1)+K(1))+K(4)*Y(4)+K(6)
C
C     P(3,1)=-K(1)*Y(3)+K(3)*Y(2)*Y(4)
C     P(3,2)=Y(2)*(+K(2)+K(2))+K(3)*Y(1)*Y(4)
C     P(3,3)=-K(1)*Y(1)-K(4)*Y(4)-K(6)
C     return
C     end
C
C  IV.  HOW TO USE THE OUTPUT OF THIS ROUTINE WITH SDRIV1,2,3  .........
C
C    A. RATE CONSTANTS
C       THE USER MUST PROVIDE VALUES FOR THE REACTION RATE CONSTANTS
C       K(1), ... , K(NREAC).  THIS ARRAY APPEARS IN THE common BLOCK
C       /SDRCOM/ IN subroutineS F AND JACOBN.  IF THESE QUANTITIES
C       CHANGE WITH TIME THE USER SHOULD CALCULATE THEIR VALUES AT THE
C       BEGINNING OF subroutine F.  OTHERWISE THESE CONSTANTS CAN BE SET
C       IN THE USER'S PROGRAM WHICH CALLS SDRIV1,2 OR 3.
C
C    B. EXTERNAL SOURCES
C       PROVISION HAS BEEN MADE FOR EXTERNAL SOURCE TERMS Q(1), ... ,
C       Q(NDE), IN THE EQUATIONS IN subroutine F.  THESE ARE
C       AUTOMATICALLY PRESET TO ZERO IN F.  FOR PROBLEMS HAVING SOURCE
C       TERMS THE USER SHOULD INCLUDE FORTRAN STATEMENTS AT THE
C       BEGINNING OF subroutine F TO CALCULATE THESE VALUES.
C
C    C. subroutine SDRIV1 SHOULD BE USED IF THE NUMBER OF DIFFERENTIAL
C       EQUATIONS, NDE, IS LESS THAN 200.  THE subroutine JACOBN IS NOT
C       CALLED BY SDRIV1 (AND MAY BE OMITTED) BECAUSE THE JACOBIAN
C       MATRIX IS APPROXIMATED BY NUMERICAL DIFFERENCING.
C       subroutine SDRIV2 SHOULD BE USED IF THE NUMBER OF DIFFERENTIAL
C       EQUATIONS, NDE, IS GREATER THAN 200.  THE subroutine JACOBN IS
C       NOT CALLED BY SDRIV2 (AND MAY BE OMITTED) BECAUSE THE JACOBIAN
C       MATRIX IS APPROXIMATED BY NUMERICAL DIFFERENCING.
C       subroutine SDRIV3 SHOULD BE USED IF THE ANALYTICAL JACOBIAN
C       subroutine, JACOBN, IS TO REPLACE THE NUMERICAL DIFFERENCING
C       APPROXIMATION.
C       WHEN CALLING SDRIV1,2,3 THE VALUE OF NDE SHOULD BE USED FOR THE
C       PARAMETER N.  WHEN USING SDRIV2,3 A VALUE OF 2 IS SUGGESTED
C       THE PARAMETER MINT.
C
C    D. SAMPLE USAGE OF subroutine F WITH SDRIV1:
C             PROGRAM CHEM
C             REAL Y(...), WORK(...)
C             common /SDRCOM/ RK(...)
C
C             open(file = 'OUTPUT', unit=1, status = 'new')
C             NREAC = ...
C             NDE = ...
C             NSPECI = ...
C             T = ...
C             DO 10 i = 1,NSPECI
C        10     Y(i) = ...                        SET INITIAL CONDITIONS
C             TOUT = T
C             MSTATE = 1
C             EPS = ...
C             LENW = ...
C             DO 20 i = 1,NREAC
C        20     RK(i) = ...                SET VALUES FOR RATE CONSTANTS
C        30   CALL SDRIV1 (NDE, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C             if(MSTATE.GT.2) STOP
C             write(1, '(...)') TOUT, (Y(i), i=1,NSPECI)
C             TOUT = ...
C             if(TOUT .LE. ...) GO TO 30
C             end(SAMPLE)
C
C  V.  REMARKS  ................................................
C
C    A. THE ERROR HANDLER ROUTINE XERROR WILL BE CALLED TO TRANSMIT A
C       DIAGNOSTIC MESSAGE TO THE USER AND TERMINATE THE RUN IF A
C       PROBABLE PROBLEM SETUP ERROR IS DETECTED, E.G. INSUFFICIENT
C       SPACE IN THE IWORK ARRAY, OR A SYMBOL IN THE REAC ARRAY IS NOT
C       IN THE SPECI ARRAY.  THE MESSAGE WILL BE WRITTEN ON THE
C       STANDARD ERROR MESSAGE FILE, WHICH SHOULD BE DETERMINED BY THE
C       USER FROM THE DOCUMENTATION FOR XERROR.
C       FOLLOWING IS A LIST OF POSSIBLE FATAL ERRORS:
C
C         NO.  MESSAGE
C         ---  -------
C          1   FROM EQFCHM:  INSUFFICIENT WORK SPACE.
C          2   FROM EQFCHM:  EXPRESSION TOO LONG TO PROCESS.
C          3   FROM EQJCBN:  INSUFFICIENT WORK SPACE.
C          4   FROM EQJCBN:  EXPRESSION TOO LONG TO PROCESS.
C          5   FROM EQDY:    SYMBOL MISSING FROM THE SPECI ARRAY.
C          6   FROM EQDY:    INSUFFICIENT WORK SPACE.
C          7   FROM EQDIFF:  INSUFFICIENT WORK SPACE.
C          8   FROM EQRFCT:  INSUFFICIENT WORK SPACE.
C          9   FROM EQFCTR:  INSUFFICIENT WORK SPACE.
C         10   FROM EQOUT:   INSUFFICIENT WORK SPACE.
C         11   FROM EQSPLT:  INSUFFICIENT WORK SPACE.
C
C    B. THE PARAMETER LINES IS THE MAXIMUM NUMBER OF LINES ALLOWED IN
C       ONE REPLACEMENT STATEMENT, AND IS SET TO A VALUE OF 20 IN THE
C       data STATEMENT BELOW.  THE OUTPUT ROUTINES F AND JACOBN MAY
C       CONTAIN EXPRESSIONS THAT ARE TOO COMPLEX TO BE PROCESSED BY SOME
C       COMPILERS.  IN THAT CASE THE VALUE OF LINES SHOULD BE REDUCED
C       SLIGHTLY, E.G. TO 15.
C
C    C. THE PARAMETER NCHAR IS THE NUMBER OF characterS IN THE FIRST
C       LINE OF A REPLACEMENT STATEMENT, AND NCHAR-1 IS THE NUMBER OF
C       characterS IN SUCCEEDING LINES.
C
C    D. OTHER ROUTINES USED: EQFCHM, EQJCBN, EQDY, EQDIFF, EQRFCT,
C       EQFCTR, EQOUT, EQSPLT, EQPRNT.
C
C***End PROLOGUE  EQNWRT
!
!      integer IWORK(LENIW)
      integer iWork(LiW)
!      character SPECI(NSPECI)*(*), REAC(1,NREAC)*(*), CWORK(LENCW)*1
      character SPECI(NSPECI)*(*), REAC(1,NREAC)*(*), CWORK(LCW)*1
      data LINES /20/, NCHAR /66/
C
      NSym = 3 + NProd
c
c      write(unit = *, fmt = *) " LenIW = ", LENIW, "NProd = ", NProd, 
c     1  "Sym =", Sym
c
      call EQFCHM (NREAC, NDE, NSPECI, NSym, REAC, SPECI, LINES,
!     8             NCHAR, LOUT, IWORK, LENIW, CWORK, LENCW)
     8  NCHAR, LOUT, IWORK, LiW, CWORK, LCW)
      call EQJCBN (NREAC, NDE, NSPECI, NSym, REAC, SPECI, LINES,
!     8             NCHAR, LOUT, IWORK, LENIW, CWORK, LENCW)
     8             NCHAR, LOUT, IWORK, LiW, CWORK, LCW)
      return
      end
!
      subroutine EqFChm (NREAC, NDE, NSPECI, NSym, REAC, SPECI, LINES,
!     8                 NCHAR, LOUT, IWORK, LENIW, CWORK, LENCW)
     8  NCHAR, LOUT, IWORK, LiW, CWORK, LCW)
!
!     This subroutine writes subroutine FChem.
!
      character BCD*92
!      integer iWork(LeniW)
      integer iWork(LiW)
      character REAC(NSym,NREAC)*(*), SPECI(NSPECI)*(*)
!      character CWORK(LENCW)*1
      character CWORK(LCW)*1
!
      IDELIM = MAX(NREAC, NSPECI) + 1
!      MAXLEN = LENCW
      MAXLEN = LCW
      IP = 1
      CALL EQSPLT (CWORK, IP, NDE, MAXLEN)
!      if(LENCW.LT.IP) GO TO 100
      if(LCW.LT.IP) GO TO 100
      CWORK(IP) = ')'
!      write(LOUT, '(6X, ''subroutine F (NDE, T, Y, YP)'' /
!     8       6X, ''REAL K, T, Y(NDE), YP(NDE), Q('', 36A1)')
!     8  (CWORK(i), i=1,IP)
      write(unit = LOut, fmt = "(6x, a29)") 
     1  "subroutine FChem(n, t, y, yp)"
c xxx
c      write(unit = *, fmt = "(6x, a29)") 
c     1  "subroutine FChem(n, t, y, yp)"
c      pause
c xxx
! *****
!      write(unit = *, fmt = *) " In EQFCHM 'xxx':  NReac =", NReac
! *****
      write(unit = LOut, fmt = "(6x,a38)") 
     1  "parameter (NyDim = 600, MxReac = 8000)"
      write(unit = LOut, fmt = "(6x, a21)") "dimension y(*), yp(*)"
      write(unit = LOut, fmt = "(6x, a59/ 5x, a53)") 
     1  "common/lc2/ A(MxReac), B(MxReac), C(MxReac), Numb(MxReac), ",
     2  "1  Cat(MxReac), ex(MxReac), TLo(MxReac), THi(MxReac) " 
!     2  "1  Cat(MxReac), ex(MxReac), bf(MxReac), cf(MxReac)"
!     2  "1  iCat(MxReac), ex(MxReac), bf(MxReac), cf(MxReac)"
c     1  "common/lc2/ A(MaxEqs), B(MaxEqs), C(MaxEqs), Numb(MaxEqs), ",
c     2  "1  icat(MaxEqs), ex(MaxEqs), bf(MaxEqs), cf(MaxEqs)"
      write(unit = LOut, fmt = "(6x, a31)") 
     1  "common/char2/ Symbol(7, MxReac)"
c     1  "common/char2/ Symbol(7, MaxEqs)"
      write(unit = LOut, fmt = "(6x, a24)") "common/stcom9/ k(MxReac)" 
c      write(unit = LOut, fmt = "(6x, a24)") "common/stcom9/ k(MaxEqs)" 
      write(unit = LOut, fmt = "(6x, a27)") 
     1  "common/stcom10/ Rat(MxReac)"
c     1  "common/stcom10/ Rat(MaxEqs)"
      write(unit = LOut, fmt = "(6x, a24)") "common/stcom11/ q(NyDim)" 
      write(unit = LOut, fmt = "(6x, a6/ a1)") "real k", "!"
!      write(unit = LOut, fmt = "(6x, a25/ a1)") 
!     1  "double precision k, y, yp", "!"
!      IP = 1
!      CALL EQSPLT (CWORK, IP, NREAC, MAXLEN)
!      if(LENCW.LT.IP) GO TO 100
!      CWORK(IP) = ')'
!      write(LOUT, '(6X, ''common /SDRCOM/ K('', 48A1)')
!     8  (CWORK(i), i=1,IP)
!      IP = 1
!      CALL EQSPLT (CWORK, IP, NDE, MAXLEN)
!      if(LENCW.LT.(IP+3)) GO TO 100
!      CWORK(IP) = ')'
!      CWORK(IP+1) = '/'
!      IP = IP + 2
!      MAXLEN = LENCW - IP
!      CALL EQSPLT (CWORK, IP, NDE, MAXLEN)
!      if(LENCW.LT.(IP+3)) GO TO 100
!      CWORK(IP) = '*'
!      CWORK(IP+1) = '0'
!      CWORK(IP+2) = '.'
!      CWORK(IP+3) = '/'
!      IP = IP + 3
!      write(LOUT, '(6X, ''data (Q(i), i=1,'', 50A1)') (CWORK(i), i=1,IP)
!      write(LOUT, '(''!'')')
      MXLENI = 0
      MXLENC = 0
      MAXWRT = (NCHAR - 1)*LINES
      do 90 j = 1, NDE
      write(LOut, '(''!     Density rate for species '', A)')SPECI(j)
!        LENGTH = LENIW
        Length = LiW
        IWORK(1) = -IDELIM
c
c      write(unit = *, fmt = *) " iWork(2)=", iWork(2), "NReac=", NReac,
c     1  "NSpeci=", NSpeci, "NSym=", NSym, "Length=", Length, "j=", j
c
        call EQDY(Reac, Speci, iWork(2), NReac, NSpeci, NSym,
     8             Length, j)
        IARROW = 2
        if(LENGTH.NE.0) THEN
          LENGTH = LENGTH + 2
          IWORK(Length) = IDELIM
!          LENFIL = LENIW - 1
          LENFIL = LiW - 1
          CALL EQRFCT (IWORK, LENFIL, LENGTH, NSPECI, IDELIM)
!          if(LENCW.LT.4) GO TO 100
          if(LCW.LT.4) GO TO 100
          CWORK(1) = '+'
          CWORK(2) = 'Y'
          CWORK(3) = 'P'
          CWORK(4) = '('
          NTAG = 5
!          MAXLEN = LENCW - 5
          MAXLEN = LCW - 5
          CALL EQSPLT (CWORK, NTAG, j, MAXLEN)
!          if(LENCW.LT.(NTAG+1)) GO TO 100
          if(LCW.LT.(NTAG+1)) GO TO 100
          CWORK(NTAG) = ')'
          CWORK(NTAG+1) = '='
!          IARROW = LENCW - (NTAG + 2) + 1
          IARROW = LCW - (NTAG + 2) + 1
          CALL EQOUT (IWORK, CWORK(NTAG+2), LENGTH, IARROW, IDELIM)
          MXLENI = MAX(MXLENI, LENGTH)
        end if
        IP = NTAG + 1 + IARROW
!        if(LENCW.LT.(IP+2)) GO TO 100
        if(LCW.LT.(IP+2)) GO TO 100
        CWORK(IP) = '+'
        CWORK(IP+1) = 'Q'
        CWORK(IP+2) = '('
        IP = IP + 3
!        MAXLEN = LENCW - IP
        MAXLEN = LCW - IP
        CALL EQSPLT (CWORK, IP, j, MAXLEN)
!        if(LENCW.LT.IP) GO TO 100
        if(LCW.LT.IP) GO TO 100
        CWORK(IP) = ')'
        MXLENC = MAX(MXLENC, IP)
        ILAST = IP
        ISTART = NTAG + 3
        ISTOP = ISTART + MAXWRT - NTAG
   20   if(ISTOP.LT.ILAST) THEN
C                                              FIND A SEQUENCE THAT WILL
C                                              FIT IN LINES OF OUTPUT.
          IC = 0
          Iend = ISTOP
          ISTOP = 0
          do 30 i = ISTART,Iend
            if(CWORK(i).EQ.')') IC = IC + 1
            if(CWORK(i).EQ.'(') IC = IC - 1
            if((IC.EQ.0) .AND. (CWORK(i).EQ.')') .AND.
     8      (CWORK(i+1).NE.'*')) ISTOP = i
   30       continue
          if(ISTOP.EQ.0) THEN
C                                 WITHIN THE LOSS TERMS, i.E. THOSE
C                                 WITH Y(j) AS A FACTOR, FIND A SEQUENCE
C                                 THAT WILL FIT IN LINES OF OUTPUT.
            Iend = Iend - 1
            do 40 i = ISTART,Iend
              if(CWORK(i).EQ.'*') GO TO 50
   40         continue
            GO TO 110
   50       IPAREN = i + 1
            IC = 1
            do 60 i = IPAREN,Iend
              if(CWORK(i).EQ.')') IC = IC + 1
              if(CWORK(i).EQ.'(') IC = IC - 1
              if((IC.EQ.0) .AND. (CWORK(i).EQ.')') .AND.
     8        (CWORK(i+1).NE.'*')) ISTOP = i
   60         continue
            if(ISTOP.EQ.0) GO TO 110
            if(ISTART.EQ.NTAG+3) THEN
              IPRNT = 3
            else
              IPRNT = 4
            end if
            CALL EQPRNT(CWORK, NTAG, ISTART, ISTOP, IPRNT, NCHAR, LOUT)
            IMOVE = ISTOP - IPAREN
            do 70 i = IPAREN,ISTART,-1
              ISUB = i + IMOVE
   70         CWORK(ISUB) = CWORK(i)
            ISTART = ISTART + IMOVE
            ISTOP = ISTART + MAXWRT - 2*NTAG
            GO TO 20
          end if
        end if
        ISTOP = MIN(ISTOP, ILAST)
        if(ISTART.EQ.NTAG+3) THEN
          IPRNT = 1
        else
          IPRNT = 2
        end if
        CALL EQPRNT(CWORK, NTAG, ISTART, ISTOP, IPRNT, NCHAR, LOUT)
        if(ISTOP.NE.ILAST) THEN
          ISTART = ISTOP + 1
          ISTOP = ISTART + MAXWRT - 2*NTAG
          GO TO 20
        end if
   90 continue
      write(unit = LOut, fmt = "(6x, a6/ 6x, a3/ a1)") 
     1  "return", "end", "!"
      return
!
  100 write(BCD, '(''EQFCHM1FE IN subroutine EQFCHM, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE CWORK ARRAY.'')')
      CALL XERROR(BCD, 91, 1, 2)
      return
!
  110 write(BCD, '(''EqFChm2FE in subroutine EQFChm, the expression'',
     8  '' for species '', I3, '' was too long to be processed.'')') j
      CALL XERROR(BCD, 92, 2, 2)
      return
      end
      
      subroutine EqJcbn(NREAC, NDE, NSPECI, NSym, REAC, SPECI, LINES,
     1                  NCHAR, LOUT, IWORK, LiW, CWORK, LCW)
!     1                  NCHAR, LOUT, IWORK, LENIW, CWORK, LENCW)
!
!     This subroutine writes subroutine Jacobn.
!
      character BCD*92
!      integer IWORK(LENIW)
      integer IWORK(LiW)
      character REAC(NSym,NREAC)*(*), SPECI(NSPECI)*(*)
!      character CWORK(LENCW)*1
      character CWORK(LCW)*1
C
      IDELIM = MAX(NREAC, NSPECI) + 1
!      MAXLEN = LENCW
      MAXLEN = LCW
      MXLENI = 0
      MXLENC = 0
      MAXWRT = (NCHAR - 1)*LINES
      IP = 1
      CALL EQSPLT (CWORK, IP, NREAC, MAXLEN)
!      if(LENCW.LT.IP) GO TO 80
      if(LCW.LT.IP) GO TO 80
      CWORK(IP) = ')'
!
      write(unit = LOut, fmt = "(6x, a47)") 
     1  "subroutine Jacobn(NDE, t, y, p, MatDim, ML, MU)"
!      write(unit = LOut, fmt = "(6x, a52)")
      write(unit = LOut, fmt = "(6x, a38)")
     1  "parameter (NyDim = 200, MxReac = 1800)"      !, NReac = 1045)" 
c     1  "parameter (NyDim = 200, MaxEqs = 1000, NReac = 820)" 
! *****
      write(unit = *, fmt = *) " In EqJcbn near start:  NReac = ", NReac
! *****
      write(unit = LOut, fmt = "(6x, a33)")
     1  "real k, P(MatDim, NDE), t, y(NDE)"
!      write(unit = LOut, fmt = "(6x, a45)")
!     1  "double precision k, P(MatDim, NDE), t, y(NDE)"
      write(unit = LOut, fmt = "(6x, 18a, 48a1)")
     1  "common /sdrcom/ k(", (CWORK(i), i = 1, iP)    !This has been used, but stcom9 used also.
      write(unit = LOut, fmt = "(a1)") "c"
      write(unit = LOut, fmt = "(a1, 5x, a25/ a1)")
     1  "c", "common /stcom9/ k(MxReac)", "!"          !Comment out second (wrong?) K
!     8  6x, ''common /STCOM9/ k('', 48A1)') (CWORK(i), i=1,IP)
      IP = 1
      CALL EQSPLT(CWORK, IP, NDE, MAXLEN)
      IP = IP - 1
      write(LOut, '(6x, ''do 10 j = 1, '', 54A1)') (CWORK(i), i=1, IP)
      write(LOut, '(8x, ''do 10 i = 1, '', 52A1)') (CWORK(i), i=1, IP)
      write(LOut, '('' 10       p(i, j) = 0.'' /  
     1  ''! p(i, j) contains'',
     2  '' the partial derivative of yp(i) with respect to y(j)'')')
      do 70 j = 1, NDE
        write(unit = LOut, fmt = "(a1)") "!"
!        LONG = LENIW
        LONG = LiW
        CALL EQDY (REAC, SPECI, IWORK, NREAC, NSPECI, NSym, LONG, j)
        if(LONG.EQ.0) GO TO 70
        do 60 L = 1,NDE
!          LIPART = LENIW - (LONG + 2) + 1
          LIPART = LiW - (LONG + 2) + 1
          CALL EQDIFF (IWORK, IWORK(LONG+2), LONG, LIPART, L)
          if(LIPART.EQ.0) GO TO 60
          LENGTH = LIPART + 2
          IWORK(LONG+1) = -IDELIM
          LSUB = LONG + LENGTH
          IWORK(LSUB) = IDELIM
!          LENFIL = LENIW - (LONG + 1) + 1
          LENFIL = LiW - (LONG + 1) + 1
          CALL EQRFCT(IWORK(LONG+1), LENFIL, LENGTH, NSPECI,
     8                IDELIM)
          MXLENI = MAX(MXLENI, LONG+LENGTH)
!          if(LENCW.LT.3) GO TO 80
          if(LCW.LT.3) GO TO 80
          CWORK(1) = '+'
          CWORK(2) = 'P'
          CWORK(3) = '('
          IP = 4
!          MAXLEN = LENCW - IP
          MAXLEN = LCW - IP
          CALL EQSPLT (CWORK, IP, j, MAXLEN)
!          if(LENCW.LT.IP) GO TO 80
          if(LCW.LT.IP) GO TO 80
          CWORK(IP) = ','
          IP = IP + 1
!          MAXLEN = LENCW - IP
          MAXLEN = LCW - IP
          CALL EQSPLT (CWORK, IP, L, MAXLEN)
!          if(LENCW.LT.(IP+1)) GO TO 80
          if(LCW.LT.(IP+1)) GO TO 80
          CWORK(IP) = ')'
          CWORK(IP+1) = '='
          NTAG = IP
!          IARROW = LENIW - (NTAG + 2) + 1
          IARROW = LiW - (NTAG + 2) + 1
          CALL EQOUT (IWORK(LONG+1), CWORK(NTAG+2), LENGTH, IARROW,
     8                IDELIM)
          MXLENC = MAX(MXLENC, NTAG + 1 + IARROW)
          IARROW = IARROW - 1
          ILAST = NTAG + 1 + IARROW
          ISTART = NTAG+3
          ISTOP = ISTART + MAXWRT - NTAG
   20     if(ISTOP.LT.ILAST) THEN
!                  Find a sequence that will fit in lines of output.
            IC = 0
            Iend = ISTOP
            ISTOP = 0
            do 30 i = ISTART,Iend
              if(CWORK(i).EQ.')') IC = IC + 1
              if(CWORK(i).EQ.'(') IC = IC - 1
              if((IC.EQ.0) .AND. (CWORK(i).EQ.')') .AND.
     8        (CWORK(i+1).NE.'*')) ISTOP = i
   30         continue
            if(iStop .eq. 0) then
              write(BCD, '(''EQJCBN4FE IN subroutine EQJCBN, THE '',
     8        ''EXPRESSION FOR P('', I3, '','', I3, '') WAS TOO LONG '',
     8        ''TO BE PROCESSED.'')') j, L
              CALL XERROR(BCD, 91, 4, 2)
              return
            end if
          end if
          ISTOP = MIN(ISTOP, ILAST)
          if(ISTART.EQ.NTAG+3) THEN
            IPRNT = 1
          else
            IPRNT = 2
          end if
          CALL EQPRNT(CWORK, NTAG, ISTART, ISTOP, IPRNT, NCHAR, LOUT)
          if(ISTOP.NE.ILAST) THEN
            ISTART = ISTOP + 1
            ISTOP = ISTART + MAXWRT - 2*NTAG
            GO TO 20
          end if
   60     continue
   70   continue
      write(LOut, '(6x, ''return'' / 6x, ''end''/ ''!'')')
      return
!
   80 write(BCD, '(''EQJCBN3FE IN subroutine EQJCBN, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE CWORK ARRAY.'')')
      CALL XERROR(BCD, 91, 3, 2)
      return
      end
!
      subroutine EQDY (REAC, SPECI, IFORM, NREAC, NSPECI, NSym,
     8                 ISPEAR, j)
C
C     IFORM(1) TO IFORM(ISPEAR) CONTAIN CODED FORMATION AND
C     REMOVAL TERMS.
C
      character BCD*93
      character REAC(NSym,NREAC)*(*), SPECI(NSPECI)*(*)
      integer IFORM(*), IRCTNT(3)
C                                           SEARCH TO FIND SPECIES j IN
C                                           POSITION IPOS OF REACTION K.
! *****
!      write(unit = *, fmt = *) NReac, NSpeci  !, iForm, NReac, NSpeci, NSym, 
 !    1  iSpear, j
!      pause
! *****
      LENLFT = ISPEAR
      ISPEAR = 0
      do 140 K=1, NREAC
        do 140 IPOS = 1, NSym
          if(REAC(IPOS, K) .NE. SPECI(j)) GO TO 140
C
C                                K-TH REACTION CONTAINS Y(j) IN POSITION
C                                IPOS.  ENCODE SYMBOLS INTO IRCTNT.
          ISET = 1
          do 30 m = 1, 3
            if(Reac(m, k) .eq. ' ') then
              IRCTNT(M) = 0
            else
              ISET = ISET + 1
              do 10 L = 1, NSPECI
                if(REAC(M, K) .EQ. SPECI(L)) then
                  GO TO 20
                else
                end if
c
! ^^^^^^^^
!      write(unit = *, fmt = *) " iSet =", iSet, ", m =", m, ", k =", k, 
!     1  ", Reac(m, k) = ", Reac(m, k), ", L =", L, ", Speci(L) = ", 
!     1  Speci(L)
  !    pause
! ********
c
   10           continue
              write(BCD, '(''.5FE REACTION '',  I5, '' contains a '',
     8        ''SYMBOL IN POSITION '', I1, '' THAT IS NOT IN THE INPU'',
     8        ''T SYMBOL TABLE.'')') K, M
      write(unit =*, fmt = "(1x, a8)") Reac(m, k)
!      stop
              CALL XERROR(BCD, 93, 5, 2)
              return
   20         IRCTNT(M) = L
            end if
   30       continue
          if(IPOS.LE.3) THEN
C                                         REMOVAL TERM.  CHECK FOR
C                                         APPEARANCE OF SPECIES j ON THE
C                                         OPPOSITE SIDE OF REACTION K
            do 40 L = 4, NSym
              if(REAC(L, K) .EQ. SPECI(j)) GO TO 50
   40         continue
            GO TO 70
C                                       TEST FOR THE FIRST APPEARANCE OF
C                                       SPECIES j ON THE SAME SIDE OF
C                                       REACTION K.
   50       if(IPOS .EQ. 1) GO TO 140
            do 60 L = 2, IPOS
              if(REAC(L - 1, K) .EQ. SPECI(j)) GO TO 70
   60         continue
            GO TO 140
C                                      K-TH REACTION IS REMOVAL REACTION
C                                      FOR Y(j).  ENCODE REMOVAL TERM
C                                      IN IFORM, STARTING AT ISPEAR.
C
   70       if(LENLFT.LT.(ISPEAR+2)) GO TO 160
            IFORM(ISPEAR + 1) = -ISET
            IFORM(ISPEAR + 2) = K
            ISPEAR = ISPEAR + 2
            do 80 II = 1, 3
              if(REAC(II, K) .NE. ' ') THEN
                ISPEAR = ISPEAR + 1
                if(LENLFT .LT. ISPEAR) GO TO 160
                IFORM(ISPEAR) = IRCTNT(II)
              end if
   80         continue
          else
C                                                         FORMATION TERM
            do 90 L = 1, 3
              if(REAC(L, K) .EQ. SPECI(j)) GO TO 100
   90         continue
            GO TO 120
  100       if(IPOS .EQ. 4) GO TO 140
            do 110 L = 5, IPOS
              if(REAC(L - 1, K) .EQ. SPECI(j)) GO TO 120
  110         continue
            GO TO 140
C                                    K-TH REACTION IS FORMATION REACTION
C                                    FOR Y(j).  ENCODE FORMATION TERM
C                                    IN IFORM, STARTING AT ISPEAR.
C
  120       if(LENLFT .LT. (ISPEAR + 2)) GO TO 160
            IFORM(ISPEAR + 1) = ISET
            IFORM(ISPEAR+2) = K
            ISPEAR = ISPEAR + 2
            do 130 II=1,3
              if(REAC(II,K).NE.' ') THEN
                ISPEAR=ISPEAR+1
                if(LENLFT.LT.ISPEAR) GO TO 160
                IFORM(ISPEAR)=IRCTNT(II)
              end if
  130         continue
          end if
  140   continue
      return
C
  160 write(BCD, '(''EQDY6FE IN subroutine EQDY, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE IWORK ARRAY.'')')
      CALL XERROR(BCD, 87, 6, 2)
      return
      end
C
      subroutine EQDIFF (INFILE, IPART, INSTOP, LENGTH, LDEN)
C
C     subroutine TAKES PORTION OF CORE FROM INFILE(1) TO AND INCLUDING
C     INFILE(INSTOP) AND CALCULATES PARTIAL DERIVATIVE OF IT WITH
C     RESPECT TO Y(LDEN). DERIVATIVE IS STORED IN IPART(1) TO
C     IPART(LENGTH).  IF DERIVATIVE IS ZERO, LENGTH IS RETURNED AS ZERO.
C
      character BCD*92
      integer INFILE(INSTOP), IPART(*)
C
      LENLFT = LENGTH
      IARROW=0
      INDEX=1
C                                         CHECK FOR end OF INPUT PORTION
   10 ISTOP=ABS(INFILE(INDEX))
      if(ISTOP.GE.2) THEN
C                                                   SCAN SERIES FOR LDEN
        ICOUNT=0
        do 20 i=INDEX+2,INDEX+ISTOP
          if(INFILE(i).EQ.LDEN) ICOUNT = ICOUNT + 1
   20     continue
        if(ICOUNT.NE.0) THEN
C                                  ONE LDEN IN SERIES, CALCULATE PARTIAL
C                                  IN IPART, STARTING AT IARROW
C
          do 60 MARK = 1,ICOUNT
            if(LENLFT.LT.(IARROW+2)) GO TO 70
            IPART(IARROW+1) = SIGN((ISTOP - 1), INFILE(INDEX))
            IPART(IARROW+2) = INFILE(INDEX+1)
            IARROW = IARROW + 2
            do 40 i=INDEX+2,INDEX+ISTOP
              if(INFILE(i).NE.LDEN) THEN
                IARROW=IARROW+1
                if(LENLFT.LT.IARROW) GO TO 70
                IPART(IARROW)=INFILE(i)
              end if
   40         continue
            if(ICOUNT.GE.2) THEN
C                                          LDEN OCCURS MORE THAN ONCE IN
C                                          SERIES, DIFFERENTIATE AGAIN
C
              if(LENLFT.LT.(IARROW + ICOUNT - 1)) GO TO 70
              do 50 LL=IARROW+1,IARROW+ICOUNT-1
   50           IPART(LL) = LDEN
              IARROW = IARROW + ICOUNT - 1
            end if
   60       continue
        end if
      end if
C                                                            RESET INDEX
      INDEX = INDEX + ISTOP + 1
      if(INDEX.LE.INSTOP) GO TO 10
      LENGTH=IARROW
      return
C
   70 write(BCD, '(''EQDIFF7FE IN subroutine EQDIFF, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE IWORK ARRAY.'')')
      CALL XERROR(BCD, 91, 7, 2)
      return
      end
C
      subroutine EQRFCT (INFILE, LENFIL, LENGTH, NSPECI, IDELIM)
C
C     subroutine COMPLETELY FACTORS INFILE(1) TO INFILE(LENGTH).
C     LENGTH IS ALTERED IN THE ROUTINE.
C
      character BCD*91
      integer INFILE(*)
      data LMIN /7/
C
      MARK = 1
   10 LMARK = MARK
   20 if((LENGTH - LMARK).LT.LMIN) return
C
C                                   SEARCH INPUT AREA TO FIND FIRST
C                                   RIGHT PARENTHESIS TO RIGHT OF LMARK
      ISTART = LMARK + 1
      do 30 INDEX = ISTART,LENGTH
        if(INFILE(INDEX).EQ.IDELIM) GO TO 40
   30   continue
      INDEX = LENGTH
   40 MARK=INDEX
C                               NOW A LEFT HAND PARENTHESIS IS AT LMARK
C                               AND A RIGHT HAND PARENTHESIS IS AT MARK.
C                               CHECK FOR ENOUGH ROOM TO FACTOR.
C
      if((MARK - LMARK).LT.LMIN) GO TO 10
      IQUIT = MARK - ISTART
C                                              FACTOR EVERYTHING BETWEEN
C                                              LMARK AND MARK
      IPFACT = LENGTH + NSPECI + 1
      if(LENFIL.LT.IPFACT) GO TO 110
      ISPRED = LENFIL - IPFACT + 1
      CALL EQFCTR (INFILE(ISTART), INFILE(LENGTH+1), INFILE(IPFACT),
     8             IQUIT, NSPECI, ISPRED, IDELIM)
      if(ISPRED.EQ.0) GO TO 10
C                                      PUT FACTORED VERSION (IWORK) INTO
C                                      INFILE BETWEEN LMARK AND MARK
      LWORK = LENGTH + NSPECI + ISPRED
      if(IQUIT.GT.ISPRED) THEN
C                                            TOO MUCH ROOM BETWEEN LMARK
C                                            AND MARK, COLLAPSE SOME
        IMOVE = IQUIT - ISPRED
        do 60 i = MARK,LWORK
          ISUB = i - IMOVE
   60     INFILE(ISUB) = INFILE(i)
        LENGTH = LENGTH - IMOVE
      else if(ISPRED.GT.IQUIT) THEN
C
C                                   NEED MORE ROOM, MOVE EVERYTHING TO
C                                   RIGHT OF AND INCLUDING MARK TO RIGHT
        IMOVE = ISPRED - IQUIT
        if(IMOVE.GT.(LENFIL - LWORK)) GO TO 110
        do 80 i = LWORK,MARK,-1
          ISUB = i + IMOVE
   80     INFILE(ISUB) = INFILE(i)
        LENGTH = LENGTH + IMOVE
      end if
C                                        DONT HAVE TO MOVE THINGS AROUND
      do 100 i=1,ISPRED
        ISUB1 = LMARK + i
        ISUB2 = LENGTH + NSPECI + i
  100   INFILE(ISUB1) = INFILE(ISUB2)
      LMARK = LMARK + 2
      GO TO 20
C
  110 write(BCD, '(''EQRFCT8FE IN subroutine EQRFCT, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE IWORK ARRAY.'')')
      CALL XERROR(BCD, 91, 8, 2)
      return
      end
C
      subroutine EQFCTR (INFILE, ICHECK, IFACT, LENGTH, NSPECI,
     8                   ISPRED, IDELIM)
C
C     subroutine FACTORS OUT LMAX, THE SPECIES OCCURRING MOST OFTEN,
C     FROM INFILE(1) TO INFILE(LENGTH). FACTORED VERSION IS RETURNED IN
C     IFACT(1) TO IFACT(ISPRED).
C
      character BCD*91
      integer INFILE(LENGTH), ICHECK(NSPECI), IFACT(*)
C     NSPECI=NO. OF SPECIES
      LENLFT = ISPRED
      do 10 i=1,NSPECI
   10   ICHECK(i)=0
C                              COUNT HOW MANY TIMES EACH SPECIES IS USED
      INDEX=1
   20 ISTOP=ABS(INFILE(INDEX))
      if(ISTOP.GE.2) THEN
        do 40 i = INDEX+2,INDEX+ISTOP
          II = INFILE(i)
          if(i.NE.INDEX+ISTOP) THEN
C                                     CHECK FOR REPEATED SPECIES IN TERM
            do 30 j = i+1,INDEX+ISTOP
              if(II.EQ.INFILE(j)) GO TO 40
   30         continue
          end if
C                                                 REPEATED SPECIES IS II
          ICHECK(II) = ICHECK(II) + 1
   40     continue
      end if
C                                           RESET INDEX TO HIT NEXT TERM
      INDEX = INDEX + ISTOP + 1
      if(INDEX.LE.LENGTH) GO TO 20
C                                 NOW THE ICHECK(L) CONTAINS THE NUMBER
C                                 OF TIMES THE SPECIES L OCCURS IN
C                                 LOCATIONS INFILE(1) TO INFILE(LENGTH).
C                                 SCAN THE ICHECK TO FIND LMAX, THE
C                                 SPECIES OCCURRING MOST OFTEN.
      LCOUNT=0
      do 50 i=1,NSPECI
        if(ICHECK(i).GT.LCOUNT) THEN
          LMAX=i
          LCOUNT=ICHECK(i)
        end if
   50   continue
      if(LCOUNT.LE.1) THEN
C                                                     CANNOT BE FACTORED
        ISPRED=0
        return
      end if
C                                         FACTOR AND ARRANGE PARENTHESIS
      if(LENLFT.LT.2) GO TO 150
      IFACT(1)=LMAX
      IFACT(2)=-IDELIM
      ISPEAR=3
      INDEX=1
   60 ISTOP=ABS(INFILE(INDEX))
      if(ISTOP.GE.2) THEN
        do 70 i=INDEX+2,INDEX+ISTOP
          if(INFILE(i).EQ.LMAX) GO TO 80
   70     continue
C                                      NO LMAX IN THIS SERIES OF SPECIES
        GO TO 100
C                                             LMAX IS IN THIS SERIES,
C                                             COLLAPSE SERIES INTO IFACT
   80   if(LENLFT.LT.(ISPEAR+1)) GO TO 150
        IFACT(ISPEAR) = SIGN((ISTOP - 1), INFILE(INDEX))
        ISPEAR=ISPEAR+1
        IFACT(ISPEAR)=INFILE(INDEX+1)
        IDROP=0
        do 90 i=INDEX+2,INDEX+ISTOP
          if((IDROP.EQ.0) .AND. (INFILE(i).EQ.LMAX)) THEN
            IDROP=1
          else
            ISPEAR=ISPEAR+1
            if(LENLFT.LT.ISPEAR) GO TO 150
            IFACT(ISPEAR)=INFILE(i)
          end if
   90     continue
        ISPEAR=ISPEAR+1
      end if
C                                                            RESET INDEX
  100 INDEX = INDEX + ISTOP + 1
      if(INDEX.LT.LENGTH) GO TO 60
C                                                  SET RIGHT PARENTHESIS
      if(LENLFT.LT.ISPEAR) GO TO 150
      IFACT(ISPEAR)=IDELIM
      INDEX=1
  110 ISTOP=ABS(INFILE(INDEX))
      if(ISTOP.GE.2) THEN
        do 120 i=INDEX+2,INDEX+ISTOP
          if(INFILE(i).EQ.LMAX) GO TO 140
  120     continue
      end if
C                                          NO LMAX IN THIS SERIES OF
C                                          SPECIES, SHIP SERIES TO IFACT
      ISPEAR = ISPEAR + 1
      if(LENLFT.LT.(ISPEAR+ISTOP)) GO TO 150
      IFACT(ISPEAR) = INFILE(INDEX)
      do 130 i=1,ISTOP
        ISUB1 = ISPEAR + i
        ISUB2 = INDEX + i
  130   IFACT(ISUB1) = INFILE(ISUB2)
      ISPEAR = ISPEAR + ISTOP
C                                                            RESET INDEX
  140 INDEX = INDEX + ISTOP + 1
      if(INDEX.LT.LENGTH) GO TO 110
      ISPRED = ISPEAR
      return
C
  150 write(BCD, '(''EQFCTR9FE IN subroutine EQFCTR, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE IWORK ARRAY.'')')
      CALL XERROR(BCD, 91, 9, 2)
      return
      end
C
      subroutine EQOUT (INFILE, OUTFIL, LENGTH, JARROW, IDELIM)
C
C     THIS subroutine DECODES THE INFORMATION ON THE INFILE ARRAY,
C     PRODUCING THE FORTRAN SOURCE ON THE OUTFIL ARRAY.
C
      character BCD*90
      character OUTFIL(*)*1
      integer INFILE(LENGTH)
C
      LENLFT = JARROW
      IARROW=0
      INDEX=1
   10 IARROW=IARROW+1
      if(INFILE(INDEX).EQ.-IDELIM .OR. INFILE(INDEX).EQ.IDELIM) THEN
        if(LENLFT.LT.IARROW) GO TO 70
        if(INFILE(INDEX).EQ.-IDELIM) THEN
          OUTFIL(IARROW) = '('
        else
          OUTFIL(IARROW) = ')'
        end if
        INDEX = INDEX + 1
        if(INDEX.GT.LENGTH) THEN
          JARROW=IARROW
          return
        end if
      else if(INFILE(INDEX + 1) .NE. -IDELIM) THEN
        if(LENLFT .LT. (IARROW + 2)) GO TO 70
        if(INFILE(INDEX) .LT. 0) THEN
          OUTFIL(IARROW) = '-'
        else
          OUTFIL(IARROW) = '+'
        end if
        OUTFIL(IARROW + 1) = 'K'
        OUTFIL(IARROW + 2) = '('
        IARROW = IARROW + 3
        IGOTO=INFILE(INDEX + 1)
        MAXLEN = LENLFT - IARROW + 1
        CALL EQSPLT (OUTFIL, IARROW, IGOTO, MAXLEN)
        if(LENLFT .LT. IARROW) GO TO 70
        OUTFIL(IARROW) = ')'
        ISTOP = ABS(INFILE(INDEX))
        if(ISTOP .NE. 1) THEN
          do 30 i = INDEX + 2, INDEX + ISTOP
            if(LENLFT .LT. (IARROW + 3)) GO TO 70
            OUTFIL(IARROW + 1) = '*'
            OUTFIL(IARROW + 2) = 'Y'
            OUTFIL(IARROW + 3) = '('
            IARROW = IARROW + 4
            IGOTO = INFILE(i)
            MAXLEN = LENLFT - IARROW + 1
            CALL EQSPLT(OUTFIL, IARROW, IGOTO, MAXLEN)
            if(LENLFT .LT. IARROW) GO TO 70
   30       OUTFIL(IARROW) = ')'
        end if
        INDEX = INDEX + ISTOP + 1
      else
         if(INFILE(INDEX - 1) .EQ. IDELIM) THEN
           if(LENLFT .LT. IARROW) GO TO 70
           OUTFIL(IARROW) = '+'
           IARROW = IARROW + 1
         end if
         if(LENLFT .LT. (IARROW + 1)) GO TO 70
         OUTFIL(IARROW) = 'Y'
         OUTFIL(IARROW+1) = '('
         IARROW=IARROW+2
         IGOTO=INFILE(INDEX)
         MAXLEN = LENLFT - IARROW + 1
         CALL EQSPLT (OUTFIL, IARROW, IGOTO, MAXLEN)
         if(LENLFT.LT.(IARROW+1)) GO TO 70
         OUTFIL(IARROW) = ')'
         OUTFIL(IARROW+1) = '*'
         IARROW=IARROW+1
         INDEX = INDEX + 1
      end if
      GO TO 10
C
   70 write(BCD, '(''EQOUT10FE IN subroutine EQOUT, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE CWORK ARRAY.'')')
      CALL XERROR(BCD, 90, 10, 2)
      return
      end
C
      subroutine EQSPLT (OUTFIL, IC, IWORD, MAXLEN)
C
C     THIS subroutine SPLITS UP A BINARY WORD INTO ITS BCD COMPONENTS.
C
      character BCD*92
      character OUTFIL(*)*1, W(10)*1
      data (W(i), i=1,10) /'0', '1', '2', '3', '4', '5', '6', '7', '8',
     8    '9'/
C
      IMAX = 1
      IJ = 1
   10 if(10*IJ.LE.IWORD) THEN
        IMAX = IMAX + 1
        IJ = 10*IJ
        GO TO 10
      end if
      if(MAXLEN.LT.IMAX) THEN
        write(BCD, '(''EQSPLT11FE IN subroutine EQSPLT, THERE WAS'',
     8  '' INSUFFICIENT SPACE ALLOCATED FOR THE CWORK ARRAY.'')')
        CALL XERROR(BCD, 92, 11, 2)
        return
      end if
      i = 1
   30 II = MOD(IWORD/IJ,10) + 1
      OUTFIL(IC)=W(II)
      IC=IC+1
      IJ=IJ/10
      i = i + 1
      if(i.LE.IMAX) GO TO 30
      return
      end
C
      subroutine EQPRNT(CWORK, NTAG, ISTART, ISTOP, IPRNT, NCHAR, LOUT)
C
C     THIS subroutine PUTS THE FORTRAN SOURCE ON THE FILE WITH LOGICAL
C     unit NUMBER LOUT.
C
      character CWORK(*)*1
      if(IPRNT.EQ.1) THEN
        if((NTAG+ISTOP-ISTART+1).LE.NCHAR) THEN
          write(LOUT, '(6X, 66A1)')
     8    (CWORK(i+1), i=1,NTAG), (CWORK(i), i=ISTART,ISTOP)
        else
          write(LOUT, '(6X, 66A1 /(5X, ''8 '', 65A1))')
     8    (CWORK(i+1), i=1,NTAG), (CWORK(i), i=ISTART,ISTOP)
        end if
      else if(IPRNT.EQ.2) THEN
        if((2*NTAG+ISTOP-ISTART+1).LE.NCHAR) THEN
          write(LOUT, '(6X, 66A1)')
     8    (CWORK(i+1), i=1,NTAG), (CWORK(i), i=ISTART,ISTOP),
     8    (CWORK(i), i=1,NTAG)
        else
          write(LOUT, '(6X, 66A1 /(5X, ''8 '', 65A1))')
     8    (CWORK(i+1), i=1,NTAG), (CWORK(i), i=ISTART,ISTOP),
     8    (CWORK(i), i=1,NTAG)
        end if
      else if(IPRNT.EQ.3) THEN
        if((NTAG+ISTOP-ISTART+2).LE.NCHAR) THEN
          write(LOUT, '(6X, 66A1)')
     8    (CWORK(i+1), i=1,NTAG), (CWORK(i), i=ISTART,ISTOP), ')'
        else
          write(LOUT, '(6X, 66A1 /(5X, ''8 '', 65A1))')
     8    (CWORK(i + 1), i = 1, NTAG), (CWORK(i), i=ISTART, ISTOP), ')'
        end if
      else if(IPRNT .EQ. 4) THEN
        if((2*NTAG + ISTOP - ISTART + 2) .LE. NCHAR) THEN
          write(LOUT, '(6X, 66A1)')
     8    (CWORK(i+1), i = 1, NTAG), (CWORK(i), i = ISTART, ISTOP), ')',
     8    (CWORK(i), i = 1, NTAG)
        else
          write(LOUT, '(6X, 66A1 /(5X, ''8 '', 65A1))')
     8    (CWORK(i+1), i = 1, NTAG), (CWORK(i), i = ISTART, ISTOP), ')',
     8    (CWORK(i), i = 1, NTAG)
        end if
      end if
      return
      end
