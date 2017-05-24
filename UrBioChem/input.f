      subroutine input(NyDim,Y,totdens,temp,rh,DiskH,tlg,dtlg,tLim,tOut,
     1     ewt,eps,nspeci,Nstate,NTask,NRoot,iError,modelname)
 
***********************************************************
*     Read the Output file from a prior run of the code and
*     use the new output parameters as an input starting guess
*     for a new run of the code.  The parameters are passed back
*     to the main code and are se[parately saved in the file NewInput.
*     A. Bayless March, 2017
*************************************************************
*     NyDim: Dimention of the Y array (# of elements) from main code
*     ylet, paren, bang (!) are used to parse out the Y values and 
*     element name from the string in Output. 
*     ylet: 'Y(', paren: ') =', and bang: '!' in Y( ###) = # !element lines
      real*8 Y(NyDim)
      real*8 ewt(NyDim)
      real*8 totdens,temp,rh,DiskH,negden,posden
      real*8 tLim,tOut,eps, xj
      CHARACTER*8 spelmnt(NyDim),element(NyDim)
      INTEGER*4 chrg(NyDim)
      integer j,nspeci,Nstate,NTask,NRoot,iError
      character*8 ylet,paren,bang,unit
      character*50 modelname,label,outdate,outtime
      CHARACTER*8 date
      CHARACTER*10 time
      CHARACTER*5 zone
      INTEGER*4 VALUES(8),ij 
 
      open (unit = 80, file = "NewInput_G1", status = "unknown")
      open (unit = 90, file = "Output_G1", status = "old")
      open (unit = 95, file = "SpeciesSN1", status = "old")

!     Initalize densities for positive = negative diagnostic check
      negden = 0.0
      posden = 0.0


      read(90,fmt="(a50)") modelname 
      read(90,fmt="(a25)") outdate !date output was made
      read(90,fmt="(a25)") outtime !time output was made
            

      read(90,fmt="(a15,E11.4,a)") label,totdens,unit
      read(90,fmt="(a6,F11.4,a1,a)") label,temp,unit
      read(90,fmt="(a4,f8.2,a3)") label,rh,unit
      read(90,fmt="(a7,f8.2,a3)") label,DiskH,unit
      read(90,fmt="(a8,i4)") label,nspeci
      read(90,fmt="(a6,f10.2)") label, tlg
      read(90,fmt="(a7,f10.2)") label, dtlg
      read(90,fmt="(a7,f10.2)") label, tLim
      read(90,fmt="(a9,i3)") label, NState
      read(90,fmt="(a7,f10.2)") label, tOut
      read(90,fmt="(a8,i2)") label, NTask
      read(90,fmt="(a8,i2)") label, NRoot
      read(90,fmt="(a6,f13.6)") label,eps
      read(90,fmt="(a9,i2)") label, iError
      do i = 1, nspeci
          read(90,fmt="(a12,f11.4)") label, ewt(i)
       enddo
      do i = 1, nspeci
         read(90,fmt="(a,i3,a3,E11.4,a2,a8)") 
     1        ylet,j,paren,Y(i),bang,element(i)
  
      enddo

******************************************************
*     Write NewInput as a record of the new parameters 
*     being used in the in the new run of the code.
*     A. Bayless March, 2017
******************************************************
      write(80,fmt="(a50)") modelname
      write(80,fmt="(a25)") outdate
      write(80,fmt="(a25)") outtime

      call date_and_time( date, time, zone, values ) 

      write(80,fmt="(a20,i4,a1,i2,a1,i2)") 
     1     'Input Date Created: ',values(1),"/",values(2),"/",values(3)
      write(80,fmt="(a20,i2,a1,i2,a1,i2)")
     1     'Input Time Created: ',values(5),":",values(6),":",values(7)

********************************************************
*     Print some diagnostic information to the file
********************************************************

!     If Total density is more than 10% of M = Y(nspeci) [last Y],
!     Then print warning and stop code
      
      if (abs(totdens/Y(nspeci) - 1.0) >  0.1) then  
         write(*,*) "  "
         write(*,*) "In input.f"
         write(*,fmt = "(a,f6.2,a)") "Warning: TotDens",
     1        abs(totdens/Y(nspeci)-1.0)*100, "% from M"
         write(*,*) "M: ",Y(nspeci)
         write(*,*) "TotDens: ",totdens
         write(*,*) "  "
         stop
      endif
c------------------------------------------------------------
!     Check that the Y elements match the Species file order

!     Read SpecRect file lines
!     Header Lines  -- don't know what everything is 
      read (95,*) modelname
      read (95,*) label
      read (95,*) jnk,jnk,jnk,jnk,jnk,jnk

!     Read Data Lines
      do i = 1, nspeci 
         read(95,fmt="(a8,i3)") spelmnt(i),chrg(i)
!     Check that the Y elements match the Species file order
!     If not, then print warning to screen and stop code
!     NOTE: THIS IS SENSITIVE TO LEADING AND TRAILING WHITESPACE
         if (spelmnt(i) /= element(i)) then
            write(*,*) "  "
            write(*,*) "In input.f"
            write(*,*) "Species in Y and Species File are NOT EQUAL "
            write(*,*) "SpecReac: ",spelmnt(i), "Y: ",element(i)
            write(*,*) "  "
            stop
         endif
 
c------------------------------------------------------------
!     Check that the density of the postive species matches
!     the density of the negative species to within 10%

         if (chrg(i) == -1) then
            negden = negden + Y(i)
         endif
         
         if (chrg(i) == 1) then
            posden = posden + Y(i)
         endif      
         
         

      enddo
      
      if (abs(posden/negden - 1.0) >  0.1) then  
         write(*,*) "  "
         write(*,*) "In input.f"
         write(80,fmt = "(a,f6.2,a)") "Warning: Positive Density",
     1        abs(posden/negden-1.0)*100, "% from Negative Density"
         write(80,*) "Positive Density: ",posden
         write(80,*) "Negative Density: ",negden
         write(80,*) "  "
         write(*,fmt = "(a,f6.2,a)") "Warning: Positive Density",
     1        abs(posden/negden-1.0)*100, "% from Negative Density"
         write(*,*) "Positive Density: ",posden
         write(*,*) "Negative Density: ",negden
         write(*,*) "  "
c     stop
      endif

c-------------------------------------------------------------
c    Calculate the new electron abundance for Y(1).
!negden includes e- abundance, subtract out e- (old Y(1)) to get negative ion 
!abundance only

      negion = negden - Y(1) !negden includes e- abundance
      Y(1) = posden - negion !new e- abundance
      write(*,*) "New ", element(1), Y(1)
c------------------------------------------
c     Check to see that the new density is within 10%
c     of solar values.  If not, print warning to screen and stop
c     Subroutine in this file below.
      call checksolar(NyDim,nspeci,element,Y)

*********************************************************
*     Print the new input parameters to a file for saving
*********************************************************

      write(80, fmt = "(a15, 1pe10.3, a9)") 
     1    "Total Density =", totdens, " cm^(-3) "
            
      write(80, fmt = "(a6, e10.3, a2)") "Temp =", Temp, "K"
      write(80,fmt="(a4,f8.2,a3)") "Rh =", rh,"AU"
      write(80,fmt="(a7,f8.2,a3)") "DiskH =",DiskH,"AU"
      write(80,fmt = "(a8, i4)") "NSpeci =", NSpeci
      write(80,fmt="(a6,f10.2)") "tlg = ", tlg
      write(80,fmt="(a7,f10.2)") "dtlg = ", dtlg
      write(80,fmt="(a7,f10.2)") "tLim = ", tLim
      write(80,fmt="(a9,i3)") "Nstate = ", NState
      write(80,fmt="(a7,f10.2)") "tOut = ", tOut
      write(80,fmt="(a8,i2)") "NTask = ", NTask
      write(80,fmt="(a8,i2)") "NRoot = ", NRoot
      write(80,fmt="(a6,f13.6)") "eps = ", eps
      write(80,fmt="(a9,i2)") "iError = ", iError

       do i = 1, nspeci
          write(80,fmt="(a4,i3,a5,f11.4)") "ewt(",i,") = ", ewt(i)
      enddo
    
      do i = 1, nspeci
         write(80,fmt="(a,i3,a3,E11.4,a1,a)") 
     1        ylet,i,paren,Y(i),bang,element(i)
      enddo
      
      close(90)
      close(80)
      return
      end

*********************************************************************
      subroutine checksolar(NyDim,nspeci,element,Y)
C     Checks to see that the new densities are still within 10%
C     of solar values.  If not print a warning and stop.
*********************************************************************
      real*8 Y(NyDim)
      real H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K
      real Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn 
      integer nspeci,jnk,chrg
      integer nH(NyDim),nHe(NyDim),nC(NyDim),nN(NyDim)
      integer nO(NyDim),nS(NyDim),nF(NyDim),nP(NyDim)
      CHARACTER*8 element(NyDim),spelmnt
      real*8 Y_H,Y_He,Y_C,Y_N,Y_O,Y_S !Sum of ion abundances
      real*8 rho_H,rho_He,rho_C,rho_N,rho_O,rho_S !Sum of ion abundances
      character*50 modelname,label
    
      open (unit = 96, file = "SpeciesSN1", status = "old")


!     Log Solar photospheric number abundances:  Asplund, M., Grevesse, 
!       N., Sauval, A. J., Scott, P.,  The chemical composition of the 
!       Sun, Annual Rev. Astron. Astrophys. 47, 481-522 (2009).
      H = 12.00
      He = 10.93
      Li = 1.05
      Be = 1.38
      B  =  2.70 
      C  =  8.43
      N  =  7.83 
      O  =  8.69 
      F  =  4.56 
      Ne =  7.93 
      Na =  6.24 
      Mg =  7.60 
      Al =  6.45 
      Si =  7.51 
      P  =  5.41 
      S  =  7.12 
      Cl =  5.50 
      Ar =  6.40 
      K  =  5.03 
      Ca =  6.34 
      Sc =  3.15 
      Ti =  4.95 
      V  =  3.93 
      Cr =  5.64 
      Mn =  5.43 
      Fe =  7.50 
      Co =  4.99 
      Ni =  6.22 
      Cu =  4.19 
      Zn =  4.56 


c     Sum up all the relevant species 
      Y_H = 0.0
      Y_He = 0.0
      Y_C = 0.0
      Y_N = 0.0
      Y_O = 0.0
      Y_S = 0.0

!     Read SpecRect file lines
!     Header Lines  -- don't know what everything is 
      read (96,*) modelname
      read (96,*) label
      read (96,*) jnk,jnk,jnk,jnk,jnk,jnk

!     Read Data Lines 
      do i = 1, nspeci 
         read(96,fmt="(a8,9i3)")spelmnt,chrg,
     1        nH(i),nHe(i),nC(i),nN(i),nO(i),nF(i),nP(i),nS(i)
      enddo
      close(96)
      
C     NOTE:  Because the sum is multiplied by the number of elements in
C     the molecule, He will not be summed with H even with a filter for "H", i.e. nHe=0.
C     O(1S) will only sum O and not read it as O and S molecule becuase nS=0, etc.

      do i = 1, nspeci
c        write(*,*) "Loop Element:",element(i)
C     Find the species that contain "H"
         int = INDEX(element(i), "H") 
         if (int > 0) then 
c           write(*,*) "Element:",element(i)
c           write(*,*) "Y(i)",Y(i)
c           write(*,*) "Old Y_H",Y_H
            Y_H = Y(i)*nH(i) + Y_H
c           write(*,*) "New Y_H", Y_H
c           write(*,*) " "
         endif

C     Find the species that contain "He"
         int = INDEX(element(i), "He")
         if (int > 0) then
c           write(*,*) "Element:",element(i)
c           write(*,*) "Y(i)",Y(i)
c           write(*,*) "Old Y_He",Y_He
            Y_He = Y(i)*nHe(i) + Y_He
c           write(*,*) "New Y_He", Y_He 
c           write(*,*) " "
        endif

C     Find the species that contain "C"
         int = INDEX(element(i), "C")
         if (int > 0) then
c             write(*,*) "Element:",element(i)
c             write(*,*) "Y(i)",Y(i)
c             write(*,*) "Old Y_C",Y_C
              Y_C = Y(i)*nC(i) + Y_C
c             write(*,*) "New Y_C", Y_C
c             write(*,*) " "
           endif


C     Find the species that contain "N"
         int = INDEX(element(i), "N")
         if (int > 0) then
c           write(*,*) "Element:",element(i)
c           write(*,*) "Y(i)",Y(i)
c           write(*,*) "Old Y_N",Y_N
            Y_N = Y(i)*nN(i) + Y_N
c           write(*,*) "Y_N", Y_N
c           write(*,*) " "
        endif

C     Find the species that contain "O"
         int = INDEX(element(i), "O")
         if (int > 0) then
c           write(*,*) "Element:",element(i)
c           write(*,*) "Y(i)",Y(i)
c           write(*,*) "Old Y_O",Y_O
            Y_O = Y(i)*nO(i) + Y_O
c           write(*,*) "Y_O", Y_O 
c           write(*,*) " "
        endif
      
C     Find the species that contain "S" -- will need other checks if Si is ever added
        int = INDEX(element(i), "S")
        if (int > 0) then
c          write(*,*) "Element:",element(i)
c          write(*,*) "Y(i)",Y(i)
c          write(*,*) "Old Y_S",Y_S
            Y_S = Y(i)*nS(i) + Y_S
c          write(*,*) "Y_S", Y_S
c          write(*,*) " "
        endif
      enddo

      write(*,fmt="(A,E13.5)") "Total H:",Y_H
      write(*,fmt="(A,E13.5)") "Total He:",Y_He
      write(*,fmt="(A,E13.5)") "Total C:",Y_C
      write(*,fmt="(A,E13.5)") "Total O:",Y_O
      write(*,fmt="(A,E13.5)") "Total N:",Y_N
      write(*,fmt="(A,E13.5)") "Total S:",Y_S
  
!     Calculate the new abundances and compare to solar
      rho_He = log10(Y_He/Y_H) + 12.0
      rho_C = log10(Y_C/Y_H) + 12.0
      rho_O = log10(Y_O/Y_H) + 12.0
      rho_N = log10(Y_N/Y_H) + 12.0
      rho_S = log10(Y_S/Y_H) + 12.0

      write(*,*) "New Helium: ", rho_He, "Solar: ",He
      write(*,*) "New Carbon: ", rho_C, "Solar: ",C
      write(*,*) "New Oxygen: ", rho_O, "Solar: ",O
      write(*,*) "New Nitrogen: ", rho_N, "Solar: ",N
      write(*,*) "New Sulfur: ", rho_S, "Solar: ",S
      
      write(80,*) "New Helium: ", rho_He, "Solar: ",He
      write(80,*) "New Carbon: ", rho_C, "Solar: ",C
      write(80,*) "New Oxygen: ", rho_O, "Solar: ",O
      write(80,*) "New Nitrogen: ", rho_N, "Solar: ",N
      write(80,*) "New Sulfur: ", rho_S, "Solar: ",S
 
      close(96)
      return
      end
