#FFLAGS=-i4 -fpp -g -traceback -fpe0 -debug extended -diag-disable 8291
FFLAGS=-i4 -fpp -g -traceback -check -fpe0 -debug extended -diag-disable 8291
#FFLAGS=-i4 -xHost -fpp -diag-disable 8291 -O 
FORTRAN=ifort

all:	UrBioChemG1

%.o:	%.f90
	$(FORTRAN) $(FFLAGS) -c -o $@ $<

%.o:	%.f
	$(FORTRAN) $(FFLAGS) -c -o $@ $<

EWS.f:	EqWrt ReactnSN1 SpeciesSN1
	./EqWrt > EqnOut

EqWrt:	EqWrt.o R1Mach.o xerror.o
	$(FORTRAN) $(FFLAGS) EqWrt.o R1Mach.o xerror.o -o EqWrt

UrBioChemG1.o: UrBioChemG1.f Y.h
	$(FORTRAN) $(FFLAGS) -c -o UrBioChemG1.o UrBioChemG1.f

UrBioChemG1: UrBioChemG1.o R1Mach.o xerror.o EWS.o input.o
	$(FORTRAN) $(FFLAGS) UrBioChemG1.o EWS.o R1Mach.o xerror.o input.o -o UrBioChemG1
	ln -f -s SRT SpecReac

clean:
	-rm -f UrBioChemG1.o R1Mach.o xerror.o EWS.f EWS.o EqWrt.o UrBioChemG1.o EqWrt \
	input.o UrBioChemG1 SRT DatFil EqnOut SpecReac Input Output Plot Results
