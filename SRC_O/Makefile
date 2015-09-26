#TARGETS= lap_sub.o geom.o integ.o f_mapping.o  rtn.o state.o dmodules.o matrix.o mesh3D.o mesh.o gmres.o blocks.o model.o model3D.o eval_sol.o io_sub.o problem.o inv_fluxes.o euler.o project.o aposter.o solve.o  f90_kind.o sparse_utils.o UMFPACK/UMFPACK/Demo/umf4_f77wrapper.o umfpack.o fem.o struct.o init.o calc.o inserting.o removing.o swapping.o moving.o rep_cbp.o angener.o marking.o adaptation.o hp_adapt.o main.o


include make.inc

# -fPIC  ... warning: variable is never used
#FFLAGS=  -fPIC -fdefault-real-8  -w

# optimize for computations
#FFLAGS=  -fPIC -fdefault-real-8  -fopenmp -O2 -w

# for debugging of the code
FFLAGS= -fPIC -fdefault-real-8 -g -fbacktrace -fbounds-check -w -Wall -finit-real=nan -finit-integer=-999999 -ffpe-trap=invalid -fno-align-commons

#FFLAGS=-Og -fPIC -fdefault-real-8 -fimplicit-none -Wall  -Wextra  -Wunused-parameter  -fwhole-file  -fcheck=all -fbacktrace

# -Wuninitialized

#-debug -traceback -check -ftrapuv

#FFLAGS2=  -fPIC -fdefault-real-8 -g -fbacktrace -fbounds-check -w -Wall -finit-real=nan -finit-integer=-999999 -ffpe-trap=invalid -fno-align-commons

#-fdump-parse-tree
#-fdump-core
#-finit-local-zero
#-falign-commons
#-fno-align-commons
#-fcheck-array-temporaries
#-fimplicit-none

#FFLAGS1=  -g -fbacktrace -fbounds-check -w
FFLAGS1=   -O2 -w

#AGMG = ./AGMG_3.2.1-aca/
#AGMGLIB = $(AGMG)libagmg.a
AGMG =
AGMGLIB =

#LIBS=-llapack -lpgplot
#LIBS=-lblas

# without UMFPACK
LIBS= $(AGMGLIB) LAPACK/lapack_Linux.a LAPACK/blas_Linux.a

# with UMFPACK
#LIBS= LAPACK/lapack_Linux.a LAPACK/blas_Linux.a UMFPACK/UMFPACK/Lib/libumfpack.a UMFPACK/AMD/Lib/libamd.a -lblas -lgfortran -lgfortranbegin  -lm

##LIBS= LAPACK/lapack_Linux.a LAPACK/blas_Linux.a UMFPACK/UMFPACK/Lib/libumfpack.a UMFPACK/AMD/Lib/libamd.a -lblas -lgfortran -lgfortranbegin -lg2c  -lm


#LIBS=-llapack -L/usr/local/pgplot -lpgplot  -lX11
#-lpgplot

#FFLAGS=  -fPIC -r8
#FXX=g95
FXX=gfortran
#FXX=ifort

all: SuperADGo dgfem2fvmx dgfem_range plotdgm  cfig Setorderx gen_tisk_texAD #fvm2jpg
#all: Adgfem AAdgfem  dgfem2fvmx dgfem_range plotdgm  cfig Setorderx gen_tisk_texAD #fvm2jpg
#Debug: Adgfem AAdgfem ADGo dgfem2fvmx dgfem_range plotdgm  cfig Setorderx gen_tisk_texAD #fvm2jpg

SuperADGo: $(TARGETS_O)
	$(FXX) $(FFLAGS) -o SuperADGo $^ $(LIBS)

#-L/usr/local/pgplot  -lpgplot
#-L/usr/local/lib -lslatec -lg2c

first: blaslib lapacklib agmglib SuperADGo

agmglib:
	( cd $(AGMG); $(MAKE) )

blaslib:
	( cd LAPACK/BLAS/SRC; $(MAKE) )

lapacklib:
	( cd LAPACK/SRC; $(MAKE) )

plotdgm: plgeom.o plotdgm.o
	$(FXX) $(FFLAGS) -o plotdgm $^

dgfem2fvmx: plgeom.o dgfem2fvmx.o
	$(FXX) $(FFLAGS) -o dgfem2fvmx $^

dgfem_range:  plgeom.o dgfem_range.o
	$(FXX) $(FFLAGS) -o dgfem_range $^

fvm2jpg: fvm2jpg.f
	gfortran -o fvm2jpg fvm2jpg.f -L/usr/local/ -lpgplot -lX11

cfig   : GGG.o cfig.o
	$(FXX) $(FFLAGS1) -O -C -o cfig cfig.o GGG.o
GGG.o : GGG.f
	$(FXX) $(FFLAGS1) -O -C -c GGG.f
cfig.o : cfig.f90
	$(FXX) $(FFLAGS1) -O -C -c cfig.f90

Setorderx : Setorderx.o paramets.o
	$(FXX) $(FFLAGS) -O -C -o Setorderx paramets.o Setorderx.o

gen_tisk_texAD : gen_tisk_texAD.o paramets.o
	$(FXX) $(FFLAGS) -O -C -o gen_tisk_texAD paramets.o gen_tisk_texAD.o

MG:
	( $(MAKE) -f Makeown uncommentmg )
	( $(MAKE) -f Makeown AdgfemUMF AAdgfemUMF )

dir_solver_dummy.o: dir_solver_dummy.f90
	$(FXX) $(FFLAGS) -c $?

dir_solver.o: dir_solver.f90
	$(FXX) $(FFLAGS) -c $?

agmg_solver_dummy.o: agmg_solver_dummy.f90
	$(FXX) $(FFLAGS) -c $?

agmg_solver.o: agmg_solver.f90
	$(FXX) $(FFLAGS) -c $?


#cleanDebug:
#	-rm -f ADGo Adgfem AAdgfem dgfem2fvmx dgfem_range fvm2jpg  Setorderx gen_tisk_texAD *.o *.mod

clean:
	-rm -f SuperADGo dgfem2fvmx dgfem_range fvm2jpg  Setorderx gen_tisk_texAD *.o *.mod

cleanSuperADGo:
	-rm -f SuperADGo *.o *.mod

%.o:%.f90
	$(FXX) $(FFLAGS) -c $?

%.o:%.f95
	$(FXX) $(FFLAGS) -c $?

%.o:%.f
	$(FXX) $(FFLAGS) -c $?

igmres.o: igmres.for
	$(FXX) $(FFLAGS) -c igmres.for

fcolor.o: fcolor.f
	$(FXX) $(FFLAGS) -c fcolor.f

