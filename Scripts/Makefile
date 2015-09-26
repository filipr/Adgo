
# -fPIC  ... warning: variable is never used
#FFLAGS=  -fPIC  -w
FFLAGS=  -fPIC  -g -fbacktrace -fbounds-check -w
#FFLAGS=  -fPIC -fdefault-real-8  -w
#FFLAGS=  -fPIC -fdefault-real-8 -O2 -w
#FFLAGS=  -fPIC -fdefault-real-8 -g -fbacktrace -fbounds-check -w

#LIBS=-llapack -lpgplot 
#LIBS= LAPACK/lapack_Linux.a LAPACK/blas_Linux.a

#FXX=g95
FXX=gfortran
#FXX=ifort
all: Setorderx extract_order_dat power pow div div_min0.1 multi setnum setnum2 moddiv gen_tisk_tex gen_tisk_texAD extr-order multiM read_icase read_icaseAD vel_1  vel_2 mach AD_gen_ini_file AD_gen_ini_file_ST gen_exact_Barenblatt re_order_dat remake_ord remake_ord_NS read_integ extract

#all: Setorderx extract_order_dat power pow div multi setnum setnum2 moddiv gen_tisk_tex gen_tisk_texAD extr-order multiM read_icase read_icaseAD vel_1  vel_2 mach AD_gen_ini_file AD_gen_ini_file_ST gen_exact_Barenblatt re_order_dat

# gen_exact_Barenblatt removed could not be translated
all: Setorderx extract_order_dat power pow div multi setnum setnum2 moddiv gen_tisk_tex gen_tisk_texAD extr-order multiM read_icase read_icaseAD vel_1  vel_2 mach AD_gen_ini_file AD_gen_ini_file_ST re_order_dat
#AD_gen_ini_fileST 



Setorderx   : Setorderx.o
	$(FXX) $(FFLAGS) -O -C -o Setorderx Setorderx.o 

read_icase   : read_icase.o
	$(FXX) $(FFLAGS) -O -C -o read_icase read_icase.o 

read_icaseAD   : read_icaseAD.o
	$(FXX) $(FFLAGS) -O -C -o read_icaseAD read_icaseAD.o 

extract_order_dat   : extract_order_dat.o
	$(FXX) $(FFLAGS) -O -C -o extract_order_dat extract_order_dat.o 

power   : power.o
	$(FXX) $(FFLAGS) -O -C -o power power.o 

pow   : pow.o
	$(FXX) $(FFLAGS) -O -C -o pow pow.o 

div   : div.o
	$(FXX) $(FFLAGS) -O -C -o div div.o 

div_min0.1   : div_min0.1.o
	$(FXX) $(FFLAGS) -O -C -o div_min0.1 div_min0.1.o 

multi   : multi.o
	$(FXX) $(FFLAGS) -O -C -o multi multi.o 

multiM   : multiM.o
	$(FXX) $(FFLAGS) -O -C -o multiM multiM.o 

setnum   : setnum.o
	$(FXX) $(FFLAGS) -O -C -o setnum setnum.o 

setnum2   : setnum2.o
	$(FXX) $(FFLAGS) -O -C -o setnum2 setnum2.o 

moddiv   : moddiv.o
	$(FXX) $(FFLAGS) -O -C -o moddiv moddiv.o 

gen_tisk_tex   : gen_tisk_tex.o
	$(FXX) $(FFLAGS) -O -C -o gen_tisk_tex gen_tisk_tex.o 

gen_tisk_texAD   : gen_tisk_texAD.o
	$(FXX) $(FFLAGS) -O -C -o gen_tisk_texAD gen_tisk_texAD.o 

extr-order   : extr-order.o
	$(FXX) $(FFLAGS) -O -C -o extr-order extr-order.o 

vel_1   : vel_1.o
	$(FXX) $(FFLAGS) -O -C -o vel_1 vel_1.o 

vel_2   : vel_2.o
	$(FXX) $(FFLAGS) -O -C -o vel_2 vel_2.o 

mach   : mach.o
	$(FXX) $(FFLAGS) -O -C -o mach mach.o 

extract   : extract.o
	$(FXX) $(FFLAGS) -O -C -o extract extract.o 

AD_gen_ini_file   : AD_gen_ini_file.o
	$(FXX) $(FFLAGS) -O -C -o AD_gen_ini_file AD_gen_ini_file.o 

AD_gen_ini_file_ST   : AD_gen_ini_file_ST.o
	$(FXX) $(FFLAGS) -O -C -o AD_gen_ini_file_ST AD_gen_ini_file_ST.o 

#AD_gen_ini_fileST   : AD_gen_ini_fileST.o
#	$(FXX) $(FFLAGS) -O -C -o AD_gen_ini_fileST AD_gen_ini_fileST.o 

re_order_dat   : re_order_dat.o
	$(FXX) $(FFLAGS) -O -C -o re_order_dat re_order_dat.o 

remake_ord   : remake_ord.o
	$(FXX) $(FFLAGS) -O -C -o remake_ord remake_ord.o 

remake_ord_NS   : remake_ord_NS.o
	$(FXX) $(FFLAGS) -O -C -o remake_ord_NS remake_ord_NS.o 

gen_exact_Barenblatt   : gen_exact_Barenblatt.o
	$(FXX) $(FFLAGS) -O -C -o gen_exact_Barenblatt gen_exact_Barenblatt.o 


read_integ   : read_integ.o
	$(FXX) $(FFLAGS) -O -C -o read_integ read_integ.o 




##$(LIBS) 

#


clean:
	-rm -f Setorderx power pow div  setnum moddiv gen_tisk_tex  gen_tisk_texAD  *.o *.mod extr-order read_icase read_icaseAD vel_1  vel_2 mach gen_exact_Barenblatt AD_gen_ini_file re_order_dat remake_ord remake_ord_NS div_min0.1, read_integ extract

%.o:%.f90
	$(FXX) $(FFLAGS) -c $? 

%.o:%.f
	$(FXX) $(FFLAGS) -c $? 


