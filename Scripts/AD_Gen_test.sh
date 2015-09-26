#!/bin/csh

## Syntax: AD_gen_ini_file i_name model Re mod par Ttim  ssTOL ctol     grid        Qk prf IPG cW Pk Tdis Tk tau  Ttol  S_est tol_mx tol_mi Lq madt alev maxI stopNS tol MG nBC M  alp
#../Scripts/AD_gen_ini_file ADtest1 NSe  0.0  0  0.0 1E+30 1E-6 1E-4  NACA.inv.r0     2  nac I  20  2  BDF 1 adapt 10.    RES  1E-4  1E-8  2.0  AMAhp  2  50   aRES 1E-1  no 2  0.5  2
#../Scripts/AD_gen_ini_file ADtest2 sca  1E+2 2  0.0 1E+00 1E-6 1e+5 square-11.str.250 1  0  I  20  3  BDF 2 fixed 0.025  RES  5E-2  5E-4  2.0  HGhp  20   4   aRES 1E-1  no 4 
#../Scripts/AD_gen_ini_file ADtest3 sca  1E+3 3  0.0 1E+30 1E-6 1e+5 square01.uns.50  1   0  I  20  4  BDF 1 fixed 1E+12  RES  2E-6  1E-5  2.0  AMAhp 10  10   aRES 1E-1  no 8 
#../Scripts/AD_gen_ini_file ADtest4 sca  1E+3 4 -1.5 1E+30 1E-6 1e+5 square01.str.250 1   0  I  20  4  BDF 1 fixed 1E+12  RES  1E-3  1E-4  2.0  HGhp  10  10   aRES 1E-1  no 4 
#../Scripts/AD_gen_ini_file ADtest5 NSe  5E+3 0  0.0 1E+30 1E-6 1E-4  NACA.NS.r0      2  nac I  20  2  BDF 1 adapt 10.    RES  1E-4  8E-4  2.0  AMAhp  8  50   aRES 1E-1  no 2  0.5  3

../Scripts/AD_gen_ini_file ADtestST sca 1E+3 1  0.0 1.0   1E-6 1e+5 square01.str.8   1   0  I  20  2 STDG 1 fixed 0.02   RES  1E-4  1E-3  0.0  HGhp   0  10  rezL2 1E-6  no 8 1.E-06
#../Scripts/AD_gen_ini_file ADtestMG sca  1E+2 4 -1.5 1E+30 0.0 1e+5 square01.str.2   1   0  I  20  5  BDF 1 fixed 1E-1   RES  1E-3  1E-4  2.0  HGhp   0   3   aRES 1E-1  lin 4  


#../SRC/AAdgfem ADtest1.ini 
#../SRC/AAdgfem ADtest2.ini 
#../SRC/AAdgfem ADtest3.ini 
