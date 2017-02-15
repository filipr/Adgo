#!/bin/csh

##                              1       2          3   4  5    6  7   8  9  10  11   12  13   14   15    16   17   18   19   20    21   22    23    24  25 
## Syntax:  gen_ini_file.sh  ini_name  grid        Qk prf dim mod par Pk IES Tk BDF iter adt alev  TOL   tol Ttim  tau  CFL  BDFt  Re tol_mx tol_mi Lq  nBC
../Scripts/gen_ini_file.sh test1 square-11.str.250 1   0   1   2  0.0  2  I  1  BDF  10  Ahp  100 1E+00 1E+5 1.00  1E-2 1E+5 0.2  1E+2  5E-5  5E-3  2.0  5 
