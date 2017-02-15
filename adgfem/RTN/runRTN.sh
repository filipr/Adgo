#!/bin/csh 

./scalar/scalar_clean.sh
./scalar.sh 125 125 1 1 1 1 0.1 1 RTN2
./mmake 
./Exe
./copy_files.sh CASE2_h125_p1_q1


#cp template.pdf template1.pdf 




# ./scalar.sh 125 2000 1 1 1 1 0.05 3 RTN2
#./mmake 
#./Exe

# CASE 42 RTNtime

#./scalar.sh 128 128 1 1 1 1 0.1 1 RTNtime
#### ./scalar.sh 32 32 1 1 0 0 0.2 2 RTN
#./mmake 
#./Exe
#cp template.pdf ledenEx7_S128_P1_Q1T01.pdf

#okular template.pdf

# Case 7 RTN
#./scalar.sh 128 4000 1 1 1 1 0.05 3 RTN

#./mmake 
#./Exe
#okular template.pdf



#cp template.pdf template_y1.pdf 

# ./scalar.sh 125 2000 1 1 0 0 0.05 3 RTN2
#./mmake 
#./Exe


#cp template.pdf template_x2.pdf

# ./scalar.sh 125 500 3 3 1 1 0.05 2 RTN2
#./mmake 
#./Exe

#cp template.pdf template_x3.pdf

# CASE 2 RTN2
# ./scalar.sh 250 4000 2 2 1 1 0.025 3 RTN2
#./mmake 
#./Exe

#./scalar.sh 125 2000 3 3 1 1 0.05 3 RTN2
#./mmake 
#./Exe




#cp template.pdf ZARI3_HTAU_CASE2_P2Q1_125.pdf
