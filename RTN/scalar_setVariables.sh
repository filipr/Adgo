#!/bin/csh

set mesh_name = square01.str
set problem = scalar

set Re = 1.0
set ttime = 0.5
set Miter = 1000
#set Miter = -1
#set nBC = 4
#set Ftau = 1E+15
set Ftau = 0.
set algeb = 1E-5
set sigma = 20.
set tau_type = "fixed"
set estim_space = "RES" 
set non_solver = "rezL2" 

set lin_tol = 1E-12


if($argv[9] == 'A') then
  set icase = 44
  set Re = 0.0
  # paramet is the size of reaction !!!!
  set paramet = 10.0   
  set nBC = "-2 1 3"
  set iord = 3
  set type_tab = "STDGM"


else if ($argv[9] == 'STDG') then 
   set icase = 42
   set paramet = 0.0
   set ttime = 1.0 
   set nBC = 4
   set iord = 3
   set Re = 1.E+3
   set type_tab = "STDGM"
   
else if ($argv[9] == 'ODE') then 
   set icase = 46
   set paramet = 100.0
   set ttime = 1.0 
   set nBC = 4
   set iord = 3
   set Re = 0.0
   set type_tab = "STDGM" 
   
else if ($argv[9] == 'RTN' ) then 
   set icase = 7 
   set paramet = 0.0 
   set ttime = 1.0 
   set nBC = 8 
   set iord = 3 
   set Re = 1.0 
   set estim_space = "RTNst"
   set type_tab = "RTNst"


else if($argv[9] == 'B') then
  set icase = 45
  set Re = 0.0
  #set paramet = 10.0
  set paramet = 0.0
  set nBC = "-2 1 3"
  set iord = 2  

  set type_tab = "AMAstO"

else if($argv[9] == 'C') then
  set icase = 49
  set Re = 10.0
  set paramet = 0.0
  set nBC = 4
  #set nBC = "-2 1 3"
  set iord = 3
  set type_tab = "RES_STDGMp"


else if($argv[9] == 'D') then
  set icase = 48
  set Re = 10.0
  set paramet = 0.0
  #set nBC = "-2 1 3"
  set nBC = 4
  set iord = 2
  set type_tab = "RES_STDGMp"

else if($argv[9] == 'F') then
  set icase = 7
  set Re = 1.0
  set paramet = 0.0
  #set nBC = "-2 1 3"
  set nBC = 4
  set iord = 3
  set tau_type = "fixed"
  
  set type_tab = "RES_STDGM"
  


else if($argv[9] == 'G') then
  set icase = 51
  #set Re = 1.
  set Re = 10.
    echo ' set Re was changed !!!'
  set paramet = 0.0
  #set nBC = "-2 1 3"
  set nBC = 4
  set iord = 3
  set type_tab = "RES_STDGMp"


else if($argv[9] == 'H') then
  set icase = 52
  set Re = 10
  set paramet = -1.5
  #set nBC = "-2 1 3"
  set nBC = 4
  set iord = 2
  set type_tab = "RES_STDGMp"
  set tau_type = "adapt"
  #set tau_type = "fixed"


else if($argv[9] == 'I') then
  set icase = 2
  set mesh_name = square-11.str
  set Re = 100.
  set paramet = 0.0
  #set nBC = "-2 1 3"
  set nBC = 4
  set iord = 2
  set type_tab = "RES_STDGM"
  #set tau_type = "adapt"
  #set tau_type = "fixed"


else
 echo 'UNKNOWN CASE'
endif


