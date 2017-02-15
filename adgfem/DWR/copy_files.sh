#! /bin/csh -f
#
# Example of a shell script to copy all files from one directory 
# to another. The input directory must not contain any subdirectories, 
# and it will not copy any so-called (hidden) dot-files.
#
##                 	 check if called properly
if ($#argv != 1) then
  echo "Usage: $0 dir1"
  echo "copies all files from current directory to another"
  goto done
endif
##                  save command line args in variables
set dir2=$1

##                  check if dir2 does not exist
if (-e $dir2) then
  echo "$dir2 already exists" ; exit 1
endif
##                  create new dir2
mkdir $dir2
if ($status != 0) goto error
##                  loop through all files in dir1
cp hpmesh* *.gnu tri* sol* template.* gnu.* order*.dat tab.txt exa-00* est-00* err-00* exaA* estA* errA* dua-00*   Gfigure.ps results  scalar.conv scalar.sol U_3D-00* U_CUT-00* U_ISO-00* mesh*0* mesh*999* U_CUT*999* U_ISO*999* tisk.pdf  GGGfig.ps errL8-00* fig1-*.ps* fig2-*.ps* *.eps.ppm Fig-00*.jpg *eps-converted-to.pdf GGGraphics.* *.tex errs*.eps error*.eps aDWR_nlErrors $dir2/

##           Labels to jump to exit OK (done) or not OK (error)
done:
echo "All files were copied to $dir2!"
exit 0
error:
exit 1


##!/bin/csh

## copy all files made by vis_tisk.sh to new folder 
#if (  $#argv != 1) then
#  echo ' Syntax: name of the new folder'
#else

#set target=$argv[1] ;

#set ddd=""; 


### Since $c is empty, this will check if the
### file exists in target.
#while [ -e "$target"$ddd/ ]; do
#   echo "$target exists"; 
##  ## If the target exists, add 1 to the value of $c
##  ## and check if a file called $target$c (for example, bar.txt1)
##  ## exists. This loop will continue until $c has a value
##  ## such that there is no file called $target$c in the directory.
##  let ddd++; 
##  target="$target"$ddd; 
#   done; 


#echo 'dsa'
### We now have everything we need, so lets copy.
##cp "$path" "$dest"/"$target"; 

##mkdir $target

##cp *.eps  $argv[1]/

## tab*.tex *.gnu tri* sol* template.* gnu.* order*.dat tab.txt exa-00* est-00* err-00* exaA* estA* errA* dua-00* figures-resC.tex figures-res.tex  Gfigure.ps results tisk.tex scalar.conv scalar.sol U_3D-00* U_CUT-00* U_ISO-00* mesh*0* mesh*999* U_CUT*999* U_ISO*999* errA00* exaA00* estA00* duaA00* tisk.pdf res.tex GGGfig.ps errL8-00* fig1-*.ps* fig2-*.ps* commands *.eps.ppm Fig-00*.jpg *eps-converted-to.pdf

#endif
