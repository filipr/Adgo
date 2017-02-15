#!/bin/csh

if($argv[9] == 'A' || $argv[9] == 'C') then
 cat << EOF >> tisk1.gnu

 set size 0.4, 0.5

 set output 'u-cut-$Cname.eps'
#p [] [:] '$dir/gnu.03' t'cut' w l, exp(-100*((x-0.75)**2 + (x-0.75)**2 )) t 'exact'
#p [] [:] '$dir/gnu.03' t'cut' w l, (0.1 + exp(5)) / (0.1 + exp(10 ))*x*x*(1-x)**2 t 'exact'
p [] [:] '$dir/gnu.03' t'cut' w l, (0.1 + exp(5)) *x*x*(1-x)**2 t 'exact'

EOF

else if( $argv[9] == 'D') then
 cat << EOF >> tisk1.gnu

 set size 0.4, 0.5

 set output 'u-cut-$Cname.eps'
#p [] [:] '$dir/gnu.03' t'cut' w l, exp(-100*((x-0.75)**2 + (x-0.75)**2 )) t 'exact'
#p [] [:] '$dir/gnu.03' t'cut' w l, (0.1 + exp(5)) / (0.1 + exp(10 ))*x*x*(1-x)**2 t 'exact'
p [] [:] '$dir/gnu.03' t'cut' w l, (1.0+ 0.5)  *x*x*(1-x)**2 t 'exact'

EOF
else if( $argv[9] == 'F') then
 cat << EOF >> tisk1.gnu

 set size 0.4, 0.5

 set output 'u-cut-$Cname.eps'
#p [] [:] '$dir/gnu.03' t'cut' w l, exp(-100*((x-0.75)**2 + (x-0.75)**2 )) t 'exact'
#p [] [:] '$dir/gnu.03' t'cut' w l, (0.1 + exp(5)) / (0.1 + exp(10 ))*x*x*(1-x)**2 t 'exact'
p [] [:] '$dir/gnu.03' t'cut' w l, exp(1+ 2*x) t 'exact'

EOF
else if( $argv[9] == 'H') then
 cat << EOF >> tisk1.gnu

 set size 0.4, 0.5

 set output 'u-cut-$Cname.eps'
#p [] [:] '$dir/gnu.03' t'cut' w l, exp(-100*((x-0.75)**2 + (x-0.75)**2 )) t 'exact'
#p [] [:] '$dir/gnu.03' t'cut' w l, (0.1 + exp(5)) / (0.1 + exp(10 ))*x*x*(1-x)**2 t 'exact'
p [] [:] '$dir/gnu.03' t'cut' w l,  2*x*x*(1-x)*(1-x)*(2*x*x)**($paramet/2)  t'exact'lw 4.0 w l

EOF
else if($argv[9] == 'B') then

 cat << EOF >> tisk1.gnu
 set size 0.4, 0.5

 set output 'u-cut-$Cname.eps'
p [] [:] '$dir/gnu.03' t'cut' w l, exp(-100*((x-0.75)**2 + (x-0.75)**2 )) t 'exact'
EOF

endif
