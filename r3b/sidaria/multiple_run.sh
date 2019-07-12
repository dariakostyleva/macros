#!/bin/bash
#  --------------------------------------------------------------------------------------------------------
#  Author: D. Kostyleva
#  This script runs MC simulation and then analysis script many times for different number of events(counts)
#  After each run the analysis script writes calculated variable (width of momentum distribution)
#  into some file
#  That's it
#  ---------------------------------------------------------------------------------------------------------

echo 'Lets begin!'
# max number of events up to which we want to go
maxcount=10
# for each number of events we run runnunb of times
runnumb=9
# below I run nested loop: for each number of counts the simulation and then analysis script are run runnumb times

for k in $(seq 2 $maxcount)
  do
  for i in $(seq 0 $runnumb)
    do 
  	  echo 'Iteration '$i
      root -l <<EOF
      .L runsim.C
      runsim($k)  
EOF
      root -l -q 'ana_width.C('$i','$runnumb','$k')'
    done
  done

