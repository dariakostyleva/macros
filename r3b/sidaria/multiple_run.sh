#!/bin/bash
# run as bash multiple_run.sh stacount maxcount runnumb
#  --------------------------------------------------------------------------------------------------------
#  Author: D. Kostyleva
#  This script runs MC simulation and then analysis script many times for different number of events(counts)
#  After each run the analysis script writes calculated variable (width of momentum distribution)
#  into some file
#  That's it
#  ---------------------------------------------------------------------------------------------------------
start=$SECONDS
# VERY IMPORTANT to delete .txt file here, because otherwise it will continue writing in it
rm mom_widths_real.txt
# number of events from which we want to start
stacount=$1
# max number of events up to which we want to go
maxcount=$2
# for each number of events we run runnumb of times
runnumb=$3
echo 'Lets begin! Start from' $stacount 'up to' $maxcount 'events'
echo 'Run sim' $runnumb 'times for each number of events'
sleep 4

if [ $stacount -lt 2 ] || [ $stacount -gt $maxcount ]
then
	echo 'Incorrect number of start counts' $stacount
fi


# below I run nested loop: for each number of counts the simulation and then analysis script are run runnumb times
echo '    '
#for (( k=$stacount ; k<=$maxcount ; k++ ))
#I want to go with larger step 
for (( k=$stacount ; k<=$maxcount ; k++ ))
  do
  for (( i=0; i<$runnumb ; i++ ))
    do 
  	  echo 'Iteration' $i 'for' $k 'counts'
      root -l <<EOF
      .L runsim_real.C
      runsim_real($k)  
EOF
	  echo '    '
      echo 'Input to script: iteration' $i 'for' $k 'counts' 
      root -l -q 'ana_width_real.C('$i','$k','$runnumb')'
    done
  echo '    '
  echo $runnumb 'iterations over' $k 'counts are over!'
  echo '    '
  done
echo 'Output in mom_widths.txt!'
duration=$(( SECONDS - start ))
echo 'It took' $duration 'sec or' $(($duration/3600)) 'hours'

#syntax to comment whole block
: <<'END'
END
