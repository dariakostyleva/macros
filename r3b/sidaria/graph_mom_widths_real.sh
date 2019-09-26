#!/bin/bash

for (( i=40 ; i<=70 ; i++))
  do
  	root -l -q 'read_mom_widths_real.C("mom_widths_real.txt",'$i')'
  done