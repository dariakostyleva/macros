#!/bin/bash

for (( i=55 ; i<=61 ; i++))
  do
  	root -l -q 'read_mom_widths_real.C("mom_widths_real.txt",'$i')'
  done