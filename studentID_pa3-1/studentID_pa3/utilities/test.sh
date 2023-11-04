#!/bin/bash

type="public"
# type="hidden"
for i in 1 2 3 4 7 8
do
    start=`date +%s.%N`

    echo "testing case $type-${i} ... "
    ./bin/cb ./inputs/$type\_case_$i.in ./outputs/$type\_case_$i.out

    end=`date +%s.%N`
    echo "$end - $start" | bc -l 
    ./pa3_checker ./inputs/$type\_case_$i.in ./outputs/$type\_case_$i.out
    
done