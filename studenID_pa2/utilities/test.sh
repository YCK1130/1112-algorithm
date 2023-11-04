for i in 12 1000 10000 100000; do
    echo "testing case ${i} ... "
    start=`date +%s.%N`
    ./bin/mps inputs/${i}.in outputs/up${i}.out 
    end=`date +%s.%N`
    echo "$end - $start" | bc -l 
done
