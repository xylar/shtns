#!/bin/bash
# script to test many sht cases

id=`hg id`
log="test_suite.log"
echo "beginning test suite for $id" > $log

for switch in "" "-oop" "-transpose" "-schmidt" "-4pi"
do
  for mode in "-quickinit" "-gauss" "-reg" "-fly"
  do
    for lmax in 1 2 3 4 11 12 13 14 31 32 33 34 121 122 123 124
    do
      for mmax in 0 1 $lmax
      do
        c="./time_SHT $lmax -mmax=$mmax $mode $switch -iter=1"
        echo $c
        echo "---" >> $log
        echo "*** $c *** " >> $log
        $c > tmp.out
        cat tmp.out | grep ERROR
        cat tmp.out >> $log
      done
    done
  done
done

# do also a huge transform :
c="./time_SHT 2047 -mres=15 -quickinit -iter=1"
echo $c
echo "---" >> $log
echo "*** $c *** " >> $log
$c > tmp.out
cat tmp.out | grep ERROR
cat tmp.out >> $log

