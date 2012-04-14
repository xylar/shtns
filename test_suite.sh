#!/bin/bash
# script to test many sht cases

id=`hg id`
log="test_suite.log"

function test1 {
    echo $1
    echo "---" >> $log
    echo "*** $1 *** " >> $log
    $1 > tmp.out
    cat tmp.out | grep ERROR
    cat tmp.out | grep -i nan
    cat tmp.out >> $log
}

echo "beginning test suite for $id" > $log

# first, do a huge transform :
test1 "./time_SHT 2047 -mres=15 -quickinit -iter=1"

# even bigger :
test1 "./time_SHT 7975 -mres=145 -quickinit -iter=1"

# without threads
test1 "./time_SHT 2047 -mres=15 -quickinit -iter=1 -nth=1"

for switch in "" "-oop" "-transpose" "-schmidt" "-4pi"
do
  for mode in "-quickinit" "-gauss" "-reg" "-fly" "-gauss -nth=1"
  do
    for lmax in 1 2 3 4 11 12 13 14 31 32 33 34 121 122 123 124
    do
      for mmax in 0 1 $lmax
      do
         test1 "./time_SHT $lmax -mmax=$mmax $mode $switch -iter=1"
      done
    done
  done
done
