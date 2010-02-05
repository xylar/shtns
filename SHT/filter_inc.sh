#!/bin/bash

inc_file='sht_generic.c'
inc_str="#include \"SHT/$inc_file\""

nl=`grep "$inc_str" $1 | wc -l`
if [ "$nl" == "0" ]; then
  cat $1
else
  sed $1 -e "s/#include .*$inc_file.*//" 

  ns=`grep "SUPARG" $1 | wc -l`
  if [ "$ns" == "0" ]; then
    sed SHT/$inc_file -e "s/SUPARG[^)]*//g" 
  else
    cat SHT/$inc_file
  fi
fi


