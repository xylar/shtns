#!/bin/bash
# this script runs xhells, logs output, store par file, and send e-mail when finished.

xs="./xshells"
xspar="xshells.par"
xsh="xshells.h"
log_file="xs.log"
log_dir="$HOME/xshells_log"

# if command line parameter is a file, then use this as programm.
if [ -f "$1" ]
then
  xs=$1
  i1=$2
  i2=$3
else
  i1=$1
  i2=$2
fi

md5file=/tmp/parfile.tmp
echo $xs > $md5file
sed $xspar -e "/^#/d" -e "/^$/d" >> $md5file
cat $xsh >> $md5file

id=`md5sum $md5file | sed -e "s/\([0-9a-f]\{3\}\).*\([0-9a-f]\{3\}\) .*/\1\2/"`
job=`cat $xspar |grep job| sed -e "s/job[ =]*//" -e "s/[\t ]*#.*//"`

# make : rebuild if parameters have changed.
make `basename $xs`

echo "" >> $log_dir/$log_file
echo "******************************************************" >> $log_file
echo "*** STARTING JOB $job :: id=$id :: $(date) ***" >> $log_file
# get header
$xs /dev/nulll | sed -e "/load/d" >> $log_file
echo "******************************************************" >> $log_file

#save .par file and .h file
cp $xspar $log_dir/par.$job.$id
cat $xsh >> $log_dir/par.$job.$id

#run job
echo "running $xs :: job $job :: id=$id"
echo "    logging to $log_dir/out.$job.$id"
$xs $i1 $i2 > $log_dir/out.$job.$id
if [ $? -eq 1 ]
then
  status='runtime error'
else
  status='all done'
fi
echo $status

echo "**** END OF JOB $job :: id=$id :: $(date) :: $status ****" >> $log_file

# send mail
send_mail_xshells "end of job $job.$id : $status" $log_dir/par.$job.$id $log_dir/out.$job.$id | telnet mailhost.ujf-grenoble.fr 25
