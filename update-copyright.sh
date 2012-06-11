#!/bin/bash
# this will update the commented copyright notice to the file given in command line

# temp file
tmp="/tmp/updcpy.tmp"

# this is the comment string
cmtstart='/*'
cmt=' * '
cmtend=' */'

if [ "$1" == "-fortran" ]; then
  cmtstart='!'
  cmt='!  '
  cmtend=''
  shift
fi

if [ "$1" == "-python" ]; then
  cmtstart="#!/usr/bin/python\n#"
  cmt='#  '
  cmtend=''
  shift
fi

if [ -f "$1" ]; then

  # is there already a copyright header ?
  n=`head -5 $1 |grep -i Copyright |wc -l`

  # display the copyright notice within a comment block
  echo -e "$cmtstart" > $tmp
  sed "s/^/$cmt/" COPYRIGHT >> $tmp
  if [ "x$cmtend" != "x" ]; then
	  echo "$cmtend" >> $tmp
  fi
  # do not forget the blank line !
  echo "" >> $tmp

  if [ "$n" != "0" ]; then
	# remove the existing copyright notice and replace it with the new one.
	sed '1,/^$/ d' $1 >> $tmp
  else
	# display the file.
	cat $1 >> $tmp
  fi

  mv $tmp $1

else
  echo "file not found"
fi
