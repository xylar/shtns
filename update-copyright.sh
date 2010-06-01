#!/bin/bash
# this will update the commented copyright notice to the file given in command line

# temp file
tmp="/tmp/updcpy.tmp"

# this is the comment string
cmtstart='/*'
cmt=' * '
cmtend=' */'

if [ "$1" == "-fortran" ]; then
  cmtstart='c'
  cmt='c  '
  cmtend=''
  shift
fi


if [ -f "$1" ]; then

  # is there already a copyright header ?
  n=`head -3 $1 |grep opyright |wc -l`

  # display the copyright notice within a comment block
  echo "$cmtstart" > $tmp
  sed "s/^/$cmt/" copyright >> $tmp
  echo "$cmtend" >> $tmp
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
