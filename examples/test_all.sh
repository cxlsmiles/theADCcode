#!/bin/bash
#Loop over all folders to make the tests

echo "Running tests over all example inputs:"
current_dir=`pwd`
for d in `ls -l | grep '^d' | awk '{print $9}'`; do
  echo "Testing $d" |  tr '[:lower:]' '[:upper:]'
  cd $d
  ./test.pl
  cd $current_dir
done
