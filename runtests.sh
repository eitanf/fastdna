#!/usr/bin/env bash

set -e

for UNROLL in 1 2
do
  rm -f fastdna fastdna.o
  make fastdna OPTFLAGS="-O0 -g -DUNROLL=$UNROLL"

  for FN in tests/A tests/AA tests/AC tests/CGAT tests/GATTACCA tests/TGCAGCATTACGCAGT tests/GCGTACGTACGTACGTACGT tests/AeCa # DNACorpus2/HoSa
  do
    for OPTS in "-t1" "-t4" "-t1 -c8" "-a" "-t1 -a" "-t1 -p" "-t1 -p -u" "-t1 -c8 -p -u" "-t9 -c8 -p -u" "-t1 -c1032 -p -u" "-t15 -c65536 -p" "-t14 -c2147483648 -p"
    do
      echo -n "Testing UNROLL=$UNROLL $OPTS on $FN"
      ./fastdna $OPTS e $FN > /dev/null
      echo -n "."
      ./fastdna $OPTS d $FN.2b > /dev/null
      echo -n "."
      cmp $FN $FN.2b.dec
      echo "."
    done
  done
done

