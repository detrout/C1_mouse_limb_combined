#!/bin/bash
cat <(python3 collapse_umis.py -c 7 -r 3 --header) \
    <(python3 collapse_umis.py -c 7 -r 4) \
    <(python3 collapse_umis.py -c 7 -r 5) \
    <(python3 collapse_umis.py -c 7 -r 6) \
    <(python3 collapse_umis.py -c 7 -r 7) \
    <(python3 collapse_umis.py -c 7 -r 8) \
    <(python3 collapse_umis.py -c 7 -r 12) \
    <(python3 collapse_umis.py -c 7 -r 13) | \
      samtools view -@ 4 -bS - -o c7_mus1_cb_ub_chr_start/10x-cluster-7_unsorted.bam
