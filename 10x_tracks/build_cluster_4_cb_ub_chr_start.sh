#!/bin/bash
cat <(python3 collapse_umis.py -c 4 -r 1 --header) \
    <(python3 collapse_umis.py -c 4 -r 3) \
    <(python3 collapse_umis.py -c 4 -r 4) \
    <(python3 collapse_umis.py -c 4 -r 5) \
    <(python3 collapse_umis.py -c 4 -r 6) \
    <(python3 collapse_umis.py -c 4 -r 8) \
    <(python3 collapse_umis.py -c 4 -r 12) \
    <(python3 collapse_umis.py -c 4 -r 13) | \
      samtools view -@ 4 -bS - -o c4_mus2_cb_ub_chr_start/10x-cluster-4_unsorted.bam
