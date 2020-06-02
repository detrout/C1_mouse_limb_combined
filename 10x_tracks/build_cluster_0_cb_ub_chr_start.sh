#!/bin/bash
cat <(python3 collapse_umis.py -c 0 -r 1 --header) \
    <(python3 collapse_umis.py -c 0 -r 3) \
    <(python3 collapse_umis.py -c 0 -r 4) \
    <(python3 collapse_umis.py -c 0 -r 5) \
    <(python3 collapse_umis.py -c 0 -r 6) \
    <(python3 collapse_umis.py -c 0 -r 7) \
    <(python3 collapse_umis.py -c 0 -r 8) \
    <(python3 collapse_umis.py -c 0 -r 12) \
    <(python3 collapse_umis.py -c 0 -r 13) | \
      samtools view -@ 4 -bS - -o c0_mesprox_cb_ub_chr_start/10x-cluster-0_unsorted.bam
