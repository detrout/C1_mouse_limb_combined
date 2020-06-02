#!/bin/bash
RUN1=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-1-encode-count-cells10000
RUN3=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-3-encode-count-cells10000
RUN4=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-4-encode-count-cells10000
RUN5=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-5-encode-count-cells10000
RUN6=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-6-encode-count-cells10000
RUN7=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-7-encode-count-cells10000
RUN8=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-8-encode-count-cells10000
RUN9=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-9-encode-count-cells10000
RUN10=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-10-encode-count-cells10000
RUN11=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-11-encode-count-cells10000
RUN12=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-12-encode-count-cells10000
RUN13=/woldlab/loxcyc/home/diane/proj/brian-2018-01-10x/Wold10x-13-encode-count-cells10000

cat <(samtools view -H $RUN1/outs/possorted_genome_bam.bam) \
    <(samtools view $RUN1/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-1.txt) \
    <(samtools view $RUN4/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-4.txt) \
    <(samtools view $RUN5/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-5.txt) \
    <(samtools view $RUN6/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-6.txt) \
    <(samtools view $RUN7/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-7.txt) \
    <(samtools view $RUN8/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-8.txt) \
    <(samtools view $RUN12/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-12.txt) \
    <(samtools view $RUN13/outs/possorted_genome_bam.bam | grep -f barcodes-10x-cluster-17-run-13.txt) | \
      samtools view -@ 4 -bS - -o 10x-cluster-17_unsorted.bam
