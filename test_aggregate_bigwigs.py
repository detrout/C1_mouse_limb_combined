from unittest import TestCase
from unittest.mock import Mock
from collections import Counter
import os
from subprocess import check_call

from make_aggregate_bigwigs import (
    get_mode,
    get_chrom_info,
    sum_alignments,
    write_contig_to_bigwig,
)


# Create bam files with
class TestAggregateBigwigs(TestCase):
    def setUp(self):
        samfiles = ['simple.sam', 'simple-no-overlap.sam']
        for samfile in samfiles:
            base, ext = os.path.splitext(samfile)
            bamfile = base + '.bam'
            baifile = bamfile + '.bai'

            if not os.path.exists(bamfile):
                check_call(['samtools', 'view', '-t', 'bam', '-o', bamfile, samfile])
            if not os.path.exists(baifile):
                check_call(['samtools', 'index', bamfile])
                
    def test_get_mode(self):
        self.assertEqual(get_mode('a.bam', 'rb'), 'rb')
        self.assertEqual(get_mode('a.bam', 'wb'), 'wb')
        self.assertEqual(get_mode('a.bam', 'rt'), 'rb')
        self.assertEqual(get_mode('a.bam', 'wt'), 'wb')
        self.assertEqual(get_mode('a.bam', 'w'), 'wb')
        self.assertEqual(get_mode('a.bam', 'wbu'), 'wbu')

        self.assertEqual(get_mode('a.sam', 'wb'), 'w')
        self.assertEqual(get_mode('a.sam', 'rb'), 'r')
        self.assertEqual(get_mode('a.sam', 'w'), 'w')
        
    def test_get_chrom_info(self):
        chrom_info = get_chrom_info('simple.sam')
        self.assertEqual(len(chrom_info), 1)
        self.assertEqual(chrom_info['chr1'], 1000)        

    def test_sum_alignments_two_reads(self):
        alignment_files = ['simple-no-overlap.bam']
        chrom_info = get_chrom_info(alignment_files[0])
        current = sum_alignments(alignment_files, chrom_info, 'chr1')

        for i in range(0, 99):
            self.assertEqual(current[i], 0)
        for i in range(99, 99+50): 
            self.assertEqual(current[i], 1)
        for i in range(99+50, 499): 
            self.assertEqual(current[i], 0)
        for i in range(499, 499+50):
            self.assertEqual(current[i], 1)

    def test_sum_alignments(self):
        for alignment_files in [ ['simple.bam'],
                                 ['simple.bam', 'simple.bam'],
                                 ['simple.bam', 'simple.bam', 'simple.bam']]:
            chrom_info = get_chrom_info(alignment_files[0])
            current = sum_alignments(alignment_files, chrom_info, 'chr1')

            for i in range(0, 99):
                self.assertEqual(current[i], 0 * len(alignment_files), i)
            for i in range(99, 124): 
                self.assertEqual(current[i], 1 * len(alignment_files), i)
            for i in range(124, 164): 
                self.assertEqual(current[i], 2 * len(alignment_files), i)
            for i in range(164, 199):
                self.assertEqual(current[i], 1 * len(alignment_files), i)
            for i in range(199, 499):
                self.assertEqual(current[i], 0 * len(alignment_files), i)
            for i in range(499, 499+25):
                self.assertEqual(current[i], 1 * len(alignment_files), i)
            for i in range(499+25, 499+100):
                self.assertEqual(current[i], 2 * len(alignment_files), i)
            for i in range(499+100, 499+125):
                self.assertEqual(current[i], 1 * len(alignment_files), i)
            for i in range(499+125, chrom_info['chr1']):
                self.assertEqual(current[i], 0 * len(alignment_files), i)

    def test_write_wiggle(self):
        current = Counter({
            500: 5,
            501: 5,
            502: 5,
            503: 5,
            504: 5,
            505: 10,
            506: 10,
            507: 10,
            508: 10,
            509: 10,
            600: 1,
            601: 1,
            602: 1,
            603: 1,
            604: 1,
            700: 2,
            701: 2,
            702: 2,
        })
        chrom_info = {'chr1': 1000}
        bigwig = Mock()
        write_contig_to_bigwig(bigwig, 'chr1', current)
        bigwig.addEntries.assert_called_with(
            ['chr1'] * 4,
            [500, 505, 600, 700],
            ends=[505, 510, 605, 703],
            values=[5, 10, 1, 2])
        
