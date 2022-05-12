#!/usr/bin/env python3
# coding: utf-8

import argparse
import logging
import re
import sys

import mappy as mp


def count_bases(rd1_file, rd2_file, map_read, five_p, three_p, map_pattern):
    primer_hit,primer_nohit,unmap,rd_n = 0,0,0,0
    edit_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0, }

    for rd1, rd2 in zip(mp.fastx_read(rd1_file), mp.fastx_read(rd2_file)):
        rd_n += 1
        if rd_n % 100000 == 0:
            logger.info(f'Processed reads: {rd_n:,}, No primer: {primer_nohit:,}, Unmapped: {unmap:,}')
        rd1, rd2 = rd1[1], rd2[1]
        if rd1.startswith(five_p) and rd2.startswith(three_p):
            primer_hit += 1
            if map_read == 5:
                edit = re.findall(map_pattern, rd1)
            else:
                edit = re.findall(map_pattern, mp.revcomp(rd2))
            if edit:
                edit_count[edit[0]] = edit_count[edit[0]] + 1
            else:
                unmap += 1
        elif rd1.startswith(three_p) and rd2.startswith(five_p):
            primer_hit += 1
            if map_read == 5:
                edit = re.findall(map_pattern, rd2)
            else:
                edit = re.findall(map_pattern, mp.revcomp(rd1))
            if edit:
                edit_count[edit[0]] = edit_count[edit[0]] + 1
            else:
                unmap += 1
        else:
            primer_nohit += 1
    logger.info('---------------------------------------')
    logger.info(f'Processed reads: {rd_n:,}, No primer: {primer_nohit:,} ({primer_nohit/rd_n:.2%}), Unmapped: {unmap:,} ({unmap/rd_n:.2%})')
    logger.info(f'Base count of editing site:')
    for base,n in edit_count.items():
        logger.info(f'{base}\t{n} ({n/sum(edit_count.values()):.2%})')


def main():
    parser = argparse.ArgumentParser(usage="Count single base editing",
                                     description="python PEcounter.py -p 62 ./ref.fa ./test_1.fq.gz ./test_2.fq.gz",)
    parser.add_argument('ref', type=str,
                        help='fasta file of reference sequence'),
    parser.add_argument('fq1', type=str, help='fastq file of read1'),
    parser.add_argument('fq2', type=str, help='fastq file of read2'),
    parser.add_argument('-p', '--edit_pos', type=int, required=True,
                        help='the position of edit site in reference sequence, 1-based', metavar=''),
    args = parser.parse_args()

    pos = args.edit_pos
    read_len = 150
    map_len = 10

    rd1_file = args.fq1
    rd2_file = args.fq2

    ref_fa = args.ref
    for ref in mp.fastx_read(ref_fa):
        ref = ref[1]
        break

    if pos < read_len-5:
        map_read = 5
    elif pos > (len(ref)-read_len+5):
        map_read = 3
    else:
        logger.error('Reads did not contain the editing site.')
        sys.exit()

    map_pattern = f'{ref[pos-6:pos-1]}([ACGT]){ref[pos:pos+5]}'
    five_p, three_p = ref[:map_len], mp.revcomp(ref[-map_len:])
    count_bases(rd1_file, rd2_file, map_read, five_p, three_p, map_pattern)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,
                        datefmt='%Y/%m/%d %H:%M:%S',
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    main()
