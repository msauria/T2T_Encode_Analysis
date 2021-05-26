#!/usr/bin/env python

import sys

import pysam
import numpy

def main():
    bam_fname, array_fname, out_fname = sys.argv[1:4]
    bam = pysam.AlignmentFile(bam_fname, 'rb')
    arrays = load_arrays(array_fname, bam)
    scores = score_arrays(bam, arrays)
    numpy.save(out_fname, scores)

def load_arrays(fname, bam):
    chroms = set([x['SN'] for x in bam.header['SQ']])
    data = []
    for line in open(fname):
        line = line.rstrip().split()
        chrom, start, end = line[:3]
        if chrom not in chroms:
            continue
        data.append((chrom, int(start), int(end)))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S30'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32)]))
    data = data[numpy.lexsort((data['start'], data['chrom']))]
    return data

def score_arrays(bam, arrays):
    scores = numpy.zeros(arrays.shape[0], dtype=numpy.int32)
    for i in range(arrays.shape[0]):
        j = 0
        s = arrays['start'][i]
        e = arrays['end'][i]
        for read in bam.fetch(bytes(arrays['chrom'][i]), s, e):
            mid = read.reference_start + read.reference_length // 2
            if mid >= s and mid < e:
                j += 1
        scores[i] = j
    return scores

if __name__ == "__main__":
    main()
