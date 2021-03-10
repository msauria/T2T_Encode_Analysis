#!/usr/bin/env python

import sys

import pysam
import numpy

def main():
    bam_fname, array_fname, out_fname = sys.argv[1:4]
    arrays = load_arrays(array_fname)
    scores = score_arrays(bam_fname, arrays)
    numpy.save(out_fname, scores)

def load_arrays(fname):
    data = []
    for line in open(fname):
        chrom, start, end, array = line.rstrip().split()[:4]
        data.append((chrom, int(start), int(end), array))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S6'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32),
                                                ('array', 'S20')]))
    return data

def score_arrays(fname, arrays):
    scores = numpy.zeros(arrays.shape[0], dtype=numpy.int32)
    bam = pysam.AlignmentFile(fname, 'rb')
    for i in range(arrays.shape[0]):
        j = 0
        for read in bam.fetch(bytes(arrays['chrom'][i]), arrays['start'][i], arrays['end'][i]):
            j += 1
        scores[i] = j
    return scores

if __name__ == "__main__":
    main()
