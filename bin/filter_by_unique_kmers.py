#!/usr/bin/env python3

import sys

import pysam
import pyBigWig
import numpy

def main():
    bam_fname, kmer_template, kmer_sizes, out_fname = sys.argv[1:5]
    bam = pysam.AlignmentFile(bam_fname, 'rb')
    sizes = [int(x) for x in kmer_sizes.split(',')]
    sizes.sort()
    sizes = numpy.array(sizes, dtype=numpy.int32)
    sizes2 = sizes[1:-1]
    DBs = {}
    for size in sizes:
        DBs[size] = pyBigWig.open(kmer_template.replace("*", str(size)))
    output = pysam.AlignmentFile(out_fname, 'wb', template=bam)
    bam_iter = bam.__iter__()
    try:
        read1 = next(bam_iter)
        while True:
            read2 = next(bam_iter)
            chrom = read1.reference_name
            if (valid_read(read1, chrom, DBs, sizes, sizes2)
                or valid_read(read2, chrom, DBs, sizes, sizes2)):
                output.write(read1)
                output.write(read2)
            read1 = next(bam_iter)
    except StopIteration:
        pass
    finally:
        del bam_iter
    output.close()
    bam.close()

def valid_read(read, chrom, DBs, sizes, sizes2):
    start = read.reference_start
    end = read.reference_end
    span = end - start
    if span < sizes[0]:
        return False
    index = numpy.searchsorted(sizes2, span, side='left')
    n = span - sizes[index] + 1
    end = start + n
    if start >= end:
        return False
    try:
        scores = DBs[sizes[index]].values(chrom, start, end)
    except:
        return False
    if 1 in scores:
        return True
    else:
        return False


if __name__ == "__main__":
    main()
