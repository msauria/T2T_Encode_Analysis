#!/usr/bin/env python3

import sys
import multiprocessing
import time

import pyBigWig
import numpy


def main():
    kmer_fname, dedup_fname, diff_fname, threads, binsize = sys.argv[1:6]
    binsize, threads = int(binsize), int(threads)
    bw_fnames1 = sys.argv[6].strip(',').split(',') # kmer
    bw_fnames2 = sys.argv[7].strip(',').split(',') # dedup
    bw1 = []
    bw2 = []
    chroms = None
    scale = numpy.zeros(2, dtype=numpy.float64)
    for fname in bw_fnames1:
        bw1.append(pyBigWig.open(fname))
    for fname in bw_fnames2:
        bw2.append(pyBigWig.open(fname))
        h = bw2[-1].header()
        scale[0] += h['nBasesCovered']
        scale[1] += h['sumData']
    scale = 5 / (scale[1] / scale[0]) / len(bw_fnames2)
    kmer_outfile = open(kmer_fname, 'w')
    dedup_outfile = open(dedup_fname, 'w')
    diff_outfile = open(diff_fname, 'w')
    chroms = list(bw1[0].chroms().items())
    chroms.sort()
    n = len(chroms)
    #for i in range(n)[::-1]:
    #    chrom, maxlen = chroms[i]
    #    if not (chrom.lstrip('chr').isdigit() or chrom == "chrX"):
    #        chroms.pop(i)
    for chrom, maxlen in chroms:
        bins = maxlen // binsize
        processes = []
        queue = multiprocessing.JoinableQueue()
        out_queue = multiprocessing.JoinableQueue()
        for i in range(threads):
            processes.append(multiprocessing.Process(
                target=load_data, args=(bw1, bw2, chrom, bins, binsize, queue, out_queue)))
            processes[-1].daemon = True
            processes[-1].start()
        for i, bw in enumerate(bw1):
            queue.put((0, i))
        for i, bw in enumerate(bw2):
            queue.put((1, i))
        for i in range(threads):
            queue.put(None)
        results = [numpy.zeros(bins, dtype=numpy.float64),
                   numpy.zeros(bins, dtype=numpy.float64)]
        finished = 0
        while finished < threads:
            temp = out_queue.get(True)
            if temp is None:
                finished += 1
            else:
                index, result = temp
                results[index] += result['data']
        prev = False
        data1, data2 = results
        print("{} {}".format(chrom, numpy.where(data2 > 0)[0].shape[0]), file=sys.stderr)
        for i in range(len(data1)):
            if data2[i] > 0:
                if not prev:
                    print("fixedStep chrom={} start={} step={}".format(
                          chrom, 1+i*binsize, binsize), file=kmer_outfile);
                    print("fixedStep chrom={} start={} step={}".format(
                          chrom, 1+i*binsize, binsize), file=dedup_outfile);
                    print("fixedStep chrom={} start={} step={}".format(
                          chrom, 1+i*binsize, binsize), file=diff_outfile);
                print(data1[i], file=kmer_outfile)
                print(data2[i], file=dedup_outfile)
                print(data2[i] - data1[i], file=diff_outfile)
                prev = True
            else:
                prev = False
    kmer_outfile.close()
    dedup_outfile.close()
    diff_outfile.close()
    for bw in bw1:
        bw.close()
    for bw in bw2:
        bw.close()

def load_data(bw1, bw2, chrom, bins, binsize, queue, out_queue):
    temp = queue.get(True)
    while temp is not None:
        index, index1 = temp
        if bins == 0:
            out_queue.put((index, {'data': numpy.zeros(0, dtype=numpy.float64)}))
        else:
            if index == 0:
                bw = bw1[index1]
            else:
                bw = bw2[index1]
            try:
                values = bw.values(chrom, 0, bins*binsize)
            except:
                out_queue.put((index, {'data': numpy.zeros(bins, dtype=numpy.float64)}))
            else:
                values = numpy.array(values)
                values[numpy.where(numpy.isnan(values))] = 0
                data = numpy.bincount(numpy.repeat(numpy.arange(bins), binsize) ,
                                      weights=values, minlength=bins)
                out_queue.put((index, {'data': data}))
        temp = queue.get(True)
    out_queue.put(None)
    return

if __name__ == "__main__":
    main()
