#!/usr/bin/env python3

import sys
import os

import numpy

def main():
    sample_fname, repeat_fname, LINE_fname, sat_fname, segdup_fname, nonrep_fname, bfolder, afolder, data_dir, genome, out_fname = sys.argv[1:12]
    samples = load_samples(sample_fname)
    repeats = load_intervals(repeat_fname)
    LINEs = load_intervals(LINE_fname)
    Sats = load_intervals(sat_fname)
    segdups = load_intervals(segdup_fname)
    nonrep = load_intervals(nonrep_fname)
    data = load_data(samples, data_dir, genome, bfolder, afolder,
                     [repeats, LINEs, Sats, segdups, nonrep])
    write_data(data, out_fname)

def load_samples(fname):
    data = []
    with open(fname) as fs:
        for line in fs:
            if line.startswith('#'):
                continue
            line = line.rstrip().split(',')
            data.append([line[0], line[1], line[5], line[6].replace(' ', '_')])
    return data

def load_intervals(fname):
    data = []
    for line in open(fname):
        line = line.rstrip().split()
        chrom, start, end = line[:3]
        data.append((chrom, int(start), int(end)))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S30'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32)]))
    data = data[numpy.lexsort((data['start'], data['chrom']))]
    return data

def load_sizes(fname):
    total = 0
    chroms = []
    for line in open(fname):
        line = line.rstrip().split()
        total += int(line[1])
        chroms.append(line[0])
    chroms = set(chroms)
    print(total)
    return total, chroms

def load_data(samples, data_dir, genome, bfolder, afolder, intervals):
    celltypes = set()
    for expID, libID, target, celltype in samples:
        if target != 'Control':
            continue
        celltypes.add((celltype, expID))
    celltypes = list(celltypes)
    groups = ['repeat', 'LINE', 'satellite', 'segdup', 'nonrep']
    data = numpy.zeros(len(celltypes) * len(groups), dtype=numpy.dtype([
        ("Celltype", "S60"), ("Group", "S20"), ("Mapped", numpy.int64),
        ("Filtered", numpy.int64), ("Lost", numpy.int64),
        ("Lost_Percent", numpy.int64), ("Size", numpy.int64)]))
    for i, interval in enumerate(intervals):
        name = groups[i]
        size = numpy.sum(interval['end'] - interval['start'])
        for j, ct in enumerate(celltypes):
            index = j * len(groups) + i
            data['Celltype'][index] = ct[0]
            data['Group'][index] = name
            temp = numpy.load("{}/{}.{}_{}_{}_counts.npy".format(
                data_dir, ct[1], genome, bfolder, name))
            data['Mapped'][index] = numpy.sum(temp)
            temp = numpy.load("{}/{}.{}_{}_{}_counts.npy".format(
                data_dir, ct[1], genome, afolder, name))
            data['Filtered'][index] = numpy.sum(temp)
            data['Size'][index] = size
    data['Lost'] = data['Mapped'] - data['Filtered']
    for i, ct in enumerate(celltypes):
        s = i * len(groups)
        e = (i + 1) * len(groups)
        data['Lost_Percent'][s:e] = data['Lost'][s:e] / numpy.sum(data['Mapped'][s:e])
    return data

def write_data(data, out_fname):
    outfile = open(out_fname, 'w')
    print("Celltype\tGroup\tMapped\tFiltered\tLost\tLost_Percent\tSize", file=outfile)
    for i in range(data.shape[0]):
        ct, g, m, f, l, p, s = data[i]
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
              ct.decode('utf8'), g.decode('utf8'), m, f, l, p, s),
              file=outfile)
    outfile.close()

if __name__ == "__main__":
    main()
