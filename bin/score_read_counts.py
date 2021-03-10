#!/usr/bin/env python3

import sys
import os

import numpy

def main():
    sample_fname, array_fname, data_dir, out_fname = sys.argv[1:5]
    samples = load_samples(sample_fname)
    arrays = load_arrays(array_fname)
    data, totals = load_data(samples, data_dir)
    score_marks(data, totals, arrays, out_fname)

def load_samples(fname):
    data = []
    with open(fname) as fs:
        for line in fs:
            if line.startswith('#'):
                continue
            line = line.rstrip().split(',')
            data.append([line[0], line[1], line[5], line[6]])
    return data

def load_arrays(fname):
    data = []
    with open(fname) as fs:
        for line in fs:
            line = line.rstrip().split()
            data.append((line[3], line[3].split('_')[0]))
    data = numpy.array(data)
    return data

def load_data(samples, data_dir):
    data = {}
    controls = {}
    totals = {}
    celltypes = set()
    for expID, libID, target, celltype in samples:
        fname = "{}/{}.chm13v1_counts.npy".format(data_dir, expID)
        if not os.path.exists(fname):
            print(expID, libID, target, celltype)
            continue
        celltypes.add(celltype)
        key = (target, celltype)
        totals.setdefault(key, 0)
        totals[key] = int(open("bt2_dedup_kmer/{}.chm13v1.stats".format(expID)).readline().split()[0]) / 2
        temp = numpy.load(fname)
        #if target == "Control":
        #    if celltype not in controls:
        #        controls[celltype] = temp
        #else:
        data.setdefault(target, {})
        if celltype not in data[target]:
            data[target][celltype] = temp

    """
    for target in data.keys():
        for celltype in data[target].keys():
            scale = totals[('Control', celltype)] / totals[(target, celltype)]
            c = numpy.copy(controls[celltype])
            d = numpy.copy(data[target][celltype]).astype(numpy.float32)
            where = numpy.where((c == 0) | (d == 0))
            c[where] = 0
            d[where] = 0
            where = numpy.where((c > 0) | (d > 0))
            d[where] /= c[where]
            d[where] *= scale
            data[target][celltype] = d
    data['Control'] = {}
    for celltype, d in controls.items():
        data['Control'][celltype] = numpy.copy(d).astype(numpy.float32)
    """
    return data, totals

def score_marks(data, totals, arrays, out_fname):
    temp = numpy.unique(arrays[:, 1])
    array_names = []
    for array in temp:
        where = numpy.where(arrays[:, 1] == array)[0]
        if where.shape[0] >= 5:
            array_names.append(array)
    marks = list(data.keys())
    marks.remove('Control')
    marks.sort()
    outfile = open(out_fname, 'w')
    #print("Cell type\tMark\tArray\tScore", file=outfile)
    print("Cell type\tMark\tArray\tControl\tTreat", file=outfile)
    for h, array in enumerate(array_names):
        where = numpy.where(arrays == array)[0]
        names = arrays[where, 0]
        for mark in marks:
            for celltype in data[mark].keys():
                d = data[mark][celltype][where]
                c = data["Control"][celltype][where]
                nonzero = numpy.where(numpy.logical_and(d > 0, c > 0))[0]
                d = d[nonzero]
                c = c[nonzero]
                n = names[nonzero]
                #d = d[numpy.where(d > 0)]
                if d.shape[0] > 0:
                    #d = d[numpy.argsort(d)]
                    #d = numpy.log2(d)
                    #for val in d:
                    for i, val in enumerate(d):
                        print("{}\t{}\t{}\t{}\t{}".format(celltype, mark, n[i], c[i], val),
                              file=outfile)
                        #print("{}\t{}\t{}\t{}".format(celltype, mark, array, val),
                        #      file=outfile)
    celltypes = list(totals.keys())
    celltypes.sort()
    for mark, ct in celltypes:
        if mark == "Control":
            continue
        print("{}\t{}\t{}\t{}\t{}".format(ct, mark, "Total", int(totals[("Control", ct)]), int(totals[(mark, ct)])),
              file=outfile)
    outfile.close()

if __name__ == "__main__":
    main()
