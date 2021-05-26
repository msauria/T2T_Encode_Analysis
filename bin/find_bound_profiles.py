#!/usr/bin/env python3

import sys
import os
import multiprocessing
import time

import pyBigWig
import numpy

def main():
    sample_fname, array_fname, size_fname, bw_dir, minsize, window, binsize, threads, out_fname = sys.argv[1:10]
    minsize, window, binsize, threads = int(minsize), int(window), int(binsize), int(threads)
    N0 = minsize // binsize
    N1 = window //  binsize
    samples = load_samples(sample_fname)
    bounds, anames, chroms = load_arrays2(array_fname, size_fname, minsize, window)
    marks = list(set([x[2] for x in samples]))

    print(marks)
    print(anames)
    M = numpy.amax(bounds[:, -1]) + 1
    queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.JoinableQueue()
    processes = []
    for i in range(threads):
        processes.append(multiprocessing.Process(
            target=find_profiles, args=(bw_dir, bounds, chroms, N0, N1, binsize, queue, out_queue)))
        processes[-1].daemon = True
        processes[-1].start()
    for i in range(len(samples)):
        if os.path.exists("{}/{}_{}.chm13v1.bw".format(bw_dir, samples[i][3], samples[i][2])):
            queue.put((samples[i][3], samples[i][2]))
    for i in range(threads):
        queue.put(None)
    results = {}
    for mark in marks:
        results[mark] = numpy.zeros((N0+N1, M, 2), dtype=numpy.float64)
    finished = 0
    while finished < threads:
        temp = out_queue.get(True)
        if temp is None:
            finished += 1
        else:
            results[temp[1]] += temp[2]
    output = open(out_fname, 'w')
    temp = []
    for m in marks:
        results[m] = results[m][:, :, 0] / numpy.maximum(1, results[m][:, :, 1])
        for a in anames:
            temp.append("{}_{}".format(m, a))
    print(",".join(temp), file=output)
    for i in range(N0+N1):
        temp = []
        for m in marks:
            temp += [str(x) for x in results[m][i, :]]
        print(",".join(temp), file=output)
    temp = [str(x) for x in numpy.bincount(bounds[:, -1])]
    print(",".join(temp), file=output)
    output.close()

def find_profiles(bw_dir, bounds, chroms, n0, n1, binsize, q, outq):
    temp = q.get(True)
    m = numpy.amax(bounds[:, -1]) + 1
    while temp is not None:
        sample, target = temp
        bw = pyBigWig.open("{}/{}_{}.chm13v1.bw".format(bw_dir, sample, target))
        scale = 1000000000. / bw.header()['sumData']
        results = numpy.zeros((n0+n1, m, 2), dtype=numpy.float64)
        for i in range(bounds.shape[0]):
            cint, start, end, window, orientation, index = bounds[i, :]
            values = numpy.array(bw.values(chroms[cint], start, end))
            bins = numpy.floor(numpy.arange(0, end - start) / (end - start) *
                               n0).astype(numpy.int32)
            valid = numpy.where(numpy.logical_not(numpy.isnan(values)))[0]
            bins = bins[valid]
            sums = numpy.bincount(bins, weights=values[valid], minlength=n0)
            valid = numpy.bincount(bins, minlength=n0)
            if orientation == -1:
                sums = sums[::-1]
                valid = valid[::-1]
            results[:n0, index, 0] += sums
            results[:n0, index, 1] += valid
            if orientation == 1:
                bins = numpy.arange(window - end)
                values = numpy.array(bw.values(chroms[cint], end, window))
                s = end
                e = window
            else:
                bins = numpy.arange(start - window)[::-1]
                values = numpy.array(bw.values(chroms[cint], window, start))[::-1]
                s = window
                e = start
            valid = numpy.where(numpy.logical_not(numpy.isnan(values)))[0]
            bins = bins[valid] // binsize
            sums = numpy.bincount(bins, weights=values[valid], minlength=n1)
            valid = numpy.bincount(bins, minlength=n1)
            results[n0:, index, 0] += sums
            results[n0:, index, 1] += valid
        results[:, :, 0] *= scale
        outq.put((sample, target, results))
        temp = q.get(True)
    outq.put(None)

def load_samples(fname):
    data = []
    with open(fname) as fs:
        for line in fs:
            if line.startswith('#'):
                continue
            line = line.rstrip().split(',')
            if line[5] not in ('H3K9me3', 'Control'):
                continue
            data.append([line[0], line[1], line[5], line[6].replace(' ', '_')])
    return data


def load_arrays(fname, size_fname, minsize, window):
    sizes = {}
    for line in open(size_fname):
        line = line.rstrip().split()
        sizes[line[0].encode('utf8')] = int(line[1])
    data = []
    for line in open(fname):
        line = line.rstrip().split()
        chrom, start, end, name = line[:4]
        name = name.split('_')[0]
        if name == 'ct' or name == 'censat':
            continue
        data.append((chrom, int(start), int(end), int(end) - int(start), name))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S6'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32),
                                                ('size', numpy.int32),
                                                ('name', 'S40')]))
    chroms = numpy.unique(data['chrom'])
    remove = []
    for chrom in chroms:
        where = numpy.where(data['chrom'] == chrom)[0]
        mid = (data['start'][where[0]] + data['end'][where[-1]]) // 2
        where1 = where[numpy.where((data['start'][where[1:]] ==
                                    data['end'][where[:-1]]) &
                                   (data['name'][where[1:]] ==
                                    data['name'][where[:-1]]))[0]]# <
                                   #window // 4)[0]]
        for i in where1[::-1]:
            data['end'][i] = data['end'][i + 1]
            #if data['end'][i] > mid:
            #    data['name'][i] = data['name'][i + 1]
            remove.append(i + 1)
    data = numpy.delete(data, remove)
    data['size'] = data['end'] - data['start']
    cint = {}
    for i, c in enumerate(chroms):
        cint[c] = i
    bounds = []
    anames = {}
    for chrom in chroms:
        cwhere = numpy.where(data['chrom'] == chrom)[0]
        where = numpy.where((data['chrom'] == chrom) &
                            (data['size'] >= minsize))[0]
        if data['start'][where[0]] >= window // 2:
            w = max(data['start'][where[0]] - window, 0)
            start = data['start'][where[0]]
            end = data['end'][where[0]]
            name = data['name'][where[0]]
            anames.setdefault(name, len(anames))
            bounds.append((cint[chrom], start, end, w, -1, anames[name]))
        start = data['start'][where[-1]]
        end = data['end'][where[-1]]
        #if where[-1] < cwhere[-1]:
        #    w = min(data['end'][where[-1]] + window, data['start'][where[-1] + 1])
        #else:
        w = data['end'][where[-1]] + window
        name = data['name'][where[-1]]
        anames.setdefault(name, len(anames))
        bounds.append((cint[chrom], start, end, w, 1, anames[name]))
    bounds = numpy.array(bounds, dtype=numpy.int32)
    anames2 = [None for x in range(len(anames))]
    for k, v in anames.items():
        anames2[v] = k.decode('utf8')
    chroms = [x.decode('utf8') for x in chroms]
    """
    output = open('temp.bounds', 'w')
    print("chrom start end orient sizeL sizeR name", file=output)
    for i in range(bounds.shape[0]):
        temp = list(bounds[i])
        mid = (temp[1] + temp[2]) // 2
        temp = [chroms[temp[0]]] + temp[1:-1] + [mid - temp[1], temp[2] - mid, anames2[temp[-1]]]
        print(" ".join([str(x) for x in temp]), file=output)
    output.close()
    """
    return bounds, anames2, chroms


def load_arrays2(fname, size_fname, minsize, window):
    data = []
    for line in open(fname):
        line = line.rstrip().split()
        chrom, start, end, name = line[:4]
        data.append((chrom, int(start), int(end), int(end) - int(start), name))
    data = numpy.array(data, dtype=numpy.dtype([('chrom', 'S6'),
                                                ('start', numpy.int32),
                                                ('end', numpy.int32),
                                                ('size', numpy.int32),
                                                ('name', 'S40')]))

    chrom_dict = {}
    arrays = []
    array_dict = {}
    bounds = []
    chroms = numpy.unique(data['chrom'])
    for i, x in enumerate(chroms):
        chrom_dict[x] = i
    for chrom in chroms:
        where = numpy.where(data['chrom'] == chrom)[0]
        for i in where:
            cname = data['name'][i].decode('utf8').rstrip(')')
            if cname.endswith('p_arm') and data['end'][i] >= window//2:
                name = data['name'][i + 1].decode('utf8').split('_')[0]
                if name not in array_dict:
                    array_dict.setdefault(name, len(arrays))
                    arrays.append(name)
                bounds.append((chrom_dict[data['chrom'][i]], data['end'][i],
                               data['end'][i] + minsize, data['end'][i] - window,
                               -1, array_dict[name]))
            elif cname.endswith('q_arm'):
                name = data['name'][i - 1].decode('utf8').split('_')[0]
                if name not in array_dict:
                    array_dict.setdefault(name, len(arrays))
                    arrays.append(name)
                bounds.append((chrom_dict[data['chrom'][i]], data['start'][i] - minsize,
                               data['start'][i], data['start'][i] + window,
                               1, array_dict[name]))
    bounds = numpy.array(bounds, dtype=numpy.int32)
    chroms = [x.decode('utf8') for x in chroms]
    return bounds, arrays, chroms

if __name__ == "__main__":
    main()
