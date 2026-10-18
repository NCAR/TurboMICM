"""Additional diagnostic species errors; the official score remains score.py."""
import csv
import math
from pathlib import Path


def read(path):
    with Path(path).open() as f:
        reader = csv.reader(f); header = next(reader)
        keyed = header[:2] == ['step', 'cell']; start = 4 if keyed else 2
        values = {}
        for i, row in enumerate(reader):
            if row:
                key = (int(row[0]), int(row[1])) if keyed else i
                values[key] = [float(v) for v in row[start:]]
        return header[start:], values


def measure(output, reference):
    names, actual = read(output); columns = {name: i for i, name in enumerate(names)}
    names, expected = read(reference)
    records = []
    for j, name in enumerate(names):
        squared = maximum = bias = absolute = 0.; count = 0
        for key, values in expected.items():
            ref = values[j]
            if ref < 1e-14:
                continue
            difference = actual[key][columns[name]]-ref
            relative = abs(difference)/ref
            squared += relative*relative; maximum = max(maximum, relative)
            absolute = max(absolute, abs(difference)); bias += difference; count += 1
        if count:
            rms = math.sqrt(squared/count)
            records.append(dict(species=name, count=count, rms_relative=rms, max_relative=maximum,
                                digits=min(16., -math.log10(rms)) if rms else 16.,
                                mean_bias=bias/count, max_absolute=absolute))
    return records
