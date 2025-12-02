#!/usr/bin/env python3
# Author: Zifo Bioinformatics (nf-core/rnasplice)
import argparse, itertools

def main(args):
    with open(args.file, "r") as handle:
        header = handle.readline()
    samples = header.split("\t")
    conditions = [sample.rsplit("_", 1)[0] for sample in samples]
    last_index = 0
    out = []
    for _, g in itertools.groupby(enumerate(conditions), lambda k: k[1]):
        l = [*g]
        out.append([last_index + 1, l[-1][0] + 1])
        last_index += len(l)
    assert len(out) == 2, "Column numbers have to be continuous, with no overlapping or missing columns between them."
    groups = ",".join([f"{start}-{end}" for start, end in out])
    print(groups, end="")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("file")
    main(p.parse_args())
