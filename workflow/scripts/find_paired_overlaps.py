#!/usr/bin/env python
import numpy as np
from numba import njit
import pandas as pd
import argparse
import sys


@njit
def intersect(a, a0, a1, b, b0, b1, dist):
    """
    check if two genomic intervals overlap
    """
    # if limit:
    #    dist = min(dist, limit * (a1 - a0 + b1 - b0))
    return (a == b) and ((a1 + dist) >= b0) and ((b1 + dist) >= a0)


def read_bed_line(args):
    for oline in args.input:
        if oline[0] == "#":
            print(oline.strip() + "\tgroup")
            continue
        line = oline.strip().split("\t")
        chr1 = line[0]
        st1 = int(line[1])
        en1 = int(line[2])
        chr2 = line[args.cols[0] - 1]
        st2 = int(line[args.cols[1] - 1])
        en2 = int(line[args.cols[2] - 1])
        yield (chr1, st1, en1, chr2, st2, en2, oline)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input",
        help="paired end bed",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
    )
    parser.add_argument(
        "-d", "--dist", help="distance allowed between intervals", type=int, default=100
    )
    parser.add_argument(
        "-c",
        "--cols",
        nargs="+",
        type=int,
        help="in order list of columns to use for second set of intervals",
        default=[11, 12, 13],
    )
    args = parser.parse_args()

    pre = None
    group = 0
    for cur in read_bed_line(args):
        if pre is None:
            pre = cur
            continue
        first_intersects = intersect(*pre[0:3], *cur[0:3], args.dist)
        second_intersects = intersect(*pre[3:6], *cur[3:6], args.dist)
        if not (first_intersects and second_intersects):
            group += 1
        # print(cur, first_intersects, second_intersects, group)
        print(f"{pre[6].strip()}\t{group}")
        pre = cur