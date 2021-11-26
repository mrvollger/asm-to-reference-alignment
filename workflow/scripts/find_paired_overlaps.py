#!/usr/bin/env python
import argparse
import sys


def overlap(a0, a1, b0, b1):
    r = min(a1, b1) - max(a0, b0)
    if r >= 0:
        return r
    return 0


def intersect(a, a0, a1, b, b0, b1, dist):
    """
    check if two genomic intervals overlap
    """
    intersects = (a == b) and ((a1 + dist) >= b0) and ((b1 + dist) >= a0)
    if intersects:
        return overlap(a0, a1, b0, b1)
    else:
        return 0


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
    parser.add_argument(
        "-n", "--no-group", action="store_true", help="do not attempt to make groups."
    )
    parser.add_argument(
        "-f",
        "--fraction",
        type=float,
        default=0.0,
        help="group things with this level of reciprocal overlap in the first region of the pair.",
    )
    args = parser.parse_args()

    pre = None
    group = 0
    for cur in read_bed_line(args):
        if pre is None:
            print(f"{cur[6].strip()}\t{group}")
            pre = cur
            continue
        first_intersects = intersect(
            pre[0], pre[1], pre[2], cur[0], cur[1], cur[2], args.dist
        )
        second_intersects = intersect(
            pre[3], pre[4], pre[5], cur[3], cur[4], cur[5], args.dist
        )
        if (
            args.fraction
            and first_intersects
            and first_intersects > args.fraction * (pre[2] - pre[1])
            and first_intersects > args.fraction * (cur[2] - cur[1])
        ):
            group += 0
        elif first_intersects and second_intersects and not args.no_group:
            group += 0
        else:
            group += 1
        # print(cur, first_intersects, second_intersects, group)
        print(f"{cur[6].strip()}\t{group}")
        pre = cur
