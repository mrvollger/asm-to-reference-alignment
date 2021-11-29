#!/usr/bin/env python
import argparse
from os import read
import sys
import pandas as pd
import pybedtools
import networkx as nx


def make_self_intersect(df, names, args):
    cols = ["chrom", "st", "en", "record_id"]
    new_header = list(cols) + [f"{col}_b" for col in cols] + ["record_overlap"]

    # add slop to the overlapping intervals
    n_df = df[names].copy()
    n_df.columns = cols
    n_df.st = n_df.st - args.distance
    n_df.en = n_df.en + args.distance
    n_df.loc[n_df.st < 0, "st"] = 0

    bed = pybedtools.BedTool().from_dataframe(n_df)
    inter = bed.intersect(
        bed, wao=True, f=args.fraction, r=args.reciprocal
    ).to_dataframe(names=new_header)
    # filter for min number of base overlaps
    inter = inter[inter.record_overlap > args.overlap]
    # make a key based on the records that are intersecting
    inter["key"] = inter.apply(
        lambda r: f"{r.record_id}_{r.record_id_b}"
        if r.record_id < r.record_id_b
        else f"{r.record_id_b}_{r.record_id}",
        axis=1,
    )
    return inter


def print_source_windows(df):
    odf = (
        df["original_source"]
        .str.extract("(.*):(\d+)-(\d+)", expand=True)
        .astype({1: "int32", 2: "int32"})
    )
    out = odf.groupby(0, as_index=False).agg({1: "min", 2: "max"})
    out.to_csv(sys.stdout, sep="\t", index=False, header=False)


def read_bed(args):
    df = pd.read_csv(args.input, sep="\t")
    df["record_id"] = df.index
    names1 = ["#reference_name", "reference_start", "reference_end", "record_id"]
    names2 = [
        "reference_name.liftover",
        "reference_start.liftover",
        "reference_end.liftover",
        "record_id",
    ]
    inter1 = make_self_intersect(df, names1, args)
    inter2 = make_self_intersect(df, names2, args)

    paired_intersect = inter1.key.isin(inter2.key)
    valid_intersect = inter1[paired_intersect]
    graph = nx.Graph()
    graph.add_edges_from(zip(valid_intersect.record_id, valid_intersect.record_id_b))

    # get ready to print
    df.drop("record_id", inplace=True, axis=1)
    df["group"] = 0
    sys.stdout.write("\t".join(df.columns) + "\n")
    for idx, group in enumerate(nx.connected_components(graph)):
        cur_df = df.loc[sorted(group)]
        cur_df["group"] = idx + 1
        if args.source_windows:
            print_source_windows(cur_df)
        else:
            cur_df.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", help="paired end bed", default="tmp.bed")
    parser.add_argument(
        "--fraction", help="fraction of overlap", type=float, default=0.01
    )
    parser.add_argument(
        "-o", "--overlap", help="minimum required overlap.", type=int, default=0
    )
    parser.add_argument(
        "-d",
        "--distance",
        help="allowed distance between overlapping segments.",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--source-windows",
        help="print the source windows that appear to be part of one gene conversion event",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--reciprocal",
        help="force overlaps to be reciprocal.",
        action="store_true",
    )
    args = parser.parse_args()
    read_bed(args)
