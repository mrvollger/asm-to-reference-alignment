#!/usr/bin/env python
import argparse
from os import read
import sys
import pandas as pd
import pybedtools


def read_bed(args):
    df = pd.read_csv(args.input, sep="\t")
    df["group"] = df.index
    new_header = (
        [f"{col}_1" for col in df.columns] + list(df.columns) + ["record_overlap"]
    )

    bed = pybedtools.BedTool().from_dataframe(df)
    inter = bed.intersect(bed, wao=True, r=True, f=args.fraction).to_dataframe(
        names=new_header
    )

    # drop the original record, just keep the its record_id
    inter["group_a"] = inter["group_1"]
    inter.drop(
        [col for col in inter.columns if col.endswith("_1")], axis=1, inplace=True
    )
    return df, inter


def find_winners(inter, args):
    winners = set()
    for group, gdf in inter.groupby("group_a"):
        gdf["record_delta"] = gdf["matches"] - gdf["matches.liftover"]
        sdf = gdf.sort_values(by=["record_delta", "matches"], ascending=False)
        if gdf.shape[0] > 1 and False:
            with pd.option_context(
                "display.max_rows",
                None,
                "display.max_columns",
                1000,
                "display.width",
                10000,
            ):
                print(sdf)
        winners.add(sdf.iloc[0]["group"])
    sys.stderr.write(
        f"Number of record reduced from {max(inter.group)} to {len(winners)}\n"
    )
    return winners
    # print(len(winners))
    # print(inter.shape)
    # print(max(inter.record_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input",
        help="paired end bed",
    )
    parser.add_argument(
        "--fraction", help="fraction of overlap", type=float, default=0.75
    )
    args = parser.parse_args()
    df, inter = read_bed(args)
    index = find_winners(inter, args)

    df.loc[index].sort_values(by="group").to_csv(sys.stdout, sep="\t", index=False)
