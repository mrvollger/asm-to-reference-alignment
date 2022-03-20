#!/usr/bin/env python
import argparse
from email.policy import strict
import sys
import pysam
from intervaltree import Interval, IntervalTree
import pprint
import logging
import math
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"


pp = pprint.PrettyPrinter(width=60)  # , compact=True)


def make_interval_tree(bed):
    # dict of [sample,hap,chromosome] -> interval tree
    super_tree = {}
    for line in open(bed).readlines():
        if line[0] == "#":
            continue
        rec = line.split("\t")
        chrm, start, end, contig = rec[0], int(rec[1]), int(rec[2]), rec[3]
        # sample, hap = contig.split("_")[:2]
        sample, hap = contig.rsplit("_", 1)
        hap = int(hap)
        # hap = 1 if hap == 2 else 2

        if (sample, hap, chrm) not in super_tree:
            logging.debug(f"Adding {sample} {hap} {chrm} to interval tree")
            super_tree[(sample, hap, chrm)] = IntervalTree()

        super_tree[(sample, hap, chrm)][start:end] = contig
    return super_tree


def is_in_hap_coverage(sample, hap, chrm, pos, hap_coverage):
    if (sample, hap, chrm) not in hap_coverage:
        logging.debug(f"'{sample}' '{hap}' '{chrm}' not in hap_coverage")
        return False
    cov = hap_coverage[(sample, hap, chrm)][pos : pos + 1]
    return len(cov) > 0


def get_cov_based_genotype_tuple(gts_tuple, sample, chrm, pos, hap_coverage):
    gts = [gts_tuple[0], gts_tuple[1]]
    for idx, gt in enumerate(gts):
        hap = idx + 1
        if gt is None and is_in_hap_coverage(sample, hap, chrm, pos, hap_coverage):
            gts[idx] = 0
    return tuple(gts)


def make_chunks_form_vcf(vcf_in, args):
    if args.chunk is None:
        return [None]
    chunk_num = args.chunk[0]
    num_chunks = args.chunk[1]
    total_bp = 0
    length_dict = {}
    for contig, rec in vcf_in.header.contigs.items():
        total_bp += rec.length
        length_dict[contig] = rec.length

    step_size = int(total_bp / num_chunks) + 1
    regions = []
    use_regions = []
    cur_pos = 0
    total_bp_covered = 0
    for contig, length in length_dict.items():
        while cur_pos < length:
            next_pos = min(cur_pos + step_size, length)
            total_bp_covered += next_pos - cur_pos
            cur_chunk = math.floor((total_bp_covered - 1) / step_size) + 1
            regions.append((contig, cur_pos, next_pos, cur_chunk))
            if cur_chunk == chunk_num:
                use_regions.append([contig, cur_pos, next_pos])
            cur_pos = next_pos
        cur_pos = 0

    #
    # check the chunking is correct
    #
    check = 0
    chunk_num_s = set()
    for contig, start, end, chunk in regions:
        check += end - start
        chunk_num_s.add(chunk)
    assert chunk_num_s == set(
        range(1, num_chunks + 1)
    ), f"Chunk numbers are not correct: {chunk_num_s}"
    assert check == total_bp

    return use_regions


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf", help="paired end bed")
    parser.add_argument("bed", help="paired end bed")
    parser.add_argument(
        "--region", "-r", help="region to process, requires index", default=None
    )
    parser.add_argument(
        "--chunk",
        "-c",
        help="Process the Nth chunk of of X total chunks [ -c N X]",
        nargs=2,
        type=int,
        default=None,
    )
    parser.add_argument("--verbose", "-v", action="count", default=1)
    parser.add_argument("--exclude", "-e", help="Chromosomes to exclude", default=[])
    parser.add_argument(
        "--strict", "-s", help="Force results to be phased", action="store_true"
    )
    parser.add_argument("--threads", "-t", type=int, default=1)
    args = parser.parse_args()
    args.verbose = 40 - (10 * args.verbose) if args.verbose > 0 else 0
    logging.basicConfig(
        level=args.verbose,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    hap_coverage = make_interval_tree(args.bed)
    # pp.pprint(hap_coverage)
    # print(is_in_hap_coverage("HG002", 1, "chr1", 100, hap_coverage))
    # print(is_in_hap_coverage("HG002", 1, "chr4", 100, hap_coverage))
    # exit()

    vcf_in = pysam.VariantFile(
        args.vcf, threads=args.threads
    )  # auto-detect input format
    regions = make_chunks_form_vcf(vcf_in, args)

    # pp.pprint(list(vcf_in.header.records))
    vcf_out = pysam.VariantFile("-", "w", header=vcf_in.header, threads=1)
    changed_gts = 0
    total_gts = 0
    none_count = 0
    un_phased = 0
    for region in regions:
        if args.region is not None:
            vcf_iter = vcf_in.fetch()
        elif region is None:
            vcf_iter = vcf_in.fetch()
        else:
            logging.info(f"Processing {region}")
            vcf_iter = vcf_in.fetch(*region)

        for idx, rec in enumerate(vcf_iter):
            if rec.chrom in args.exclude:
                continue
            for sample in rec.samples:
                gts = rec.samples[sample]["GT"]
                total_gts += 1

                # if not strict then we can check for ambiguous genotypes that are fake
                if not args.strict and "GAP2" in rec.filter and gts[1] is None:
                    rec.samples[sample].phased = True

                # strict mode
                # one genotypes is missing and it is not phased
                if (
                    not rec.samples[sample].phased
                    and args.strict
                    and gts.count(None) == 1
                ):
                    rec.samples[sample]["GT"] = (None, None)
                    un_phased += 1
                elif None in gts:
                    none_count += 1
                    new_gt = get_cov_based_genotype_tuple(
                        gts, sample, rec.chrom, rec.pos, hap_coverage
                    )
                    if new_gt[0] != gts[0] or new_gt[1] != gts[1]:
                        rec.samples[sample]["GT"] = new_gt
                        logging.debug(f"Updated to: {new_gt}")
                        changed_gts += 1
                    # always need phased outputs
                    rec.samples[sample].phased = True
            vcf_out.write(rec)
            logging.debug(f"{idx+1} variants proccessed")
    logging.info(
        "{:.2%} of genotypes changed: {:,} of {:,}. {:,} total missing GTs. {:,} un-phased variants.".format(
            changed_gts / (total_gts + 0.00001),
            changed_gts,
            total_gts,
            none_count,
            un_phased,
        )
    )
