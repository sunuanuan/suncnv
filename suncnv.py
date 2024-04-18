#!/usr/bin/env python3

import argparse
import gzip
import random
import sys
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

now = lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def get_target_regions(bed: Path) -> list:
    target_regions = []
    if bed.suffix == ".gz":
        f = gzip.open(bed, "rt")
    else:
        f = open(bed)
    for line in f:
        row = line.strip().split("\t")
        if len(row) == 3:
            chrom, start, end = row
            gene = "."
        else:
            chrom, start, end, gene, *_ = row
        target_regions.append((chrom, int(start), int(end), gene))
    f.close()
    return target_regions


def count_reads(bam: Path, regions: list) -> list:
    counts = []
    with pysam.AlignmentFile(bam, threads=4) as bamfile:
        for chrom, start, end, gene in regions:
            read_names = set()
            read_names_uniq = set()
            fragments = fragments_uniq = 0
            for alignment in bamfile.fetch(region=f"{chrom}:{start}-{end}"):
                if not (alignment.is_secondary or alignment.is_duplicate or alignment.is_unmapped):
                    read_name = alignment.query_name
                    if read_name not in read_names:
                        fragments += 1
                        read_names.add(read_name)
                    if alignment.mapping_quality >= 40 and read_name not in read_names_uniq:
                        fragments_uniq += 1
                        read_names_uniq.add(read_name)
            counts.append((chrom, start, end, gene, fragments, fragments_uniq))
    return counts


def norm_counts(counts: dict) -> pd.DataFrame:
    df = pd.DataFrame(counts, columns=["chrom", "start", "end", "gene", "fragment", "ufragment"])
    df["length"] = df["end"] - df["start"]
    df["FR"] = df["fragment"] / df["length"]
    df["FR"] = df["FR"] / df["FR"].median()
    df["FR"] = df["FR"].round(4)
    df["UFR"] = df["ufragment"] / df["length"]
    df["UFR"] = df["UFR"] / df["UFR"].median()
    df["UFR"] = df["UFR"].round(4)
    return df


def count_bam(input_bam: Path, out_count: Path, bed: Path) -> pd.DataFrame:
    print(f"{now()} INFO: counting reads of {input_bam}", flush=True)
    regions = get_target_regions(bed)
    counts = count_reads(input_bam, regions)
    df = norm_counts(counts)
    df.to_csv(out_count, sep="\t", index=False)
    print(f"{now()} INFO: count {input_bam} -> {out_count}", flush=True)
    return df


def build_reference(count_array: list, nearest_samples: list) -> pd.DataFrame:
    print(f"{now()} INFO: building reference")
    nearest_counts = []
    for i, _ in nearest_samples:
        nearest_counts.append(count_array[i])
    reference_fr = pd.DataFrame(nearest_counts[0])[["FR"]]
    reference_fr.columns = ["FR_1"]
    reference_ufr = pd.DataFrame(nearest_counts[0])[["UFR"]]
    reference_ufr.columns = ["UFR_1"]
    for i, count in enumerate(nearest_counts[1:], start=2):
        df_fr = pd.DataFrame(count)[["FR"]]
        df_fr.columns = [f"FR_{i}"]
        df_ufr = pd.DataFrame(count)[["UFR"]]
        df_ufr.columns = [f"UFR_{i}"]
        reference_fr = pd.concat([reference_fr, df_fr], axis=1)
        reference_ufr = pd.concat([reference_ufr, df_ufr], axis=1)
    reference = pd.DataFrame()
    reference[["FR_q25", "FR_q50", "FR_q75"]] = reference_fr.quantile([0.25, 0.50, 0.75], axis=1).round(4).T
    reference["FR_std"] = reference_fr.std(axis=1).round(4)
    reference[["UFR_q25", "UFR_q50", "UFR_q75"]] = reference_ufr.quantile([0.25, 0.50, 0.75], axis=1).round(4).T
    reference["UFR_std"] = reference_ufr.std(axis=1).round(4)
    print(f"{now()} INFO: reference built")
    return reference


def get_input_count(file: Path, outdir: Path, bed: Path, overwrite: bool, only_count: bool = False) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    if file.suffix == ".bam":
        if not bed:
            sys.exit(f"{now()} ERROR: bed file [-b/--bed] is required for .bam input")
        out_count = outdir / f"{file.stem}.count"
        if overwrite or not out_count.exists():
            input_count = count_bam(input_bam=file, out_count=out_count, bed=bed)
        else:
            print(f"{now()} INFO: {out_count} exists, skipping count", flush=True)
            input_count = pd.read_csv(out_count, sep="\t") if not only_count else None
    else:
        input_count = pd.read_csv(file, sep="\t")
    return input_count


def guess_sex(count: pd.DataFrame) -> tuple:
    x_median = round(count["FR"][count["chrom"] == "chrX"].median(), 4)
    y_median = round(count["FR"][count["chrom"] == "chrY"].median(), 4)
    x_count = round(x_median / 0.5)
    y_count = round(y_median / 0.7)
    sex = "x" * x_count + "y" * y_count
    return sex, x_median, y_median


def get_nearest_samples(reference: list, count_array: list, input_count: pd.DataFrame, input_sex: str, number: int):
    near_count = 0
    nearest_samples = []
    for i, count in enumerate(count_array):
        sex, *_ = guess_sex(count)
        if sex == input_sex:
            corr = round(input_count["FR"].corr(count["FR"]), 4)
            if corr > 0.99:
                print(f"{now()} WARN: correlation with {reference[i].name} is {corr} > 0.99, assume the same sample, skipping")
            else:
                if corr > 0.90:
                    nearest_samples.append((i, corr))
                if corr > 0.95:
                    near_count += 1
                if near_count == number:
                    remaining = len(reference) - i
                    print(f"{now()} INFO: already found {number}/{i} samples with correlation > 0.95, skipping the remaining {remaining} samples")
                    break
    nearest_samples = sorted(nearest_samples, key=lambda x: x[1], reverse=True)
    sample_count = len(nearest_samples)
    if sample_count < number:
        print(f"{now()} WARN: number of nearest samples is {sample_count} < {number}")
    nearest_samples = nearest_samples[:number]
    median_corr = nearest_samples[min(number, sample_count) // 2][1]
    print(f"{now()} INFO: median correlation is {median_corr}")
    nearest_files = ", ".join(sorted([reference[i].name for i, _ in nearest_samples]))
    print(f"{now()} INFO: nearest samples are {nearest_files}")
    return nearest_samples


def reset_input_sex(input_count: pd.DataFrame) -> str:
    sex, x_median, y_median = guess_sex(input_count)
    if sex not in {"xx", "xy"}:
        print(f"{now()} WARN: sex karyotype of input is {sex}, x_median={x_median}, y_median={y_median}, reset to xx")
        sex = "xx"
    else:
        print(f"{now()} INFO: sex karyotype of input is {sex}, x_median={x_median}, y_median={y_median}")
    return sex


def calc_ratio(input_count: pd.DataFrame, reference_count: pd.DataFrame) -> pd.DataFrame:
    df = pd.merge(input_count, reference_count, right_index=True, left_index=True)
    df["rFR_q25"] = ((df["FR"] + 0.05) / (df["FR_q25"] + 0.05)).round(4)
    df["rFR_q50"] = ((df["FR"] + 0.05) / (df["FR_q50"] + 0.05)).round(4)
    df["rFR_q75"] = ((df["FR"] + 0.05) / (df["FR_q75"] + 0.05)).round(4)
    df["rUFR_q25"] = ((df["UFR"] + 0.05) / (df["UFR_q25"] + 0.05)).round(4)
    df["rUFR_q50"] = ((df["UFR"] + 0.05) / (df["UFR_q50"] + 0.05)).round(4)
    df["rUFR_q75"] = ((df["UFR"] + 0.05) / (df["UFR_q75"] + 0.05)).round(4)
    return df


def get_status(ratio: pd.DataFrame) -> pd.DataFrame:

    def _get_status(row, uniq=False):
        fr = row["UFR"] if uniq else row["FR"]
        fr_q25 = row["UFR_q25"] if uniq else row["FR_q25"]
        fr_q50 = row["UFR_q50"] if uniq else row["FR_q50"]
        fr_q75 = row["UFR_q75"] if uniq else row["FR_q75"]
        rfq_q25 = row["rUFR_q25"] if uniq else row["rFR_q25"]
        rfq_q50 = row["rUFR_q50"] if uniq else row["rFR_q50"]
        rfr_q75 = row["rUFR_q75"] if uniq else row["rFR_q75"]
        fr_std = row["UFR_std"] if uniq else row["FR_std"]
        status = ""
        if fr_q50 < 0.05:
            status = "????"
        if fr_q50 < 0.10:
            status = "???"
        elif fr_q50 < 0.15:
            status = "??"
        elif fr_q50 < 0.20:
            status = "?"
        else:
            status = ""
        if rfq_q50 < 0.70:
            status += "-"
            if rfq_q50 < 0.50:
                status += "-"
            if rfq_q25 < 0.70:
                status += "-"
            if rfq_q25 < 0.50:
                status += "-"
            if fr < fr_q50 - fr_std * 3:
                status += "*"
            if fr < fr_q25 - (fr_q75 - fr_q25) * 1.5:
                status += "*"
        elif rfq_q50 > 1.30:
            status += "+"
            if rfq_q50 > 1.50:
                status += "+"
            if rfr_q75 > 1.30:
                status += "+"
            if rfr_q75 > 1.50:
                status += "+"
            if fr > fr_q50 + fr_std * 3:
                status += "*"
            if fr > fr_q75 + (fr_q75 - fr_q25) * 1.5:
                status += "*"
        else:
            status += "."
        return status

    _get_ustatus = lambda row: _get_status(row, uniq=True)
    ratio["status"] = ratio.apply(_get_status, axis=1)
    ratio["ustatus"] = ratio.apply(_get_ustatus, axis=1)
    columns = ["chrom", "start", "end", "gene", "rFR_q50", "rUFR_q50", "status", "ustatus"]
    return ratio[columns]


def find_seeds(array: np.array, uarray: np.array) -> list:
    seeds = []
    for i in range(len(array)):
        status, ustatus = array[i], uarray[i]
        if len(status) >= 3 and ("?" in ustatus or len(ustatus) >= 3):
            seeds.append(i)
    return seeds


def adjacent_seeds(seeds: list, max_gap_count: int) -> list:
    segments = []
    start = 0
    end = start + 1
    while end < len(seeds):
        if seeds[end] - seeds[end - 1] > max_gap_count + 1:
            if end - start > 1:
                segments.append((start, end - 1))
            start = end
            end = start + 1
        else:
            end += 1
    return segments


def filter_segment(status_drop: pd.DataFrame, seeds: list, segments: list, min_bin_count: int) -> list:
    filtered_segments = []
    for start, end in segments:
        segment = status_drop.iloc[seeds[start]:seeds[end] + 1]
        if len(segment) >= min_bin_count:
            ratio = segment["rFR_q50"].median().__round__(4)
            uratio = segment["rUFR_q50"].median().__round__(4)
            pos_chrom = segment["chrom"].iloc[0]
            pos_start = segment["start"].min()
            pos_end = segment["end"].max()
            pos_length = pos_end - pos_start + 1
            genes = ",".join(segment["gene"].unique().tolist())
            if ratio < 0.7:
                filtered_segments.append((pos_chrom, pos_start, pos_end, pos_length, ratio, uratio, "DEL", genes))
            elif ratio > 1.3:
                filtered_segments.append((pos_chrom, pos_start, pos_end, pos_length, ratio, uratio, "DUP", genes))
            else:
                pass
    return filtered_segments


def make_segments(status: pd.DataFrame, min_bin_count: int, max_gap_count: int, output: Path):
    print(f"{now()} INFO: segment by status")
    status_drop = status[~status["status"].str.startswith("?")]
    ratio = 100 - len(status_drop) / len(status) * 100
    print(f"{now()} INFO: keep {len(status_drop)} of {len(status)} bins, drop {ratio:.2f}% bins")
    array = status_drop["status"].to_numpy()
    uarray = status_drop["ustatus"].to_numpy()
    seeds = find_seeds(array, uarray)
    print(f"{now()} INFO: {len(seeds)} seeds found")
    if seeds:
        segments = adjacent_seeds(seeds, max_gap_count)
        print(f"{now()} INFO: {len(segments)} segments found")
        if segments:
            filtered_segments = filter_segment(status_drop, seeds, segments, min_bin_count)
            print(f"{now()} INFO: {len(filtered_segments)} segments kept")
            with output.open("w") as f:
                header = ["chrom", "start", "end", "length", "ratio", "uratio", "type", "genes"]
                print("\t".join(header), file=f)
                for row in filtered_segments:
                    print("\t".join(map(str, row)), file=f)


def plot():
    pass


def main(input: Path, outdir: Path, bed: Path, reference: list, number: int, threads: int, overwrite: bool, min_bin_count: int, max_gap_count: int):
    if input:
        if input.suffix in {".bam", ".count", ".status"}:
            outdir.mkdir(parents=True, exist_ok=True)
            if input.suffix != ".status":
                if input.suffix == ".count" and not reference:
                    sys.exit(f"{now()} ERROR: reference is required if input is .count")
                input_count = get_input_count(input, outdir, bed, overwrite)
                sex = reset_input_sex(input_count)
                if reference:
                    random.shuffle(reference)
                    with Pool(threads) as p:
                        count_array = p.starmap(get_input_count, [(file, outdir, bed, overwrite) for file in reference])
                nearest_samples = get_nearest_samples(reference, count_array, input_count, sex, number)
                reference = build_reference(count_array, nearest_samples)
                ratio = calc_ratio(input_count, reference)
                status = get_status(ratio)
                status.to_csv(outdir / (input.stem + ".status"), sep="\t", index=False)
            else:
                status = pd.read_csv(input, sep="\t")
            make_segments(status, min_bin_count, max_gap_count, outdir / (input.stem + ".cnv"))
        else:
            sys.exit(f"{now()} ERROR: input file suffix must be .bam or .count or .status")
    else:
        if reference:
            print(f"{now()} INFO: count bam for {len(reference)} samples")
            random.shuffle(reference)
            with Pool(threads) as p:
                p.starmap(get_input_count, [(file, outdir, bed, overwrite, True) for file in reference])
        else:
            sys.exit(f"{now()} ERROR: reference is required if input is not provided")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, help="input .bam or .count or .status file")
    parser.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    parser.add_argument("-b", "--bed", type=Path, help="target .bed or .bed.gz file, required for .bam input")
    parser.add_argument("-r", "--reference", type=Path, nargs="+", help="reference .bam or .count files")
    parser.add_argument("-n", "--number", type=int, default=20, help=f"number of samples to build reference, default: 20")
    parser.add_argument("-t", "--threads", type=int, default=8, help=f"number of threads to use, default: 8")
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing .count files, default: False")
    parser.add_argument("--min-bin-count", type=int, default=3, help="minimum bin count to keep, default: 3")
    parser.add_argument("--max-gap-count", type=int, default=3, help="maximum gap count to keep, default: 3")
    args = parser.parse_args()
    main(**vars(args))
