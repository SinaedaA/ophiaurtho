import pandas as pd
import glob
import re
import argparse
import os
import warnings

def main():
    indir, outfile = parsing_arguments()
    all_benchmarks = glob.glob(os.path.join(indir, "*", "*.tsv"))
    bench_dfs = [pd.read_csv(tsv, sep = "\t") for tsv in all_benchmarks]
    bench_concat = pd.concat(bench_dfs)
    bench_summary = bench_concat[["rule_name", "wildcards", "threads", "cpu_usage", "cpu_time", "mean_load", "s", "h:m:s", "threads", "params"]]
    bench_summary.to_csv(outfile, sep = "\t", index = False)

def parsing_arguments():
    """
    Parses the arguments from the command line.
    """
    parser = argparse.ArgumentParser(description="Combine benchmark files from different rules into a single summary table.")
    parser.add_argument("--benchmark_dir", type=str, required=True, help="Directory containing benchmark files.")
    parser.add_argument("--output_file", type=str, required=True, help="Output file to save the combined benchmarks.")

    args = parser.parse_args()
    benchmark_dir = os.path.abspath(args.benchmark_dir)
    outfile = os.path.abspath(args.output_file)

    return benchmark_dir, outfile

if __name__ == '__main__':
    main()