import gzip
import pandas as pd
from Bio import SeqIO
import argparse
import threading
import os

def read_clustering_results(cluster_file):
    df = pd.read_csv(cluster_file, sep='\t')
    return df

def read_fastq_file(fastq_file):
    with gzip.open(fastq_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fastq"))
    return records

def write_clustered_reads(cluster_results, fastq_records, output_prefix, num_clusters=6):
    # creat dict
    clustered_reads = {i: [] for i in range(num_clusters)}

    # reads name into list
    read_dict = {record.id: record for record in fastq_records}
    for _, row in cluster_results.iterrows():
        read_name = row['ReadName']
        cluster_id = row['Cluster']
        if read_name in read_dict:
            clustered_reads[cluster_id].append(read_dict[read_name])

    # clusters into file
    for cluster_id in range(num_clusters):
        output_file = f"{output_prefix}_cluster{cluster_id:03d}.fq.gz"
        with gzip.open(output_file, "wt") as handle:
            SeqIO.write(clustered_reads[cluster_id], handle, "fastq")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster reads and write to separate files")
    parser.add_argument('-i', '--input', type=str, required=True, help='Input clustering results file')
    parser.add_argument('-s', '--fastq-seq', type=str, required=True, help='Input fastq.gz file with reads')
    parser.add_argument('-o', '--output-prefix', type=str, required=True, help='Output prefix for clustered fastq.gz files')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('-c', '--cluster-num', type=int, default=10, help='Number of clusters to be separated')

    args = parser.parse_args()

    # set threads limit
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)

    # read cluster files and sequence files
    cluster_results = read_clustering_results(args.input)
    fastq_records = read_fastq_file(args.fastq_seq)

    # output cluster-reads sequence files
    write_clustered_reads(cluster_results, fastq_records, args.output_prefix, num_clusters=args.cluster_num)
