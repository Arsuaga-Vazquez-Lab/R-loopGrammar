import region_handler
import bed
import fasta
import random

from typing import *

FASTA_FILE = "pFC53.fa"
PADDING = 12
WIDTH = 3
BED_FILE = "pFC53_SUPERCOILED.bed"
GENE_START = 1
GENE_END = 1929


def make_divisible_by_width(genomic_region: bed.GenomicRegion, width: int):
    genomic_segment_length = genomic_region.end - genomic_region.start
    remainder = genomic_segment_length % width
    if remainder > width // 2:
        genomic_region.start += genomic_segment_length % -width
    else:
        genomic_region.start += remainder


def create_training_set(rloops, size):
    return random.sample(rloops, size)


def main():
    with open(FASTA_FILE, 'r') as fasta_file:
        fasta_reader = fasta.Reader(fasta_file)
        plasmid_sequence = next(fasta_reader).sequence

    with open(BED_FILE, 'r') as bed_file:
        bed_reader = bed.Reader(bed_file)
        rloops = list(bed_reader)

    training_set = create_training_set(rloops, 63)

    for genomic_region in training_set:
        make_divisible_by_width(genomic_region, WIDTH)

    with open('training_set.bed', 'w') as bed_file:
        bed_writer = bed.Writer(bed_file)
        bed_writer.write_genomic_regions(training_set)

    region_summaries = region_handler.generate_region_summaries(
        plasmid_sequence, training_set, WIDTH, PADDING, GENE_START, GENE_END
    )

    region_summaries_thresholed = region_handler.threshold_shannon(
        region_summaries)

    region_handler.write_region_summaries(region_summaries, 'output.xlsx')
    region_handler.threshold_region_summaries_shannon(
        region_summaries_thresholed, 'output_threshold.xlsx')


if __name__ == "__main__":
    main()
