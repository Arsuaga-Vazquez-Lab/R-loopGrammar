from typing import *
from typing.io import *

import dataclasses
import enum

"""
with open(...) as bed_file:
	bed_reader = bed.Reader(bed_file)
	for genomic_region in bed_reader:
		...
"""


class StrandOrientation(enum.Enum):
    POSITIVE = "+"
    NEGATIVE = "-"
    NONE = "."

    def __str__(self):
        return self.value


@dataclasses.dataclass
class GenomicRegion:
    chrom: str
    start: int
    end: int
    name: Optional[str] = None
    score: Optional[int] = None
    strand: Optional[StrandOrientation] = None
    thick_start: Optional[int] = None
    thick_end: Optional[int] = None
    item_rgb: Optional[str] = None
    block_count: Optional[int] = None
    block_sizes: Optional[List[int]] = None
    block_starts: Optional[List[int]] = None


class Reader:
    def __init__(self, bed_file: TextIO):
        self.bed_file = bed_file

    def __iter__(self):
        return self

    def __next__(self):
        line = self.bed_file.readline()

        if not line:
            raise StopIteration

        separated_line = line.split()
        segment_iter = iter(separated_line)

        try:
            arguments = []

            # mandatory arguments
            chrom = next(segment_iter)
            start = int(next(segment_iter))
            end = int(next(segment_iter))

            arguments.extend([chrom, start, end])

            name = next(segment_iter)
            arguments.append(name)

            score = int(next(segment_iter))
            arguments.append(score)

            strand = next(segment_iter)

            if strand == "+":
                orientation = StrandOrientation.POSITIVE
            elif strand == "-":
                orientation = StrandOrientation.NEGATIVE
            elif strand == ".":
                orientation = StrandOrientation.NONE
            else:
                orientation = None

            arguments.append(orientation)

            thick_start = int(next(segment_iter))
            arguments.append(thick_start)

            thick_end = int(next(segment_iter))
            arguments.append(thick_end)

            item_rgb = next(segment_iter)
            arguments.append(item_rgb)

            block_count = int(next(segment_iter))
            arguments.append(block_count)

            block_sizes = list(map(int, next(segment_iter).split(",")))
            arguments.append(block_sizes)

            block_starts = list(map(int, next(segment_iter).split(",")))
            arguments.append(block_starts)
        except StopIteration:
            pass
        finally:
            return GenomicRegion(*arguments)


class Writer:
    def __init__(self, bed_file: TextIO):
        self.bed_file = bed_file

    def write_genomic_region(self, genomic_region: GenomicRegion):
        genomic_region_arguments = list(
            map(str, filter(lambda x: x != None, dataclasses.astuple(genomic_region)))
        )

        self.bed_file.write("\t".join(genomic_region_arguments) + "\n")

    def write_genomic_regions(self, genomic_regions: List[GenomicRegion]):
        for genomic_region in genomic_regions:
            self.write_genomic_region(genomic_region)
