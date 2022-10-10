import bed
import dataclasses
import enum

import region_extractor as reex

from typing import *


class ParsingBlockLocation(enum.Enum):
    BEFORE = enum.auto()
    INSIDE = enum.auto()
    AFTER = enum.auto()


@dataclasses.dataclass
class TaggedRegionWeightEntry:
    entry: reex.RegionWeightEntry
    region: int


@dataclasses.dataclass
class TaggedOccurenceCount:
    count: int
    location: ParsingBlockLocation


@dataclasses.dataclass
class RloopParsingBlocks:
    before: List[str]
    inside: List[str]
    after: List[str]


class Dictionary:
    def __init__(self, sequence: str, gene_position: Tuple[int, int], ntuple_size: int):
        self.sequence = sequence
        self.gene_position = gene_position
        self.ntuple_size = ntuple_size
        self.rloops_parsing_blocks: List[RloopParsingBlocks] = []

    def build(
        self,
        training_set: List[bed.GenomicRegion],
        region_weight_summaries: List[reex.RegionWeightSummary],
        thresholded_region_summaries: List[reex.RegionWeightSummary],
    ):
        def find_maximum_weighted_regions(ntuple: str):
            weights_found: List[TaggedRegionWeightEntry] = []
            for region_weight_summary in thresholded_region_summaries:
                entry = next(
                    (x for x in region_weight_summary.entries if x.ntuple == ntuple),
                    None,
                )
                if entry:
                    weights_found.append(
                        TaggedRegionWeightEntry(entry, region_weight_summary.region)
                    )

            weights_found_sorted = sorted(
                weights_found, key=lambda x: x.entry.weight, reverse=True
            )

            max_weights = filter(
                lambda x: x.entry.weight == weights_found_sorted[0].entry.weight,
                weights_found_sorted,
            )

            return list(max_weights)

        def find_maximum_occurences(ntuple: str):
            before_count = 0
            inside_count = 0
            after_count = 0

            for rloop_parsing_blocks in self.rloops_parsing_blocks:
                before_count += rloop_parsing_blocks.before.count(ntuple)
                inside_count += rloop_parsing_blocks.inside.count(ntuple)
                after_count += rloop_parsing_blocks.after.count(ntuple)

            occurences = [
                TaggedOccurenceCount(before_count, ParsingBlockLocation.BEFORE),
                TaggedOccurenceCount(inside_count, ParsingBlockLocation.INSIDE),
                TaggedOccurenceCount(after_count, ParsingBlockLocation.AFTER),
            ]

            sorted_occurences = sorted(occurences, key=lambda x: x.count, reverse=True)
            maximum_occurences = filter(
                lambda x: x.count == sorted_occurences[0].count, sorted_occurences
            )

            return list(maximum_occurences)

        gene_sequence = self.sequence[self.gene_position[0] : self.gene_position[1]]

        for rloop in training_set:
            self.rloops_parsing_blocks.append(
                self.__split_into_parsing_blocks(rloop, gene_sequence)
            )

        for rloop_parsing_blocks in self.rloops_parsing_blocks:
            for before_subblock in rloop_parsing_blocks.before:
                determine_letter(before_subblock, ParsingBlockLocation.BEFORE)

    def determine_letter(ntuple: str, location: ParsingBlockLocation):
        maximum_weighted_regions = find_maximum_weighted_regions(before_subblock)
        if len(maximum_weighted_regions) == 0:
            maximum_occurences = find_maximum_occurences(before_subblock)
            print(before_subblock, maximum_occurences)
        else:
            print(before_subblock, maximum_weighted_regions)

    def __split_into_parsing_blocks(self, rloop: bed.GenomicRegion, gene_sequence: str):

        rloop_location_adjusted = (
            rloop.start - self.gene_position[0],
            rloop.end - self.gene_position[1],
        )

        parsing_block_before = gene_sequence[: rloop_location_adjusted[0]]
        parsing_block_inside = gene_sequence[
            rloop_location_adjusted[0] : rloop_location_adjusted[1]
        ]
        parsing_block_inside_reverse = parsing_block_inside[::-1]
        parsing_block_after = gene_sequence[rloop_location_adjusted[1] :]

        parsing_block_before_subblocks = [
            parsing_block_before[i : i + self.ntuple_size]
            for i in range(0, len(parsing_block_before), self.ntuple_size)
        ]

        parsing_block_before_inside_subblocks = [
            parsing_block_inside[i : i + self.ntuple_size]
            for i in range(0, len(parsing_block_inside), self.ntuple_size)
        ]

        parsing_block_before_after_reverse_subblocks = [
            parsing_block_inside_reverse[i : i + self.ntuple_size][::-1]
            for i in range(0, len(parsing_block_inside_reverse), self.ntuple_size)
        ]

        return RloopParsingBlocks(
            parsing_block_before_subblocks,
            parsing_block_before_inside_subblocks,
            parsing_block_before_after_reverse_subblocks,
        )
