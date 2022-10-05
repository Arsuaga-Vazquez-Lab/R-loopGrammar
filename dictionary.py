import bed
import dataclasses

import region_extractor as reex

from typing import *


@dataclasses.dataclass
class RloopParsingBlocks:
    before: List[str]
    inside: List[str]
    after_reversed: List[str]


class Dictionary:
    def __init__(self, sequence: str, gene_position: Tuple[int, int], ntuple_size: int):
        self.sequence = sequence
        self.gene_position = gene_position
        self.ntuple_size = ntuple_size
        self.rloop_parsing_blocks: List[RloopParsingBlocks] = []

    def build(
        self,
        training_set: List[bed.GenomicRegion],
        region_weight_summaries: List[reex.RegionWeightSummary],
        thresholded_region_summaries: List[reex.RegionWeightSummary],
    ):
        gene_sequence = self.sequence[self.gene_position[0] : self.gene_position[1]]

        for rloop in training_set:
            self.rloop_parsing_blocks.append(
                self.__split_into_parsing_blocks(rloop, gene_sequence)
            )

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
