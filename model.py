import collections
import dataclasses
import math
import re

import bed
import dictionary
import region_extractor as reex

from typing import *


class Model:
    def __init__(self, ntuple_size: int = 4, padding_size: int = 13):
        self.ntuple_size = ntuple_size
        self.padding_size = padding_size

    def train(
        self,
        sequence: str,
        training_set: List[bed.GenomicRegion],
        gene_position: Optional[Tuple[int, int]] = None,
    ):
        gene_position = gene_position or (0, len(sequence))

        for rloop in training_set:
            self.__make_divisible_by_width(rloop)

        extractor = reex.RegionExtractor(self.ntuple_size, self.padding_size)

        self.region_weight_summaries = extractor.generate_region_weight_summaries(
            sequence, gene_position, training_set
        )

        self.thresholded_region_weight_summaries = (
            reex.RegionExtractor.threshold_region_summaries_shannon(
                self.region_weight_summaries
            )
        )

        self.dictionary = dictionary.Dictionary(
            sequence, gene_position, self.ntuple_size
        )
        self.dictionary.build(
            training_set,
            self.region_weight_summaries,
            self.thresholded_region_weight_summaries,
        )
        # translate training set
        # determine probabilities
        pass

    def __make_divisible_by_width(self, rloop: bed.GenomicRegion):
        genomic_segment_length = abs(rloop.end - rloop.start)
        remainder = genomic_segment_length % self.ntuple_size
        if remainder > self.ntuple_size // 2:
            rloop.start += genomic_segment_length % -self.ntuple_size
        else:
            rloop.start += remainder

    def union(self, other: "Model") -> "Model":
        pass

    def predict(self, gene_sequence: str):
        pass
