import collections
import dataclasses
import math
import re

import bed

from typing import *


@dataclasses.dataclass
class RegionWeightEntry:
    ntuple: str
    region_count: int
    nucleotide_count: Dict[str, int]
    gene_count: int
    weight: float


@dataclasses.dataclass
class RegionWeightSummary:
    region: int
    entries: List[RegionWeightEntry]


class RegionExtractor:
    def __init__(self, ntuple_size: int, padding_size: int):
        self.ntuple_size = ntuple_size
        self.padding_size = padding_size

    def generate_region_weight_summaries(
        self,
        sequence: str,
        gene_position: Tuple[int, int],
        rloops: List[bed.GenomicRegion],
    ) -> List[RegionWeightSummary]:
        def increment_occurences(occurences_dict, regions):
            for window in regions:
                occurences_dict[window] += 1

        regions_occurences: List[DefaultDict[str, int]] = [
            collections.defaultdict(lambda: 0) for _ in range(4)
        ]

        rloop_count = 0
        for regions in self.__extract_all_regions(sequence, gene_position, rloops):
            rloop_count += 1

            for i in range(len(regions)):
                increment_occurences(regions_occurences[i], regions[i])

        region_summaries = [RegionWeightSummary(i, []) for i in range(0, 4)]
        for i, region_occurences in enumerate(regions_occurences):
            for ntuple, region_count in region_occurences.items():
                a_count = ntuple.count("A")
                t_count = ntuple.count("T")
                c_count = ntuple.count("C")
                g_count = ntuple.count("G")

                nucleotide_count = {
                    "A": a_count,
                    "T": t_count,
                    "C": c_count,
                    "G": g_count,
                }

                gene = sequence[gene_position[0] : gene_position[1]]
                gene_count = len(re.findall(f"(?={ntuple})", gene))
                weight = region_count / (gene_count * rloop_count)

                entry = RegionWeightEntry(
                    ntuple, region_count, nucleotide_count, gene_count, weight
                )
                region_summaries[i].entries.append(entry)

        for region_summary in region_summaries:
            region_summary.entries = sorted(
                region_summary.entries, key=lambda e: e.weight, reverse=True
            )

        return region_summaries

    def __extract_all_regions(
        self,
        sequence: str,
        gene_position: Tuple[int, int],
        rloops: List[bed.GenomicRegion],
    ) -> Iterator[Tuple[List[str], List[str], List[str], List[str]]]:
        for rloop in rloops:
            yield self.__extract_regions(sequence, rloop)

    @staticmethod
    def threshold_region_summaries_shannon(
        region_summaries: List[RegionWeightSummary],
    ) -> List[RegionWeightSummary]:
        thresholded_region_summaries = [
            RegionWeightSummary(i, []) for i in range(0, len(region_summaries))
        ]

        for i, region_summary in enumerate(region_summaries):
            count = 1
            expected_surprise = 0.0
            max_weight = 0.0
            prev_entropy = -math.inf

            for entry in region_summary.entries:
                weight = entry.weight
                if count == 1 or max_weight == 0:
                    max_weight = weight

                p = weight / max_weight
                surprise = -p * math.log(p, 10)
                expected_surprise += surprise
                entropy = expected_surprise / count

                if entropy >= prev_entropy:
                    thresholded_region_summaries[i].entries.append(entry)
                    prev_entropy = entropy
                else:
                    break

                count += 1

        return thresholded_region_summaries

    def __divide_regions(self, sequence: str, index: int):
        before_region = sequence[index - self.ntuple_size - self.padding_size : index]
        after_region = sequence[index : index + self.ntuple_size + self.padding_size]

        before_region_sliding_windows = [
            before_region[i : i + self.ntuple_size]
            for i in range(len(before_region) - self.ntuple_size + 1)
        ]

        after_region_sliding_windows = [
            after_region[i : i + self.ntuple_size]
            for i in range(len(after_region) - self.ntuple_size + 1)
        ]

        return (before_region_sliding_windows, after_region_sliding_windows)

    def __extract_regions(self, sequence: str, rloop: bed.GenomicRegion):
        before_start, after_start = self.__divide_regions(sequence, rloop.start)

        before_end, after_end = self.__divide_regions(sequence, rloop.end)

        return (before_start, after_start, before_end, after_end)
