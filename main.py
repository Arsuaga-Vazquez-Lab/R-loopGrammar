import bed
import fasta
import model

import random

if __name__ == "__main__":
    with open("pFC53.fa", "r") as fasta_file:
        fasta_reader = fasta.Reader(fasta_file)
        pfc53_plasmid_sequence = next(fasta_reader).sequence

    with open("pFC53_GYRASECR.bed", "r") as bed_file:
        bed_reader = bed.Reader(bed_file)
        rloop_locations = list(bed_reader)

    pfc53_training_set = random.sample(rloop_locations, 10)

    pfc53_model = model.Model(ntuple_size=4, padding_size=13)
    pfc53_model.train(
        pfc53_plasmid_sequence, pfc53_training_set, gene_position=(80, 1829)
    )

    print(pfc53_model.dictionary.rloop_parsing_blocks)

    print(pfc53_model.region_weight_summaries)

    """
    pfc53_model = model.Model(ntuple_size=4, padding_size=13)
    pfc53_model.train(pfc53_sequence, pfc53_training_set)

    combined_model = pfc53_model | pfc8_model
    combined_model = pfc53_model | pfc8_model
    combined_model = pfc53_model | pfc8_model
    
    combined_model = pfc53_model.union(pfc8_model)

    combined_model.predict(...)
    """
