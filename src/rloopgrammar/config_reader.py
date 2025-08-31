import configparser
import dataclasses


@dataclasses.dataclass
class Plasmid:
    name: str
    gene_start: int
    gene_end: int
    fasta_file: str
    bed_file: str


def read_plasmids(ini_filepath):
    config = configparser.ConfigParser()
    config.read(ini_filepath)

    plasmids = []
    plasmid_names = config.sections()

    for plasmid_name in plasmid_names:
        gene_start = int(config[plasmid_name]["GeneStart"])
        gene_end = int(config[plasmid_name]["GeneEnd"])
        fasta_file = config[plasmid_name]["FastaFile"]
        bed_file = config[plasmid_name]["BEDFile"]

        plasmids.append(
            Plasmid(plasmid_name, gene_start, gene_end, fasta_file, bed_file)
        )

    return plasmids
