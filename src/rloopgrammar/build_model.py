import multiprocessing
import sys
import os
import pathlib
import dataclasses
import configparser
import math
import random
import logging
import argparse
import shutil

from typing import *

NUMBER_OF_PROCESSES = 10

import rloopgrammar.model.regions_extractor as region_extractor
import rloopgrammar.model.regions_threshold as region_threshold
import rloopgrammar.model.training_set as training_set
import rloopgrammar.model.grammar_dict as grammar_dict
import rloopgrammar.model.grammar_word as grammar_word
import rloopgrammar.model.grammar_training as grammar_training

from rloopgrammar.config_reader import read_plasmids
from rloopgrammar.config_reader import Plasmid

CONFIG_MODEL_PARAMETER_NAME = "Model Parameters"

PROGRAM_NAME = pathlib.Path(sys.argv[0]).parts[-1][:-4]
DESCRIPTION = "Build a collection of R-loop grammar models."


def build_output_folder_name(output_folder_name, **kwargs):
    for key, value in kwargs.items():
        output_folder_name = output_folder_name.replace(f"%{key}", str(value))

    return pathlib.Path(output_folder_name)


class SupressOutput:
    def __init__(self):
        pass

    def __enter__(self):
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, type, value, traceback):
        sys.stdout = sys.__stdout__


@dataclasses.dataclass
class ModelParameters:
    parent_folder: pathlib.Path
    run_number: int
    plasmid: Plasmid
    window_length: int
    padding_length: int
    training_set_percent: float
    seed_file: Optional[pathlib.Path]
    training_set_file: Optional[pathlib.Path]


def build_model(mp: ModelParameters) -> None:
    run_folder = (
        mp.parent_folder
        / f"Model_{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}"
    )

    try:
        os.mkdir(run_folder)
    except FileExistsError:
        pass

    logging.basicConfig(
        filename=mp.parent_folder / "model_build_log.txt",
        format=f"%(asctime)s [{mp.run_number}] - %(message)s",
        level=logging.DEBUG,
    )
    logger = logging.getLogger("r-loop_grammar")

    logger.info(dataclasses.asdict(mp))
    logger.info(f"Extracting critical regions.")

    if mp.seed_file:
        with open(mp.seed_file, "rb") as seed_file_handle:
            seed_bytes = seed_file_handle.read()
    else:
        seed_bytes = os.urandom(32)

    int.from_bytes(seed_bytes, "little")
    seed = seed_bytes
    random.seed(seed)

    seed_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_SEED.b"
    )

    with open(seed_filename, "wb") as seed_file:
        seed_file.write(seed)

    bed_extra_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}.bed_extra.bed"
    )
    weight_xlsx_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_weight.xlsx"
    )
    weight_shannon_entropy_xlsx_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_weight_shannon.xlsx"
    )
    bed_extra_training_set_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}.bed_extra_training-set.bed"
    )
    dict_shannon_xlsx_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_DICT_SHANNON.xlsx"
    )
    dict_shannon_json_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_DICT_SHANNON.xlsx.json"
    )
    training_set_words_filename = str(
        run_folder
        / f"{mp.plasmid.name}_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_training-set_WORDS_SHANNON.txt"
    )
    probabilities_filename = str(
        run_folder
        / f"{mp.plasmid.name}_SHANNON_p{mp.padding_length}_w{mp.window_length}_{mp.run_number}_probabilities.json"
    )

    with open(mp.plasmid.bed_file, "r") as bed_file_fd:
        bed_file_length = len(bed_file_fd.readlines())

    training_set_size = math.ceil(bed_file_length * (mp.training_set_percent / 100.0))
    logger.info("Training set size: {training_set_size}.")

    if not mp.training_set_file:
        with SupressOutput():
            region_extractor.RegionsExtractor.extract_regions(
                mp.plasmid.fasta_file,
                mp.plasmid.bed_file,
                mp.plasmid.gene_start,
                mp.plasmid.gene_end,
                window_length=mp.window_length,
                num_regions=4,
                padding=mp.padding_length,
                bed_extra=True,
                bed_extra_output=bed_extra_filename,
                create_weights=False,
            )

        logger.info("Creating training set.")
        training_set.TrainingSet.training_set(
            bed_extra_filename, training_set_size, bed_extra_training_set_filename
        )
    else:
        logger.info("Duplicating training set.")
        shutil.copyfile(mp.training_set_file, bed_extra_training_set_filename)

    with SupressOutput():
        region_extractor.RegionsExtractor.extract_regions(
            mp.plasmid.fasta_file,
            bed_extra_training_set_filename,
            mp.plasmid.gene_start,
            mp.plasmid.gene_end,
            window_length=mp.window_length,
            out_pref=weight_xlsx_filename,
            num_regions=4,
            padding=mp.padding_length,
            bed_extra=False,
            create_weights=True,
        )

    logger.info("Thresholding critical regions.")

    with SupressOutput():
        region_threshold.RegionsThreshold.extract_regions(
            weight_xlsx_filename,
            weight_shannon_entropy_xlsx_filename,
            True,  # Shannon Entropy
        )

    logger.info("Creating dictionary.")

    with SupressOutput():
        grammar_dict.GrammarDict.extract_regions(
            mp.plasmid.fasta_file,
            bed_extra_training_set_filename,
            weight_xlsx_filename,
            mp.plasmid.gene_start,
            mp.plasmid.gene_end,
            mp.window_length,
            weight_shannon_entropy_xlsx_filename,
            dict_shannon_xlsx_filename,
        )

    logger.info("Extracting training set words.")

    with SupressOutput():
        grammar_word.GrammarWord.extract_word(
            mp.plasmid.fasta_file,
            bed_extra_training_set_filename,
            dict_shannon_json_filename,
            mp.plasmid.gene_start,
            mp.plasmid.gene_end,
            mp.window_length,
            training_set_words_filename,
        )

    logger.info("Finding probabilities.")
    grammar_training.GrammarTraining.find_probabilities(
        training_set_words_filename, mp.window_length, probabilities_filename
    )


parser = argparse.ArgumentParser(prog=PROGRAM_NAME, description=DESCRIPTION)
parser.add_argument("output_folder")
parser.add_argument('-i', '--ini_file', type=str)
parser.add_argument("-c", "--count", type=int, default=10)
parser.add_argument("-p", "--paddings", type=int, nargs="+", default=[13])
parser.add_argument("-k", "--width", type=int, default=4)
parser.add_argument("--plasmid", type=str)
parser.add_argument("-sp", "--sampling_precent", type=int, default=10)
parser.add_argument(
    "-d",
    "--duplicate",
    help="Duplicate the seed from another set of runs, this will override the count.",
)


def main() -> None:
    args = parser.parse_args()
    print(args)

    plasmids = read_plasmids(args.ini_file)

    window_length = args.width
    number_of_models = args.count
    paddings = args.paddings
    model_plasmid_names = [args.plasmid]
    training_set_percent = args.sampling_precent

    if args.duplicate:
        duplicate_collection_model_folders = [x[1] for x in os.walk(args.duplicate)][0]
        number_of_models = len(duplicate_collection_model_folders)

        duplicate_collection_model_folders_sorted = sorted(
            duplicate_collection_model_folders, key=lambda x: int(x.split("_")[-1])
        )

        get_model_files = lambda x: next(os.walk(pathlib.Path(args.duplicate) / x))[2]
        model_files_find = lambda y, z: list(filter(lambda x: y in x, z))[0]
        duplicate_training_set_instead = False

        try:
            seed_files = [
                pathlib.Path(args.duplicate)
                / x
                / model_files_find("SEED", get_model_files(x))
                for x in duplicate_collection_model_folders_sorted
            ]
            get_run_seed_file = lambda x: seed_files[x]
            get_training_set_file = lambda x: None
            print("Copying seed.")

        except IndexError:
            training_set_files = [
                pathlib.Path(args.duplicate)
                / x
                / model_files_find("bed_extra_training-set.bed", get_model_files(x))
                for x in duplicate_collection_model_folders_sorted
            ]
            get_run_seed_file = lambda x: None
            get_training_set_file = lambda x: training_set_files[x]
            print("Copying training set.")

    else:
        get_run_seed_file = lambda x: None
        get_training_set_file = lambda x: None

    model_plasmids: list[Plasmid] = [
        list(filter(lambda x: x.name == k, plasmids))[0] for k in model_plasmid_names
    ]

    for plasmid in model_plasmids:
        for padding in paddings:
            parent_folder = build_output_folder_name(
                args.output_folder,
                plasmid=plasmid.name,
                padding=padding,
                width=window_length,
                number_of_models=number_of_models,
            )
            os.mkdir(parent_folder)

            run_config = configparser.ConfigParser()
            run_config[CONFIG_MODEL_PARAMETER_NAME] = {}
            run_config[CONFIG_MODEL_PARAMETER_NAME]["Plasmid"] = plasmid.name
            run_config[CONFIG_MODEL_PARAMETER_NAME]["WindowLength"] = str(window_length)
            run_config[CONFIG_MODEL_PARAMETER_NAME]["Padding"] = str(padding)
            run_config[CONFIG_MODEL_PARAMETER_NAME]["NumberOfModels"] = str(
                number_of_models
            )

            with open(parent_folder / "model_settings.ini", "w") as configfile:
                run_config.write(configfile)

        runs = [
            ModelParameters(
                parent_folder,
                run_number,
                plasmid,
                window_length,
                padding,
                training_set_percent,
                get_run_seed_file(run_number),
                get_training_set_file(run_number),
            )
            for run_number in range(number_of_models)
            for padding in paddings
        ]

        with multiprocessing.Pool(NUMBER_OF_PROCESSES) as pool:
            pool.map(build_model, runs)


if __name__ == "__main__":
    main()
