## Synopsis
R-loops are transient three-stranded nucleic acids that form during transcription when the nascent RNA hybridizes to the template DNA, freeing the DNA non-template strand. There is growing evidence that R-loops play important roles in physiological processes such as control of gene expression, while also contributing to chromosomal instability and disease. It is known that R-loop formation is influenced by both the sequence and the topology of the DNA substrate, but many questions remain about how R-loops form and the 3-dimensional structures that they adopt. 

Here we represent an R-loop as a word in a formal grammar called the _R-loop grammar_ and predict R-loop formation. We train the R-loop grammar on experimental data obtained by single-molecule R-loop footprinting and sequencing (SMRF-seq). Despite not containing explicit topological information, the R-loop grammar accurately predicts R-loop formation on plasmids with varying starting topologies and outperforms previous methods in R-loop prediction. 

<!---
Might need a reference for the abstract.
-->

## Experimental data set:
[Fasta & BED files](https://github.com/Arsuaga-Vazquez-Lab/R-loopGrammar/releases/download/v0.0.1-alpha/fasta_bed.zip)

## Dependencies:
- [Python3](https://www.python.org/downloads/)
- [matplotlib](https://pypi.org/project/matplotlib/)
- [openpyxl](https://pypi.org/project/openpyxl/)

_________________

## Installation
You can install directly from GitHub using the command below,
```sh
pip install git+https://github.com/Arsuaga-Vazquez-Lab/R-loopGrammar.git
```
Or you can clone the repo and install from the repository itself.
```sh
git clone https://github.com/Arsuaga-Vazquez-Lab/R-loopGrammar.git
cd R-loopGrammar
pip install .
```
Both of these methods will install the depedencies themselves.

## Usage

After installation, you'll have access to the following programs,
- `rloop-grammar-build-model`
- `rloop-grammar-union-models`
- `rloop-grammar-predict`

These programs depend on a configuration file of plasmids/genes named `plasmids.ini`.

This file will contain the information regarding the plasmids that each program operates on, 
e.g.
```
[Plasmid1]
GeneStart = 80
GeneEnd = 1829
FastaFile = Plasmid1.fa
BEDFile = Plasmid1.bed

[Plasmid2]
GeneStart = 80
GeneEnd = 1512
FastaFile = Plasmid2.fa
BEDFile = Plasmid2.bed
```
This file defines the start and end of the gene inside the plasmid (0-based indexing, start inclusive and end exclusive), the fasta file location, and the BED file containing the experimental R-loop data.

The BED files are of the form,
```
plasmid1_rloop1 85 125
plasmid1_rloop2 67 342
...
```
These will contain the experimental results where the first index and last index are where an R-loop has been experimentally found inside the given plasmid.

### Data splitting
-----

0. Split data as desired into full training data set defined as its own `BEDFile` and the holdout set.

Included in `training_data.zip` is `make_training.py`. This is a script that will separate out any BED files in the working directory into both a training and test set. The seed that was used in generation of the `model_data.zip` is provided.
```
python make_training.py
```
We split the BED file into three folds,
```
+--------+--------+--------+
| Fold 1 | Fold 2 | Fold 3 |
+--------+--------+--------+
```
The folds are then used to build three different training and test sets as follows,
```
		+--------+--------+
Training 1	| Fold 1 | Fold 2 |
		+--------+--------+
		+--------+--------+
Training 2	| Fold 1 | Fold 3 |
		+--------+--------+
		+--------+--------+
Training 3	| Fold 2 | Fold 3 |
		+--------+--------+
```
The `model_data.zip` is built off the Training 1 dataset, the other two are used for validation of results.

### Building the model
-----

1. This command provides a single `plasmid` dictionary with corresponding weights for the `sp` of `BEDFile` using `k` and `p` as k-mer size and padding. `c` is the number of runs for an ensemble. It could also be used to build a model with probabilities if one chooses to use only this dictionary.
```sh
rloop-grammar-build-model Collection_Plasmid1_runs_30 -c 30 -p 13 -k 4 --plasmids Plasmid1 -sp 10
```
* `--plasmid` The plasmid which to build the model from, as defined in `plasmids.ini`.
* `-sp` The sampling percent from the `BEDFile` to use for the accomponing BED file for each plasmid.
* `-c` The desired number of models for the ensemble.
* `-k` The k-mer size.
* `-p` The padding used for the sliding windows in the critical regions.
* `-d` (Optional) To duplicate a run utilizing the same seed or training set, use this option to select a model to copy from; this will override `-c`.

2. This command generates the ensemble of `c` models on the union of `Plasmid1` and `Plasmid2`. To avoid ambiguities we use a `stochastic` or `deterministic` union of the dictionaries.
```sh
rloop-grammar-union-models UnionCollection_Plasmid1_Plasmid2 -m stochastic -i Collection_Plasmid1_runs_30 Collection_Plasmid2_runs_30
```
* `-i` The two input plasmids which are used to build the union model.
* `-m` The method used to take the union, the current supported options are `stochastic` and `deterministic` (default is `stochastic`).

### Generating predictions
-----

3. This command generates an ensemble of `c` predictions on a `plasmid` for each model in the ensemble.
```sh
rloop-grammar-predict UnionCollection_Plasmid1_Plasmid2_predict_on_Plasmid3 -i UnionCollection_Plasmid1_Plasmid2 --plasmid Plasmid3
```
* `-i` The input model used to build the prediction.
* `-plasmid` The plasmid used to predict upon. 

4. Take average of the ensemble of `c` predictions and then graph.
```sh
rloop-grammar-graph-prediction UnionCollection_Plasmid1_Plasmid2_predict_on_Plasmid3 -n Prediction_Plasmid3
```
* `-n` The name displayed on the graph.

## Reproduce model data

If you would like to reproduce the model data found [here](https://github.com/Arsuaga-Vazquez-Lab/R-loopGrammar/releases/tag/v0.0.1-alpha), download the zip file.
We will be using the random seed for each run to recreate the data utilizing the `-d` duplicate option when building the model.

First we need to patition the bed-files utilizing the [training data](https://github.com/Arsuaga-Vazquez-Lab/R-loopGrammar/releases/download/v0.0.1-alpha/training-data.zip) seed.
The `make_training.py` script will partition the BED files into three sepearate sets utilizing the seed file. We will be using the first partition to build our model.

The `plasmids.ini` will expect the BED files to be located in a folder called `bed-files`, also with the FASTA files in `fasta-files`.
You can find the `plasmids.ini` file also in the `training-data.zip`.

The following Makefile can be used to then recreate the model.
```
all: models union predict

models:
	rloop-grammar-build-model new_Collection_pFC53_SUPERCOILEDCR_runs_30 -c 30 -p 13 -w 4 --plasmids pFC53_SUPERCOILEDCR_training -tp 10 -d original_Collection_pFC53_SUPERCOILEDCR_runs_30
	rloop-grammar-build-model new_Collection_pFC8_SUPERCOILEDCR_runs_30 -c 30 -p 13 -w 4 --plasmids pFC8_SUPERCOILEDCR_training -tp 10 -d original_Collection_pFC8_SUPERCOILEDCR_runs_30
	rloop-grammar-build-model new_Collection_pFC53_GYRASECR_runs_30 -c 30 -p 13 -w 4 --plasmids pFC53_GYRASECR_training -tp 10 -d original_Collection_pFC53_GYRASECR_runs_30
	rloop-grammar-build-model new_Collection_pFC8_GYRASECR_runs_30 -c 30 -p 13 -w 4 --plasmids pFC8_GYRASECR_training -tp 10 -d original_Collection_pFC8_GYRASECR_runs_30
	rloop-grammar-build-model new_Collection_pFC53_LINEARIZED_runs_30 -c 30 -p 13 -w 4 --plasmids pFC53_LINEARIZED_training -tp 10 -d original_Collection_pFC53_LINEARIZED_runs_30
	rloop-grammar-build-model new_Collection_pFC8_LINEARIZED_runs_30 -c 30 -p 13 -w 4 --plasmids pFC8_LINEARIZED_training -tp 10 -d original_Collection_pFC8_LINEARIZED_runs_30

union:
	rloop-grammar-union-models new_UnionCollection_SUPERCOILEDCR_runs_30 -m stochastic -i new_Collection_pFC53_SUPERCOILEDCR_runs_30 new_Collection_pFC8_SUPERCOILEDCR_runs_30
	rloop-grammar-union-models new_UnionCollection_GYRASECR_runs_30 -m stochastic -i new_Collection_pFC53_GYRASECR_runs_30 new_Collection_pFC8_GYRASECR_runs_30
	rloop-grammar-union-models new_UnionCollection_LINEARIZED_runs_30 -m stochastic -i new_Collection_pFC53_LINEARIZED_runs_30 new_Collection_pFC8_LINEARIZED_runs_30

predict:
	rloop-grammar-predict new_UnionCollection_SUPERCOILEDCR_predict_on_pFC53 -i new_UnionCollection_SUPERCOILEDCR_runs_30 --plasmids pFC53_SUPERCOILEDCR
	rloop-grammar-predict new_UnionCollection_SUPERCOILEDCR_predict_on_pFC8 -i new_UnionCollection_SUPERCOILEDCR_runs_30 --plasmids pFC8_SUPERCOILEDCR
	rloop-grammar-predict new_UnionCollection_GYRASECR_predict_on_pFC53 -i new_UnionCollection_GYRASECR_runs_30 --plasmids pFC53_GYRASECR
	rloop-grammar-predict new_UnionCollection_GYRASECR_predict_on_pFC8 -i new_UnionCollection_GYRASECR_runs_30 --plasmids pFC8_GYRASECR
	rloop-grammar-predict new_UnionCollection_LINEARIZED_predict_on_pFC53 -i new_UnionCollection_LINEARIZED_runs_30 --plasmids pFC53_LINEARIZED
	rloop-grammar-predict new_UnionCollection_LINEARIZED_predict_on_pFC8 -i new_UnionCollection_LINEARIZED_runs_30 --plasmids pFC8_LINEARIZED

.PHONY: models union predict
```
```
make all
```

## References

1. Stolz R, Sulthana S, Hartono SR, Malig M, Benham CJ, Chedin F. Interplay between DNA sequence and negative superhelicity drives R-loop structures. Proc Natl Acad Sci U S A. 2019 Mar 26;116(13):6260-6269. doi: 10.1073/pnas.1819476116. Epub 2019 Mar 8. PMID: 30850542; PMCID: PMC6442632.
2. Jonoska N, Obatake N, Poznanovik S, Price C, Riehl M, Vazquez M. (2021). Modeling RNA:DNA Hybrids with Formal Grammars. 10.1007/978-3-030-57129-0_3.
3. Ferrari M, Poznanović S, Riehl M, Lusk J, Hartono S, Gonzalez-Isunza G, Chedin F, Vazquez M, Jonoska N (2025). The R-loop grammar predicts R-loop formation under different topological constraints. PLOS Computational Biology, 21(8), e1013376.

## Acknowledgements

* The entire team acknowledges support from the National Science Foundation and the National Institutes of Health, DMS/NIGMS awards #2054347 and #2054321. 
* Mariel Vazquez acknowledges support from NSF grants DMS #1716987 and DMS #1817156
* Svetlana Poznanovic acknowledges support by Simons Foundation gift MP-TSM-00002798 and NSF DMS #1815832
* Margherita Maria Ferrari acknowledges the support of the Natural Sciences and Engineering Research Council of Canada (NSERC) [funding reference numbers DGECR-2023-00131, RGPIN-2023-04722] and the University of Manitoba (research start-up funds) and NJ acknowledges support in part by NSF CCF #2107267 and the W.M. Keck Foundation.
* Frederic Chedin and Stella Hartono were also supported in part by NIH grant R35 GM139549. 
* Margherita Maria Ferrari, Natasa Jonoska, Svetlana Poznanovic, Manda Riehl and Mariel Vazquez thank the Institute of Pure and Applied Mathematics and the Association for Women in Mathematics for seeding this research.
