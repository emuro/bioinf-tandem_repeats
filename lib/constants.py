# python3
# ################################################################## #
# constants.py (C) July-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: tandem repeats
#
# Purpose: Constant to be used in the project
#
# Important:
# Change: DIVISION. There rest adapt to this selection
# ################################################################## #
from datetime import date


DIVISION = "protist"  # "viruses", "bacteria", "protist", "fungi", "plants", "metazoa", "vertebrates"
# ################################################################################################### #

DB = "uniprot"   # ensembl for vertebrates
                        # ensemblgenomes for bacteria,fungi,protist,plants,metazoa and viruses

BOOL_CHECK_SOME_SPECIES = 1  # for debugging
ANNOTATION_NAMES_TO_CHECK = ["Homo_sapiens",
                             ""]  # ["Homo_sapiens", "Danio_rerio"] ["Homo_sapiens"]


BASE_PATH = "/Volumes/Wes/"  # HOME_PATH; Wes; BIRD

PROJECT = "tandem_repeats"

OUT_DATA_LOCAL_PATH_ROOT = BASE_PATH + "results/" + PROJECT + "/"
OUTPUT_INPUT_FILES_PATH = OUT_DATA_LOCAL_PATH_ROOT + "outputInputFiles/"  # or PYCHARM_PROJECT_PATH

DATA_PATH_NAME  = "data/uncompressed/"
DATA_PATH_NAME_COMPRESSED  = "data/compressed/"

# TANDEM_REPEATS DATA
TANDEM_REPEATS_DATA_PATH_NAME = BASE_PATH + DATA_PATH_NAME + PROJECT + "/" + "pre-computed_results_78_proteomes" + "/"

# MODEL_ORGANISMS
SOME_STATISTICS_PATH_NAME = "analysis/some_statistics/stat_description/"
REFERENCE_PROTEOMES_PATH_NAME = "some_tables/reference_proteomes/"
REFERENCE_PROTEOMES_FILE = "reference_proteomes_table_28.5.2021.tsv"

# TAXID
ENSEMBL_TAXID_FILE_NAME = "some_tables/species_Ensembl_taxid/species_Ensembl.tsv"

# LOG FILES
############
LOG_STATISTICS_REFERENCE_PROTEOME_FILE = BASE_PATH + \
    "results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/" + \
    "stat_description.protein.uniprot_reference_proteome.tsv"

LOG_SOME_STATISTICS_TAXID_MERGED_FILE = BASE_PATH + \
    "results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/taxid_merged/" + \
    "stat_description.taxid_merged.ensembl_and_ref_proteome.tsv"

# /Volumes/Wes/data/uncompressed/tandem_repeats/pre-computed_results_78_proteomes