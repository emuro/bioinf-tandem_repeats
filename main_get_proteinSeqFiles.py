from lib import constants as c
from lib import EM_biopython_extras as bioEx
import os
from os import listdir
import pandas as pd
import sys
import gzip
from Bio import SeqIO
import pathlib

if 0: # This are files I am going to use as input
    print(c.LOG_STATISTICS_REFERENCE_PROTEOME_FILE)
    print(c.LOG_SOME_STATISTICS_TAXID_MERGED_FILE)
    print(c.TANDEM_REPEATS_DATA_PATH_NAME)

# Proteome:
# reference proteome annotations (all species)
df_proteomes = pd.read_csv(c.LOG_STATISTICS_REFERENCE_PROTEOME_FILE, sep="\t")
if 0:
    print(c.LOG_STATISTICS_REFERENCE_PROTEOME_FILE)
    pd.options.display.max_columns=None
    print(df_proteomes.head(9))
    print(df_proteomes.columns.to_list())
    print(df_proteomes.shape)
    sys.exit()


# For each species.
species_files = []
for file in listdir(c.TANDEM_REPEATS_DATA_PATH_NAME):
    if file.startswith("UP"):
        species_files.append(file)

    if file == "UP000005640_9606.csv": # Only human at the moment
        uniprotID = file.split("_")[0]
        taxID = file.split("_")[1]
        taxID = taxID.replace(".csv", "")
        print(uniprotID + " " + taxID)
    else:
        continue

    # Proteins of the proteome of this organism
    # reference proteome
    #
    # get files
    fa_proteome = df_proteomes[df_proteomes['proteome_id'] == "UP000005640"]
    aux_proteome_file = fa_proteome['uniprot_fasta_file'].tolist()[0]
    proteome_file = os.path.basename(aux_proteome_file)
    proteome_path = aux_proteome_file.split(proteome_file)[0]
    proteome_file_with_path = c.BASE_PATH + c.DATA_PATH_NAME_COMPRESSED + proteome_path + proteome_file
    #
    protein_path = c.OUT_DATA_LOCAL_PATH_ROOT + proteome_path + "proteins/"
    print(proteome_file_with_path)
    print(protein_path)

    # get all the proteome and save each protein sequence (a *.fa file per protein)
    #
    new_dict_row={"db": [],"UniprotID": [],"EntryName": [],"ProteinName": [],
                  "OS": [],"OX": [],"GN": [],"PE": [],"SV": [],"length": []}
    with gzip.open(proteome_file_with_path, "rt") as handle:
        for record in SeqIO.parse(handle,"fasta"):
            db, UniprotID, EntryName, ProteinName, OS, OX, GN, PE, SV = bioEx.parse_uniprot_header(record)
            seq = record.seq
            if 0:
                print(
                    "\ndb:(%s)\tUniprotID:(%s)\tEntryName:(%s)\tProteinName:(%s)" % (db, UniprotID, EntryName, ProteinName))
                print("OS:(%s)\tOX:(%s)\tGN:(%s)\tPE:(%s)\tSV:(%s)" % (OS, OX, GN, PE, SV))
                print("seq:(%s)" % seq)
            out_protein_file = protein_path + UniprotID + "." + db + ".fa"
            if 0:
                print(protein_path)
                print(out_protein_file)
                sys.exit()
            #
            # save the protein sequence file (*.fa)
            pathlib.Path(protein_path).mkdir(parents=True, exist_ok=True)
            f = open(out_protein_file, 'w') # lengths
            f.writelines(seq)
            f.close()

print(len(species_files))
print(species_files)
