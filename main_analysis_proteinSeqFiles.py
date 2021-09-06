from lib import constants as c
from lib import EM_biopython_extras as bioEx
import os
from os import listdir
import glob
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

        #
        # capture the TR annotations of this organism
        df_TR_of_species = pd.read_csv(c.TANDEM_REPEATS_DATA_PATH_NAME + file, sep=",")
        df_TR_of_species = df_TR_of_species.applymap(lambda x: str(x).strip())
        df_TR_of_species.columns = ['uniprotID', 'sp_or_tr', 'TR_kind', 'aa_start', 'aa_end', 'TR_seq', 'score', 'p_val']
        if 0:
            if 1:
                cond_sp = (df_TR_of_species["sp_or_tr"] == "sp")
                df_TR_of_species = df_TR_of_species[cond_sp]
            else:
                cond_tr = (df_TR_of_species["sp_or_tr"] == "tr")
                df_TR_of_species = df_TR_of_species[cond_tr]
        df_TR_of_species.sort_values(by=["uniprotID"], inplace=True)
        prots_of_df = df_TR_of_species["uniprotID"].unique().tolist()
        if 0:
            print(prots_of_df)
            print(len(prots_of_df))
            print(type(prots_of_df))
        if 0:
            print(c.TANDEM_REPEATS_DATA_PATH_NAME + file)
            pd.options.display.max_columns=None
            print(df_TR_of_species.head(9))
            print(df_TR_of_species.columns.to_list())
            #print(df_TR_of_species.loc[0, ].to_list())
            print(df_TR_of_species.shape)
            sys.exit()
    else:
        continue

    # Proteins of the proteome of this organism
    # reference proteome

    # get files
    fa_proteome = df_proteomes[df_proteomes['proteome_id'] == "UP000005640"]
    aux_proteome_file = fa_proteome['uniprot_fasta_file'].tolist()[0]
    proteome_path = aux_proteome_file.split( os.path.basename(aux_proteome_file) )[0]
    #
    protein_path = c.OUT_DATA_LOCAL_PATH_ROOT + proteome_path + "proteins/"
    print(protein_path)

    os.chdir(protein_path)
    proteins_fa = glob.glob("*.fa")
    proteins_fa.sort()
    if 1:
        #print(proteins_fa)
        print("There are " + str(len(proteins_fa)) + " proteins in " + uniprotID + "/" + taxID)

    df_prot = pd.DataFrame( columns=["uniprotID", "sp_or_tr", "bytes", "gzip", "gzip_ratio", "has_TR"])
    for f in proteins_fa:
        prot = f.split(".")[0]
        db = f.split(".")[1]
        prot_size = os.path.getsize(protein_path+f)  # bytes
        gzip_f = f + ".gz"
        gzip_prot_size = os.path.getsize(protein_path+gzip_f) # bytes
        compress_ratio = gzip_prot_size/prot_size
        has_TR = "no_TR"
        if prot in prots_of_df:
            has_TR = "TR"
        if 0:
            print(prot+ "/"+db+ "/"+str(prot_size)+"/"+str(gzip_prot_size)+"/"+str(compress_ratio)+"/"+str(has_TR))
        #
        df_prot.loc[len(df_prot.index)] = [prot, db, prot_size, gzip_prot_size, compress_ratio, has_TR]
    if 0:
        pd.options.display.max_columns=None
        print(df_prot.head(9))
        print(df_prot.columns.to_list())
        print(df_prot.shape)
        sys.exit()
    df_prot.to_csv(c.OUT_DATA_LOCAL_PATH_ROOT + proteome_path + uniprotID + ".tsv", sep="\t", index=False)

if 0:
    print(len(species_files))
    print(species_files)
