import pandas as pd
import sys
#import plotnine as p9
from plotnine import *
from scipy import stats
if 0:
    print(ggplot.__doc__)


# proteins
#
###stat_file = "/Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/proteins/"
###stat_file += "stat_description.protein.uniprot_reference_proteome.tsv"
stat_file = "/Volumes/Wes/results/geneLength/outputInputFiles/analysis/some_statistics/stat_description/taxid_merged/stat_description.taxid_merged.ensembl_and_ref_proteome.tsv"

df_prot_stat_desc = pd.read_csv(stat_file, sep="\t")
if 0:
    print(stat_file)
    pd.options.display.max_columns = None
    print(df_prot_stat_desc.head(9))
    print(df_prot_stat_desc.columns.to_list())
    print(df_prot_stat_desc.shape)
    sys.exit()

# Tandem Repeats (TR)
#
TR_file = "/Volumes/Wes/data/uncompressed/tandem_repeats/rep2_TR.tsv"
df_TR = pd.read_csv(TR_file, sep="\t")
df_TR = df_TR.rename(columns={'taxid': 'tax_id', 'species': 'species_TR', 'proteins': 'proteins_TR'})
df_TR["TR_normalized"] = df_TR["TR"]/df_TR["proteins_TR"]
if 0:
    print(TR_file)
    pd.options.display.max_columns = None
    print(df_TR.head(9))
    print(df_TR.columns.to_list())
    print(df_TR.shape)
    sys.exit()

# filter by tax_id
df_prot_stat_desc = df_prot_stat_desc[df_prot_stat_desc['tax_id'].isin(df_TR["tax_id"].to_list())]
if 0:
    pd.options.display.max_columns = None
    print(df_prot_stat_desc.head(9))
    print(df_prot_stat_desc.columns.to_list())
    print(df_prot_stat_desc.shape)
    sys.exit()


df_merged = pd.merge(df_prot_stat_desc, df_TR, on="tax_id", how='inner')
if 0:
    pd.options.display.max_columns = None
    print(df_merged.head(9))
    print(df_merged.columns.to_list())
    print(df_merged.shape)
    sys.exit()

if 1:
    df_prot_stat_desc_4MA=df_prot_stat_desc[['tax_id','genes_mean','prots_mean']]
    if 0:
        pd.options.display.max_columns=None
        print(df_prot_stat_desc_4MA.head(9))
        print(df_prot_stat_desc_4MA.columns.to_list())
        print(df_prot_stat_desc_4MA.shape)
        sys.exit()
    df_merged_4MA=pd.merge(df_TR,df_prot_stat_desc_4MA,on="tax_id",how='inner')
    df_merged_4MA=df_merged_4MA[['tax_id','species_TR','proteins_TR','TR','genes_mean','prots_mean']]
    df_merged_4MA.to_csv("/Volumes/Wes/data/uncompressed/tandem_repeats/rep2_TR_4MA.tsv", sep="\t", index=False)
    if 1:
        pd.options.display.max_columns = None
        print(df_merged_4MA.head(9))
        print(df_merged_4MA.columns.to_list())
        print(df_merged_4MA.shape)
        sys.exit()
    sys.exit()

df2plot = df_merged.copy()
#
p = (ggplot(df2plot, aes("prots_mean", "TR", color="merged_division_superregnum")) + geom_point(size=1.75)
     + labs(title="Tandem Repeats")
     + xlab("average protein length")
    ) +  theme(legend_position=(0.35, 0.75), legend_key_size=3, legend_background=element_rect(fill='grey', alpha=0.01))
if 1:
    print(p)


#
p = (ggplot(df2plot, aes("genes_mean", "TR", color="merged_division_superregnum")) + geom_point(size=1.75)
     + labs(title="Tandem Repeats")
     + xlab("average gene length")
     + scale_x_log10(breaks=[10 ** power for power in range(6)],
                     limits=[min(df2plot["genes_mean"].to_list())/2, 2*max(df2plot["genes_mean"].to_list())]
                   )
    ) +  theme(legend_position=(0.35,0.75), legend_key_size=3, legend_background=element_rect(fill='grey', alpha=0.01))
if 1:
    print(p)



