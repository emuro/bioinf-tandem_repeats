import pandas as pd
import sys
#import plotnine as p9
from plotnine import *
from scipy import stats



# proteome
#
ratio_file ="/Volumes/Wes/results/tandem_repeats/ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640.tsv"

df_prot_ratios = pd.read_csv(ratio_file,sep="\t")
df_prot_ratios.sort_values(by=["bytes"], inplace=True)
if 0:
    print(ratio_file)
    pd.options.display.max_columns = None
    print(df_prot_ratios.head(9))
    print(df_prot_ratios.columns.to_list())
    print(df_prot_ratios.shape)
    sys.exit()



df2plot = df_prot_ratios.copy()
#
cond_sp = df_prot_ratios["sp_or_tr"] == "sp"
df2plot = df2plot[cond_sp]



#
p = (ggplot(df2plot, aes("bytes", "gzip_ratio", color="has_TR")) + geom_point(size=0.2)
     + labs(title="Human (compressing protein seq files)")
     + xlim([800,900])
     + ylim([0.5, 0.625])
     + xlab("bytes (protein seq file)")
     + ylab("(ratio: gzip_file/file)")
    ) +  theme(legend_position=(0.7,0.8), legend_key_size=3, legend_background=element_rect(fill='grey', alpha=0.01))
if 1:
    print(p)




