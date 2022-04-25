import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import anndata as ad
import scanpy as sc

adata = ad.read_csv(f"data/iPSC_DESeq_JUL2015_all_merged_ensemble_DataLevel3.csv").T
metadata = pd.read_excel(f"data/RNA-seq iPSC Sample Experimental Metadata_20160911.xlsx")
adata.obs = metadata.set_index('Sample_ID')

sc.pl.highest_expr_genes(adata, n_top=20, )

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)

save_path = f"{BASE_DIR}/nbdata/iPSC_DESeq_JUL2015_all_merged_ensemble_DataLevel3_scpp.h5"
adata.write(save_path)