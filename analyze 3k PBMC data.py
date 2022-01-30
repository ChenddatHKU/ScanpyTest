import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


results_file = '/Users/chendd/data/write/pbmc3k.h5ad'  # the file that will store the analysis results
adata = sc.read_10x_mtx(
    '/Users/chendd/data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file and '.tsv' file. # where to get these formated files?
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
# print(adata)  # I added print() function, since I think this code was run on jupyter notebook
# rows of the anndata object are cell barcode, column name are genes.
# how to view different layer of the data set and their index?
# sc.pl.highest_expr_genes(adata, n_top=20)
#basic filtering # before filtering, how to check the QC of the count matrix?
#optional: plot qc matrixs
# qcmatrix = sc.pp.calculate_qc_metrics(adata)
# plt.subplot(2, 2, 1)
# qcmatrix[0]['total_counts'].plot.hist(bins = 100)
# plt.title('total_counts')

# plt.subplot(2, 2, 2)
# qcmatrix[0]['n_genes_by_counts'].plot.hist(bins = 100)
# plt.title('n_genes_by_counts')

# plt.subplot(2, 2, 3)
# qcmatrix[0]['pct_counts_in_top_50_genes'].plot.hist(bins = 100)
# plt.title('pct_counts_in_top_50_genes')
# plt.show()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)
# #why store the .mtx file into h5 file format, is this the default operation in scRNN-Seq data analysis?
# sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


##mainly, scripts below filter cells based on n_gens and pct_count_mt.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)  #Total-count normalize (library-size correct) the data matrix  to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.log1p(adata)  #Logarithmize the data

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)  # Identify highly-variable genes.
# sc.pl.highly_variable_genes(adata)
adata.raw = adata  # Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) #Regress out (mostly) unwanted sources of variation.
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)

print(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
adata.write(results_file)



# print(adata)
# print(adata.obs)
# print(adata.var)






