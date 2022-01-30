import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = '/Users/chendd/write/pbmc3k.h5ad'  # the file that will store the analysis results
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

print(adata[adata.obs.n_genes_by_counts < 2500, :])

# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

