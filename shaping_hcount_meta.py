import pandas as pd
import scanpy as sc
import os
import sys
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# path to the file
path = sys.argv[1]
path_meta = sys.argv[2]
base_path = sys.argv[3]

# Check if files exist
if not os.path.exists(path):
    raise FileNotFoundError(f"The file {path} does not exist.")
if not os.path.exists(path_meta):
    raise FileNotFoundError(f"The file {path_meta} does not exist.")

# Load data
temp = pd.read_table(path).set_index("ensembl_gene_id")
temp.columns = temp.columns.str.replace('-', '_')
temp_meta = pd.read_table(path_meta)
temp_meta = temp_meta.iloc[:-1, :]
temp_meta.columns = ['barcode_sample', 'cell_type']
temp_meta['barcode_sample'] = temp_meta['barcode_sample'].str.replace('-', '_')
temp.to_csv(os.path.join(base_path, "normalised_log_counts.tsv"), sep="\t")
temp_meta.to_csv(os.path.join(base_path, "metadata.tsv"), sep="\t", index=None)

# Paths for CellphoneDB
cpdb_file_path = '/home/saito_h/sc_RNAseq/CrosstalkeR/v5.0.0/cellphonedb.zip'
meta_file_path = os.path.join(base_path, 'metadata.tsv')
counts_file_path_tsv = os.path.join(base_path, 'normalised_log_counts.tsv')
counts_file_path_h5ad = os.path.join(base_path, 'normalised_log_counts.h5ad')
out_path = base_path

# Load meta data
meta_data = temp_meta

# Load counts data
counts_data = temp

# Verify indexes and columns
counts_data.index.name = 'Gene'
if not all(meta_data['barcode_sample'].isin(counts_data.columns)):
    raise ValueError("Some barcodes in meta_data are not present in counts_data columns")

meta_barcodes = set(meta_data['barcode_sample'])
counts_barcodes = set(counts_data.columns)

missing_in_meta = counts_barcodes - meta_barcodes
missing_in_counts = meta_barcodes - counts_barcodes

if missing_in_meta:
    print(f"Missing in meta: {missing_in_meta}")
if missing_in_counts:
    print(f"Missing in counts: {missing_in_counts}")

if counts_data.sum().sum() == 0:
    raise ValueError("All counts are zero. Please check your counts data.")

meta_data = meta_data.set_index('barcode_sample')
counts_data = counts_data.loc[:, meta_data.index]

print(meta_data.head())
print(counts_data.head())

adata = sc.AnnData(X=counts_data.T, obs=meta_data)
print(adata.obs.head())

adata.write(counts_file_path_h5ad)

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_file_path,
    meta_file_path=meta_file_path,
    counts_file_path=counts_file_path_h5ad,
    counts_data='ensembl',
    active_tfs_file_path=None,
    microenvs_file_path=None,
    score_interactions=True,
    iterations=1000,
    threshold=0.1,
    threads=5,
    debug_seed=42,
    result_precision=3,
    pvalue=0.05,
    subsampling=False,
    subsampling_log=False,
    subsampling_num_pc=100,
    subsampling_num_cells=1000,
    separator='|',
    debug=False,
    output_path=out_path,
    output_suffix=None
)

