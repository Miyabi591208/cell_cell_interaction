# load packages -----------------------------------------------------------
# BiocManagerパッケージがインストールされていない場合はインストール
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# biomaRtパッケージのインストール
BiocManager::install("biomaRt")
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
require(Seurat)

# define function -----------------------------------------------------------
shaping_fun <- function(obj, name_, out_path) {
  # shaping
  alldata <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  alldata@active.assay <- "RNA"
  allgenes <- rownames(alldata)
  matrix1 <- as.data.frame(alldata@assays$RNA@data)
  matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]]) != 0,]
  
  ### If you are using a mouse data, then its needed to convert the gene names to human orthologs
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
  
  # Error handling for biomaRt
  tryCatch({
    genesV2 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                     filters = "hgnc_symbol", 
                     values = rownames(alldata@assays$RNA@data), 
                     mart = human, uniqueRows = T)
  }, error = function(e) {
    message("Error in biomaRt query: ", e$message)
    stop("Failed to retrieve gene annotations.")
  })
  
  # matrix1 の結合とフィルタリング
  matrix1$hgnc_symbol <- rownames(matrix1)
  matrix1_joined <- merge(matrix1, genesV2, by = "hgnc_symbol", all.x = TRUE, all.y = FALSE)
  matrix1_joined <- matrix1_joined[!is.na(matrix1_joined$ensembl_gene_id), ]
  matrix1_joined <- matrix1_joined %>% distinct(hgnc_symbol, .keep_all = TRUE)
  rownames(matrix1_joined) <- matrix1_joined$hgnc_symbol
  matrix1_joined <- matrix1_joined[, colnames(matrix1_joined) != "hgnc_symbol"]
  
  # meta.data の cell_abbr を factor から numeric に変換
  # alldata@meta.data$cell_abbr <- as.numeric(as.factor(alldata@meta.data$cell_abbr))
  metadata <- data.frame(cells = colnames(matrix1), cluster = alldata@meta.data[colnames(matrix1), "cell_abbr"])
  
  # 出力の保存
  write.table(matrix1_joined, file.path(out_path, paste0(name_, "_filtered_hcount.txt")), row.names = TRUE, sep = "\t")
  write.table(metadata, file.path(out_path, paste0(name_, "_filtered_meta.txt")), row.names = FALSE, sep = "\t")
}

# main -----------------------------------------------------------
# obtain argument
args <- commandArgs(trailingOnly = TRUE)
obj_path <- args[1]
proj_name <- args[2]
out_path <- args[3]
# obj_path <- "/Users/hitoshisaito/Bioinfo_Project/scRNAseq_GAPFREE2/Data/GF2_N_annotation.rds"
# proj_name <- "GAPFREE2"
# out_path <- "/Users/hitoshisaito/Bioinfo_Project/scRNAseq_GAPFREE2/script/Output"

# load seurat object
seu_obj <- readRDS(obj_path)
seu_obj <- UpdateSeuratObject(object = seu_obj)

# CNV分類追加
seu_obj@meta.data$CNV <- case_when(
  seu_obj@meta.data$case == "NOR1" ~ "Nor",
  seu_obj@meta.data$case == "NOR2" ~ "Nor",  
  seu_obj@meta.data$case == "NOR3" ~ "Nor",
  seu_obj@meta.data$case == "PRE1" ~ "Mid",
  seu_obj@meta.data$case == "PRE2" ~ "High",
  seu_obj@meta.data$case == "PRE3" ~ "High",
  seu_obj@meta.data$case == "PRE4" ~ "Low",
  seu_obj@meta.data$case == "PRE5" ~ "High",
  seu_obj@meta.data$case == "POST1" ~ "Low",
  seu_obj@meta.data$case == "POST2" ~ "Mid",
  seu_obj@meta.data$case == "POST3" ~ "Mid"
)

seu_obj@meta.data <- seu_obj@meta.data %>% 
  mutate(CNV2 = case_when(CNV == "Mid" ~ "Low",
                          CNV == "Nor" ~ "Low",
                          TRUE ~ as.character(CNV)),
         CNV3 = case_when(CNV == "Mid" ~ "Low",
                          TRUE ~ as.character(CNV))
  )

# Categoryの追加
seu_obj@meta.data$Category <- case_when(
  seu_obj@meta.data$case == "NOR1" ~ "NOR",
  seu_obj@meta.data$case == "NOR2" ~ "NOR",  
  seu_obj@meta.data$case == "NOR3" ~ "NOR",
  seu_obj@meta.data$case == "PRE1" ~ "PRE",
  seu_obj@meta.data$case == "PRE2" ~ "PRE",
  seu_obj@meta.data$case == "PRE3" ~ "PRE",
  seu_obj@meta.data$case == "PRE4" ~ "PRE",
  seu_obj@meta.data$case == "PRE5" ~ "PRE",
  seu_obj@meta.data$case == "POST1" ~ "POST",
  seu_obj@meta.data$case == "POST2" ~ "POST",
  seu_obj@meta.data$case == "POST3" ~ "POST"
)

seu_obj@meta.data$cell_abbr <- seu_obj@active.ident
# cell_abbr列をfactor型からcharacter型に変換
seu_obj@meta.data[["cell_abbr"]] <- as.character(seu_obj@meta.data[["cell_abbr"]])

print(levels(Idents(seu_obj)))
# 必要な idents のサブセットを取る
seu_obj <- subset(x = seu_obj, idents = c("T cell","Neutrophil","Epithelial cell","B/plasma cell","Macrophage",
                                         "Endothelial cell","Dendritic cell","Mast cell","Fibroblast","NK cell"))
# seu_obj <- subset(x = seu_obj, idents = c("Macrophage", "Epithelial cell"))

# run 
shaping_fun(seu_obj, proj_name, out_path)
