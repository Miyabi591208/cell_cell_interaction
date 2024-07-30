library(Seurat)
library(ggplot2)
library(BiocManager)
library(dplyr)
library(viridis)
library(data.table)
library(ggrepel)
library(patchwork)

# Annotation Dataを割り当てる
dat <- readRDS("/home/saito_h/scRNAseq_GAPFREE2/Data/GF2_N_annotation.rds")
dat <- UpdateSeuratObject(object = dat)

# CNV分類追加
dat@meta.data$CNV <- case_when(dat@meta.data$case == "NOR1" ~ "Nor",
                               dat@meta.data$case == "NOR2" ~ "Nor",  
                               dat@meta.data$case == "NOR3" ~ "Nor",
                               dat@meta.data$case == "PRE1" ~ "Mid",
                               dat@meta.data$case == "PRE2" ~ "High",
                               dat@meta.data$case == "PRE3" ~ "High",
                               dat@meta.data$case == "PRE4" ~ "Low",
                               dat@meta.data$case == "PRE5" ~ "High",
                               dat@meta.data$case == "POST1" ~ "Low",
                               dat@meta.data$case == "POST2" ~ "Mid",
                               dat@meta.data$case == "POST3" ~ "Mid"
)


dat@meta.data <- dat@meta.data %>% 
  mutate(CNV2 = case_when(CNV == "Mid" ~ "Low",
                          CNV == "Nor" ~ "Low",
                          TRUE ~ as.character(CNV)),
         CNV3 = case_when(CNV == "Mid" ~ "Low",
                          TRUE ~ as.character(CNV)),
  )

# Categoryの追加
dat@meta.data$Category <- case_when(dat@meta.data$case == "NOR1" ~ "NOR",
                                    dat@meta.data$case == "NOR2" ~ "NOR",  
                                    dat@meta.data$case == "NOR3" ~ "NOR",
                                    dat@meta.data$case == "PRE1" ~ "PRE",
                                    dat@meta.data$case == "PRE2" ~ "PRE",
                                    dat@meta.data$case == "PRE3" ~ "PRE",
                                    dat@meta.data$case == "PRE4" ~ "PRE",
                                    dat@meta.data$case == "PRE5" ~ "PRE",
                                    dat@meta.data$case == "POST1" ~ "POST",
                                    dat@meta.data$case == "POST2" ~ "POST",
                                    dat@meta.data$case == "POST3" ~ "POST"
)



print(head(dat@meta.data))

dat_cnv_high <- subset(dat, subset = CNV2 == "High")
dat_cnv_low <- subset(dat, subset = CNV2 == "Low")

# 必要に応じてサブセット化したデータを保存
saveRDS(dat_cnv_high, file = "/home/saito_h/script/cnv_high.rds")
saveRDS(dat_cnv_low, file = "/home/saito_h/script/cnv_low.rds")