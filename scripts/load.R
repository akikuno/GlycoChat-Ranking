library(Seurat)
library(tidyverse)

ls()

load("data/pdac_ctype.RData")
load("data/pdac_siglec.RData")

str(total5)

total4


rownames(total5) %>% head()
colnames(total5) %>% head()


rownames(total4) %>% head()
colnames(total4) %>% head()

total4@assays$glycan


# RNAアッセイのカウント行列を取得
rna_counts <- total5@assays$RNA@layers$counts
names(total5@assays$RNA@layers)
dim(rna_counts)
print(rna_counts[1:5, 1:5]) # 最初の5行と5列を表示
