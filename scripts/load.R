library(Seurat)
library(tidyverse)

load("data/pdac_ctype.RData")
load("data/pdac_siglec.RData")

str(total5)
glimpse(total5)
levels(total5@active.ident) # Cell typeの可視化

glycan_data <- total5@assays$glycan
rna_data <- total5@assays$RNA
cell_types <- total5@active.ident

# 発現マトリックスを取得
glycan_counts <- total5[["glycan"]]$counts

glycan_data_normalized <- GetAssayData(total5, assay = "glycan", slot = "data")  # 正規化データ
glycan_data_scaled <- GetAssayData(total5, assay = "glycan", slot = "scale.data")  # スケーリング済みデータ

tibble(glycan_data)

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



# 特異性を計算する関数 - 1つの細胞種での高発現を重視
calculate_specificity <- function(values_by_celltype) {
  # 基本の変動係数
  cv <- sd(values_by_celltype, na.rm = TRUE) / mean(values_by_celltype, na.rm = TRUE)
  
  # 最大値と次に高い値の比率を計算
  sorted_values <- sort(values_by_celltype, decreasing = TRUE, na.last = TRUE)
  if (length(sorted_values) >= 2 && !is.na(sorted_values[2]) && sorted_values[2] > 0) {
    max_to_second_ratio <- sorted_values[1] / sorted_values[2]
  } else {
    max_to_second_ratio <- ifelse(is.na(sorted_values[1]) || sorted_values[1] == 0, 
                                 1, Inf)
  }
  
  # 最大値と平均値の比率
  max_to_mean_ratio <- sorted_values[1] / mean(values_by_celltype, na.rm = TRUE)
  
  # 複合スコア：単一細胞種特異性を重視
  # max_to_second_ratioが高いほど、1つの細胞種に特異的
  specificity_score <- cv * (0.3 + 0.7 * min(max_to_second_ratio / 5, 1))
  
  return(specificity_score)
}

# がん細胞と免疫細胞の特異性スコアを計算する関数 - 単一細胞種特異性を重視
calculate_cancer_immune_specificity <- function(values_by_celltype, cancer_cells, immune_cells) {
  # がん細胞と免疫細胞の両方が存在するか確認
  if (!all(c(immune_cells, cancer_cells) %in% names(values_by_celltype))) {
    return(calculate_specificity(values_by_celltype))
  }
  
  # がん細胞の中での特異性を計算
  cancer_values <- values_by_celltype[cancer_cells]
  cancer_spec_score <- calculate_specificity(cancer_values)
  
  # 免疫細胞の中での特異性を計算
  immune_values <- values_by_celltype[immune_cells]
  immune_spec_score <- calculate_specificity(immune_values)
  
  # がん細胞で最も高い細胞種とその値
  if (length(cancer_values) > 0) {
    max_cancer_cell <- names(cancer_values)[which.max(cancer_values)]
    max_cancer_value <- max(cancer_values, na.rm = TRUE)
  } else {
    max_cancer_cell <- NULL
    max_cancer_value <- 0
  }
  
  # 免疫細胞で最も高い細胞種とその値
  if (length(immune_values) > 0) {
    max_immune_cell <- names(immune_values)[which.max(immune_values)]
    max_immune_value <- max(immune_values, na.rm = TRUE)
  } else {
    max_immune_cell <- NULL
    max_immune_value <- 0
  }
  
  # 単一細胞種特異性評価
  # 1. がん細胞または免疫細胞の中で1つの細胞種だけが高発現
  if (max_cancer_value > 0 && max_immune_value > 0) {
    max_to_second_overall <- max(max_cancer_value, max_immune_value) / 
                            min(max_cancer_value, max_immune_value)
  } else if (max_cancer_value > 0) {
    # 免疫細胞での発現がない場合
    max_to_second_overall <- Inf
  } else if (max_immune_value > 0) {
    # がん細胞での発現がない場合
    max_to_second_overall <- Inf
  } else {
    max_to_second_overall <- 1
  }
  
  # 最終的な特異性スコア
  # がん細胞内特異性、免疫細胞内特異性、全体での単一細胞種特異性を組み合わせる
  # 単一細胞種特異性に高い重みづけ
  combined_spec_score <- (0.2 * cancer_spec_score + 
                         0.2 * immune_spec_score + 
                         0.6 * min(max_to_second_overall / 3, 1))
  
  return(combined_spec_score)
}

# 単一細胞種での高発現を特定する関数
identify_top_cell_type <- function(values_by_celltype) {
  # 最大値を持つ細胞種
  max_cell_type <- names(values_by_celltype)[which.max(values_by_celltype)]
  max_value <- max(values_by_celltype, na.rm = TRUE)
  
  # 二番目に高い値
  sorted_values <- sort(values_by_celltype, decreasing = TRUE, na.last = TRUE)
  if (length(sorted_values) >= 2 && !is.na(sorted_values[2])) {
    second_max <- sorted_values[2]
    max_to_second_ratio <- max_value / second_max
  } else {
    second_max <- NA
    max_to_second_ratio <- Inf
  }
  
  # 平均値との比較
  mean_value <- mean(values_by_celltype, na.rm = TRUE)
  max_to_mean_ratio <- max_value / mean_value
  
  # 結果を返す
  return(list(
    top_cell_type = max_cell_type,
    max_value = max_value,
    max_to_second_ratio = max_to_second_ratio,
    max_to_mean_ratio = max_to_mean_ratio
  ))
}

# メイン関数のscore_glycansを修正して、単一細胞種特異性の情報を追加
score_glycans <- function(seurat_obj, 
                         weight_cell_spec_rna = 9, 
                         weight_expr_level_rna = 7, 
                         weight_cell_spec_glycan = 6, 
                         weight_binding_level = 7) {
  
  # 細胞型情報を取得
  cell_types <- as.character(Idents(seurat_obj))
  
  # データを取得
  rna_data <- seurat_obj[["RNA"]]$data
  glycan_data <- seurat_obj[["glycan"]]$data
  glycan_list <- rownames(glycan_data)
  
  # 結果を格納するデータフレーム
  results <- data.frame(
    glycan = glycan_list, 
    rna_spec_score = 0, 
    rna_expr_score = 0,
    glycan_spec_score = 0,
    glycan_binding_score = 0,
    total_score = 0,
    top_cell_type_glycan = "",
    max_to_second_ratio_glycan = 0,
    top_cell_type_rna = "",
    max_to_second_ratio_rna = 0,
    stringsAsFactors = FALSE
  )
  
  # 細胞型の分類（実際のデータに合わせて調整）
  immune_cells <- c("M1-like TAM", "M2-like TAM", "Dendritic cells") # 例
  cancer_cells <- c("Basal-Like", "Intermediate", "Classical") # 例
  
  # 最大発現値を計算（正規化）
  max_expr_value <- max(rna_data, na.rm = TRUE)
  
  # 各glycanに対してスコアを計算
  for (g in 1:length(glycan_list)) {
    glycan <- glycan_list[g]
    
    # glycanのバインディングデータを取得
    glycan_binding <- glycan_data[glycan, ]
    
    # 細胞型ごとのglycanバインディングを計算
    glycan_by_celltype <- calculate_mean_by_celltype(glycan_binding, cell_types)
    
    # 最も高いバインディングを示す細胞種を特定
    top_glycan_info <- identify_top_cell_type(glycan_by_celltype)
    
    # glycanスコアを計算
    glycan_scores <- calculate_glycan_scores(
      glycan_binding, cell_types, cancer_cells, immune_cells,
      weight_cell_spec_glycan, weight_binding_level
    )
    
    # 対応するRNA遺伝子の発現スコアを計算
    gene_for_glycan <- glycan # 実際のマッピングに置き換え
    top_rna_info <- list(top_cell_type = NA, max_to_second_ratio = NA)
    
    # 対応する遺伝子がRNA assayに存在するか確認
    if (gene_for_glycan %in% rownames(rna_data)) {
      gene_expr <- rna_data[gene_for_glycan, ]
      
      # 細胞型ごとの遺伝子発現を計算
      expr_by_celltype <- calculate_mean_by_celltype(gene_expr, cell_types)
      
      # 最も高い発現を示す細胞種を特定
      top_rna_info <- identify_top_cell_type(expr_by_celltype)
      
      rna_scores <- calculate_rna_scores(
        gene_expr, cell_types, weight_cell_spec_rna, weight_expr_level_rna, max_expr_value
      )
      
      results$rna_spec_score[g] <- rna_scores$spec_score
      results$rna_expr_score[g] <- rna_scores$expr_score
    }
    
    # glycanスコアを格納
    results$glycan_spec_score[g] <- glycan_scores$spec_score
    results$glycan_binding_score[g] <- glycan_scores$binding_score
    
    # 単一細胞種特異性の情報を格納
    results$top_cell_type_glycan[g] <- top_glycan_info$top_cell_type
    results$max_to_second_ratio_glycan[g] <- top_glycan_info$max_to_second_ratio
    results$top_cell_type_rna[g] <- ifelse(is.na(top_rna_info$top_cell_type), "", top_rna_info$top_cell_type)
    results$max_to_second_ratio_rna[g] <- ifelse(is.na(top_rna_info$max_to_second_ratio), 0, top_rna_info$max_to_second_ratio)
    
    # 総合スコア計算 - 単一細胞種特異性を重視
    single_cell_specificity_bonus <- min(top_glycan_info$max_to_second_ratio / 3, 1) * 2
    results$total_score[g] <- results$rna_spec_score[g] + 
                              results$rna_expr_score[g] + 
                              results$glycan_spec_score[g] + 
                              results$glycan_binding_score[g] +
                              single_cell_specificity_bonus
  }
  
  # 総合スコアで降順にソート
  results <- results %>% arrange(desc(total_score))
  
  return(results)
}

# 拡張された結果の可視化
plot_top_glycans_extended <- function(glycan_ranking, n = 10) {
  top_n_glycans <- head(glycan_ranking, n)
  
  # 基本スコアの可視化
  p1 <- plot_top_glycans(top_n_glycans)
  
  # 細胞種特異性情報の表
  specificity_data <- top_n_glycans %>%
    select(glycan, top_cell_type_glycan, max_to_second_ratio_glycan, 
           top_cell_type_rna, max_to_second_ratio_rna)
  
  # 表をgrid形式で表示
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    grid_table <- gridExtra::tableGrob(specificity_data)
    gridExtra::grid.arrange(p1, grid_table, nrow = 2, 
                           heights = c(2, 1))
  } else {
    print(p1)
    print(specificity_data)
  }
}
