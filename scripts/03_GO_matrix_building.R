# 03_GO_matrix_building.R
# ------------------------
# 03_GO_matrix_building.R
# This script defines a function build_GO_matrix() that constructs binary gene-by-GO-term matrices for GO enrichment results from clusterProfiler::enrichGO. 
# It supports GO categories including Biological Process (BP), Cellular Component (CC), and Molecular Function (MF). The output can be used for downstream scoring (e.g., GOscore), heatmaps, or clustering analyses.
# Requires: clusterProfiler enrichment result and gene expression data rownames
# Author: complement1205

# ------------------------
# Load required package
library(stringr)

# ------------------------
# Function: build_GO_matrix
# ------------------------
# Inputs:
#   ego_result:   enrichGO result object from clusterProfiler (ego_BP, ego_CC, ego_MF, etc.)
#   expr_data:    expression matrix or data frame with gene symbols as rownames
# Output:
#   A binary matrix (data.frame) with genes as rows and GO terms as columns (0/1)


build_GO_matrix <- function(ego_result, expr_data) {
  # Step 1: extract gene lists per GO term
  gene_list <- ego_result$geneID %>% str_split("/")
  
  # Step 2: create empty matrix matching gene rownames
  df <- matrix(ncol = 1, nrow = nrow(expr_data))
  rownames(df) <- rownames(expr_data)
  df <- as.data.frame(df)
  
  # Step 3: populate matrix with 1/0 values
  for (i in seq_along(ego_result$ID)) {
    df[, i + 1] <- ifelse(rownames(expr_data) %in% gene_list[[i]], "1", "0")
  }
  
  # Step 4: clean up first column and set column names
  df <- df[, -1]
  colnames(df) <- ego_result$ID
  
  return(df)
}

# ------------------------
# Example usage (uncomment and modify paths as needed)
# ------------------------

# expr_data <- read.csv("data/1205_allsymbol_rawcount.csv", row.names = 1)

# load GO enrichment files
# load("results/ego_BP.RData")  # should contain object 'ego_BP'
# load("results/ego_CC.RData")  # should contain object 'ego_CC'
# load("results/ego_MF.RData")  # should contain object 'ego_MF'
# load("results/all_GO.RData")  # should contain object 'all_GO'


# df_bp <- build_GO_matrix(ego_BP, expr_data)
# df_cc <- build_GO_matrix(ego_CC, expr_data)
# df_mf <- build_GO_matrix(ego_MF, expr_data)
# df_all <- build_GO_matrix(all_GO, expr_data)

# write.csv(df_bp, "output/GO_matrix_BP.csv")
# write.csv(df_cc, "output/GO_matrix_CC.csv")
# write.csv(df_mf, "output/GO_matrix_MF.csv")
# write.csv(df_all, "output/GO_matrix_all.csv")