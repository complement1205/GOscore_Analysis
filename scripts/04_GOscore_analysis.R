# 04_GOscore_analysis.R
# ------------------------------------------------------------------------------
# Title    : GOscore calculation for GO term-wise expression comparison
# Author   : Complement1205
# Purpose  : For each GO term (column of binary matrix), perform t-tests comparing
#            W0 expression vs other time points across genes belonging to that GO term.
# 
# Input:
#   - go_matrix: data.frame or matrix, rows = gene symbols, columns = GO terms (1/0)
#   - expr_matrix: data.frame or matrix, rows = gene symbols, columns = sample expression (e.g., W0~W11)
#       => Row names of both must be gene symbols!
#
# Output:
#   - A data.frame with GO term and p-values of W0 vs W2/W5/W8/W11 comparisons
#
# Dependencies: reshape2
# ------------------------------------------------------------------------------


run_GOscore <- function(go_matrix, expr_data, 
                        verbose = FALSE, 
                        filter_na = TRUE) {
  library(reshape2)
  
  go_matrix[] <- lapply(go_matrix, function(col) as.numeric(as.character(col)))
  
  # Sample column names
  sample_names <- c("W0", "W2", "W5", "W8", "W11")
  lst <- data.frame()
  
  for (i in 1:ncol(go_matrix)) {
    if (verbose) cat("Processing", i, "/", ncol(go_matrix), ":", colnames(go_matrix)[i], "\n")
    
    genes_in_go <- rownames(go_matrix)[go_matrix[, i] == 1]
    go_expr <- expr_data[genes_in_go, sample_names]
    if (nrow(go_expr) == 0) next
    
    # Long format transformation
    go_expr$Gene <- rownames(go_expr)
    df <- melt(go_expr, id.vars = "Gene")
    colnames(df) <- c("Gene", "Group", "Value")
    df$Group <- factor(df$Group, levels = sample_names)
    
    # Perform t-tests
    t1 <- tryCatch(t.test(df$Value[df$Group == "W0"], df$Value[df$Group == "W2"])$p.value, error = function(e) NA)
    t2 <- tryCatch(t.test(df$Value[df$Group == "W0"], df$Value[df$Group == "W5"])$p.value, error = function(e) NA)
    t3 <- tryCatch(t.test(df$Value[df$Group == "W0"], df$Value[df$Group == "W8"])$p.value, error = function(e) NA)
    t4 <- tryCatch(t.test(df$Value[df$Group == "W0"], df$Value[df$Group == "W11"])$p.value, error = function(e) NA)
    
    lst <- rbind(lst, data.frame(
      GO = colnames(go_matrix)[i],
      W0_W2 = t1, W0_W5 = t2, W0_W8 = t3, W0_W11 = t4
    ))
  }
  
  # Optional: filter NA results
  if (filter_na) {
    lst <- lst[complete.cases(lst), ]
  }
  
  return(lst)
}

# ------------------------------------------------------------------------------
# Example usage:
# rownames(expr_data) <- expr_data$X; expr_data <- expr_data[, -1]
# rownames(dfbp) <- dfbp$X; dfbp <- dfbp[, -1]
# bp_result <- run_GOscore(dfbp, expr_data, verbose = TRUE)
# ------------------------------------------------------------------------------

# Save GOscore results for downstream analysis and visualization
# write.csv(bp_result, "results/GOscore_BP_results.csv", row.names = FALSE)
# write.csv(cc_result, "results/GOscore_CC_results.csv", row.names = FALSE)
# write.csv(mf_result, "results/GOscore_MF_results.csv", row.names = FALSE)

cat("GOscore results saved. You may use Excel or other software for custom visualization (e.g., -log10(p.value) ranking and scatterplots).\n")

# ------------------------------------------------------------------------------
# Downstream analysis (not included in this script):
#   - Use GOscore results as input for KOBAS-i to perform secondary enrichment,
#     reducing noise and enhancing biological relevance.
#   - Use cirFunMap to cluster enriched terms and visualize temporal functional landscapes.
#   - Optional manual annotation of clusters (Level 2) with QuickGO.
#
# These steps were applied in our study to reconstruct the CAC progression landscape.
# ------------------------------------------------------------------------------