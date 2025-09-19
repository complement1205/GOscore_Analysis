## 01_data_preparation.R
# This script loads a gene expression matrix with gene symbols as the first column.
# The matrix should be a CSV file with the following format:
# Gene,W0,W2,W5,W8,W11
# Gnai3,2168,2103,3671,2863,2229
# Pbsn,1,8,2,0,0
# ...

load_expression_matrix <- function(filepath) {
  expr <- read.csv(filepath, stringsAsFactors = FALSE)
  rownames(expr) <- expr[[1]]           # Set gene symbols as rownames
  expr <- expr[, -1]                    # Remove the gene column
  return(expr)
}

convert_to_entrez <- function(expr) {
  library(clusterProfiler)
  library(org.Mm.eg.db)
  
  gene_symbols <- rownames(expr)
  gene_mapping <- bitr(gene_symbols,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Mm.eg.db)
  
  entrez_ids <- gene_mapping$ENTREZID
  return(entrez_ids)
}
# Example usage:
# expr_data <- load_expression_matrix("data/your_matrix.csv")
# entrez_ids <- convert_to_entrez(expr_data)
