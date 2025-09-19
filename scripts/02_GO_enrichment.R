# 02_GO_enrichment.R

# Perform GO enrichment analysis for BP, CC, and MF categories
# This script performs Gene Ontology (GO) enrichment analysis using clusterProfiler, based on a gene list derived from input expression data.
#
# It runs enrichGO separately for:
#   - Biological Process (BP)
#   - Cellular Component (CC)
#   - Molecular Function (MF)
#
# Gene symbols from the input file (data/example_counts.csv) are converted to ENTREZ IDs
# before enrichment. Results for each GO category are filtered (minimum gene count = 3)
# and saved separately, as well as combined into one output file.
#
# Output:
#   - results/GO_enrichment_combined.csv
#   - results/GO_enrichment_BP.csv
#   - results/GO_enrichment_CC.csv
#   - results/GO_enrichment_MF.csv


# ----------------------
# GO Enrichment Analysis
# ----------------------

# Function wrapper for enrichGO by ontology
run_GO_enrichment <- function(ont) {
  
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(stringr)
  
  ego <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Mm.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  df <- ego@result
  df <- df[df$Count >= 3, ]  # Optional filter: remove terms with low gene count
  df$Category <- ont
  return(df)
}

# Example usage:
# Run enrichments
#ego_BP <- run_GO_enrichment("BP")
#ego_CC <- run_GO_enrichment("CC")
#ego_MF <- run_GO_enrichment("MF")

# ----------------------
# Save combined and separate results
# ----------------------

# Combine all results
#all_GO <- bind_rows(ego_BP, ego_CC, ego_MF)

# Write output
#write.csv(all_GO, "results/GO_enrichment_combined.csv", row.names = FALSE)
#write.csv(ego_BP, "results/GO_enrichment_BP.csv", row.names = FALSE)
#write.csv(ego_CC, "results/GO_enrichment_CC.csv", row.names = FALSE)
#write.csv(ego_MF, "results/GO_enrichment_MF.csv", row.names = FALSE)

cat("GO enrichment for BP/CC/MF completed. Results saved in /results.\n")
