# Tested with R version 4.5.0

required_packages <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "dplyr",
  "reshape2",
  "stringr",
  "tidyr",

)

# Automatically install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
