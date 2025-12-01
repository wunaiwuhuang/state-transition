################################################################################
# Install Required Packages for GSEA Analysis
# Run this script once before using the GSEA pipeline

# I use conda to install, conda install -c conda-forge and conda install -c bioconda command are useful
# please use conda env r43 to do next analysis
################################################################################

cat("========== Installing Required Packages ==========\n\n")

# Function to check and install packages
install_if_missing <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        if (bioc) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
            }
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else {
            install.packages(pkg)
        }
        cat(pkg, "installed successfully!\n\n")
    } else {
        cat(pkg, "already installed\n")
    }
}

# ========== Install CRAN packages ==========
cat(">>> Installing CRAN packages...\n")
cran_packages <- c("tidyverse", "ggplot2")

for (pkg in cran_packages) {
    install_if_missing(pkg, bioc = FALSE)
}

# ========== Install Bioconductor packages ==========
cat("\n>>> Installing Bioconductor packages...\n")
bioc_packages <- c(
    "clusterProfiler",      # Main enrichment analysis package
    "org.Mm.eg.db",         # Mouse genome annotation
    "enrichplot",           # Visualization
    "DOSE",                 # Disease Ontology
    "rtracklayer",          # GTF file reading
    "msigdbr"               # MSigDB gene sets
)

for (pkg in bioc_packages) {
    install_if_missing(pkg, bioc = TRUE)
}

# ========== Verify installation ==========
cat("\n========== Verifying Installation ==========\n")

required_packages <- c(
    "clusterProfiler",
    "org.Mm.eg.db",
    "enrichplot",
    "DOSE",
    "tidyverse",
    "rtracklayer",
    "msigdbr"
)

all_installed <- TRUE
for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat("✓", pkg, "is ready\n")
    } else {
        cat("✗", pkg, "installation failed\n")
        all_installed <- FALSE
    }
}

if (all_installed) {
    cat("\n========== All packages installed successfully! ==========\n")
    cat("You can now run the GSEA analysis pipeline.\n")
} else {
    cat("\n========== Some packages failed to install ==========\n")
    cat("Please check the error messages above and try installing them manually.\n")
}

# ========== Test basic functionality ==========
cat("\n========== Testing Basic Functionality ==========\n")

tryCatch({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(msigdbr)
    
    # Test msigdbr
    cat("Testing MSigDB access...\n")
    test_msigdb <- msigdbr(species = "Mus musculus", category = "H")
    cat("  MSigDB Hallmark gene sets:", length(unique(test_msigdb$gs_name)), "gene sets\n")
    
    # Test org.Mm.eg.db
    cat("Testing mouse annotation database...\n")
    test_genes <- c("ENSMUSG00000069805", "ENSMUSG00000020651")
    test_convert <- bitr(test_genes, 
                        fromType = "ENSEMBL",
                        toType = c("ENTREZID", "SYMBOL"),
                        OrgDb = org.Mm.eg.db)
    cat("  Converted", nrow(test_convert), "test genes successfully\n")
    
    cat("\n✓ All tests passed! Ready for GSEA analysis.\n")
    
}, error = function(e) {
    cat("\n✗ Test failed:", e$message, "\n")
    cat("Please check package installation.\n")
})

cat("\n========== Setup Complete ==========\n")