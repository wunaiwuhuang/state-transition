################################################################################
# GSEA Analysis Usage Example
# This script demonstrates how to use the GSEA pipeline
# please use conda env r43 to do next analysis
################################################################################

# Load the GSEA functions
source("01.pathway_GSEAandORA_script.r")

# Set up paths
base_dir <- "/disk2/cai113/data/stateTrans"
gtf_file <- file.path(base_dir, "02.alignment/STAR_ref/mouse_GRCm39/Mus_musculus.GRCm39.110.gtf")
deg_dir <- file.path(base_dir, "06.DEG_pathway/DEG_results")
gsea_dir <- file.path(base_dir, "07.pathway_analysis/GSEA_results")


################################################################################
# Example 1: Single Contrast GSEA (All DEGs, default settings)
################################################################################

gsea_c1_vs_c1star <- run_comprehensive_gsea(
    deg_file = file.path(deg_dir, "c1_vs_c1_star/DEG_significant.txt"),
    gtf_file = gtf_file,
    gene_filter = "all",        # Use all DEGs (up + down)
    methods = c("GO", "KEGG", "MSigDB"),
    go_ontologies = c("BP", "CC", "MF"),
    msigdb_categories = c("H", "C2", "C5"),
    gsea_method = "GSEA",       # Use GSEA (not ORA)
    output_dir = file.path(gsea_dir, "c1_vs_c1_star"),
    prefix = "GSEA"
)


################################################################################
# Example 2: Analyze Only Up-regulated Genes
################################################################################

gsea_c2_vs_c1_up <- run_comprehensive_gsea(
    deg_file = file.path(deg_dir, "c2_vs_c1/DEG_significant.txt"),
    gtf_file = gtf_file,
    gene_filter = "up",         # Only up-regulated genes
    padj_cutoff = 0.01,         # More stringent cutoff
    lfc_cutoff = 1.5,           # Require |log2FC| > 1.5
    output_dir = file.path(gsea_dir, "c2_vs_c1_up"),
    prefix = "GSEA_up"
)


################################################################################
# Example 3: Analyze Only Down-regulated Genes
################################################################################

gsea_c2_vs_c1_down <- run_comprehensive_gsea(
    deg_file = file.path(deg_dir, "c2_vs_c1/DEG_significant.txt"),
    gtf_file = gtf_file,
    gene_filter = "down",       # Only down-regulated genes
    output_dir = file.path(gsea_dir, "c2_vs_c1_down"),
    prefix = "GSEA_down"
)


################################################################################
# Example 4: Use ORA (Over-Representation Analysis) Instead of GSEA
################################################################################

ora_c3_vs_c2 <- run_comprehensive_gsea(
    deg_file = file.path(deg_dir, "c3_vs_c2/DEG_significant.txt"),
    gtf_file = gtf_file,
    gsea_method = "ORA",        # Use ORA instead of GSEA
    enrich_pvalue = 0.05,
    output_dir = file.path(gsea_dir, "c3_vs_c2_ORA"),
    prefix = "ORA"
)


################################################################################
# Example 5: Custom MSigDB Categories and GO Ontologies
################################################################################

gsea_custom <- run_comprehensive_gsea(
    deg_file = file.path(deg_dir, "c3_vs_c1_star/DEG_significant.txt"),
    gtf_file = gtf_file,
    methods = c("GO", "MSigDB"),        # Only GO and MSigDB
    go_ontologies = c("BP"),             # Only Biological Process
    msigdb_categories = c("H", "C5"),    # Hallmark + GO gene sets
    min_geneset_size = 15,               # Require at least 15 genes
    max_geneset_size = 300,              # Maximum 300 genes
    enrich_pvalue = 0.01,                # More stringent p-value
    output_dir = file.path(gsea_dir, "c3_vs_c1_star_custom")
)


################################################################################
# Example 6: BATCH ANALYSIS - All Contrasts at Once (RECOMMENDED!)
################################################################################

# Define all contrasts
contrast_list <- list(
    c("cp", "c1", "c1_star"),
    c("cp", "c2", "c1_star"),
    c("cp", "c3", "c1_star"),
    c("cp", "c2", "c1"),
    c("cp", "c3", "c1"),
    c("cp", "c3", "c2")
)

# Run batch GSEA analysis
batch_results_all <- batch_gsea_analysis(
    deg_base_dir = deg_dir,
    contrast_list = contrast_list,
    gtf_file = gtf_file,
    output_base_dir = file.path(gsea_dir, "batch_all"),
    gene_filter = "all",
    methods = c("GO", "KEGG", "MSigDB"),
    go_ontologies = c("BP", "CC", "MF"),
    msigdb_categories = c("H", "C2", "C5")
)


################################################################################
# Example 7: Batch Analysis - Only Up-regulated Genes
################################################################################

batch_results_up <- batch_gsea_analysis(
    deg_base_dir = deg_dir,
    contrast_list = contrast_list,
    gtf_file = gtf_file,
    output_base_dir = file.path(gsea_dir, "batch_up"),
    gene_filter = "up",
    padj_cutoff = 0.05,
    lfc_cutoff = 1
)


################################################################################
# Example 8: Batch Analysis - Only Down-regulated Genes
################################################################################

batch_results_down <- batch_gsea_analysis(
    deg_base_dir = deg_dir,
    contrast_list = contrast_list,
    gtf_file = gtf_file,
    output_base_dir = file.path(gsea_dir, "batch_down"),
    gene_filter = "down",
    padj_cutoff = 0.05,
    lfc_cutoff = 1
)


################################################################################
# Accessing Results
################################################################################

# Access specific results from batch analysis
c1_vs_c1star_results <- batch_results_all$c1_vs_c1_star

# View GO BP results
if (!is.null(c1_vs_c1star_results$GO_BP)) {
    print(head(c1_vs_c1star_results$GO_BP, 10))
}

# View KEGG results
if (!is.null(c1_vs_c1star_results$KEGG)) {
    print(head(c1_vs_c1star_results$KEGG, 10))
}

# View annotated gene list
print(head(c1_vs_c1star_results$annotated_genes))


################################################################################
# Quick Parameter Reference
################################################################################

# gene_filter options:
#   - "all": Use all DEGs (both up and down)
#   - "up": Only up-regulated genes
#   - "down": Only down-regulated genes

# methods options:
#   - c("GO", "KEGG", "MSigDB"): Run all three
#   - c("GO", "KEGG"): Skip MSigDB
#   - c("GO"): Only GO enrichment

# go_ontologies options:
#   - c("BP", "CC", "MF"): All three GO aspects
#   - c("BP"): Only Biological Process
#   - c("BP", "MF"): BP and Molecular Function

# msigdb_categories options:
#   - "H": Hallmark gene sets (50 gene sets)
#   - "C1": Positional gene sets
#   - "C2": Curated gene sets (includes KEGG, Reactome, etc.)
#   - "C3": Regulatory target gene sets
#   - "C5": Ontology gene sets (GO)
#   - "C6": Oncogenic signature gene sets
#   - "C7": Immunologic signature gene sets
#   - "C8": Cell type signature gene sets

# gsea_method options:
#   - "GSEA": Gene Set Enrichment Analysis (uses ranked list)
#   - "ORA": Over-Representation Analysis (uses gene set only)

# Cutoff parameters:
#   - padj_cutoff: Filter genes by adjusted p-value (default: 0.05)
#   - lfc_cutoff: Filter genes by |log2FC| (default: 0)
#   - enrich_pvalue: P-value cutoff for enrichment (default: 0.05)
#   - enrich_qvalue: Q-value cutoff for enrichment (default: 0.2)
#   - min_geneset_size: Minimum genes in pathway (default: 10)
#   - max_geneset_size: Maximum genes in pathway (default: 500)