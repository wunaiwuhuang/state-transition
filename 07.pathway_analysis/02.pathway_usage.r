# please use conda env r43 to do next analysis
# Load the GSEA functions
source("01.pathway_GSEAandORA_script.r")

# Set up paths
base_dir <- "/disk2/cai113/data/stateTrans"
gtf_file <- file.path(base_dir, "02.alignment/STAR_ref/mouse_GRCm39/Mus_musculus.GRCm39.110.gtf")
deg_dir <- file.path(base_dir, "06.DEG_analysis/DEG_results")
gsea_dir <- file.path(base_dir, "07.pathway_analysis/GSEA_results")

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
    methods = c("GO", "MSigDB"),
    go_ontologies = c("BP"),
    msigdb_categories = c("H", "C2", "C7", "C8")
)
rm(batch_results_all)
gc()

batch_results_up <- batch_gsea_analysis(
    deg_base_dir = deg_dir,
    contrast_list = contrast_list,
    gtf_file = gtf_file,
    output_base_dir = file.path(gsea_dir, "batch_up"),
    gene_filter = "up",
    methods = c("GO", "MSigDB"),
    go_ontologies = c("BP"),
    msigdb_categories = c("H", "C2", "C7", "C8")
)
rm(batch_results_up)
gc()

batch_results_down <- batch_gsea_analysis(
    deg_base_dir = deg_dir,
    contrast_list = contrast_list,
    gtf_file = gtf_file,
    output_base_dir = file.path(gsea_dir, "batch_down"),
    gene_filter = "down",
    methods = c("GO", "MSigDB"),
    go_ontologies = c("BP"),
    msigdb_categories = c("H", "C2", "C7", "C8")
)
rm(batch_results_down)
gc()