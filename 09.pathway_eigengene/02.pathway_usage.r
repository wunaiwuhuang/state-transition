# please use conda env r43 to do next analysis

BASE_DIR <- "/disk2/cai113/data/stateTrans"
load(file.path(BASE_DIR, "08.eigengene", "driver_gene_selection_results.rdata"))
source(file.path(BASE_DIR, "09.pathway_eigengene", "01.pathway_script.r"))

all <- batch_driver_genes_gsea(
    driver_genes_list = results$driver_genes,
    output_base_dir = file.path(BASE_DIR, "09.pathway_eigengene", "GSEA_results"),
    min_geneset_size = 3,
    max_geneset_size = 100,
    methods = c("GO", "MSigDB"),
    go_ontologies = c("BP"),
    msigdb_categories = c("H", "C2", "C7", "C8")
)
rm(all)
gc()