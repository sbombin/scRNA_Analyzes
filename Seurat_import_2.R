DGE_folder <- "/Users/sbombin/Desktop/analysis_mm38.comb/Ms802_BEC_sham/DGE_filtered/"
mat <- readMM(paste0(DGE_folder, file="DGE.mtx"))
cell_meta <- read.delim(paste0(DGE_folder, "cell_metadata.csv"),
                        stringsAsFactor = FALSE, sep = ",")
genes <- read.delim(paste0(DGE_folder, "all_genes.csv"),
                    stringsAsFactor = FALSE, sep = ",")

cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")
rownames(cell_meta) <- cell_meta$bc_wells
genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")

# Setting column and rownames to expression matrix
colnames(mat) <- genes$gene_name
rownames(mat) <- rownames(cell_meta)
mat_t <- t(mat)

# Remove empty rownames, if they exist
mat_t <- mat_t[(rownames(mat_t) != ""),]
pbmc <- CreateSeuratObject(mat_t, min.genes = 100, min.cells = 2, meta.data = cell_meta)
