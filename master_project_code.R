library(GEOquery)
library(Seurat)
library(Matrix)
# library(matrixStats)
library(tools)
library(openxlsx)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(future)
# library(tibble)
library(data.table)
plan(sequential)


# Data load

# This function will be used to collapse gene names into one row: 
# e.g. GENE1, GENE1.1, GENE1.2, GENE1.1.2.1 will be collapsed to GENE1 and the median value will be taken for each cell.

collapse_by_median <- function(mat) {
  base_names <- sub("\\..*", "", rownames(mat))
  u_base_names <- unique(base_names)
  
  # If already unique, return original
  if (length(u_base_names) == nrow(mat)) return(mat)
  
  # Convert to data.table
  dt <- as.data.table(as.matrix(mat))
  # Add a new column named gene to dt, with the values from base_names
  dt[, gene := base_names]
  
  # Group by the gene column and compute the median of each column
  dt_median <- dt[, lapply(.SD, median), by = gene]
  
  # Convert back to matrix
  result <- as.matrix(dt_median[, -1])
  rownames(result) <- dt_median$gene
  colnames(result) <- colnames(mat)
  
  # Convert back to sparse
  result <- Matrix::Matrix(result, sparse = TRUE)
  
  return(result)
}


# Load GSE161529:
  

# # Make symbol gene names unique
# features161529 <- read.delim('GSE161529/GSE161529_features.tsv', header = FALSE)
# features161529$V2 <- make.unique(features161529$V2)
# write.table(features161529, file = "GSE161529/GSE161529_features_unique.tsv", 
#             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Function that reformats sample names
clean_sample_name_161529 <- function(name) {
  parts <- strsplit(name, "-", fixed = TRUE)[[1]]
  
  if (length(parts) == 3 && parts[1] == "TN" && parts[2] == "B1") {
    # TN-B1-XYZnnnn to TN-B1-nnnn
    cleaned_last <- sub("^[A-Za-z]+", "", parts[3])
    return(paste(parts[1], parts[2], cleaned_last, sep = "-"))
    
  } else if (length(parts) >= 2 && parts[1] == "TN") {
    # TN-XYZnnnn[-suffix] to TN-nnnn[-suffix]
    cleaned_second <- sub("^[A-Za-z]+", "", parts[2])
    suffix <- if (length(parts) > 2) paste(parts[3:length(parts)], collapse = "-") else NULL
    result <- paste("TN", cleaned_second, sep = "-")
    if (!is.null(suffix)) result <- paste(result, suffix, sep = "-")
    return(result)
    
  } else {
    # Change other types
    cleaned_second <- sub("^[A-Za-z]+", "", parts[2])
    suffix <- if (length(parts) > 2) paste(parts[3:length(parts)], collapse = "-") else NULL
    result <- paste(parts[1], cleaned_second, sep = "-")
    if (!is.null(suffix)) result <- paste(result, suffix, sep = "-")
    return(result)
  }
}

features161529 <- 'GSE161529/GSE161529_features_unique.tsv'
matrix_files161529 <- list.files("GSE161529/suppl", pattern = "*.mtx$", full.names = TRUE)
barcode_files161529 <- list.files("GSE161529/suppl", pattern = "barcodes\\.tsv$", full.names = TRUE)

# Sort to maintain consistency
matrix_files161529 <- sort(matrix_files161529)
barcode_files161529 <- sort(barcode_files161529)

# Clean sample names
sample_names161529 <- sub("^[^_]+_([^-]+-.*)-matrix\\.mtx$", "\\1", basename(matrix_files161529))
sample_names161529 <- sapply(sample_names161529, clean_sample_name_161529, USE.NAMES = FALSE)
sample_names161529[sample_names161529 == "ER-0040"] <- "ER-0040-T"

# Create matrix for each sample
create_matrix_161529 <- function(mtx_file, barcode_file, sample_name) {
  mat <- ReadMtx(
    mtx = mtx_file,
    features = features161529,
    cells = barcode_file,
    feature.column = 2
  )
  colnames(mat) <- paste(sample_name, colnames(mat), sep = "_")
  mat <- collapse_by_median(mat)
  print(sample_name)
  return(mat)
}

count_list_161529 <- mapply(
  create_matrix_161529,
  mtx_file = matrix_files161529,
  barcode_file = barcode_files161529,
  sample_name = sample_names161529,
  SIMPLIFY = FALSE
)

# # Check that gene names are the same and in the same order
# ref_rownames161529 <- rownames(count_list_161529[[1]])
# all(sapply(count_list_161529, function(mat) {
#   identical(rownames(mat), ref_rownames161529)
# }))
# # TRUE

# Combine into a single matrix and make a Seurat object 
counts_161529 <- do.call(cbind, count_list_161529)
GSE161529 <- CreateSeuratObject(counts = counts_161529, project = "GSE161529")
GSE161529$sample_id <- sapply(strsplit(colnames(GSE161529), "_"), `[`, 1)
GSE161529$dataset_id <- "GSE161529"
rm(count_list_161529)
rm(counts_161529)


# Load GSE176078:
  

archive_paths176078 <- list.files("GSE176078/GSE176078_RAW/", pattern = "*.tar.gz$", full.names = TRUE)

excluded_samples176078 <- c("CID3963", "CID4066", "CID4398", "CID4513", "CID4523")
exclude_pattern176078 <- paste(excluded_samples176078, collapse = "|")

archive_paths176078 <- archive_paths176078[!grepl(exclude_pattern176078, archive_paths176078)]

# Initialize list to store matrices
count_list_176078 <- list()

for (archive in archive_paths176078) {
  # Extract sample name
  sample_name <- file_path_sans_ext(basename(archive))
  sample_name <- file_path_sans_ext(sample_name)
  sample_name <- sub("^.*?_", "", sample_name)
  
  # Extract to temporary directory
  extract_dir <- file.path(tempdir(), sample_name)
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  untar(archive, exdir = extract_dir)
  
  # Locate matrix files
  inner_dir <- list.dirs(extract_dir, full.names = TRUE, recursive = FALSE)[1]
  matrix_file   <- file.path(inner_dir, "count_matrix_sparse.mtx")
  genes_file    <- file.path(inner_dir, "count_matrix_genes.tsv")
  barcodes_file <- file.path(inner_dir, "count_matrix_barcodes.tsv")
  
  # Read matrix
  mat <- ReadMtx(
    mtx = matrix_file,
    features = genes_file,
    cells = barcodes_file,
    feature.column = 1
  )
  
  mat <- collapse_by_median(mat)
  print(sample_name)
  
  # Store matrix in list
  count_list_176078[[sample_name]] <- mat
  
  # Optional cleanup
  unlink(extract_dir, recursive = TRUE)
}

# # Check that gene names are the same and in the same order
# ref_rownames176078 <- rownames(count_list_176078[[1]])
# all(sapply(count_list_176078, function(mat) {
#   identical(rownames(mat), ref_rownames176078)
# }))
# # TRUE

# Make a single matrix
counts_176078 <- do.call(cbind, count_list_176078)
GSE176078 <- CreateSeuratObject(counts = counts_176078, project = "GSE176078")
GSE176078$sample_id <- sapply(strsplit(colnames(GSE176078), "_"), `[`, 1)
GSE176078$dataset_id <- "GSE176078"
rm(count_list_176078)
rm(counts_176078)
gc()



# Load GSE167036:
  

# Load count matrix
counts167036 <- readMM("GSE167036/GSE167036_cell_counts_matrix.mtx")
features167036 <- read.csv("GSE167036/GSE167036_features.csv")
barcodes167036 <- read.csv("GSE167036/GSE167036_barcodes.csv")

# Reformat cell IDs (replace _ with -)
barcodes167036$cell_id <- gsub("_", "-", barcodes167036$cell_id)

# Assign row and column names
rownames(counts167036) <- features167036$x
colnames(counts167036) <- barcodes167036$cell_id

# Read metadata
metadata167036 <- read.csv("GSE167036/GSE167036_meta_all.csv")
metadata167036$X <- NULL

# Replace _ with - in metadata$cell_id
metadata167036$cell_id <- gsub("_", "-", metadata167036$cell_id)
rownames(metadata167036) <- metadata167036$cell_id

# Exclude TCR cells
tcr167036 <- read.csv("GSE167036/GSE167036_meta_tcr.csv")
tcr167036$cell_id <- gsub("_", "-", tcr167036$cell_id)  # match format
cells_to_exclude <- tcr167036$cell_id

metadata167036 <- metadata167036[!(rownames(metadata167036) %in% cells_to_exclude), ]
counts167036 <- counts167036[, !(colnames(counts167036) %in% cells_to_exclude)]

# Check that cell names are the same and in the same order in counts and meta
stopifnot(all(rownames(metadata167036) == colnames(counts167036)))

# Create sample_id
metadata167036$sample_id <- paste(metadata167036$patient_id, metadata167036$sample_type, sep = "-")
metadata167036$sample_id <- gsub("[ _]", "-", metadata167036$sample_id)

# Update cell IDs to: <sample_id>_<cell_id>
new_cell_ids <- paste(metadata167036$sample_id, metadata167036$cell_id, sep = "_")

# Apply new cell IDs
colnames(counts167036) <- new_cell_ids
rownames(metadata167036) <- new_cell_ids

# Filter metadata to keep only selected columns
metadata167036 <- metadata167036[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample_id")]

# Final check
stopifnot(all(rownames(metadata167036) == colnames(counts167036)))

counts167036 <- collapse_by_median(counts167036)

# Create Seurat object
GSE167036 <- CreateSeuratObject(
  counts = counts167036,
  project = "GSE167036"
)
GSE167036$sample_id <- sapply(strsplit(colnames(GSE167036), "_"), `[`, 1)
GSE167036$dataset_id <- "GSE167036"
rm(counts167036)


# Load Lambrecht_PD1:
  

lampd1 <- readRDS('Lambrecht_PD1/1863-counts_cells_cohort1.rds')
# Remove on-treatment samples
lampd1 <- lampd1[, grepl("Pre", colnames(lampd1))]

# Reformat cell names to distinguish the sample_id
rename_id_lampd1 <- function(id) {
  # Find position of the last underscore
  last_underscore_pos <- max(gregexpr("_", id)[[1]])
  
  # Split into two parts: sample_id and cell_id
  sample_id_part <- substr(id, 1, last_underscore_pos - 1)
  cell_id_part <- substr(id, last_underscore_pos + 1, nchar(id))
  
  # Replace underscores with dashes in sample part
  sample_id_fixed <- gsub("_", "-", sample_id_part)
  
  # Combine with underscore separator
  paste0(sample_id_fixed, "_", cell_id_part)
}

colnames(lampd1) <- sapply(colnames(lampd1), rename_id_lampd1, USE.NAMES = FALSE)

lampd1 <- collapse_by_median(lampd1)

# Create Seurat object
LambrechtPD1 <- CreateSeuratObject(
  counts = lampd1,
  project = "LambrechtPD1"
)

# Add sample_id to meta
LambrechtPD1@meta.data$sample_id <- LambrechtPD1@meta.data$orig.ident
LambrechtPD1$dataset_id <- "LambrechtPD1"
rm(lampd1)



# Load GSE114727 samples BC09, BC10, BC11:
  

data_dir114727 <- "GSE114727/BC091011/"

# List all matrix files and sample names
matrix_files114727 <- list.files(data_dir114727, pattern = "_matrix.mtx.gz$", full.names = TRUE)
sample_names114727 <- str_replace(basename(matrix_files114727), "_matrix.mtx.gz", "")

read_sample_matrix_114727 <- function(sample) {
  matrix_path <- file.path(data_dir114727, paste0(sample, "_matrix.mtx.gz"))
  genes_path <- file.path(data_dir114727, paste0(sample, "_genes.tsv.gz"))
  barcodes_path <- file.path(data_dir114727, paste0(sample, "_barcodes.tsv.gz"))
  
  # Read files
  expression_matrix <- readMM(matrix_path)
  genes <- read.delim(genes_path, header = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE)
  
  # Assign row/col names
  rownames(expression_matrix) <- make.unique(genes$V2)
  colnames(expression_matrix) <- barcodes$V1
  colnames(expression_matrix) <- paste(sample, colnames(expression_matrix), sep = "_")
  
  # Collapse duplicate genes
  expression_matrix <- collapse_by_median(expression_matrix)
  
  return(expression_matrix)
}

count_list_114727 <- lapply(sample_names114727, read_sample_matrix_114727)

# # Check that genes are the same and in the same oder in all matrice
# ref_rownames114727 <- rownames(count_list_114727[[1]])
# all(sapply(count_list_114727, function(mat) {
#   identical(rownames(mat), ref_rownames114727)
# }))
# # TRUE

# Combine all samples into a single sparse matrix
counts_114727 <- do.call(cbind, count_list_114727)

# Change cell names for consistency
clean_cell_name_114727 <- function(name) {
  # Split on first underscore, discard the first part
  name_no_prefix <- sub("^[^_]+_", "", name)
  # Replace first underscore after that with dash
  sub("([^_]+)_", "\\1-", name_no_prefix)
}

colnames(counts_114727) <- unname(sapply(colnames(counts_114727), clean_cell_name_114727))

# Create Seurat object
GSE114727 <- CreateSeuratObject(
  counts = counts_114727,
  project = "GSE114727"
)

# Add sample_id metadata by extracting prefix before first underscore in cell names
GSE114727$sample_id <- sapply(strsplit(colnames(GSE114727), "_"), `[`, 1)
GSE114727$dataset_id <- "GSE114727"
rm(count_list_114727)
rm(counts_114727)



# Load GSE110686

data_dir110686 <- "GSE110686"

# List matrix files and extract sample names
matrix_files110686 <- list.files(data_dir110686, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
sample_names110686 <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files110686))

# Sort to ensure consistent order
matrix_files110686 <- sort(matrix_files110686)
sample_names110686 <- sort(sample_names110686)

# Define paths to genes and barcodes
genes_paths <- file.path(data_dir110686, paste0(sample_names110686, "_genes.tsv.gz"))
barcodes_paths <- file.path(data_dir110686, paste0(sample_names110686, "_barcodes.tsv.gz"))

# Read count matrix
read_sample_matrix_110686 <- function(mtx_file, genes_file, barcodes_file, sample_name) {
  mat <- readMM(mtx_file)
  genes <- read.table(genes_file, header = FALSE, sep = "\t")
  barcodes <- read.table(barcodes_file, header = FALSE, sep = "\t")$V1
  
  # Assign names
  rownames(mat) <- make.unique(genes$V2)
  colnames(mat) <- paste0(sample_name, "_", barcodes)
  
  # Collapse duplicate gene rows
  mat <- collapse_by_median(mat)
  
  return(mat)
}

count_list_110686 <- mapply(
  read_sample_matrix_110686,
  mtx_file = matrix_files110686,
  genes_file = genes_paths,
  barcodes_file = barcodes_paths,
  sample_name = sample_names110686,
  SIMPLIFY = FALSE
)

# # Check that gene names are the same and in the same order
# all(rownames(count_list_110686[[1]]) == rownames(count_list_110686[[2]]))
# # [1] TRUE

# Combine all into a single matrix
counts_110686 <- do.call(cbind, count_list_110686)

colnames(counts_110686) <- sub("^GSM[0-9]+_", "", colnames(counts_110686))


# Create Seurat v5 object
GSE110686 <- CreateSeuratObject(
  counts = counts_110686,
  project = "GSE110686"
)

# Add sample_id metadata
GSE110686$sample_id <- sapply(strsplit(colnames(GSE110686), "_"), `[`, 1)
GSE110686$dataset_id <- "GSE110686"
rm(counts_110686)
rm(count_list_110686)



# Load GSE148673:
  

files148673 <- list.files('GSE148673/counts/', full.names = TRUE)

# Initialize a list to store matrices
expr_list148673 <- list()

for (file in files148673) {
  # Extract sample name from filename
  sample_name <- sub(".*_([^_]+)\\.txt\\.gz$", "\\1", basename(file))
  
  # Read and clean the expression matrix
  expr <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  expr <- expr[3:nrow(expr), , drop = FALSE]
  expr[] <- lapply(expr, as.numeric)
  expr <- Matrix(as.matrix(expr), sparse = TRUE)
  
  # Prefix cell barcodes
  colnames(expr) <- paste(sample_name, colnames(expr), sep = "_")
  
  expr <- collapse_by_median(expr)
  
  # Store the matrix
  expr_list148673[[sample_name]] <- expr
}

# Combine matrices using full join on genes (rownames)
all_genes148673 <- Reduce(union, lapply(expr_list148673, rownames))

# Fill in missing genes with zeros
fill_missing_genes_148673 <- function(mat) {
  missing_genes <- setdiff(all_genes148673, rownames(mat))
  if (length(missing_genes) > 0) {
    filler <- matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                     dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, filler)
  }
  # Reorder for match
  mat[all_genes148673, , drop = FALSE]
}

expr_list148673_filled <- lapply(expr_list148673, fill_missing_genes_148673)

# ref_rownames148673 <- rownames(expr_list148673_filled[[1]])
# all(sapply(expr_list148673_filled, function(mat) {
#   identical(rownames(mat), ref_rownames148673)
# }))
# # TRUE
# 
counts_148673 <- do.call(cbind, expr_list148673_filled)

# # Create Seurat v5 object
GSE148673 <- CreateSeuratObject(counts = counts_148673, project = "GSE148673")
GSE148673$sample_id <- sapply(strsplit(colnames(GSE148673), "_"), `[`, 1)
GSE148673$dataset_id <- "GSE148673"
rm(counts_148673)
rm(expr_list148673)
rm(expr_list148673_filled)




# Make a single list of Seurat objects:
  
seurat_list_all <- list(GSE161529, GSE176078, GSE167036, LambrechtPD1, 
                        GSE114727, GSE110686, GSE148673)
names(seurat_list_all) <- c('GSE161529', 'GSE176078', 'GSE167036', 'LambrechtPD1', 
                            'GSE114727', 'GSE110686', 'GSE148673')
saveRDS(seurat_list_all, 'GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673.rds')
rm(GSE161529)
rm(GSE176078)
rm(GSE167036)
rm(LambrechtPD1)
rm(GSE114727)
rm(GSE110686)
rm(GSE148673)
gc()



# Add metadata


# seurat_list_all <- readRDS('GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673.rds')
meta_full <- read.csv('meta.csv', sep=';', row.names = 1)
rownames(meta_full) <- gsub("_", "-", rownames(meta_full))
meta_full$Dataset_ID <- NULL



add_patient_metadata <- function(seu, metadata_df) {
  meta <- seu@meta.data
  meta$.__cell_barcode__ <- rownames(meta)
  
  metadata_df$sample_id <- rownames(metadata_df)
  
  merged_meta <- merge(meta, metadata_df, by = "sample_id", all.x = TRUE, sort = FALSE)
  
  rownames(merged_meta) <- merged_meta$.__cell_barcode__
  merged_meta$.__cell_barcode__ <- NULL
  
  seu@meta.data <- merged_meta
  return(seu)
}

seurat_list_all <- lapply(seurat_list_all, add_patient_metadata, meta_full)


# QC

# For each Seurat object, filtering out cells with < 100 and > 2500 genes (> 7000 genes for GSE148673)
# and with > 10% of mitochondrial content:
  

filter_min <- function(seu) {
  # Filter out cells with < 100 genes expressed
  seu <- subset(seu, subset = nFeature_RNA >= 100)
  return(seu)
}

filter_mito <- function(seu) {
  # Add percent.mt if missing
  if (!"percent.mt" %in% colnames(seu[[]])) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  }
  # Filter out cells with > 10% mitochondrial content
  seu <- subset(seu, subset = percent.mt <= 10)
  return(seu)
}

filter_max <- function(seu) {
  # Filter out cells with > 2500 genes expressed
  dataset_id <- seu@project.name
  if (dataset_id == "GSE148673") {
    seu <- subset(seu, subset = nFeature_RNA <= 7000)
  } else {
    seu <- subset(seu, subset = nFeature_RNA <= 2500)
  }
  return(seu)
}

filter_2sd <- function(seu) {
  # Add log10_nCount_RNA if missing
  if (!"log10_nCount_RNA" %in% colnames(seu[[]])) {
    seu[["log10_nCount_RNA"]] <- log10(seu$nCount_RNA + 1)
  }
  # Filter out cells within 2 SDs of log10_nCount_RNA
  log10_vals <- seu$log10_nCount_RNA
  upper_bound <- mean(log10_vals) + 2 * sd(log10_vals)
  lower_bound <- mean(log10_vals) - 2 * sd(log10_vals)
  
  seu <- subset(seu, subset = log10_nCount_RNA >= lower_bound & log10_nCount_RNA <= upper_bound)
  
  return(seu)
}



seurat_list_min <- lapply(seurat_list_all, filter_min)
seurat_list_mito <- lapply(seurat_list_min, filter_mito)
seurat_list_max <- lapply(seurat_list_mito, filter_max)
seurat_list_2sd <- lapply(seurat_list_max, filter_2sd)



saveRDS(seurat_list_2sd, 
        "GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673_list_qc.rds")


## QC stats

# Stats by dataset:
  
get_cell_counts <- function(seu, colname) {
  df <- data.frame(dataset_name = seu@project.name)
  df[[colname]] <- ncol(seu)
  return(df)
}

cell_count_init <- do.call(rbind, lapply(seurat_list_all, get_cell_counts, colname = 'initial_cells'))
cell_count_min <- do.call(rbind, lapply(seurat_list_min, get_cell_counts, colname = 'cells_min_qc'))
cell_count_mito <- do.call(rbind, lapply(seurat_list_mito, get_cell_counts, colname = 'cells_mito_qc'))
cell_count_max <- do.call(rbind, lapply(seurat_list_max, get_cell_counts, colname = 'cells_max_qc'))
cell_count_2sd <- do.call(rbind, lapply(seurat_list_2sd, get_cell_counts, colname = 'cells_2sd_qc'))
cell_count_after_qc <- do.call(rbind, lapply(seurat_list_2sd, get_cell_counts, colname = 'cells_after_qc'))


cell_count_list = list(cell_count_init, cell_count_after_qc, cell_count_min, cell_count_mito,
                       cell_count_max, cell_count_2sd)
cell_count_merge <- Reduce(function(x, y) merge(x, y, by = "dataset_name", all = TRUE), cell_count_list)


# QC plots per dataset:
  
cell_count_steps <- cell_count_merge[, -3]

# Calculate percentage for each QC step
for (i in 3:ncol(cell_count_steps)) {
  cell_count_steps[[i]] <- round((cell_count_steps[[i]] / cell_count_steps[[2]]) * 100, 2)
}
cell_count_steps$initial_cells <- 100

# Prepare the data for plotting
cell_count_long <- pivot_longer(cell_count_steps,
                                cols = -dataset_name,
                                names_to = "qc_step",
                                values_to = "percentage")

# Preserve original QC step order
qc_step_levels <- colnames(cell_count_steps)[-1]  # exclude dataset_name
cell_count_long$qc_step <- factor(cell_count_long$qc_step, levels = qc_step_levels)

ggplot(cell_count_long, aes(x = qc_step, y = percentage, fill = qc_step)) +
  geom_col() +
  facet_wrap(~ dataset_name, scales = "free_y") +
  labs(x = "QC Step", y = "Percentage of Cells",
       title = "Percentage of Cells after QC Steps") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10)) + guides(fill = "none")



# Number of samples:
  
get_sample_counts <- function(seu) {
  data.frame(
    dataset_name = seu@project.name,
    sample_number = length(unique(seu$sample_id))
  )
}
sample_num <- do.call(rbind, lapply(seurat_list_all, get_sample_counts))
stat_list <- list(sample_num, cell_count_merge)
cell_sample_count <- Reduce(function(x, y) merge(x, y, by = "dataset_name", all = TRUE), stat_list)

write.xlsx(cell_sample_count, file = "QC.xlsx", sheetName = "dataset_untreated", rowNames = FALSE)


# Stats by sample:
  
# Function to get cell counts per sample for a single Seurat object
count_cells_per_sample_single <- function(seu, add_dataset_id = TRUE, cell_count_col = "cell_number") {
  meta <- seu@meta.data
  counts <- as.data.frame(table(meta$sample_id))
  colnames(counts) <- c("sample_id", cell_count_col)
  
  if (add_dataset_id) {
    counts$dataset_id <- seu@project.name
    counts <- counts[, c("sample_id", "dataset_id", cell_count_col)]
  }
  
  return(counts)
}

count_cells_per_sample <- function(seurat_list, cell_count_col, add_dataset_id = TRUE) {
  counts_list <- lapply(seurat_list, count_cells_per_sample_single, cell_count_col = cell_count_col, 
                        add_dataset_id = add_dataset_id)
  combined_counts <- do.call(rbind, counts_list)
  rownames(combined_counts) <- NULL
  return(combined_counts)
}

qc_per_sample_init <- count_cells_per_sample(seurat_list_all, "cells_initial", add_dataset_id = TRUE)
qc_per_sample_min <- count_cells_per_sample(seurat_list_min, "cells_after_min_qc", add_dataset_id = F)
qc_per_sample_mito <- count_cells_per_sample(seurat_list_mito, "cells_after_mito_qc", add_dataset_id = F)
qc_per_sample_max <- count_cells_per_sample(seurat_list_max, "cells_after_max_qc", add_dataset_id = F)
qc_per_sample_2sd <- count_cells_per_sample(seurat_list_2sd, "cells_after_2sd_qc", add_dataset_id = F)

qc_per_sample_stats <- list(qc_per_sample_init, qc_per_sample_min, qc_per_sample_mito,
                            qc_per_sample_max, qc_per_sample_2sd)
qc_sample_table <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), qc_per_sample_stats)
qc_sample_table <- qc_sample_table[order(qc_sample_table$dataset_id), ]

write.xlsx(qc_sample_table, file = "QC.xlsx", sheetName = "sample_untreated", rowNames = FALSE)



# Shared features

seurat_list_all <- readRDS('GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673_list_qc.rds')

get_single_gene_stats <- function(seu, shared_genes) {
  dataset_name <- seu@project.name
  gene_names <- rownames(seu)
  total_genes <- length(gene_names)
  shared_count <- sum(gene_names %in% shared_genes)
  shared_pct <- round((shared_count / total_genes) * 100, 2)
  
  data.frame(
    dataset = dataset_name,
    total_features = total_genes,
    shared_features = shared_count,
    shared_percentage = shared_pct
  )
}

get_shared_gene_stats <- function(seurat_list) {
  # Get the list of gene names from each object
  gene_lists <- lapply(seurat_list, rownames)
  
  # Identify shared genes across all datasets
  shared_genes <- Reduce(intersect, gene_lists)
  
  # Apply the helper function to each Seurat object
  stats_list <- lapply(seurat_list, get_single_gene_stats, shared_genes = shared_genes)
  
  # Combine into a single data.frame
  do.call(rbind, stats_list)
}
shared_gene_table <- get_shared_gene_stats(seurat_list_all)

shared_features <- Reduce(intersect, lapply(seurat_list_all, rownames))

write.table(shared_features, file = "shared_features.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)



# Merge


seurat_list_qc <- readRDS('GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673_list_qc.rds')
merged_all <- merge(
  x = seurat_list_qc[[1]],
  y = seurat_list_qc[-1],
  add.cell.ids = names(seurat_list_qc),
  project = "GSE161529_GSE176078_GSE167036_LambrechtPD1_GSE114727_GSE110686_GSE148673"
)

rm(seurat_list_qc)
gc()

# Process the merged object
merged_all <- NormalizeData(merged_all, verbose = FALSE)
merged_all <- ScaleData(merged_all, verbose = FALSE)
merged_all <- FindVariableFeatures(merged_all)
merged_all <- RunPCA(merged_all, reduction.name = "pca", verbose = FALSE)

saveRDS(merged_all, file = "merged_161529_176078_167036_LamPD1_114727_110686_GSE148673_datasets.rds")



# Integration

# RPCA integration


# merged_all <- readRDS('merged_161529_176078_167036_LamPD1_114727_110686_GSE148673_datasets.rds')
integrated_all_rpca <- IntegrateLayers(
  object = merged_all,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "rpca_integrated",
  verbose = T
)

rm(merged_all)
gc()

saveRDS(integrated_all_rpca, file = "integrated161529_176078_167036_LamPD1_114727_110686_GSE148673_rpca.rds")


# Harmony integration


merged_all <- readRDS('merged_161529_176078_167036_LamPD1_114727_110686_GSE148673_datasets.rds')
integrated_all_harmony <- IntegrateLayers(
  object = merged_all,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = T
)

rm(merged_all)
gc()

saveRDS(integrated_all_harmony, file = "integrated161529_176078_167036_LamPD1_114727_110686_GSE148673_harmony.rds")



# UMAP Merged

merged_all <- readRDS('merged_161529_176078_167036_LamPD1_114727_110686_GSE148673.rds')
# Calculate the percentage of variance explained by each PC
pct_merged <- merged_all@reductions$pca@stdev / sum(merged_all@reductions$pca@stdev) * 100
# Calculate the cumulative variance explained across PCs
cumu_merged <- cumsum(pct_merged)
# Find the first PC where: cumulative variance > 90% 
# and individual PC contributes less than 5% variance
co1_merged <- which(cumu_merged > 90 & pct_merged < 5)[1] 
# Find the "elbow point" where the drop in explained variance between consecutive 
# PCs is still meaningful (> 0.1% drop)
co2_merged <- sort(which((pct_merged[1:length(pct_merged) - 1] - 
                            pct_merged[2:length(pct_merged)]) > 0.1), 
                   decreasing = T)[1] + 1 
# Take the more conservative estimate with fewer PCs between co1 and co2
co3_merged <- min(co1_merged, co2_merged) # 16


# Elbow plot merged
  
# Variance explained df
elbow_df_merged <- data.frame(
  PC = 1:length(pct_merged),
  Variance = pct_merged
)

ggplot(elbow_df_merged, aes(x = PC, y = Variance)) +
  geom_point(color = "black") +
  labs(
    title = "Elbow plot for the merged object",
    x = "PC", y = "Variance Explained") + theme_minimal()


# Find neighbors, find clusters and UMAP before batch analysis:
  
merged_all <- FindNeighbors(object = merged_all, dims = 1:co3_merged, verbose = F)
merged_all <- FindClusters(object = merged_all, verbose = F, resolution=0.2)
merged_all <- RunUMAP(object = merged_all, dims = 1:co3_merged, verbose = F)



# Number of dims for UMAP - RPCA


integrated_all_rpca <- readRDS('integrated161529_176078_167036_LamPD1_114727_110686_GSE148673_rpca.rds')

# Extract the PC matrix from the integrated reduction
pc_matrix_rpca <- Embeddings(integrated_all_rpca, reduction = "rpca_integrated")

# Calculate standard deviations
pc_stdev_rpca <- apply(pc_matrix_rpca, 2, sd)

# Compute percentage variance
pct_var_rpca <- pc_stdev_rpca^2 / sum(pc_stdev_rpca^2) * 100

# Determine optimal dims
cumu_var_rpca <- cumsum(pct_var_rpca)
co1_rpca <- which(cumu_var_rpca > 90 & pct_var_rpca < 5)[1]
co2_rpca <- sort(which((pct_var_rpca[1:length(pct_var_rpca)-1] - 
                          pct_var_rpca[2:length(pct_var_rpca)]) > 0.1),
                 decreasing = TRUE)[1] + 1
co3_rpca <- min(co1_rpca, co2_rpca) # 33



# Elbow plot:
  

# # Not working because no data in integrated_all_rpca[["rpca_integrated"]]@stdev
# ElbowPlot(integrated_all_rpca, reduction = "rpca_integrated", ndims = 50) +
#   ggtitle("Elbow plot for RPCA integration")

# Variance explained df
elbow_df_rpca <- data.frame(
  PC = 1:length(pct_var_rpca),
  Variance = pct_var_rpca
)

# # SD df
# elbow_df_rpca <- data.frame(
#   PC = 1:length(pc_stdev_rpca),
#   StdDev = pc_stdev_rpca
# )

ggplot(elbow_df_rpca, aes(x = PC, y = Variance)) +
  geom_point(color = "black") +
  labs(
    title = "Elbow plot for the RPCA integrated object",
    x = "PC", y = "Variance Explained") + theme_minimal()



# Number of dims for UMAP - Harmony


integrated_all_harmony <- readRDS('integrated161529_176078_167036_LamPD1_114727_110686_GSE148673_harmony.rds')

# Extract the PC matrix from the integrated reduction
pc_matrix_harmony <- Embeddings(integrated_all_harmony, reduction = "harmony")

# Calculate standard deviations
pc_stdev_harmony <- apply(pc_matrix_harmony, 2, sd)

# Compute percentage variance
pct_var_harmony <- pc_stdev_harmony^2 / sum(pc_stdev_harmony^2) * 100

# Determine optimal dims
cumu_var_harmony <- cumsum(pct_var_harmony)
co1_harmony <- which(cumu_var_harmony > 90 & pct_var_harmony < 5)[1]
co2_harmony <- sort(which((pct_var_harmony[1:length(pct_var_harmony)-1] - 
                             pct_var_harmony[2:length(pct_var_harmony)]) > 0.1),
                    decreasing = TRUE)[1] + 1
co3_harmony <- min(co1_harmony, co2_harmony) # 33


elbow_df_harmony <- data.frame(
  PC = 1:length(pct_var_harmony),
  Variance = pct_var_harmony
)

ggplot(elbow_df_harmony, aes(x = PC, y = Variance)) +
  geom_point(color = "black") +
  labs(
    title = "Elbow plot for the Harmony integrated object",
    x = "PC", y = "Variance Explained") + theme_minimal()



# Correlation between RPCA and Harmony

# Correlation between SDs
cor(pc_stdev_rpca, pc_stdev_harmony)



# UMAP RPCA and Harmony


integrated_all_rpca <- FindNeighbors(integrated_all_rpca, dims = 1:co3_rpca, 
                                     verbose = F, reduction = "rpca_integrated")
integrated_all_rpca <- FindClusters(integrated_all_rpca, verbose = F, resolution = 0.2)
integrated_all_rpca <- RunUMAP(integrated_all_rpca, dims = 1:co3_rpca, 
                               verbose = F, reduction = "rpca_integrated")



integrated_all_harmony <- FindNeighbors(integrated_all_harmony, dims = 1:co3_harmony, 
                                        verbose = F, reduction = "harmony")
integrated_all_harmony <- FindClusters(integrated_all_harmony, verbose = F, resolution = 0.2)
integrated_all_harmony <- RunUMAP(integrated_all_harmony, dims = 1:co3_harmony, 
                                  verbose = F, reduction = "harmony")




# Results

## Batch effect

### By datasset


DimPlot(object = merged_all, reduction = "umap", group.by = "dataset_id") +
  ggtitle('Batches by dataset - merged')



DimPlot(object = integrated_all_rpca, reduction = "umap", group.by = "dataset_id") +
  ggtitle('Batches by dataset - RPCA integrated')



DimPlot(object = integrated_all_harmony, reduction = "umap", group.by = "dataset_id") +
  ggtitle('Batches by dataset - Harmony integrated')



### By Seurat cluster


DimPlot(object = merged_all, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle('Batches by Seurat cluster - merged')



DimPlot(object = integrated_all_rpca, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle('Batches by Seurat cluster - RPCA integrated')



DimPlot(object = integrated_all_harmony, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle('Batches by Seurat cluster - Harmony integrated')


### By Sample Type


DimPlot(object = merged_all, reduction = "umap", group.by = "Sample_type")+
  ggtitle("Batches by sample type - merged")



DimPlot(object = integrated_all_rpca, reduction = "umap", group.by = "Sample_type")+
  ggtitle("Batches by sample type - RPCA integrated")



DimPlot(object = integrated_all_harmony, reduction = "umap", group.by = "Sample_type")+
  ggtitle("Batches by sample type - Harmony integrated")


## Cluster composition

### Samples by cluster heatmap


cluster_sample_merged <- prop.table(table(merged_all@meta.data$seurat_clusters, 
                                          merged_all@meta.data$sample_id), margin=1) * 100

pheatmap(cluster_sample_merged,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_col = 4,
         angle_col = 315,
         # display_numbers = TRUE,
         main = "Cluster Composition by Sample - Merged")




cluster_sample_rpca <- prop.table(table(integrated_all_rpca@meta.data$seurat_clusters, 
                                        integrated_all_rpca@meta.data$sample_id), margin=1) * 100

pheatmap(cluster_sample_rpca,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_col = 4,
         angle_col = 315,
         # display_numbers = TRUE,
         main = "Cluster Composition by Sample - RPCA")




cluster_sample_harmony <- prop.table(table(integrated_all_harmony@meta.data$seurat_clusters, 
                                           integrated_all_harmony@meta.data$sample_id), margin=1) * 100

pheatmap(cluster_sample_harmony,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_col = 4,
         angle_col = 315,
         # display_numbers = TRUE,
         main = "Cluster Composition by Sample - Harmony")



### Cells from datasets by cluster


cluster_dataset_merged <- prop.table(table(merged_all@meta.data$seurat_clusters, 
                                           merged_all@meta.data$dataset_id), margin=1) * 100

pheatmap(t(cluster_dataset_merged),
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 0,
         # display_numbers = TRUE,
         main = "Cluster Composition by Dataset Heatmap - Merged")



cluster_dataset_merged_df <- as.data.frame(cluster_dataset_merged)

# Rename columns for clarity
colnames(cluster_dataset_merged_df) <- c("Cluster", "Dataset", "Percentage")

# Plot
ggplot(cluster_dataset_merged_df, aes(x = Cluster, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity") +
  ylab("% of cells in cluster") +
  xlab("Cluster") +
  ggtitle("Cluster Composition by Dataset - Merged") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )




cluster_dataset_rpca <- prop.table(table(integrated_all_rpca@meta.data$seurat_clusters, 
                                         integrated_all_rpca@meta.data$dataset_id), margin=1) * 100

pheatmap(t(cluster_dataset_rpca),
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 0,
         # display_numbers = TRUE,
         main = "Cluster Composition by Dataset Heatmap - RPCA")



cluster_dataset_rpca_df <- as.data.frame(cluster_dataset_rpca)

# Rename columns for clarity
colnames(cluster_dataset_rpca_df) <- c("Cluster", "Dataset", "Percentage")

# Plot
ggplot(cluster_dataset_rpca_df, aes(x = Cluster, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity") +
  ylab("% of cells in cluster") +
  xlab("Cluster") +
  ggtitle("Cluster Composition by Dataset - RPCA") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )




cluster_dataset_harmony <- prop.table(table(integrated_all_harmony@meta.data$seurat_clusters, 
                                            integrated_all_harmony@meta.data$dataset_id), margin=1) * 100

pheatmap(t(cluster_dataset_harmony),
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 0,
         # display_numbers = TRUE,
         main = "Cluster Composition by Dataset Heatmap - Harmony")



cluster_dataset_harmony_df <- as.data.frame(cluster_dataset_harmony)

# Rename columns for clarity
colnames(cluster_dataset_harmony_df) <- c("Cluster", "Dataset", "Percentage")

# Plot
ggplot(cluster_dataset_harmony_df, aes(x = Cluster, y = Percentage, fill = Dataset)) +
  geom_bar(stat = "identity") +
  ylab("% of cells in cluster") +
  xlab("Cluster") +
  ggtitle("Cluster Composition by Dataset - Harmony") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


## Features

### CD45


VlnPlot(merged_all, features = "PTPRC", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("Clusters encriched with CD45 - merged")



VlnPlot(integrated_all_rpca, features = "PTPRC", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("CD45 Expression Levels in Clusters - PRCA")



VlnPlot(integrated_all_harmony, features = "PTPRC", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("CD45 Expression Levels in Clusters - Harmony")


### MALAT1


VlnPlot(merged_all, features = "MALAT1", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("Clusters encriched with MALAT1 - merged")



VlnPlot(integrated_all_rpca, features = "MALAT1", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("MALAT1 Expression Levels in Clusters - RPCA")



VlnPlot(integrated_all_harmony, features = "MALAT1", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("MALAT1 Expression Levels in Clusters - Harmony")


## Boxplots

### Percent.mt per sample by cluster


meta_merged <- merged_all@meta.data
meta_merged$seurat_clusters <- as.factor(meta_merged$seurat_clusters)

ggplot(meta_merged, aes(x = seurat_clusters, y = percent.mt)) +
  geom_boxplot(fill = "darkgreen") + 
  stat_compare_means(method = "anova", 
                     label.y = max(meta_merged$percent.mt, na.rm = TRUE),
                     label.x = length(unique(meta_merged$seurat_clusters)) / 2 + 0.5) +
  theme_minimal() +
  labs(title = "Mito Content Means by Clusters with ANOVA - Merged",
       x = "Cluster", y = "Percent MT")



meta_integrated_rpca <- integrated_all_rpca@meta.data
meta_integrated_rpca$seurat_clusters <- as.factor(meta_integrated_rpca$seurat_clusters)

ggplot(meta_integrated_rpca, aes(x = seurat_clusters, y = percent.mt)) +
  geom_boxplot(fill = "darkgreen") + 
  stat_compare_means(method = "anova", 
                     label.y = max(meta_integrated_rpca$percent.mt, na.rm = TRUE),
                     label.x = length(unique(meta_integrated_rpca$seurat_clusters)) / 2 + 0.5) +
  theme_minimal() +
  labs(title = "Mito Content Means by Clusters with ANOVA - RPCA",
       x = "Cluster", y = "Percent MT")



meta_integrated_harmony <- integrated_all_harmony@meta.data
meta_integrated_harmony$seurat_clusters <- as.factor(meta_integrated_harmony$seurat_clusters)

ggplot(meta_integrated_harmony, aes(x = seurat_clusters, y = percent.mt)) +
  geom_boxplot(fill = "darkgreen") + 
  stat_compare_means(method = "anova", 
                     label.y = max(meta_integrated_harmony$percent.mt, na.rm = TRUE),
                     label.x = length(unique(meta_integrated_harmony$seurat_clusters)) / 2 + 0.5) +
  theme_minimal() +
  labs(title = "Mito Content Means by Clusters with ANOVA - Harmony",
       x = "Cluster", y = "Percent MT")


### Features per sample


ggplot(meta_merged, aes(x = seurat_clusters, y = nFeature_RNA)) +
  geom_boxplot(fill = "steelblue") +
  stat_compare_means(
    method = "anova",
    label.y = max(meta_merged$nFeature_RNA, na.rm = TRUE) * 1.05,
    label.x = length(unique(meta_merged$seurat_clusters)) / 2 + 0.5
  ) +
  theme_minimal() +
  labs(
    title = "nFeature Means by Clusters with ANOVA - Merged",
    x = "Cluster",
    y = "Number of Features"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(meta_integrated_rpca, aes(x = seurat_clusters, y = nFeature_RNA)) +
  geom_boxplot(fill = "steelblue") +
  stat_compare_means(
    method = "anova",
    label.y = max(meta_integrated_rpca$nFeature_RNA, na.rm = TRUE) * 1.05,
    label.x = length(unique(meta_integrated_rpca$seurat_clusters)) / 2 + 0.5
  ) +
  theme_minimal() +
  labs(
    title = "nFeature Means by Clusters with ANOVA - RPCA",
    x = "Cluster",
    y = "Number of Features"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(meta_integrated_harmony, aes(x = seurat_clusters, y = nFeature_RNA)) +
  geom_boxplot(fill = "steelblue") +
  stat_compare_means(
    method = "anova",
    label.y = max(meta_integrated_harmony$nFeature_RNA, na.rm = TRUE) * 1.05,
    label.x = length(unique(meta_integrated_harmony$seurat_clusters)) / 2 + 0.5
  ) +
  theme_minimal() +
  labs(
    title = "nFeature Means by Clusters with ANOVA - Harmony",
    x = "Cluster",
    y = "Number of Features"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Cell cycle

### Obtain cell cycle info


# Seurat lists of S phase and G2/M phase genes
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes



assign_cell_cycle_phase <- function(norm_data, s.genes, g2m.genes) {
  # Ensure genes are present in the data
  s.genes.use <- intersect(rownames(norm_data), s.genes)
  g2m.genes.use <- intersect(rownames(norm_data), g2m.genes)
  
  if (length(s.genes.use) == 0 || length(g2m.genes.use) == 0) {
    stop("No overlap between gene lists and expression matrix rows.")
  }
  
  # Compute average scores
  s_score <- Matrix::colMeans(norm_data[s.genes.use, , drop = FALSE])
  g2m_score <- Matrix::colMeans(norm_data[g2m.genes.use, , drop = FALSE])
  
  # Assign phases
  phase <- ifelse(s_score > g2m_score, "G1S", "G2M")
  
  # Return scores and assigned phases as a data frame
  return(data.frame(S.Score = s_score, G2M.Score = g2m_score, Phase = phase))
}


### Cell cycle plots


# Assign cell phases for the merged data
norm_data_merge <- merged_all[["RNA"]]$data
cell_cycle_merge <- assign_cell_cycle_phase(norm_data_merge, s.genes, g2m.genes)

merged_all$S.Score <- cell_cycle_merge$S.Score
merged_all$G2M.Score <- cell_cycle_merge$G2M.Score
merged_all$Phase <- cell_cycle_merge$Phase



cell_cycle_long <- merged_all@meta.data %>%
  select(seurat_clusters, S.Score, G2M.Score) %>%
  pivot_longer(cols = c(S.Score, G2M.Score), 
               names_to = "CellCycleScore", 
               values_to = "Score") %>%
  mutate(
    CellCycleScore = factor(CellCycleScore, 
                            levels = c("S.Score", "G2M.Score"),
                            labels = c("S Score", "G2M Score")),
    seurat_clusters = factor(seurat_clusters)
  )

# # Not normally distributed
# ggplot(cell_cycle_long, aes(x = Score)) +
#   geom_histogram(aes(y = ..density..), bins = 50, fill = "skyblue", color = "black") +
#   geom_density(color = "red", size = 1) +
#   facet_wrap(~ CellCycleScore, scales = "free") +
#   theme_minimal() +
#   labs(title = "Distribution of Cell Cycle Scores - Merged", x = "Score", y = "Density")

cell_cycle_filtered <- cell_cycle_long %>%
  group_by(seurat_clusters, CellCycleScore) %>%
  filter(!is.na(Score)) %>%
  mutate(
    Q1 = quantile(Score, 0.25),
    Q3 = quantile(Score, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(Score >= lower_bound & Score <= upper_bound) %>%
  ungroup()

ggplot(cell_cycle_filtered, aes(x = seurat_clusters, y = Score, fill = CellCycleScore)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_compare_means(method = "anova", 
                     aes(group = CellCycleScore)) +
  # label.y = max(cell_cycle_long$Score, na.rm = TRUE) * 1.05, 
  # label.y = 0.29, label.params = list(angle = 90, hjust = 0.5)) +
  theme_minimal() + 
  labs(
    title = "Cell Cycle ANOVA - Merged",
    x = "Cluster",
    y = "Cell Cycle Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Assign cell phases for the Harmony integrated data
norm_data_rpca <- integrated_all_rpca[["RNA"]]$data
cell_cycle_rpca <- assign_cell_cycle_phase(norm_data_rpca, s.genes, g2m.genes)

integrated_all_rpca$S.Score <- cell_cycle_rpca$S.Score
integrated_all_rpca$G2M.Score <- cell_cycle_rpca$G2M.Score
integrated_all_rpca$Phase <- cell_cycle_rpca$Phase



cell_cycle_long_rpca <- integrated_all_rpca@meta.data %>%
  select(seurat_clusters, S.Score, G2M.Score) %>%
  pivot_longer(cols = c(S.Score, G2M.Score), 
               names_to = "CellCycleScore", 
               values_to = "Score") %>%
  mutate(
    CellCycleScore = factor(CellCycleScore, 
                            levels = c("S.Score", "G2M.Score"),
                            labels = c("S Score", "G2M Score")),
    seurat_clusters = factor(seurat_clusters)
  )

cell_cycle_filtered_rpca <- cell_cycle_long_rpca %>%
  group_by(seurat_clusters, CellCycleScore) %>%
  filter(!is.na(Score)) %>%
  mutate(
    Q1 = quantile(Score, 0.25),
    Q3 = quantile(Score, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(Score >= lower_bound & Score <= upper_bound) %>%
  ungroup()


ggplot(cell_cycle_filtered_rpca, aes(x = seurat_clusters, y = Score, fill = CellCycleScore)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_compare_means(method = "anova", 
                     aes(group = CellCycleScore)) +
  theme_minimal() + 
  labs(
    title = "Cell Cycle ANOVA - RPCA",
    x = "Cluster",
    y = "Cell Cycle Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Assign cell phases for the Harmony integrated data
norm_data_harmony <- integrated_all_harmony[["RNA"]]$data
cell_cycle_harmony <- assign_cell_cycle_phase(norm_data_harmony, s.genes, g2m.genes)

integrated_all_harmony$S.Score <- cell_cycle_harmony$S.Score
integrated_all_harmony$G2M.Score <- cell_cycle_harmony$G2M.Score
integrated_all_harmony$Phase <- cell_cycle_harmony$Phase



cell_cycle_long_harmony <- integrated_all_harmony@meta.data %>%
  select(seurat_clusters, S.Score, G2M.Score) %>%
  pivot_longer(cols = c(S.Score, G2M.Score), 
               names_to = "CellCycleScore", 
               values_to = "Score") %>%
  mutate(
    CellCycleScore = factor(CellCycleScore, 
                            levels = c("S.Score", "G2M.Score"),
                            labels = c("S Score", "G2M Score")),
    seurat_clusters = factor(seurat_clusters)
  )

cell_cycle_filtered_harmony <- cell_cycle_long_harmony %>%
  group_by(seurat_clusters, CellCycleScore) %>%
  filter(!is.na(Score)) %>%
  mutate(
    Q1 = quantile(Score, 0.25),
    Q3 = quantile(Score, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(Score >= lower_bound & Score <= upper_bound) %>%
  ungroup()


ggplot(cell_cycle_filtered_harmony, aes(x = seurat_clusters, y = Score, fill = CellCycleScore)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  stat_compare_means(method = "anova", 
                     aes(group = CellCycleScore)) +
  theme_minimal() + 
  labs(
    title = "Cell Cycle ANOVA - Harmony",
    x = "Cluster",
    y = "Cell Cycle Score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



