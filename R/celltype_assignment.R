# Load required libraries
#devtools::install_github("JinmiaoChenLab/Rphenograph")

library(reticulate)
library(anndata)
library(tidyverse)
library(Rphenograph)
library(uwot)
library(pheatmap)
library(viridis)
library(ggplot2)
library(Seurat)

python_dir="/Users/romansankowski/anaconda3/bin/python"
use_python(python_dir, required = T)
sys = import("sys")

## load data
pth <- file.path("data","adata_objects")
fls <- list.files(pth)

if (!file.exists(file.path("results","merged_codex_object_meningioma.rds"))) {
  lst <- map(fls, function(x) {
    adata <- anndata::read_h5ad(file.path(pth, x))
    
    counts <- as.data.frame(adata$X)
    
    counts <- counts[ counts$DAPI > quantile(counts$DAPI , 0.05 ) , ]
    counts <- t(as.matrix(counts))
    
    seurat_obj <- CreateSeuratObject(counts = counts) %>% 
      AddMetaData(x, "core") %>% 
      AddMetaData(adata$obs$area, "area")
    
    coords <- adata$obsm$spatial[as.data.frame(adata$X)$DAPI > quantile(as.data.frame(adata$X)$DAPI , 0.05 ) , ] %>%
      as.data.frame()
    
    rownames(coords) <- colnames(seurat_obj)
    
    seurat_obj@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = coords)
    
    seurat_obj
    
  })
  
  obj <- merge(lst[[1]],lst[2:134])
  
  obj <- obj %>% JoinLayers()
  
  names(lst) <- fls
  ## save object
  saveRDS(obj, file.path("results","merged_codex_object_meningioma.rds"))
  saveRDS(lst, file.path("results","single_codex_object_list_meningioma.rds"))
} else {
  ## load data
  pth <- file.path("data","adata_objects")
  fls <- list.files(pth)
  
  obj <- readRDS(file.path("results","merged_codex_object_meningioma.rds"))
  data <- readRDS(file.path("results","single_codex_object_list_meningioma.rds"))
  
  names(data) <- fls
}

# ===== 1. LOAD AND PREPARE DATA =====

# Assuming your data is in a CSV with columns:
# - cell_id, X, Y (spatial coordinates)
# - marker columns (protein expression values)

# Extract marker expression matrix
# Adjust column names based on your specific markers
marker_cols <- c("CD45", "CD3e", "CD4", "CD8", "CD20", "CD79a", "CD68", "CD163", "CD206",
                 "SMA", "CD31","CD34","Ki67", "Collagen IV", "Claudin-5","Pax5","Iba1",
                 "GFAP", "Podoplanin", "FOXP3", "NeuN", "Neurofilament", "CD38", "Olig2", "CSF1R", "P2RY12")

# Keep only markers present in your data
marker_cols <- marker_cols[marker_cols %in% rownames(data[[1]])]

expr_matrix <- map(data, function(x) {
  as.data.frame(t(x[["RNA"]]$counts[marker_cols,]))
})

# ===== 2. DATA TRANSFORMATION =====

# Arcsinh transformation (standard for CyTOF/CODEX)
# Cofactor of 5 is typical, adjust if needed
cofactor <- 5
expr_transformed <- map(expr_matrix, function(x) {
  asinh(x / cofactor)
})
  
# ===== 3. UNSUPERVISED CLUSTERING =====

# Phenograph clustering
set.seed(42)
if (!file.exists(file.path("results", "results_phenograph_results_codex.rds"))) {
  phenograph_result <- map(expr_transformed, function(x) {
    tryCatch({
    Rphenograph(x, k = 30)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
  
  saveRDS(phenograph_result, file.path("results", "results_phenograph_results_codex.rds"))
} else {
  phenograph_result <- readRDS(file.path("results", "results_phenograph_results_codex.rds"))
}

## cluster the data
clusters <- map(phenograph_result, function(x) {
  tryCatch({
  membership(x[[2]])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# Add clusters to data
data <- map(1:length(data), function(x) {
  tryCatch({
     cluster <- as.factor(clusters[[x]])
    names(cluster) <- colnames(data[[x]])
    data[[x]] <- data[[x]] %>% 
      AddMetaData(cluster, "cluster")
    
    data[[x]]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})


# ===== 4. DIMENSIONALITY REDUCTION =====

# UMAP for visualization
set.seed(42)
if (!file.exists(file.path("results", "umap_results_codex.rds"))) {
  umap_result <- map(expr_transformed, function(x) {
    tryCatch({
      umap(x, n_neighbors = 15, min_dist = 0.1)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      })
  saveRDS(umap_result, file.path("results", "umap_results_codex.rds"))
} else {
  umap_result <- readRDS(file.path("results", "umap_results_codex.rds"))
}

## add umap coordinates
data <- map(1:length(data), function(x) {
  tryCatch({
    .df <- data.frame(UMAP1=umap_result[[x]][,1],
                      UMAP2=umap_result[[x]][,2])
    rownames(.df) <- colnames(data[[x]])
    
    data[[x]] <- data[[x]] %>% 
      AddMetaData(.df)
    
    data[[x]]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# Visualize clusters
dir.create(file.path("plots","umap","phenograph"))

walk(names(data), function(x) {
  tryCatch({
    plt <- ggplot(data[[x]][[]], aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme_minimal() +
    labs(title = "Cell Clusters in UMAP Space") +
    theme(legend.position = "right")
  
    print(plt)
  ggsave(file.path("plots","umap","phenograph",paste0(x,"_umap_clusters.pdf")), width = 10, height = 8)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# ===== 5. EXAMINE MARKER EXPRESSION BY CLUSTER =====

# Calculate mean expression per cluster
dir.create(file.path("plots","heatmaps","phenograph"))
    
walk(1:length(data), function(x) {
        tryCatch({  
          cluster_means <- as.data.frame(expr_transformed[[x]]) %>% 
            bind_cols(data.frame(cluster=data[[x]]$cluster)) %>% 
            #select(cluster, all_of(marker_cols)) %>%
            group_by(cluster) %>%
            summarise(across(everything(), mean)) %>%
            column_to_rownames("cluster")
          
        # Heatmap of marker expression
        plt <- pheatmap(t(cluster_means),
                 scale = "row",
                 color = viridis(100),
                 cluster_cols = TRUE,
                 cluster_rows = FALSE,
                 main = "Mean Marker Expression by Cluster",
                 filename = "cluster_heatmap.png",
                 width = 12, height = 8)
        
          pdf(file.path("plots","heatmaps","phenograph", paste0(gsub("20250919_MGN_TMA_Akoya_core_|\\.h5ad","",names(data)[x]), "_cluster_mean_expr_heatmap.pdf")))
        print(plt)
        dev.off()
        
        cluster_means
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# ===== 6. CELL TYPE ANNOTATION =====

# Define cell types based on marker expression
# These are meningioma-specific annotations

## celltype layer 1
annotate_cell_type_1 <- function(df) {
  df %>%
    mutate(cell_type_1 = case_when(
     
      # T cells
      CD45 > 3 ~ "Immune",
      
      # Endothelial/ Mural cells
      CD31 > 2 | CD34 > 2 | SMA > 2 ~ "Endo/Mural",
      
      # CNS Cells
      NeuN > 2 | GFAP > 2 | Olig2 > 2 ~ "CNS",
      
      # Proliferating cells
      Ki67 > 2 ~ "Proliferating",
      
      TRUE ~ "Tumor"
    ))
}

## layer 2
annotate_cell_type_2 <- function(df) {
  df %>%
    mutate(cell_type_2 = case_when(
      
      # T cells
      cell_type_1 == "Immune" & CD3e > 2 & CD4 > 2 ~ "CD4_T",
      cell_type_1 == "Immune" & CD3e > 2 & CD8 > 2 ~ "CD8_T",
      cell_type_1 == "Immune" & CD3e > 2 ~ "other_T",
      
      # B cells
      cell_type_1 == "Immune" & CD20 > 1 & Pax5 > 1 ~ "B/Plasma",
      
      # Macrophages
      cell_type_1 == "Immune" & Iba1 > 2 & CD68 > 2.5  ~ "MSR1_TAM",
      cell_type_1 == "Immune" & Iba1 > 2 & P2RY12 < 2.5 ~ "P2RY12lo_TAM",
      cell_type_1 == "Immune" & Iba1 > 2 & P2RY12 > 2.5 ~ "P2RY12hi_TAM/Mg",
      cell_type_1 == "Immune" & Iba1 > 2 ~ "DC/Mono",
      cell_type_1 == "Immune" ~ "Other_Immune",
      
      # Endothelial/Mural cells
      cell_type_1 == "Endo/Mural" & SMA > 2 ~ "Mural",
      cell_type_1 == "Endo/Mural" ~ "Endothelial",
      
      # CNS
      cell_type_1 == "CNS" & NeuN > 2  ~ "Neuron",
      cell_type_1 == "CNS" & GFAP > 2  ~ "Astro",
      cell_type_1 == "CNS" & Olig2 > 2  ~ "Oligo",
      
      # Proliferating cells
      cell_type_1 == "Proliferating" ~ "Proliferating",
       
      TRUE ~ "Tumor"
    ))
}

# Apply annotation
annot <- map(names(data), function(x) {
  tryCatch({  
  .df <- expr_transformed[[x]]
  cell_type_1 <- annotate_cell_type_1(.df)
  cell_type_2 <- annotate_cell_type_2(cell_type_1)
  
  cell_type_2
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

names(annot) <- fls

## inspect
annot2 <- bind_rows(annot, .id = "sample")
annot2$sample <- gsub("20250919_MGN_TMA_Akoya_core_|\\.h5ad","",annot2$sample)

# Review cell type distribution
table(annot2$cell_type_2)
df <- as.data.frame(table(annot2$cell_type_2, annot2$sample))

## add annotations to the data object
data <- map(names(data), function(x) {
  tryCatch({
  data[[x]] %>% 
    AddMetaData(annot[[x]][,c("cell_type_1","cell_type_2")])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

names(data) <- fls

## optional:
# ===== 7. REFINE ANNOTATIONS USING CLUSTER INFORMATION =====

# Check which clusters map to which cell types
cluster_celltype <- cbind(annot2, bind_rows(clusters)) %>%
  group_by(cluster, cell_type) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(cluster, desc(proportion))

print(cluster_celltype)

# Manually assign cell types to clusters based on dominant marker expression
# Adjust these based on your heatmap results
cluster_annotations <- c(
  "1" = "Meningioma_tumor",
  "2" = "CD4_T_cell",
  "3" = "CD8_T_cell",
  "4" = "Macrophage",
  "5" = "B_cell",
  "6" = "Endothelial",
  "7" = "Fibroblast"
  # Add more based on your number of clusters
)

data$cell_type_refined <- cluster_annotations[as.character(data$cluster)]

# ===== 8. VISUALIZE RESULTS =====

# UMAP colored by cell type
dir.create(file.path("plots","umap","celltypes"))

walk(names(data), function(x) {
  tryCatch({
  plt <- ggplot(data[[x]][[]], aes(x = UMAP1, y = UMAP2, color = cell_type_2)) +
    geom_point(size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Cell Types in UMAP Space", color = "Cell Type") +
    theme(legend.position = "right")
  
  print(plt)

  ggsave(file.path("plots","umap","celltypes", paste0(gsub("20250919_MGN_TMA_Akoya_core_|\\.h5ad","",x), "_celltypes.pdf")))
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# Spatial visualization
ggplot(data, aes(x = X, y = Y, color = cell_type_refined)) +
  geom_point(size = 0.3, alpha = 0.8) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Spatial Distribution of Cell Types", color = "Cell Type") +
  theme(legend.position = "right")

ggsave("spatial_celltypes.png", width = 12, height = 10, dpi = 300)

# ===== 9. QUALITY METRICS =====

if (F) {# Calculate proportion of classified cells
classification_rate <- sum(data$cell_type != "Unclassified") / nrow(data)
print(paste("Classification rate:", round(classification_rate * 100, 2), "%"))

# Check marker specificity
marker_specificity <- data %>%
  group_by(cell_type_refined) %>%
  summarise(across(all_of(marker_cols), mean))

print(marker_specificity)
}
# ===== 10. SAVE RESULTS =====

saveRDS(data, file.path("results","codex_MGN_TMA_celltypes_classified.rds"))
