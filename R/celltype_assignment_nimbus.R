# Load required libraries

## Nimbus cell type assignment
library(tidyverse)
library(Rphenograph)
library(uwot)
library(pheatmap)
library(viridis)
library(Seurat)
library(ggthemes)

if (!file.exists(file.path("results","merged_nimbus_codex_object_meningioma.rds"))) {
  tma_batch_1 <- read.csv(file.path("data","nimbus","nimbus_cell_table_with_coordinates.csv"))
  tma_batch_2 <- read.csv(file.path("data","nimbus","nimbus_cell_table_with_coordinates_batch2.csv"))
  
  tma_merged <- rbind(tma_batch_1, tma_batch_2)
  tma_merged$fov <- gsub("20251009_MGN_TMA_Akoya_core_","",tma_merged$fov)
  
  lst <- split(tma_merged, tma_merged$fov)
  nms <- map_vec(lst, function(x) unique(x$fov))
  names(lst) <- nms
  
  lst <- map(lst, function(x) {
    counts <- x[, 3:70]
    rownames(counts) <- x$label
    
    meta <- x[, c(2,71:78)]
    rownames(meta) <- x$label
    
    ## remove dim cells
    counts <- counts[ counts$DAPI > quantile(counts$DAPI , 0.05 ) , ]
    meta <- meta[rownames(counts),]
    
    seurat_obj <- CreateSeuratObject(counts = t(counts)) %>% 
      NormalizeData()
    
    VariableFeatures(seurat_obj) <- rownames(seurat_obj)[-1]
    
    rownames(meta) <- colnames(seurat_obj)
    
    seurat_obj <- seurat_obj %>% 
      AddMetaData(metadata = meta) %>% 
      ScaleData()
    
    coords <- data.frame(x=meta$centroid_x,
                         y=meta$centroid_y)
    
    rownames(coords) <- colnames(seurat_obj)
    
    seurat_obj@images$image =  new(
      Class = 'SlideSeq',
      assay = "Spatial",
      key = "image_",
      coordinates = coords)
    
    seurat_obj
    
  })
  
  obj <- merge(lst[[1]],lst[2:134], add.cell.ids = nms)
  
  obj <- obj %>% JoinLayers()
  
  saveRDS(nms, file.path("results","codex_object_names.rds"))
  ## save object
  saveRDS(obj, file.path("results","merged_nimbus_codex_object_meningioma.rds"))
  saveRDS(lst, file.path("results","single_nimbus_codex_object_list_meningioma.rds"))
} else {
  ## load data
  pth <- file.path("data","adata_objects")
  
  obj <- readRDS(file.path("results","merged_nimbus_codex_object_meningioma.rds"))
  data <- readRDS(file.path("results","single_nimbus_codex_object_list_meningioma.rds"))
  nms <- readRDS(file.path("results","codex_object_names.rds"))
}

# ===== 1. LOAD AND PREPARE DATA =====

# Assuming your data is in a CSV with columns:
# - cell_id, X, Y (spatial coordinates)
# - marker columns (protein expression values)

# Extract marker expression matrix
# Adjust column names based on your specific markers
marker_cols <- c("CD45", "CD3e", "CD4", "CD8", "CD20", "CD79a", "CD68", "CD163", "CD206","SMA", "CD31","CD34","Ki67", "Collagen-IV", "Claudin-5","Pax5","Iba1","GFAP", "Podoplanin", "FOXP3", "NeuN", "Neurofilament", "CD38", "Olig2", "CSF1R", "P2RY12","TMEM119","Map2","GNLY","MPO","AXL")

# Keep only markers present in your data
marker_cols <- marker_cols[marker_cols %in% rownames(data[[1]])]

expr_matrix <- map(data, function(x) {
  as.data.frame(t(x[["RNA"]]$scale.data[marker_cols,]))
})

# ===== 3. UNSUPERVISED CLUSTERING =====

# Phenograph clustering
set.seed(42)
if (!file.exists(file.path("results", "results_nimbus_phenograph_results_codex.rds"))) {
  phenograph_result <- map(expr_matrix, function(x) {
    tryCatch({
    Rphenograph(x, k = 30)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    })
  
  saveRDS(phenograph_result, file.path("results", "results_nimbus_phenograph_results_codex.rds"))
} else {
  phenograph_result <- readRDS(file.path("results", "results_nimbus_phenograph_results_codex.rds"))
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

## add object names
names(data) <- nms

# ===== 4. DIMENSIONALITY REDUCTION =====

# UMAP for visualization
set.seed(42)
if (!file.exists(file.path("results", "umap_nimbus_results_codex.rds"))) {
  umap_result <- map(expr_matrix, function(x) {
    tryCatch({
      umap(x, n_neighbors = 15, min_dist = 0.1)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      })
  saveRDS(umap_result, file.path("results", "umap_nimbus_results_codex.rds"))
} else {
  umap_result <- readRDS(file.path("results", "umap_nimbus_results_codex.rds"))
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
dir.create(file.path("plots","umap","phenograph_nimbus"))

names(data) <- nms

walk(names(data), function(x) {
  tryCatch({
    plt <- ggplot(data[[x]][[]], aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme_minimal() +
    labs(title = "Cell Clusters in UMAP Space") +
    theme(legend.position = "right")
  
    print(plt)
  ggsave(file.path("plots","umap","phenograph_nimbus",paste0(x,"_umap_clusters.pdf")), width = 10, height = 8)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# ===== 5. EXAMINE MARKER EXPRESSION BY CLUSTER =====

# Calculate mean expression per cluster
dir.create(file.path("plots","heatmaps","phenograph_nimbus"))
    
walk(1:length(data), function(x) {
        tryCatch({  
          cluster_means <- as.data.frame(expr_matrix[[x]]) %>% 
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
        
          pdf(file.path("plots","heatmaps","phenograph_nimbus", paste0(gsub("20250919_MGN_TMA_Akoya_core_|\\.h5ad","",names(data)[x]), "_cluster_mean_expr_heatmap.pdf")))
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
     
      # CD45 low Myeloid cells
      Iba1 > 1 ~ "Myeloid",
      
      # T cells
      CD45 > 2 ~ "Immune",
      
      # Endothelial/ Mural cells
      CD31 > 2 | CD34 > 2 | SMA > 2 ~ "Endo/Mural",
      
      # CNS Cells
      NeuN > 4 | GFAP > 4 | Olig2 > 4 ~ "Other",
      
      # Proliferating cells
      #Ki67 > .05 ~ "Proliferating",
      
      TRUE ~ "Tumor"
    ))
}

## layer 2
annotate_cell_type_2 <- function(df) {
  df %>%
    mutate(cell_type_2 = case_when(
      
      # T cells
      cell_type_1 == "Immune" & CD3e > 1 & CD4 > 1 ~ "CD4_T",
      cell_type_1 == "Immune" & CD3e > 1 & CD8 > 1 ~ "CD8_T",
      cell_type_1 == "Immune" & CD3e > 1 ~ "other_T",
      cell_type_1 == "Immune" & MPO > 4 ~ "Granulo",
      cell_type_1 == "Immune" & CD163 > 1  ~ "Trans_Mono",
      cell_type_1 == "Immune" & CD206 > 1  ~ "P2RY12lo_TAM",
      cell_type_1 == "Immune" ~ "Other_Immune",
      
      # B cells
      cell_type_1 == "Immune" & CD20 > 1 | cell_type_1 == "Immune" & CD79a  > 1 | cell_type_1 == "Immune" & Pax5 > 1 ~ "B_Plasma",
      
      # Macrophages
      cell_type_1 == "Myeloid" & MPO > 4 ~ "Granulo",
      cell_type_1 == "Myeloid" & Ki67 > 3  ~ "Prolif_TAM",
      cell_type_1 == "Myeloid" & P2RY12 < 1 ~ "P2RY12lo_TAM",
      cell_type_1 == "Myeloid" & P2RY12 > 1  ~ "P2RY12hi_TAM",
      cell_type_1 == "Myeloid" & CD206 > 1  ~ "P2RY12lo_TAM",
      cell_type_1 == "Myeloid" & CD163 > 1  ~ "Trans_Mono",
      cell_type_1 == "Myeloid" ~ "Mono_DC",
      
      # Endothelial/Mural cells
      cell_type_1 == "Endo/Mural" & SMA > 1 ~ "Mural",
      cell_type_1 == "Endo/Mural" ~ "Endothelial",
      
      # CNS
      cell_type_1 == "Other" ~ "Other",
      
      # Proliferating cells
      #cell_type_1 == "Proliferating" ~ "Proliferating",
       
      TRUE ~ "Tumor"
    ))
}

names(data) <- nms

# Apply annotation
annot <- map(names(data), function(x) {
  tryCatch({  
  .df <- expr_matrix[[x]]
  cell_type_1 <- annotate_cell_type_1(.df)
  cell_type_2 <- annotate_cell_type_2(cell_type_1)
  
  cell_type_2
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

names(annot) <- nms

## inspect
annot2 <- bind_rows(annot, .id = "sample")

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

names(data) <- nms

# ===== 8. VISUALIZE RESULTS =====

## define cell type colors
my_colors <- readRDS(file.path("/Users/romansankowski/Library/CloudStorage/GoogleDrive-romansankowski@gmail.com/Other computers/My MacBook Pro/single_cell_analysis/Visium_HD_meningioma_TMA/results/celltype_colors.rds"))
names(my_colors)[5] <- "Prolif_TAM"
names(my_colors)[10] <- "B_Plasma"
names(my_colors)[8] <- "Tumor"
names(my_colors)[9] <- "Endothelial"
names(my_colors)[20] <- "Mural"
names(my_colors)[6] <- "Mono_DC"
names(my_colors)[4] <- "CD4_T"
names(my_colors)[13] <- "CD8_T"
names(my_colors)[14] <- "other_T"
names(my_colors)[15] <- "Neuroectoderm"
names(my_colors)[16] <- "Proliferating"
names(my_colors)[16] <- "Other_Immune"


my_colors

# ===== EXAMINE UMAP COLORED BY CELL TYPE =====
# UMAP colored by cell type
dir.create(file.path("plots","umap","celltypes_nimbus"))

walk(names(data), function(x) {
  tryCatch({
  plt <- ggplot(data[[x]][[]], aes(x = UMAP1, y = UMAP2, color = cell_type_2)) +
    geom_point(size = 1) + #, alpha = 0.7
    theme_minimal() +
    labs(title = "Cell Types in UMAP Space", color = "Cell Type") +
    theme(legend.position = "right") +
    scale_color_manual(values = my_colors[unique(data[[x]][[]]$cell_type_2)]) +
    theme_void()
  
  print(plt)

  ggsave(file.path("plots","umap","celltypes_nimbus", paste0(x, "_celltypes.pdf")))
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

 # ===== EXAMINE MARKER EXPRESSION BY CELL TYPE =====

# Calculate mean expression per cluster
dir.create(file.path("plots","heatmaps","celltypes_nimbus"))

walk(1:length(data), function(x) {
  tryCatch({  
    cluster_means <- as.data.frame(expr_matrix[[x]]) %>% 
      bind_cols(data.frame(cluster=data[[x]]$cell_type_2)) %>% 
      #select(cluster, all_of(marker_cols)) %>%
      group_by(cluster) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("cluster")
    
    # Heatmap of marker expression
    plt <- pheatmap(t(cluster_means),
                    scale = "row",
                    color = viridis(100),
                    cluster_cols = TRUE,
                    cluster_rows = TRUE,
                    main = "Mean Marker Expression by Cluster",
                    filename = "cluster_heatmap.png",
                    width = 12, height = 8,
                    cellheight=12, cellwidth = 12,
                    border_color = NA
                    )
    
    pdf(file.path("plots","heatmaps","celltypes_nimbus", paste0(names(data)[x], "_cluster_mean_expr_heatmap.pdf")))
    print(plt)
    dev.off()
    
    cluster_means
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

# ===== SAVE RESULTS =====

saveRDS(data, file.path("results","codex_MGN_TMA_celltypes_classified.rds"))

## add diagnosis metadata
meta <- read.csv(file.path("data","tma_3_24_metadata.csv"))
colnames(meta)[1] <- "sample"

## add cell types to seurat object 
rnames <- rownames(annot2)

rnames <-  gsub("\\..*", "", rnames)
rnames <- paste(annot2$sample,rnames,sep="_")

## sanity check
sum(rnames%in% colnames(obj)) == ncol(obj)
identical(rnames, colnames(obj))

rownames(annot2) <- rnames
annot2 <- annot2 %>% 
  left_join(meta)

## add celltype annotation to seurat object
obj <- obj %>% 
  AddMetaData(annot2[,31:37])

##save object
saveRDS(obj, file.path("results","merged_nimbus_with_celltypes_codex_object_meningioma.rds"))
