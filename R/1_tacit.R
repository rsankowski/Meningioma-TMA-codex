library(reticulate)
library(anndata)
library(tidyverse)

source(file.path("R","TACIT_function.R"))

python_dir="/Users/romansankowski/anaconda3/bin/python"
use_python(python_dir, required = T)
sys = import("sys")


if (!file.exists(file.path("results","merged_codex_object_meningioma.rds"))) {
  ## load data
  pth <- file.path("data","adata_objects")
  fls <- list.files(pth)
  
  lst <- map(fls, function(x) {
    adata <- anndata::read_h5ad(file.path(pth, x))
    
    counts <- as.data.frame(adata$X)
    
    counts <- counts[ counts$DAPI > quantile(counts$DAPI , 0.05 ) , ]
    counts <- t(as.matrix(counts))
    
    seurat_obj <- CreateSeuratObject(counts = counts) %>% 
      AddMetaData(x, "core")
    
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
  obj <- readRDS(file.path("results","merged_codex_object_meningioma.rds"))
  lst <- readRDS(file.path("results","single_codex_object_list_meningioma.rds"))
}

## load cell reference 
CELLxFEATURE=t(obj[["RNA"]]$counts)
if (!file.exists(file.path("data","TYPExMARKER_transposed.csv"))) {
  TYPExMARKER=read.csv(file.path("data", "cellxmarker_adj.csv"), row.names = 1) %>% t %>% as.data.frame()
  TYPExMARKER$cell_type <- rownames(TYPExMARKER)
  write.csv(TYPExMARKER, file.path("data","TYPExMARKER.csv")) # I manually adjust this file in excel
} else {
  TYPExMARKER=read.csv(file.path("data", "TYPExMARKER_transposed.csv")) 
  }

##make the same rownames
colnames(CELLxFEATURE) <- colnames(TYPExMARKER)[-1]
identical(colnames(CELLxFEATURE),colnames(TYPExMARKER)[-1] )

## set global limits
options(future.globals.maxSize = 1.0 * 1e20)

## run TACIT
TACIT=TACIT(CELLxFEATURE,TYPExMARKER,r=10,p=10)

print(TACIT)

## align with concatenated macrophage reference
TYPExMARKER2=read.csv(file.path("data", "TYPExMARKER_concat_pos_markers_only.csv")) 

## adjust counts
colnames(CELLxFEATURE) <- gsub("\\.|-| ","_", colnames(CELLxFEATURE))
CELLxFEATURE <- CELLxFEATURE[,colnames(TYPExMARKER2)[-1]]

## run TACIT
TACIT2=TACIT(CELLxFEATURE,TYPExMARKER2,r=10,p=10)

table(TACIT2$TACIT)

saveRDS(TACIT2, file.path("results","TACIT_on_the_whole_object.rds"))

## TACIT list of individual TMA cores
tma_celltypes <- map(names(lst), function(x) {
  tryCatch({
    if (grepl("G2", x)) { ## G2 is a core with CNS tissue only
      CELLxFEATURE <- t(lst[[x]][["RNA"]]$counts)
      colnames(CELLxFEATURE) <- gsub("\\.|-| |/","_", colnames(CELLxFEATURE))
      colnames(CELLxFEATURE) <- gsub("Map_2","Map2", colnames(CELLxFEATURE))
      CELLxFEATURE <- CELLxFEATURE[,colnames(TYPExMARKER2)[-1]]
      
      TACIT=TACIT(CELLxFEATURE,TYPExMARKER2[TYPExMARKER2$cell_type != "Meningioma",],r=10,p=10)
      as.data.frame(TACIT)
    } else {
      CELLxFEATURE <- t(lst[[x]][["RNA"]]$counts)
      colnames(CELLxFEATURE) <- gsub("\\.|-| |/","_", colnames(CELLxFEATURE))
      colnames(CELLxFEATURE) <- gsub("Map_2","Map2", colnames(CELLxFEATURE))
      CELLxFEATURE <- CELLxFEATURE[,colnames(TYPExMARKER2)[-1]]
      
      TACIT=TACIT(CELLxFEATURE,TYPExMARKER2[TYPExMARKER2$cell_type != "Astro_Oligo_Neurons",],r=10,p=10)
      as.data.frame(TACIT)
    }
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})
