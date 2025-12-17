library(reticulate)
library(anndata)
library(tidyverse)

python_dir="/Users/romansankowski/anaconda3/bin/python"
use_python(python_dir, required = T)
sys = import("sys")



if (!file.exists(file.path("results","xmin_ymin_codex_object_meningioma.rds"))) {
  pth <- file.path("data","adata_objects")
  fls <- list.files(pth)
  
  lst <- map(fls, function(x) {
    tryCatch({ 
    adata <- anndata::read_h5ad(file.path(pth, x))
    
    coords <- adata$obsm$spatial %>%
      as.data.frame()
    
    rownames(coords) <- adata$obs_names
    colnames(coords) <- c("x", "y")
    
    .df <- data.frame(xmin=min(coords$x), ymin=min(coords$y))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  })
  
  names(lst) <- fls
  xy_df <- lst %>% 
    bind_rows(.id="fov")
  
  xy_df$fov <- gsub("(20251009_MGN_TMA_Akoya_core_|\\.h5ad)","", xy_df$fov)
  
  ## save object
  saveRDS(xy_df, file.path("results","xmin_ymin_codex_object_meningioma.rds"))
} else {
  ## load data
  xy_df <- readRDS(file.path("results","xmin_ymin_codex_object_meningioma.rds"))
  
}

