library(tidyverse)
library(ggthemes)
library(viridis)
library(patchwork)
library(tessera)
library(Seurat)

obj <- readRDS(file.path("results","merged_nimbus_with_celltypes_codex_object_meningioma.rds"))
data <- readRDS(file.path("results","codex_MGN_TMA_celltypes_classified.rds"))

## set parameters
fig.size <- function(h, w) {
  options(repr.plot.height = h, repr.plot.width = w)
}

verbose = TRUE
show_plots = TRUE

###### STEP 0 ######
npcs = 5
## Graph pruning
prune_thresh_quantile = 0.95
prune_min_cells = 5


###### STEP 1: GRADIENTS ######
smooth_distance = c('none', 'euclidean', 'projected', 'constant')[3] 
smooth_similarity = c('none', 'euclidean', 'projected', 'constant')[3] 


###### STEP 2: DMT ######
## ... no options


###### STEP 3: AGGREGATION ######
max_npts = 50
min_npts = 5
alpha = 1 ## 0.2 = conservative merging, 2 = liberal merging 

## subset data
obj2 <- obj[,obj$Histology =="Meningioma" & obj$Fixation_protocol=="FFPE"]

## subset counts
counts <- obj2[["RNA"]]$counts

## get coordinates
## calculate bounding boxes
bboxes <- data.frame(sample=unique(obj2$fov),
                     y=gsub("-.*","",unique(obj2$fov)),
                     x=gsub(".*-","",unique(obj2$fov)))

bboxes$y <- as.numeric(as.factor(bboxes$y))

## add x intercept
bboxes$x_int <- (as.numeric(bboxes$x) - 2) * 4250
bboxes$y_int <- (bboxes$y - 1) * 3680

coords <- map(names(obj2@images), function(x) {
  tryCatch({  
  data.frame(obj2@images[[x]]@coordinates, image=x)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

coords <- coords %>% 
  bind_rows()
coords2 <- coords %>% 
  mutate(fov=gsub("_.*","",coords$cells)) %>% 
  distinct(fov, image)

rownames(bboxes) <- coords2$image

## re-run coord_extraction with intercepts
coords <- map(names(obj2@images), function(x) {
  tryCatch({  
    .df = obj2@images[[x]]@coordinates
    
    .df$x = .df$x + bboxes[x,]$x_int
    .df$y = .df$y + bboxes[x,]$y_int
    
    .df
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

coords <- coords %>% 
  bind_rows()

meta_data <- data.frame(X= coords[colnames(obj2),]$x,
                        Y= coords[colnames(obj2),]$y,
                        type=obj2$cell_type_2,
                        fov=obj2$fov)

## plot metadata
fig.size(8, 8)
ggplot() + 
  geom_point(data = meta_data, aes(X, Y, color = type), size=.5) + 
  theme_void() + 
  scale_color_tableau(palette = "Tableau 20") + 
  coord_sf(expand = FALSE) + 
  NULL

## run tessera
meta_vars_include = c('type')

## split the data to not get error message for the graph

# Split spatially and process separately
set.seed(79104)
n_chunks <- 50
meta_data$spatial_chunk <- cut(meta_data$X, breaks = n_chunks, labels = FALSE)

results_list <- list()

for (i in 1:n_chunks) {
  chunk_idx <- which(meta_data$spatial_chunk == i)
  
  res_chunk = GetTiles(
    X = meta_data[chunk_idx,]$X, 
    Y = meta_data[chunk_idx,]$Y, 
    counts = counts[,chunk_idx], 
    meta_data = meta_data[chunk_idx,], 
    meta_vars_include = meta_vars_include,
  )
  
  results_list[[i]] <- res_chunk
}

# Combine results afterward

res = GetTiles(
  X = meta_data$X, 
  Y = meta_data$Y, 
  counts = counts, 
  meta_data = meta_data, 
  meta_vars_include = meta_vars_include,
)

dmt = res$dmt
aggs = res$aggs








## Prepare data
meta_vars_include = c('type')
dmt = init_data(meta_data$X,meta_data$Y, counts, meta_data, meta_vars_include)
dmt = prune_graph(dmt, thresh_quantile = prune_thresh_quantile, mincells = prune_min_cells) 

## add triangles
dmt = add_exterior_triangles(dmt)

fig.size(15, 15)
# fig.size(20, 20)

if (show_plots) {    
  ggplot() + 
    
    # ## big data 
    # geom_point(data = dmt$pts, aes(X, Y), shape = '.', alpha = .1) + 
    # geom_segment(data = dmt$edges[boundary == TRUE], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'red', lwd = .2) + 
    
    ## small data 
    geom_segment(data = dmt$edges, aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'black', lwd = .2) + 
    geom_segment(data = dmt$edges[boundary == TRUE], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'red', lwd = .2) + 
    geom_point(data = dmt$pts, aes(X, Y), size = .5) + 
    
    theme_void(base_size = 20) + 
    coord_cartesian(expand = FALSE) + 
    labs(title = 'Pruned adjacency graph') + 
    NULL
}

## pca
dmt$udv_cells = do_pca(dmt$counts, npcs)

# Step 1: compute gradients on all data structures
field = compute_gradients(dmt, smooth_distance, smooth_similarity)
field = compress_gradients_svd(field)

if (show_plots) {    
  
  len_plot_constant = .8
  fig.size(15, 15)
  # fig.size(20, 20)
  ggplot() + 
    geom_point(data = dmt$pts, aes(X, Y), size = .5) + 
    geom_segment(data = dmt$edges[boundary == TRUE, ], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'red') + 
    
    ## Triangle Gradients
    geom_segment(
      data = data.table(dmt$tris, field$tris_svd), 
      aes(
        x=X-len_plot_constant*(len_grad+len_ortho)*dx_ortho, 
        y=Y-len_plot_constant*(len_grad+len_ortho)*dy_ortho, 
        xend=X+len_plot_constant*(len_grad+len_ortho)*dx_ortho, 
        yend=Y+len_plot_constant*(len_grad+len_ortho)*dy_ortho
      ), 
      linewidth = .4, alpha = 1, 
      color = 'blue'
    ) + 
    
    theme_void() + 
    coord_fixed(expand = FALSE) + 
    NULL
}

# Step 2: DMT

## compute f
dmt = dmt_set_f(dmt, field)

if (show_plots) {    
  ntri = max(which(dmt$tris$external == FALSE))
  i = Matrix::t(dmt$tri_to_pt[1:ntri, ])@i+1
  plt_df = data.table(
    X = dmt$pts$X[i],
    Y = dmt$pts$Y[i],
    f = rep(dmt$tris$f[1:ntri], each = 3)
  )[
    , id := rep(1:ntri, each = 3)
  ][]
  
  
  fig.size(15, 15)
  ggplot() + 
    geom_polygon(data = plt_df, aes(X, Y, group = id, fill = f, color = f)) + 
    theme_void() + 
    coord_fixed(expand = FALSE) + 
    scale_fill_viridis() + 
    scale_color_viridis() + 
    NULL
}

## forests
dmt$prim = do_primary_forest(dmt)
dmt$dual = do_dual_forest(dmt)

if (show_plots) {    
  fig.size(15, 15)
  ggplot() +     
    ## primary forest
    geom_point(data = dmt$tris[dmt$dual$maxima, ], aes(X, Y), color = 'blue', size = 2) + 
    geom_segment(data = dmt$dual$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'blue') + 
    
    ## primary forest
    geom_point(data = dmt$pts[dmt$prim$minima, ], aes(X, Y), color = 'red', size = 2) + 
    geom_segment(data = dmt$prim$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'red') + 
    
    theme_void() + 
    coord_cartesian(expand = FALSE) + 
    NULL
}


## extract epaths
dmt$e_sep = dmt_get_separatrices(dmt)

if (show_plots) {    
  fig.size(15, 15)
  ggplot() + 
    
    geom_segment(data = dmt$edges[dmt$e_sep, ], aes(x = x0_tri, y = y0_tri, xend = x1_tri, yend = y1_tri), lwd = 1, color = 'blue') + 
    geom_segment(data = dmt$edges[boundary == TRUE], aes(x = x0_pt, y = y0_pt, xend = x1_pt, yend = y1_pt), color = 'blue', lwd = 1) + 
    
    ## primary forest
    geom_point(data = dmt$pts[dmt$prim$minima, ], aes(X, Y), color = 'red', size = 2) + 
    geom_segment(data = dmt$prim$edges, aes(x=x0, y=y0, xend=x1, yend=y1), color = 'red') + 
    
    theme_void() + 
    coord_cartesian(expand = FALSE) + 
    NULL
}

## extract tiles
dmt = dmt_assign_tiles(dmt)
aggs = dmt_init_tiles(dmt)

if (show_plots) {    
  set.seed(2)
  fig.size(15, 20)
  ggplot() + 
    geom_sf(data = aggs$meta_data$shape) + 
    # geom_point(data = dmt$pts, aes(X, Y, color = factor(agg_id, sample(nrow(aggs$meta_data)))), size = 1) + 
    # scale_color_tableau() + 
    theme_void() + 
    coord_sf(expand = FALSE) + 
    # coord_cartesian(expand = FALSE) + 
    guides(color = 'none') + 
    
    NULL 
}

if (show_plots) {    
  set.seed(2)
  fig.size(15, 20)
  ggplot() + 
    geom_sf(data = aggs$meta_data$shape) + 
    geom_point(data = dmt$pts, aes(X, Y, color = type), size=.5) + 
    theme_void() + 
    coord_sf(expand = FALSE) + 
    scale_color_manual(values = colors_fig) + 
    guides(color = 'none') + 
    NULL 
}

# Step 3: Aggregation

## Merge main
aggs = init_scores(aggs, agg_mode=2, alpha=alpha, max_npts=max_npts)
aggs = merge_aggs(aggs, agg_mode=2, max_npts=max_npts)
dmt = update_dmt_aggid(dmt, aggs)
aggs = update_agg_shapes(dmt, aggs)

## Merge small outliers
aggs = init_scores(aggs, agg_mode=3, alpha=alpha, min_npts=min_npts)
aggs = merge_aggs(aggs, agg_mode=3, min_npts=min_npts)
dmt = update_dmt_aggid(dmt, aggs)
aggs = update_agg_shapes(dmt, aggs)

## Final tiles 
if (show_plots) { 
  purrr::map(1:5, function(i) {
    ggplot(cbind(aggs$meta_data, val=aggs$pcs[, i])) + 
      geom_sf(aes(geometry = shape, fill = val)) + 
      theme_void(base_size = 16) + 
      coord_sf(expand = FALSE) + 
      scale_fill_gradient2_tableau() + 
      guides(color = 'none') + 
      labs(title = paste0('PC', i)) + 
      NULL 
  }) %>% 
    purrr::reduce(`|`)
}

# Results

## Aggregates 
## The primary output is the tiles. Each tile has a row in the meta_data table: 
##  - npts denotes the number of cells in the tile. 

head(aggs$meta_data)


## We also have pooled gene counts, for differential gene expression analysis. 
aggs$counts[1:5, 1:5]

## And we have PCA embeddings for the tiles. 
head(aggs$pcs)


