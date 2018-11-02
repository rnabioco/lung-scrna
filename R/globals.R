# globals shared across markdown docs

library(Seurat)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Matrix)
library(matrixStats)
library(purrr)
library(R.utils)
library(viridis)
library(ComplexHeatmap)
library(doParallel)
library(itertools)
library(gridExtra)
library(here)
#### Paths ####

project_dir <- here() 
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
docs_dir <- file.path(project_dir, "docs")
db_dir <- file.path(project_dir, "dbases")

# vector of figure paths
figs_dir <-  file.path(results_dir, "Figures") %>%
  dir(pattern = "Figure_[1-4]$",
      include.dirs = TRUE,
      full.names = T)


##### Functions ####

#' When writing out excel workbooks using openxlsx::write.xlsx()
#' this function will set the class attributes for a column, which
#' enforces a column type in the resulting xlsx file. 
#' Useful for avoid gene names being clobbered to dates and 
#' setting scientific number formatting

set_xlsx_class <- function(df, col, xlsx_class){
  for(i in seq_along(col)){
    class(df[[col[i]]]) <- xlsx_class
  }
  df
}

####
#' set up defaults for tsne plot
plot_tsne <- function(seurat_obj,
                      ident = "ident",
                      pt.size = 0.5,
                      pt.alpha = 1,
                      label_text = FALSE,
                      label.size = 6,
                      label.color = "grey",
                      .cols = NULL,
                      legend_names = NULL,
                      cell_filter = NULL,
                      legend_order = NULL,
                      ...){
  
  
  tsne_dat <- seurat_obj@dr$tsne@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell")
  
  xlimits <- c(min(tsne_dat$tSNE_1),
               max(tsne_dat$tSNE_1))
  ylimits <- c(min(tsne_dat$tSNE_2),
               max(tsne_dat$tSNE_2))
  
  if (!is.null(cell_filter)){
    seurat_obj <- SubsetData(seurat_obj, cells.use = cell_filter)
  }
  
  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")
  tsne_dat <- left_join(mdata, tsne_dat, by = "cell")
  
  p <- ggplot(tsne_dat, 
              aes(tSNE_1, tSNE_2)) +
    ylim(ylimits) +
    xlim(xlimits)
  
  if (is.null(ident)){
    if (is.null(.cols)){
      p <- p + geom_point(
        size = pt.size,
        alpha = pt.alpha) 
    } else {
      p <- p + geom_point(
        color = .cols[1],
        size = pt.size,
        alpha = pt.alpha) 
    }
  } else {
    p <- p + geom_point(aes_string(color = ident),
                        size = pt.size,
                        alpha = pt.alpha)
  }
  
  p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(legend.title = element_blank())
  
  if (label_text) {
    tsne_mean_dat <- tsne_dat %>% 
      group_by_at(vars(one_of(ident))) %>% 
      summarize(mean_dim_1 = mean(tSNE_1), 
                mean_dim_2 = mean(tSNE_2))
    
    p <- p + 
      geom_text(data = tsne_mean_dat, 
                aes_string(x = "mean_dim_1",
                           y = "mean_dim_2",
                           label = ident),
                size = label.size, color = label.color)
  }
  
  ## handle colors
  if (is.null(.cols)){
    if (is.null(legend_names)){
      p <- p + scale_color_viridis(discrete = T,
                                   direction = -1)
    } else {
      p <- p + scale_color_viridis(labels = legend_names,
                                   discrete = T,
                                   direction = -1)
    }
  } else {
    if (is.null(legend_names) & is.null(legend_order)){
      p <- p + scale_color_manual(values = .cols)
    } else if(is.null(legend_order)) {
      p <- p + scale_color_manual(labels = legend_names, values = .cols)
    } else {
      if(is.null(legend_names)) legend_names <- legend_order
      p <- p + scale_color_manual(values = .cols, 
                                  breaks = legend_order,
                                  labels = legend_names)
    }
  }
  p
}


### feature

discrete_palette_default <- c(brewer.pal(12, "Paired"),
                              brewer.pal(9, "Set1"),
                              brewer.pal(8, "Set2"),
                              brewer.pal(8, "Dark2"))
                              
#' Plot cells in reduced dimensionality 2D space 
#' 
#' @description Cells can be colored by gene or feature in meta.data dataframe
#' 
#' @param seurat_obj object of class Seurat 
#' @param feature feature to plot, either gene name or column in seurat_obj@meta.data
#' @param plot_dat supplemental data.frame containing feature to plot. 
#' Must have a column named cell that contains matching colnames in seurat_obj@data
#' @param pt_size size of points produced by geom_point
#' @param pt_alpha alpha value for points plotted by geom_point
#' @param label_text if TRUE display feature labels on plot
#' @param label_size size of label text
#' @param label_color color of label text
#' @param .cols vector of colors to use for plot. 
#' @param cell_filter character vector of cell names to include in plot
#' @param palette_type color palette type to use (either viridis, brewer, or cloupe)
#' defaults to using cellranger loupe-like colors
#' @param col_pal palette name to use if palette_type is brewer
#' @param max_y maximum feature value to set scale to. Defaults to max of the feature
#' @param legend_title string to supply for title for the legend
#' @param embedding dimensionality reduction to extract from seurat_obj. Can be any
#' dr method present in seurat_obj@dr (e.g. umap, pca, tsne). defaults to tsne
#' 
plot_feature <- function(seurat_obj,
                         feature = NULL,
                         plot_dat = NULL,
                         pt_size = 0.001,
                         pt_alpha = 1,
                         label_text = FALSE,
                         label_size = 6,
                         label_color = "grey",
                         .cols = NULL,
                         cell_filter = NULL,
                         palette_type = "cloupe",
                         col_pal = "Reds",
                         max_y = NULL,
                         legend_title = NULL,
                         embedding = "tsne"){
  
  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")
  
  if(!embedding %in% names(seurat_obj@dr)){
    stop(paste0(embedding, " not found in seurat object"))
  }
  
  embed_dat <- seurat_obj@dr[[embedding]]@cell.embeddings %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell")
  
  embed_cols <- colnames(embed_dat)
  xcol <- embed_cols[2]
  ycol <- embed_cols[3]
  
  embed_dat <- left_join(mdata, embed_dat, by = "cell")
  
  if (!is.null(cell_filter)){
    embed_dat <- dplyr::filter(embed_dat,
                               cell %in% cell_filter)
  }
  
  meta_data_col <- feature %in% colnames(embed_dat)
  
  if (!is.null(feature) & !meta_data_col) {
    feature_dat <- FetchData(seurat_obj, feature) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column("cell")
    embed_dat <- left_join(embed_dat, feature_dat, by = "cell")
  }
  
  if (!is.null(plot_dat)){
    embed_dat <- left_join(embed_dat, plot_dat, by = "cell")
  }
  
  color_aes_str <- feature
  
  color_aes_str_q <- quo(color_aes_str)
  embed_dat <- embed_dat %>% arrange_at(.vars = color_aes_str)
  
  p <- ggplot(embed_dat, 
              aes_string(xcol, ycol)) +
    geom_point(aes_string(color = color_aes_str),
               size = pt_size,
               alpha = pt_alpha)
  
  ## discrete or continuous data?
  if (typeof(embed_dat[[feature]]) %in% c(
    "character",
    "logical"
  ) | is.factor(embed_dat[[feature]])) {
    discrete <- T
  } else {
    discrete <- F
  }
  
  ## increase legend size
  if (discrete) {
    p <- p + guides(colour = guide_legend(override.aes = list(size = 4))) +
      theme(legend.title = element_blank())
  }
  
  if (label_text) {
    if(discrete) {
    tsne_mean_dat <- embed_dat %>% 
      group_by_at(vars(one_of(feature))) %>% 
      summarize(med_dim_1 = median(tSNE_1), 
                med_dim_2 = median(tSNE_2))
    
    p <- p + 
      geom_text(data = tsne_mean_dat, 
                aes_string(x = "med_dim_1",
                           y = "med_dim_2",
                           label = feature),
                size = label_size, 
                color = label_color)
    } else {
      warning("label_text not compatible with continuous features")
    }
  }
  
  ## handle legend limit 
  if (is.null(max_y) & !discrete) {
    max_y <- c(0, max(embed_dat[[color_aes_str]]))
  } else if (discrete & is.null(max_y)){
    max_y <- c(NA, NA)
  } 
  
  # loupe-like colors
  cols <- rev(brewer.pal(11, "RdGy")[c(1:5, 7)])
  
  #handle legend name
  if(is.null(legend_title)) legend_title <- color_aes_str
  
  ## handle zero expression
  if (!all(is.na(max_y)) && all(max_y == c(0, 0))){
    p <- p + scale_color_gradient(low = cols[1], high = cols[1], name = legend_title)
    return(p)
  }
  
  ## handle colors
  if (is.null(.cols) && !discrete){
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis(discrete = F,
                                   direction = -1,
                                   option = col_pal,
                                   limits = max_y, name = legend_title)
    } else if (palette_type == "brewer") {
      p <- p + scale_color_distiller(limits = max_y,
                                     palette = col_pal,
                                     direction = 1, name = legend_title)
    } else if (palette_type == "cloupe") {
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = cols, name = legend_title)
    }
  } else if (!is.null(.cols) && !discrete){
      p <- p + scale_color_gradientn(limits = max_y,
                                     colors = .cols, name = legend_title)
  } else {
    
    if(!is.null(.cols)) {
      # use colors provided
      p <- p + scale_color_manual(
        values = .cols,
        name = legend_title
      )
    } else {
      p <- p + scale_color_manual(
        values = discrete_palette_default,
        name = legend_title
      )
    }
  } 
  p
}


### violin 

#' set up defaults for violin plot
plot_violin <- function(df, .x, .y, 
                        .fill = "",
                        .size = 0.50,
                        .width = 1, 
                        .scale = "width",
                        .alpha = 1,
                        cols = scale_fill_viridis(discrete = T),
                        single_col = NULL,
                        jitter = F,
                        rotate_x_text = TRUE){
  
  p <- ggplot(df, aes_string(x = .x, y = .y))
  
  if (jitter){
    p <- p  + geom_jitter(size = 0.1, alpha = 0.2, color = "black")
  }
  
  if (!is.null(single_col)){
    p <- p + 
      geom_violin(size = .size,
                  scale = .scale,
                  fill = single_col,
                  alpha = .alpha)
  } else {
    p <- p +
      geom_violin(aes_string(fill = .fill), 
                  size = .size,
                  scale = .scale,
                  alpha = .alpha) +
      cols 
  }
  
  if(rotate_x_text){
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                            hjust = 1, 
                                            vjust = 0.5))
  }
  p <- p + theme(legend.title = element_blank()) 
  p
}


#' set up defaults for violin plot
plot_box <- function(df, .x, .y, 
                     .fill = "",
                     .size = 0.50,
                     .width = 1, 
                     ...){
  
  ggplot(df, aes_string(x = .x, y = .y)) +
    geom_boxplot(aes_string(fill = .fill), 
                 coef = 1e6
    ) +
    scale_fill_viridis(discrete = T) +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5),
          legend.title = element_blank()) 
}



#### 

plot_expt_tsnes <- function(seurat_obj = atcells, 
                            expt1_cells = expt1,
                            expt2_cells = expt2,
                            expt1_name = "Experiment 1",
                            expt2_name = "Experiment 2",
                            out_name = "all",
                            tsne_cols = brewer.pal(11, "Paired"),
                            outdir = file.path("fig2", "tsne_per_expt")) {
  
  paired_tsne_theme  <- theme(axis.title.x = element_blank(),
                              axis.text = element_text(size = 10),
                              axis.title = element_text(size = 10),
                              legend.pos = "top")
  
  limit_setting <- lims(x = c(min(seurat_obj@dr$tsne@cell.embeddings[, 1]), 
                              max(seurat_obj@dr$tsne@cell.embeddings[, 1])),
                        y = c(min(seurat_obj@dr$tsne@cell.embeddings[, 2]),
                              max(seurat_obj@dr$tsne@cell.embeddings[, 2])))
  
  tsne_expts <- plot_tsne(seurat_obj, 
                          ident = "expt",
                          pt.alpha = 0.5,
                          .cols = tsne_cols, 
                          cell_filter = c(expt1_cells, expt2_cells),
                          legend_names = c(expt1_name,
                                           expt2_name)) +
    paired_tsne_theme
  
  
  tsne_expt1 <- plot_tsne(seurat_obj, 
                          ident = "expt",
                          pt.alpha = 0.5,
                          .cols = tsne_cols[1], 
                          cell_filter = expt1_cells,
                          legend_names = expt1_name) +
    limit_setting +
    paired_tsne_theme
  
  
  tsne_expt2 <- plot_tsne(seurat_obj, 
                          ident = "expt",
                          pt.alpha = 0.5,
                          .cols = tsne_cols[2], 
                          cell_filter = expt2_cells,
                          legend_names = expt2_name) +
    limit_setting +
    paired_tsne_theme
  
  plt <- plot_grid(tsne_expt1,
                   tsne_expt2,
                   tsne_expts, 
                   nrow = 1, 
                   align = 'hv',
                   axis = 'l')
  
  save_plot(file.path(outdir, paste0("tsne_", out_name, "_expts.pdf")), 
            plt,
            ncol = 3, 
            base_height = 3,
            base_aspect_ratio = 1.1)
  
  save_plot(file.path(outdir, paste0("tsne_", out_name, "_expt1.pdf")), 
            tsne_expt1,
            base_height = 3,
            base_aspect_ratio = 1.1)
  
  save_plot(file.path(outdir, paste0("tsne_", out_name, "_expt2.pdf")), 
            tsne_expt2,
            base_height = 3,
            base_aspect_ratio = 1.1)
  
  save_plot(file.path(outdir, paste0("tsne_", out_name, "_expt1+2.pdf")), 
            tsne_expts,
            base_height = 3,
            base_aspect_ratio = 1.1)
  
}


plot_expt_tech_tsnes <- function(seurat_obj = atcells, 
                                 expt1_cells = expt1,
                                 expt2_cells = expt2,
                                 expt1_name = paste("Technical replicate", 
                                                    c("1", "2")),
                                 expt2_name = paste("Technical replicate", 
                                                    c("1", "2")),
                                 out_name = "all",
                                 # tsne_cols = brewer.pal(9, "Paired")[c(6, 5, 2, 1)],
                                 tsne_cols = brewer.pal(4, "Set1"),
                                 combined_cols = tsne_cols,
                                 outdir = file.path("fig2", 
                                                    "tsne_per_expt",
                                                    "technical_replicates"),
                                 save_plots = T) {
  
  paired_tsne_theme  <- theme(axis.title.x = element_blank(),
                              axis.text = element_text(size = 10),
                              axis.title = element_text(size = 10)) 
  
  limit_setting <- lims(x = c(min(seurat_obj@dr$tsne@cell.embeddings[, 1]), 
                              max(seurat_obj@dr$tsne@cell.embeddings[, 1])),
                        y = c(min(seurat_obj@dr$tsne@cell.embeddings[, 2]),
                              max(seurat_obj@dr$tsne@cell.embeddings[, 2])))
  
  tsne_expts <- plot_tsne(seurat_obj, 
                          ident = "experiment_and_tech_rep",
                          pt.alpha = 0.5,
                          .cols = combined_cols, 
                          cell_filter = c(expt1_cells, expt2_cells)) +
    paired_tsne_theme +
    guides(color = guide_legend(nrow = 4 , byrow=TRUE,
                                label.position = "left",
                                override.aes = list(size = 5)))
  
  if (save_plots){
    tsne_expts <- tsne_expts  + theme(legend.pos = "none")
  } else{
    tsne_expts <- tsne_expts + theme(legend.pos = "top")
  }
  
  tsne_expt1 <- plot_tsne(seurat_obj, 
                          ident = "technical_replicate",
                          pt.alpha = 0.5,
                          .cols = tsne_cols[1:2], 
                          cell_filter = expt1_cells,
                          legend_names = expt1_name) +
    limit_setting +
    paired_tsne_theme +
    theme(legend.pos = "none")
  
  
  tsne_expt2 <- plot_tsne(seurat_obj, 
                          ident = "technical_replicate",
                          pt.alpha = 0.5,
                          .cols = tsne_cols[3:4], 
                          cell_filter = expt2_cells,
                          legend_names = expt2_name) +
    limit_setting +
    paired_tsne_theme +
    theme(legend.pos = "none")
  
  plt <- plot_grid(tsne_expt1,
                   tsne_expt2,
                   tsne_expts, 
                   nrow = 1, 
                   align = 'hv',
                   axis = 'l')
  
  if(save_plots){
    save_plot(file.path(outdir, paste0("tsne_", out_name, "_expts.pdf")), 
              plt,
              ncol = 3, 
              base_height = 3,
              base_aspect_ratio = 1.1)
  } else{
    return(tsne_expts)
  }
  
}


#### 


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


tsne_by_sample_type <- function(seurat_obj, .gene, cell_sets, 
                                top_label = NULL, 
                                limit_gene_exp = NULL, 
                                resize_tsne = FALSE,
                                ...){
  
  if (is.null(limit_gene_exp )){
    limit_gene_exp <- c(0, max(FetchData(seurat_obj, 
                                         .gene, 
                                         cells.use = unique(unlist(cell_sets)))))
  }
  #if(limit_gene_exp[2] == 0) return()
  
  if(resize_tsne){
    all_cells <- unique(unlist(cell_sets)) 
    sub_obj <- SubsetData(seurat_obj, cells.use = all_cells)
  } else{
    sub_obj <- seurat_obj
  }
  
  limit_setting <- lims(x = c(min(sub_obj@dr$tsne@cell.embeddings[, 1]), 
                              max(sub_obj@dr$tsne@cell.embeddings[, 1])),
                        y = c(min(sub_obj@dr$tsne@cell.embeddings[, 2]),
                              max(sub_obj@dr$tsne@cell.embeddings[, 2])))
  
  plt_theme <- theme(axis.line = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     plot.background = element_rect(colour = "grey", 
                                                    linetype = "solid", size = 0.5))
  
  if(is.null(top_label)){ top_label <- rep("", length(cell_sets))}
  
  plts = list()
  for (i in seq_along(1:length(cell_sets))){
    plts[[i]] <- plot_feature(seurat_obj, feature = .gene, 
                              cell_filter = cell_sets[[i]],
                              pt_alpha = 1, max_y = limit_gene_exp, ...) + 
      limit_setting +
      labs(title = top_label[i]) +
      theme(legend.position = "none") +
      plt_theme
    
    if(i == length(cell_sets)){
      tmp <- plot_feature(seurat_obj, feature = .gene, 
                          cell_filter = cell_sets[[i]],
                          pt_alpha = 1, max_y = limit_gene_exp, ...) + 
        limit_setting +
        labs(title = top_label[i]) +
        plt_theme
      plts[[i + 1]] <- g_legend(tmp)
    }
  }
  
  plot_grid(plotlist = plts, nrow = 1,
            rel_widths = c(rep(1, length(cell_sets)), .33))
}


plot_grouped_features <- function(seurat_obj, 
                                  gene_to_plot = "Krt5",
                                  grp_1_cells = non_type2,
                                  grp_2_cells = type2,
                                  grp_1_title = "Naive Non-ATII",
                                  grp_2_title = "Naive ATII",
                                  all_title = "All Cells",
                                  plt_theme = paired_tsne_theme,
                                  out_dir = "fig2") {
  
  limit_gene_exp <- c(0, max(FetchData(seurat_obj, 
                                       gene_to_plot, 
                                       cells.use = c(grp_1_cells,
                                                     grp_2_cells))))
  
  limit_setting <- lims(x = c(min(seurat_obj@dr$tsne@cell.embeddings[, 1]), 
                              max(seurat_obj@dr$tsne@cell.embeddings[, 1])),
                        y = c(min(seurat_obj@dr$tsne@cell.embeddings[, 2]),
                              max(seurat_obj@dr$tsne@cell.embeddings[, 2])))
  
  tsne_all <- plot_feature(seurat_obj, 
                           gene = gene_to_plot,
                           cell_filter = NULL,
                           max_y = limit_gene_exp) + 
    limit_setting + 
    labs(title = all_title) +
    plt_theme
  
  tsne_grp1 <- plot_feature(seurat_obj, 
                            gene = gene_to_plot,
                            cell_filter = grp_1_cells,
                            max_y = limit_gene_exp) + 
    limit_setting + 
    labs(title = grp_1_title) +
    plt_theme
  
  tsne_grp2 <- plot_feature(seurat_obj, 
                            gene = gene_to_plot,
                            cell_filter = grp_2_cells,
                            max_y = limit_gene_exp) + 
    limit_setting + 
    labs(title = grp_2_title) +
    plt_theme
  
  plt <- plot_grid(tsne_grp1,
                   tsne_grp2, 
                   tsne_all,
                   nrow = 1, 
                   align = 'hv',
                   axis = 'l')
  
  grp_1_title <- str_replace_all(grp_1_title, "( |-)", "_")
  grp_2_title <- str_replace_all(grp_2_title, "( |-)", "_")
  
  save_plot(file.path(out_dir, paste0(gene_to_plot, "_tsne_all_plots.pdf")), 
            plt,
            ncol = 3, 
            base_height = 3,
            base_aspect_ratio = 1.25)
  
  save_plot(file.path(out_dir, paste0(gene_to_plot, "_tsne_all.pdf")), 
            tsne_all,
            base_height = 3,
            base_aspect_ratio = 1.25)
  
  save_plot(file.path(out_dir, paste0(gene_to_plot, "_tsne_", grp_1_title, ".pdf")), 
            tsne_grp1,
            base_height = 3,
            base_aspect_ratio = 1.25)
  
  save_plot(file.path(out_dir, paste0(gene_to_plot, "_tsne_", grp_2_title, ".pdf")), 
            tsne_grp2,
            base_height = 3,
            base_aspect_ratio = 1.25)
}


plot_violin_tsne <- function(seurat_obj, 
                             .ident = "labeled_clusters",
                             .gene = "Hopx",
                             cell_filter = NULL){
  
  gene_dat <- FetchData(seurat_obj, .gene) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell")
  
  if (!is.null(cell_filter)){
    gene_dat <- dplyr::filter(gene_dat,
                              cell %in% cell_filter)
  }
  
  mdata <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")
  
  plot_dat <- left_join(gene_dat, mdata, by = "cell")
  vln_plt <- plot_violin(plot_dat, 
                         .x = .ident,
                         .y = .gene,
                         .fill = .ident,
                         .scale = "width",
                         jitter = T,
                         .alpha = 0.5,
                         single_col = brewer.pal(11, "RdGy")[7]) +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.pos = "none")
  
  tsne_plt <- plot_feature(seurat_obj, 
                           gene = .gene,
                           cell_filter = cell_filter,
                           pt.alpha = 1) +
    guides(color = guide_colorbar(barwidth = 0.5, barheight = 5)) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  empty_plot <- ggplot(tsne_plt$data) +
    geom_blank()
  
  plt <- plot_grid(tsne_plt,
                   empty_plot,
                   ncol = 1,
                   rel_heights = c(4, 1))
  
  plt <- plot_grid(plt,
                   vln_plt,
                   nrow = 1, 
                   axis = 'h',
                   align = 'bl',
                   rel_widths = c(1.5, 1),
                   rel_heights = c(1.5, 1))
  plt
}


plot_violin_gene <- function(seurat_obj, 
                             .ident = "labeled_clusters",
                             .genes = "Hopx"){
  
  gene_dat <- FetchData(seurat_obj, .genes) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell")
  
  mdata <- seurat_obj@meta.data %>% tibble::rownames_to_column("cell")
  
  plot_dat <- left_join(mdata, gene_dat, by = "cell")
  
  plts <- map(.genes, 
              ~plot_violin(plot_dat, 
                           .x = .ident,
                           .y = .x,
                           .fill = .ident,
                           .scale = "width",
                           jitter = T,
                           .alpha = 0.5,
                           single_col = brewer.pal(11, "RdGy")[7]) +
                theme(axis.title.x = element_blank(),
                      axis.text = element_text(size = 10),
                      axis.title = element_text(size = 10),
                      legend.pos = "none")
  )
  plot_grid(plotlist = plts)
}


violin_tsne_cell_split <- function(gene, outdir,
                                   cells_to_plot = NULL){
  output <- file.path("fig3", "all_cells_violin_tsnes", 
                      outdir, str_c(gene, ".pdf"))
  
  plts <- map(cells_to_plot, 
              ~plot_violin_tsne(atcells, 
                                .gene = gene,
                                cell_filter = .x))
  
  plt <- plot_grid(plotlist = plts, 
                   nrow = 3,
                   rel_heights = 4,
                   rel_widths = 6,
                   labels = "AUTO")
  
  save_plot(output, plt,
            base_height = 4, 
            base_aspect_ratio = 1.5,
            nrow = 3)
}

