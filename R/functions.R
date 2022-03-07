# Load required libraries
suppressPackageStartupMessages({
  library(JADE)
  library(Hmisc)
  library(corrplot)
  library(pheatmap)
  library(reshape2)
  library(scales) #! create my own scales::rescale function?
  library(ggpubr) 
  library(ggplotify)
  library(patchwork)
})
################################################################################

# Customizable corrplot functions (modified from https://www.khstats.com/blog/corr-plots/corr-plots/)
cors <- function(df, cor.stat) {
  M <- Hmisc::rcorr(as.matrix(df), type = cor.stat)
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

formatted_cors <- function(df, cor.stat){
  cors(df, cor.stat) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}


################################################################################
#' Title
#'@description 
#'
#'@details
#'
#'
#'@param param
#'
#'@return
#'
#'@examples


Probes_to_whatever <- function(df, # dataset with ProbeIDs as rownames
                               annotation_table, # Bioconductor annotation table ie. AnnotationDbi::select(hgu219.db, probes, c("SYMBOL", "ENSEMBL", "GENENAME"))
                               new_IDs # Desired IDs ("SYMBOL", "ENSEMBL", "GENENAME")
)
  #TBI: choose function of aggregation, choose starting format, choose 
{
  final_df <- merge(df,annotation_table, by.x=0, by.y="PROBEID")
  
  # Stats of conversion
  non_agg <- nrow(final_df)
  non_agg_uniq <- length(unique(final_df[new_IDs]))
  non_agg_nas <- sum(is.na(final_df[new_IDs]))
  non_agg_nonas <- non_agg-non_agg_nas
  
  # Conversion
  final_df <- aggregate(final_df,
                        final_df[new_IDs],
                        FUN = mean) %>% 
    column_to_rownames(new_IDs) %>% 
    subset(select = -c(Row.names, SYMBOL, GENENAME, ENSEMBL))
  
  agg <- nrow(final_df)
  
  print(paste(100*(non_agg-agg)/non_agg,
              "% of the originally merged df have been agregated/removed"))
  
  print(paste(100*(non_agg_nas)/non_agg,
              "% of the originally merged df have been removed due no", new_IDs))
  
  print(paste(100*(non_agg_nonas - agg)/non_agg,
              "% of the non NAs df have been aggregated"))
  
  return(final_df)
}
################################################################################

################################################################################
#' Iterative ICA
#'@description
#' Run ICA over a specific range of components
#'@details
#' `Range_ICA` returns a list of lists with the A (sample weights)
#'  and S (gene weights) matrices of the specified range of components with 
#'  100 as max iterations.
#'
#'@param df expression dataframe with sample names in columns
#' and genes/probes in rows 
#' 
#'@param range.comp a range in the format n:n to specify 
#'the range of number of components to analyze
#'
#'@return
#'
#'@examples

Range_ICA <- function(df, range.comp) {
  A_mats <- list()
  S_mats <- list()
  samples <- colnames(df)
  
  for (n.comp in range.comp) {
    jade_result <- JADE(df,
                        n.comp = n.comp, maxiter = 100
    )
    l_name <- paste(c("nc", n.comp), collapse = "")
    ic_names <- paste("IC", 1:n.comp, sep = ".")
    
    A_mats[[l_name]] <- jade_result[["A"]]
    colnames(A_mats[[l_name]]) <- ic_names
    rownames(A_mats[[l_name]]) <- samples
    
    S_mats[[l_name]] <- jade_result[["S"]]
    colnames(S_mats[[l_name]]) <- ic_names
  }
  
  return(list("A" = A_mats, "S" = S_mats, "samples" = samples))
}
################################################################################

################################################################################
#'Find the best number of components
#'
#'@description
#' Analyze the ICA with different number of components to extract the best
#' related to a specific factors
#'
#'@details
#' `Best_nc` takes the `Range_ICA` output and analyze it throughout
#'the different range of components with respect to one or more than one variable.
#'
#'If the target variable is continuous spearman R is used to evaluate its association to the components
#'
#'Else, the ratio between the each component IQR and the mean of each component splitted by group IQRs
#' is used for this purpose. This is needed given that A mtrix weights dont usually fullfill normality and homocaedasdicity assumptions
#'
#'
#'
#'@param icas_list output of `Range_ICA`
#'@param range.comp a range in the format n:n to specify 
#'the range of number of components
#'@param metadata a tibble containing the variable/s to 
#'associate and the id_column. Not restricted to only these.
#'@param metadata_id a string indicating the name of the column 
#'containing the sample identifier used when running `Range_ICA`
#'@param vars a string indicating the variable for which you want
#' to maximize the association. For discrete variables more than one 
#' variable can be provided in form of a vector
#'@param is.categorical boolean indicating wether the var provided is categorical.
#'
#'@return
#'
#'@examples
#'
#'TBD 
#'- indpendence of range.ncomp

Best_nc <- function(icas_list,
                          range.comp,
                          metadata,
                          metadata_id = 0,
                          vars,
                          is.categorical = FALSE) {

  to_plot <- tibble(nc = character(),
                    var = character(),
                    IC = character(),
                    R = numeric(),
                    Pval = numeric())
  
  correlations <- list()
  for (cv in vars){
    correlations[[cv]] <- list()
    correlations[[cv]][["Pvals"]] <- list() #! store only Rs?
    correlations[[cv]][["Rs"]] <- list()
  }

  
  for (n.comp in range.comp) {
    nc_str <- paste("nc", n.comp, sep ="")
    ica.nc <- icas_list$A[[nc_str]] %>% as_tibble(rownames = "samples")
    test_vars <- dplyr::rename(metadata, samples = all_of(metadata_id)) %>% inner_join(ica.nc, by = "samples")
    
    if (is.categorical == FALSE){
      continuous_var <- c(names(ica.nc), vars)
      test_cont <- dplyr::select(test_vars, all_of(continuous_var)) %>%
        drop_na() %>%
        column_to_rownames("samples")
      
      res2 <- rcorr(x=as.matrix(test_cont), type = "spearman")
        
      res2$r <- res2$r[vars, names(ica.nc)[-1],drop=F]
      res2$P <- res2$P[vars, names(ica.nc)[-1],drop=F]
      
      for (cv in vars) { 
        correlations[[cv]]$Rs[[nc_str]] <- res2$r[cv,]
        correlations[[cv]]$Pvals[[nc_str]] <- res2$P[cv,]
      }
    } else if (is.categorical == TRUE){
      discrete_var <- c(names(ica.nc), vars)
      test_disc <- dplyr::select(test_vars, all_of(discrete_var)) %>%
        drop_na()
    
      disp0 <- test_disc %>% dplyr::select(names(ica.nc)) %>% summarise(across(where(is.double), IQR))
      dispgs <-  test_disc %>% dplyr::group_by(get(vars)) %>%
        summarise(across(where(is.double), IQR)) %>% #get group variances
        summarise(across(where(is.double), mean)) # compute mean var groups
      
      disp_ratio <- disp0/dispgs
      
      for (dv in vars) { 
        correlations[[dv]]$Rs[[nc_str]] <- as.numeric(disp_ratio) %>% set_names(.,names(disp_ratio))
        correlations[[dv]]$Pvals[[nc_str]] <- NULL
      }
    }

    # Check the best IC and save it as a representative of the nc
    for (cv in vars){
      best_c <- sort(abs(correlations[[cv]]$Rs[[nc_str]]),decreasing = T)[1] %>% names()
      to_plot %>% add_row(nc = nc_str,
                          var = cv,
                          IC = best_c,
                          R = abs(correlations[[cv]]$Rs[[nc_str]][[best_c]]),
                          Pval = correlations[[cv]]$Pvals[[nc_str]][[best_c]]) -> to_plot
    }
    }
  
  order_bars <- group_by(to_plot, nc) %>%
    summarise(meanR = mean(R)) %>%
    arrange(desc(meanR)) %>%
    dplyr::select(nc)
  
  to_plot$nc <- to_plot$nc %>%
    factor(order_bars$nc)
  
  to_plot %>% ggplot(aes(x=nc, y=R, fill=var)) +
    geom_bar(stat='identity', position=position_dodge()) +
    geom_text(aes(x = nc, y = 0.1, label = IC), angle = 90, position = position_dodge(width = 0.9)) +
    theme_minimal() +
    ggtitle(paste0("ICA Space correlation with variable/s: ", paste(vars, sep = ", ")))
}
################################################################################

################################################################################
#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples

plot_sample_weights <- function(A_mat, annotations, plotfile = NA){
  stopifnot(rownames(A_mat) == rownames(annotations))
  if (!is.na(plotfile)){pdf(file=paste(plotfile, ".pdf", sep=""))}
  
  # Sample weights
  for (ann in colnames(annotations)){
    
    rug_aes <- annotations[[ann]]
    rug_name <- ann
    comps_plots <- lapply(colnames(A_mat), function(ic){
      p <- 
        ggplot(data.frame(A_mat)) +
        aes_string(ic) +
        geom_density() + 
        geom_rug(aes(color = rug_aes), length = unit(0.1, "npc")) +
        labs(color = ann) +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())
      
      if(is.integer(rug_aes)) {
        p <- p  +
          scale_color_discrete()
        
      } else if(is.double(rug_aes)){
        p <- p  +
          scico::scale_color_scico(palette = "berlin")
      } 
      
      p
      
    })
    ggarrange(plotlist = comps_plots, common.legend = T,  legend = "bottom") %>% annotate_figure(
      top = text_grob(ann, color = "black", face = "bold", size = 14), ) %>% print()
  }
  if (!is.na(plotfile)) {dev.off()}
}
################################################################################

################################################################################
#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples

ICA_explorator <- function(ica,
                             df,
                             df_cont = NA,
                             df_disc = NA,
                             df_id = 0, # 0 for rownmaes, otherwise put the colname
                             plotfile = NA, # Define the path and prefix of the file to be created if wnated
                             interest_IC = FALSE) {
  
  A <- as_tibble(ica[["A"]], rownames = df_id)
  corr_values <- A %>% inner_join(df, by = df_id)
  returning_list <- list()
  # Continuous
  if (!any(is.na(df_cont))){
    
  corr_cont <- corr_values %>% dplyr::select(names(A), all_of(df_cont))
  
  returning_list[["cont"]] <- corr_cont %>%
    column_to_rownames(df_id) %>%
    formatted_cors("spearman") %>% #! adjust for conditional scaling
    filter(measure1 %in% names(A), measure2 %in% names(df)) %>% #not square corr
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation", title= "ICA Sample weights correlations", 
         subtitle="Only significant correlation coefficients shown (95% I.C.)") +
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) + 
    ggpubr::rotate_x_text(angle = 90)
  
  if (!is.na(plotfile)){
    ggsave(coerrelogram, paste0(plotfile, "cont.pdf"), dpi = "print", width = 210, height = 150, units = "mm", scale = 1.25) #! adjust scale and implement auto scale
    } else {
    #print(coerrelogram)
    }
  }


  # Categorical
  if (!any(is.na(df_disc))){
    
    corr_disc <- corr_values %>%
      dplyr::select(names(A), all_of(df_disc)) %>%
      mutate_at(names(A)[-1], ~(rescale(., to = c(-1, 1)) %>% as.vector)) # standarization is needed for visualization
    
    plot_list <- list()
    legend_list <- list()
    for (var in df_disc){
      to_plot <- corr_disc %>%
        dplyr::select(names(A), all_of(var)) %>%
        pivot_longer(names(A)[-1], names_to = 'component', values_to = 'weight') %>%
        mutate_at(vars(component), ~(factor(., levels=names(A)[-1])))
      
      plot_list[[var]] <- to_plot %>% ggplot() + aes_string(color = var, x = var, y = "weight") +
        geom_point() +
        facet_grid(. ~ component)  +
        labs(x = NULL, y = NULL, color = var) +
        rremove("axis.text") + 
        rremove("ticks")
    }
    returning_list[["disc"]] <- wrap_plots(plot_list, ncol = 1, guides = "collect")
  }
  return(returning_list)
}
################################################################################

################################################################################

#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples
Sampleweights_indepth <- function(ica,
                                  interest_IC = FALSE,
                                  df,
                                  df_id = 0,
                                  var,
                                  test = NULL) {
  A <- as_tibble(ica[["A"]], rownames = df_id)

  to_test <- inner_join(A, df, by = df_id)
  base_gg <- ggplot(data = to_test) +
    aes_string(x = var, y = interest_IC)


  # check assumptions
  ## variable type
  if (summarise_all(to_test, class) %>% dplyr::select(all_of(var)) == "numeric") {
    vartype <- "continuous"
  } else {
    vartype <- "discrete"
  }

  ## normality
  n.pv <- shapiro.test(to_test[, interest_IC, drop = T])$p.value
  if (n.pv < 0.05) {
    norm <- FALSE
  } else {
    norm <- TRUE
  }

  ## equal variance (fligner test as its better with non normality)
  v.pv <- fligner.test(get(interest_IC) ~ get(var), data = to_test)$p.value
  if (v.pv < 0.05) {
    eqvar <- FALSE
  } else {
    eqvar <- TRUE
  }


  cat("Assumption check:\n")
  cat(paste0("- var interpreted as: ", vartype, "\n"))
  cat(paste0("- Normality: ", norm, "\n"))
  cat(paste0("- Homoscedasticity: ", eqvar, "\n"))

  if (vartype == "continuous" & !is.null(test)) {
    to_test_df <- to_test %>% dplyr::select(df_id, all_of(interest_IC), all_of(var)) %>% column_to_rownames(df_id)
    corr <- rcorr(x = data.matrix(to_test_df), type = test)
    base_gg <- base_gg + ggtitle(paste0(interest_IC, " vs. ", var),
                                 subtitle = paste0(test, " R: ", signif(corr$r[var,interest_IC], 2),", p-value = ", signif(corr$P[var,interest_IC], 2)))
    
  } else if (vartype == "discrete" & !is.null(test)) {
    if (test == "ttest") {
      stats <- t.test(get(interest_IC) ~ get(var), data = to_test, var.equal = eqvar)
    } else if (test == "wilcox") {
      stats <- wilcox.test(get(interest_IC) ~ get(var), data = to_test)
    } else if (test == "anova") {
      stats <- aov(get(interest_IC) ~ get(var), data = to_test)
      stats$method <- "Analysis of variance"
    } else if (test == "welch") {
      stats <- oneway.test(get(interest_IC) ~ get(var), data = to_test)
    } else if (test == "kruskal") {
      stats <- kruskal.test(get(interest_IC) ~ get(var), data = to_test)
    }

    base_gg <- base_gg + ggtitle(paste0(interest_IC, " vs. ", var), subtitle = paste0(stats$method,": p-value = ", signif(stats$p.value, 2)))
  }
  return(base_gg)
}
################################################################################

################################################################################
#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples

Corr_by_ct <- function(A, decon){
  A <- as.data.frame(A)
  dir.create("02_Output/deconv_by_ct/")
  
  # segregate the deconvolution (by cell type)
  colnames(decon) %>% str_split(., "_", n = Inf, simplify = TRUE) -> col_interpreter
  col_interpreter[,1] %>% unique(.) -> unique_methods
  decon_by_methods <- list()
  ci_by_methods <- list()
  for (dm in unique_methods){
    decon %>% .[,which(col_interpreter[,1] == dm)] -> decon_by_methods[[dm]]
    col_interpreter %>% .[which(.[,1] == dm),] -> ci_by_methods[[dm]]
    
  }
  
  
  decon_by_ct <- list()
  cell_types <- c("General_metrics", "B", "T", "NK", "Ms", "Neu", "Fibroblast", "Others")
  for (ct in cell_types){
    # Creation of DF
    decon_by_ct[[ct]] <- data.frame(row.names = rownames(decon))
    
    for (dm in unique_methods){
      
      decon_m <- decon_by_methods[[dm]]
      ci_m <- ci_by_methods[[dm]]
      decon_by_ct.temp <- ""
      
      # get the specific value
      if (ct == "General_metrics" & dm =="XCELL"){
        decon_m %>% 
          .[,which(ci_m[,2] %in% c("ImmuneScore",
                                   "StromaScore",
                                   "MicroenvironmentScore")), drop=FALSE] -> decon_by_ct.temp
      } # OK
      
      if (ct == "B"){
        if (dm %in% c("Quantiseq", "MCP")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("B", "B-cells")), drop=FALSE] -> decon_by_ct.temp
        }
        
        if (dm == "XCELL"){
          decon_m %>% 
            .[,which((paste(ci_m[,3], ci_m[,4], sep = "") %in%
                        c("B-cells", "memoryB-cells")) | (ci_m[,2] %in% c("B", "B-cells"))), drop=FALSE] -> decon_by_ct.temp
          
        }
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>% 
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       (ci_m[,3] == "B")), drop=FALSE] -> decon_by_ct.temp
          print(colnames(decon_by_ct.temp))
        }
        
      } # OK
      
      if (ct == "T"){
        if (dm %in% c("Quantiseq", "MCP")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("T", "CD8", "Cytotoxic")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
        if (dm %in% c("XCELL")){
          
          decon_m %>% 
            .[,grepl( "T", paste(ci_m[,2], ci_m[,3], ci_m[,4], sep = "")), drop=FALSE] -> decon_by_ct.temp
          
          
        }
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>% 
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       (ci_m[,3] %in% c("CD4", "CD8"))), drop=FALSE] -> decon_by_ct.temp
          
        }
        
      } # OK
      
      if (ct == "NK"){
        if (dm %in% c("Quantiseq", "MCP", "XCELL")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("NK", "NKT")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>% 
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       (ci_m[,3] == "NK")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
      } # OK
      
      if (ct == "Ms"){
        if (dm %in% c("Quantiseq", "MCP")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("Macrophage", "Monocyte",
                                     "Monocytic", "Myeloid")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
        if (dm %in% c("XCELL")){
          
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("aDC", "cDC", "DC", "iDC", "pDC" ,"Macrophages", "Monocytes",
                                     "Myeloid")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>% 
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       (ci_m[,3] %in% c("Mono", "M0", "M1", "M2"))), drop=FALSE] -> decon_by_ct.temp
          
        }
        
      } # OK
      
      if (ct == "Neu"){
        if (dm %in% c("Quantiseq", "MCP", "XCELL")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("Neutrophil", "Neutrophils")), drop=FALSE] -> decon_by_ct.temp
        }
        
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>% 
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       (ci_m[,3] == "Neu")), drop=FALSE] -> decon_by_ct.temp
          
        }
        
      } # OK
      
      if (ct == "Fibroblast"){
        if (dm %in% c("MCP", "XCELL")){
          decon_m %>% 
            .[,which(ci_m[,2] %in% c("Fibroblasts")), drop=FALSE] -> decon_by_ct.temp
        }
        
      } #OK
      
      if (ct == "Others"){
        used_cols <- unlist(lapply(X = decon_by_ct,
                                   FUN = colnames),use.names = FALSE)
        
        if (dm %in% c("Epidish", "DeconRNASeq")){
          decon_m %>%
            .[,which((ci_m[,2] == "BPRNACan3DProMet") &
                       !(colnames(.) %in% used_cols)), drop=FALSE] -> decon_by_ct.temp
        }
        else {
          decon_m %>%
            .[,which(!(colnames(.) %in% used_cols)), drop=FALSE] -> decon_by_ct.temp
        }
        
      }
      if (typeof(decon_by_ct.temp) == "double"){
        decon_by_ct[[ct]][,colnames(decon_by_ct.temp)] <- decon_by_ct.temp
      }
    }
    
    corr_values <- merge(A, decon_by_ct[[ct]], by.x=0, by.y=0)
    row.names(corr_values) <- corr_values[,"Row.names"]
    corr_values <- subset(corr_values, select = -c(Row.names))
    
    res2 <- rcorr(x=as.matrix(corr_values), type = "spearman")
    
    res2$r <- res2$r[colnames(A), colnames(decon_by_ct[[ct]])]
    res2$P <- res2$P[colnames(A), colnames(decon_by_ct[[ct]])]
    title <- paste(ct, "(Spearman correlation p < 0.05)")
    
    # pdf(paste(c("02_Output/deconv_by_ct/",
    #             elected_ncomp, "_", ct, ".pdf"), collapse = ""))
    
    corrplot(res2$r, type = "full", order="original", method = "circle",
             p.mat = res2$P, sig.level = 0.05, insig = "blank", is.corr = T,
             main = title , mar=c(0,0,1,0)) # http://stackoverflow.com/a/14754408/54964
    # dev.off()
  }
  return(decon_by_ct)
}
################################################################################

################################################################################
#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples

Explore_GW <- function(S_sym, S_ensembl, elected_ncomp, interest_IC){
  
  # new
  scran <- read_rds("../CellTypes_Markers_All_cells.rds")
  
  #old
  # scran <- read_rds("../OLD_marker_list_By_levels.rds")
  # scran <- scran$all_pops
  
  
  gedepirRes <- gedepir::enrich(S_sym,scran, ICAbased = T)
  
  ct_enrich_tmp <- tibble(.rows = length(scran), )
  ct_enrich_tmp$celltype <- names(scran)
  ct_enrich_tmp %>% column_to_rownames("celltype") -> ct_enrich_tmp
  
  do.call(cbind, lapply(1:elected_ncomp, function(n) {
    print(n)
    ic = paste("IC.", n)
    cts <- gedepirRes[[n]][["pathway"]]
    es <- gedepirRes[[n]][["ES"]]
    if(is_null(es)){
      ct_enrich_tmp[, ic] = NA
    }else{
      ct_enrich_tmp[cts, ic] = es
    }
    print(c(cts, es))
    ct_enrich_tmp
  })) -> ct_enrich
  colmat <- colorRampPalette(c("white", "black"))
  data.matrix(ct_enrich) %>% t() %>% corrplot(is.corr = F, method = "color", col = colmat(200))
  
  ## GSVA on the Gene weights
  ### gene set preparation
  all_genesets <- msigdbr("Homo sapiens")
  all_genesets %>% filter(gs_cat %in% c("H", "C2", "C7")) -> use_genesets
  msigdbr_list = split(x = use_genesets$human_ensembl_gene, f = use_genesets$gs_name)
  msigdb_descriptions <- use_genesets[c("gs_name", "gs_description")] %>%
    unique() %>% column_to_rownames("gs_name")

  gsvaRes <- gsva(data.matrix(S_ensembl), msigdbr_list, min.sz = 15)
  gsvaRes[order(gsvaRes[,interest_IC]),]

  gsvaRes %>% as.data.frame() %>%  mutate(the_rank = rank(-get(interest_IC), ties.method = "random")) %>%
    dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% arrange(the_rank) %>%
    dplyr::select(!c(the_rank)) -> gsvaTop


  gsvaTop %>% pheatmap(cluster_rows = F) # to have a pheatmap
  gsvaTop  %>% rownames() %>% msigdb_descriptions[.,] %>%
    cbind(gsvaTop[interest_IC], .) %>% rownames_to_column("msigdb") %>% write_tsv("02_Output/gsvaRes_Puleo.tsv") #to write
  
  gsvaTop[,interest_IC, drop=F] %>% rownames_to_column("pathway") -> gsvaTop_intICtoplot
  gsvaTop_intICtoplot$sign <- "pos"
  gsvaTop_intICtoplot[gsvaTop_intICtoplot[,interest_IC] < 0, "sign"] <- "neg" 
  
  
  ggplot(gsvaTop_intICtoplot[gsvaTop_intICtoplot$sign == "pos",], aes(reorder(str_sub(pathway, end = 50), get(interest_IC)),get(interest_IC), color = sign))  +
    geom_point() + coord_flip() + scale_color_discrete(type = c("blue")) + theme(legend.position = "none", plot.title = element_text(lineheight=.8, face="bold")) +
    labs(title = paste(interest_IC, "25 top gene sets"),y = "ES", x = NULL)
  
  ggplot(gsvaTop_intICtoplot[gsvaTop_intICtoplot$sign == "neg",], aes(reorder(str_sub(pathway, end = 50), get(interest_IC)),get(interest_IC), color = sign))  +
    geom_point() + coord_flip() + scale_color_discrete(type = c("red")) + theme(legend.position = "none", plot.title = element_text(lineheight=.8, face="bold")) +
    labs(title = paste(interest_IC, "25 bottom gene sets"),y = "ES", x = NULL)
  
  return(gsvaTop)
}
################################################################################

################################################################################
#'Title
#'
#'@description
#'
#'@details
#'
#'@param param
#'
#'@return
#'
#'@examples

PlotGeneWeights <- function(ica, expression, df_id, n_genes, column_annotation = NA, interest_IC){
  return_plots <- list()
  
  S <- as_tibble(ica[["S"]], rownames = df_id)
  S_sort <- arrange(S, get(interest_IC))
  S_mneg <- head(S_sort, n_genes)
  S_mpos <- tail(S_sort, n_genes)
  
  return_plots[["densplot"]] <- ggplot(S) +
    aes_string(y = interest_IC) +
    geom_density(alpha=.5, fill="#AED3FA") +
      geom_hline(yintercept = 0, colour = "black") +
      geom_hline(yintercept = max(S_mneg[interest_IC]), colour = "#7DB0DD", linetype = "dashed") + 
      geom_hline(yintercept = min(S_mpos[interest_IC]), colour = "#EAAFBB", linetype = "dashed") + 
      theme_classic()

    tibble(name = S_mneg[[df_id]], class = "most negative") %>% add_row(name = S_mpos[[df_id]], class = "most possitive") %>% column_to_rownames("name") -> annot_row 
    
    return_plots[["heatmap"]] <- dplyr::filter(expression, get(df_id) %in% rownames(annot_row)) %>%
      dplyr::arrange(match(get(df_id), rownames(annot_row))) %>%
      column_to_rownames(df_id) %>% 
      pheatmap(scale = "row", annotation_row = annot_row, cluster_rows = F, 
               annotation_col = column_annotation) %>% as.ggplot()
    dev.off()
    
    return(return_plots)
  }
################################################################################