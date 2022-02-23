# Load required libraries
suppressPackageStartupMessages({
  library(JADE)
  library(Hmisc)
  library(corrplot)
  library(pheatmap)
  library(reshape2)
})
################################################################################


# documentation to be done. Newer version in  CAFs Alexia analysis

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



#' JADE ICA of a specific range of components
#'@description
#' `Range_ICA` returns a list of lists with the A (sample weights)
#'  and S (gene weights) matrices of the specified range of components with 
#'  100 as max iterations.
#'
#'@param df expression dataframe with sample names in columns
#' and genes/probes in rows 
#' 
#'@param range.comp a range in the format n:n to specify 
#'the range of number of components to analyze

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




#' Analyze the ICA with different number of components to extract the best
#' related to a specific factors
#'@description
#' `Best_nc` takes the `Range_ICA` output and analyze it throughout
#'the different range of components.
#'
#'prints a list of spearman correlations of the desired metadata factor
#'with the IC that best correlate with it, from the ICAs with the different number of components 
#'
#'@param icas_list output of `Range_ICA`
#'@param range.comp a range in the format n:n to specify 
#'the range of number of components
#'
#'TBD 
#'- indpendence of range.ncomp
#'- input of discrete variables (ttest, anova)
#'- .tsv output to better choose
#'

Best_nc <- function(icas_list,
                          range.comp,
                          metadata,
                          metadata_id = 0,
                          cont_vars,
                          disc_vars = NULL) {
  
  correlations <- list()
  for (cv in cont_vars){
    correlations[[cv]] <- list()
    correlations[[cv]][["Pvals"]] <- list()
    correlations[[cv]][["Rs"]] <- list()
    
  }

  
  for (n.comp in range.comp) {
    nc_str <- paste("nc", n.comp, sep ="")
    ica.nc <- icas_list$A[[nc_str]]
    
    stopifnot(rownames(ica.nc) == icas_list$samples)
    test_vars <- merge(ica.nc, metadata, by.x=0,by.y=metadata_id) %>% column_to_rownames("Row.names")
    continuous_var <- c(colnames(ica.nc), cont_vars)
    
    if (length(cont_vars) > 0){
      
      test_cont <- test_vars[,continuous_var]
        
      test_cont <- drop_na(test_cont)
      res2 <- rcorr(x=as.matrix(test_cont), type = "spearman")
        
      res2$r <- res2$r[cont_vars, colnames(ica.nc)]
      res2$P <- res2$P[cont_vars, colnames(ica.nc)]
      
      for (cv in cont_vars) { 
        correlations[[cv]]$Rs[[nc_str]] <- res2$r[cv,]
        correlations[[cv]]$Pvals[[nc_str]] <- res2$P[cv,]
      }
    }
  }
  ## Check the best
  to_plot <- tibble(nc = character(),
                    var = character(),
                    IC = character(),
                    R = numeric(),
                    Pval = numeric())
  for (n.comp in range.comp) {
    nc_str <- paste("nc", n.comp, sep ="")
    for (cv in cont_vars){
      best_c <- sort(abs(correlations[[cv]]$Rs[[nc_str]]),decreasing = T)[1] %>% names()
      to_plot %>% add_row(nc = nc_str,
             var = cv,
             IC = best_c,
             R = abs(correlations[[cv]]$Rs[[nc_str]][[best_c]]),
             Pval = correlations[[cv]]$Pvals[[nc_str]][[best_c]]) -> to_plot
    }
  }
  
  to_plot %>% group_by(nc) %>% summarise(meanR = mean(R)) %>% arrange(desc(meanR)) %>% dplyr::select(nc) -> order_bars
  to_plot$nc %>% factor(order_bars$nc) -> to_plot$nc
  to_plot %>%  ggplot(aes(x=nc, y=R, fill=var)) +
    geom_bar(stat='identity', position=position_dodge()) +
    geom_text(aes(x = nc, y = 0.1, label = IC), angle = 90, position = position_dodge(width = 0.9)) +
    theme_minimal() +
    ggtitle("ICA Space correlation with continuous variables")
}




##############################################################################
#block of functions to analyze ICA analysis!
################################################################################

plot_sample_weights <- function(A_mat, annotations, analysis_name){
  stopifnot(rownames(A_mat) == rownames(annotations))
  pdf(file=paste("02_Output/", analysis_name, ".pdf", sep=""))
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
  dev.off()
}


flattenCorrMatrix <- function(cormat, pmat) { #V
  # ut <- upper.tri(cormat)
  ut <- matrix(TRUE, nrow(cormat), ncol(cormat)) # mod to work with trimmed matrix
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    # column = rownames(cormat)[col(cormat)[ut]],
    column = colnames(cormat)[col(cormat)[ut]],  # mod to work with trimmed matrix
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

ICA_explorator <- function(ica,
                             df,
                             df_cont,
                             df_disc,
                             A_id = 0,
                             df_id = 0, # 0 for rownmaes, otherwise put the colname
                             plot_dir, analysis_name, # if this is the case no subdir will be created
                             df_cont_special = FALSE,
                             df_disc_special = FALSE,
                             interest_IC = FALSE) {
  
  S <- as.data.frame(ica[["S"]])
  A <- as.data.frame(ica[["A"]])
  #dir.create(plot_dir)
  # General correlation
  corr_values <- merge(A, df, by.x=A_id, by.y=df_id)
  row.names(corr_values) <- corr_values[,"Row.names"]
  corr_values <- subset(corr_values, select = -c(Row.names))
  continuous_var <- c(colnames(A), df_cont)
  discrete_var <- c(colnames(A), df_disc)
  
  ## Continuous
  corr_cont <- corr_values[,continuous_var]
  #corr_cont <- na.exclude(corr_cont) # no NAs allowed in rcorr
  res2 <- rcorr(x=as.matrix(corr_cont), type = "spearman")

  res2$r <- res2$r[colnames(A), df_cont]
  res2$P <- res2$P[colnames(A), df_cont]
  title <- "Spearman correlation p < 0.05"
  
  # if( length(df_cont) > 30){
  #   # melt and order by greatest absolute correlations
  #   flat_res <- flattenCorrMatrix(res2$r, res2$P)
  #   flat_res_o <- flat_res[order(abs(flat_res$cor), decreasing = TRUE),]
  #   if (interest_IC != FALSE){
  #     flat_res_o <- flat_res_o[flat_res_o$column %in% interest_IC,]
  #   }
  #   
  #   # get and filter by greatest correlations
  #   cell_type_order <- unique(flat_res_o$row)
  #   
  #   res2$r <- res2$r[cell_type_order[1:30],]
  #   res2$P <- res2$P[cell_type_order[1:30],]
  # }
  
  p_dir <- paste(plot_dir, analysis_name, sep = "/")
  dir.create(p_dir)
  # pdf(paste(c(p_dir, "corrplot.pdf"), collapse = "/")) # problem with sizes, for decon is super big
  corrplot(res2$r, type = "full", order="original", method = "circle",
           p.mat = res2$P, sig.level = 0.05, insig = "blank", is.corr = T)#,
           #main = title, mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964
  # dev.off()
  
  ## Discrete
  corr_disc <- corr_values[,discrete_var]
  
  if (interest_IC == FALSE){
    corr_disc.m <- melt(corr_disc, measure.vars = paste("IC", 1:elected_ncomp, sep = "."))
  } else {
    corr_disc.m <- melt(corr_disc, measure.vars = interest_IC)
  }
  ### boxplots (+ stat significance) # STATS ONLY WORK WHEN interest_IC == TRUE!!!!!!
  base_gg <- ggplot(data = corr_disc.m, aes(x=variable, y=value))
  
  for (i in df_disc){
    groups <- df[,i,drop = T] %>% unique()
    n_groups <- groups %>% length()
    
    # assumptions for t-test/1ANOVA
    ## normality
    for (g in groups) {
      n.pv = with(corr_disc.m, shapiro.test(value[get(i) == g]))$p.value
      if (n.pv < 0.05){norm = FALSE; break} else {norm = TRUE}
    }

    ## equal variance (fligner test as its better with non normality)
    v.pv = fligner.test(value ~ get(i), data = corr_disc.m)$p.value
    if (v.pv < 0.05){eqvar = FALSE} else {eqvar = TRUE}

    print(c(i,g,norm,eqvar))
    # test choice
    ## 2 groups
    if (n_groups < 3){
      if (all.equal(norm, eqvar, TRUE) == TRUE){
        stats <- t.test(value ~ get(i), data = corr_disc.m, var.equal = TRUE)
      } else if(eqvar == TRUE) {
        stats <- wilcox.test(value ~ get(i), data = corr_disc.m)
      } else {
        stats <- kruskal.test(value ~ get(i), data = corr_disc.m) # no need for equal variance nor normality
      }
    }

    ## 3 or more groups
    else {
      if (all.equal(norm, eqvar, TRUE) == TRUE){
        stats <- aov(value ~ get(i), data = corr_disc.m)

      } else if(norm == TRUE) {
        stats <- oneway.test(value ~ get(i), data = corr_disc.m)

      } else {
        stats <- kruskal.test(value ~ get(i), data = corr_disc.m)
      }
    }

    
    (base_gg + geom_violin(aes(fill = get(i)), position=position_dodge(0.8), width=0.5) + 
        geom_boxplot(aes(color = get(i)), position=position_dodge(0.8), width=0.1) +
        labs(fill = "values") +
        theme_classic() + 
        #ggtitle(toupper(i)) +
        ggtitle(toupper(i), subtitle = paste(stats$method,signif(stats$p.value, 2))) +
        ylab("weight") +
        xlab("component")) %>% print()
  }
  
  
  # Detailed samplemaps
  # for (ann in colnames(df)[-1]) {
  #   print(ann)
  #   rug_aes <- df[,ann]
  #   rug_name <- ann
  #   lapply(colnames(A), function(ic) {
  #     print(ic)
  #     p <-
  #       ggplot(A) +
  #       aes_string(ic) +
  #       labs(title = paste(ic, ann, sep = "_")) +
  #       theme_classic() +
  #       theme(
  #         axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank()
  #       )
  # 
  #     if (ann %in% df_cont) {
  #       p <- p + geom_density() +
  #         geom_rug(aes(color = rug_aes)) +
  #         scale_colour_gradientn(colours = terrain.colors(10), na.value = NA) +
  #         guides(color = guide_colourbar())
  #     }
  # 
  #     if (ann %in% df_disc) {
  #       p <- p + geom_density() +
  #         geom_rug(aes(color = rug_aes)) +
  #         theme(axis.title.y = element_blank())
  #     }
  #     p_dir <- paste(plot_dir, analysis_name, ic , "", sep = "/")
  #     dir.create(p_dir)
  #     print(p_dir)
  #     ggsave(
  #       filename = paste(p_dir, ann, ".pdf", sep = ""),
  #       plot = p, width = 10, height = 5
  #     )
  #   })
  # }
}






# Cell type segregated Deconvolution
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


# Gene weights

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

##  enrichmentfun: enrichmentfun function to apply enrichiment analysis on data refering to a set of markers
##  authors : Yasmina K, Yuna B and Magali R
##  arguments : x: dataframe with nr: nrow (genes) and nc : ncolumn (componenent, cells, patients...), markerscelltypes: a list of cell markers
##  example of use : MyEnrichExp = enrichmentfun(mtx, SetOfMarker)


enrichmentfun = function(x, markerscelltypes){
    fgseaResall <- apply(x, 2, function(x,mcell=markerscelltypes){
        genes=sort(x,decreasing=T)
        fgseaRes = fgseaSimple(mcell, genes, minSize=1, maxSize=2000, nperm=1000, scoreType = "pos")
        fgseaResDF = data.frame(fgseaRes)
        fgseaResDF. = fgseaResDF[,c("pval", "padj", "ES")]
        rownames(fgseaResDF.) = fgseaResDF$pathway
        fgseaResDF.$pvalpositiveES = ifelse(fgseaResDF.[,"ES"]>0, fgseaResDF.$pval, 1)
        fgseaResDF.$padjpositiveES = ifelse(fgseaResDF.[,"ES"]>0, fgseaResDF.$padj, 1)
        fgseaResDF.list <- setNames(split(fgseaResDF., seq(nrow(fgseaResDF.))), rownames(fgseaResDF.))
        fgseaResDFtab = do.call(cbind, fgseaResDF.list)
        return(fgseaResDFtab)
    })
}
