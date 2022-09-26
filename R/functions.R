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

################################################################################
#' Iterative ICA
#'@description
#' `Range_ICA` runs ICA over the specified number of components in order to
#'  later explore them.
#' 
#'@usage 
#'\code{Range_ICA(df, range.comp = 2:20)}
#'
#'@param expression tibble with sample names in columns and genes/probes in
#' rows. Should be mean centered and fairly filtered.  
#' 
#'@param df_id name of the column containing the gene ids
#' 
#'@param range.comp a range in the format n:n to specify the range of number of 
#'components to analyze. Not recomended to use more than 20 components. 
#'Specially with small datasets.
#'
#'@details
#' `Range_ICA` uses JADE algorithm with a maximum of 100 iterations to 
#' converge as ICA implementation.
#'@return
#' `Range_ICA` returns a list of lists with the A (sample weights)
#'  and S (gene weights) matrices of the specified range.comp
#'
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'@TBD
#'do last conditional more elegant and change name to Run_ICA

Range_ICA <- function(expression, df_id, range.comp = 2:20) {
  
  df <- column_to_rownames(expression, df_id)
  
  A_mats <- list()
  S_mats <- list()
  samples <- colnames(df)
  
  for (n.comp in range.comp) {
    jade_result <- tryCatch(
      {
        JADE(df, n.comp = n.comp, maxiter = 100)
      },
      error = function(cond) {
        message(paste("JADE with ", n.comp, "did not converge, skipping iter"))
        return(NA)
      }
    )
    if (jade_result[1] %>% is.na()) {
      next
    }

    l_name <- paste(c("nc", n.comp), collapse = "")
    ic_names <- paste("IC", 1:n.comp, sep = ".")
    
    A_mats[[l_name]] <- jade_result[["A"]]
    colnames(A_mats[[l_name]]) <- ic_names
    rownames(A_mats[[l_name]]) <- samples
    
    S_mats[[l_name]] <- jade_result[["S"]]
    colnames(S_mats[[l_name]]) <- ic_names
  }
  
  if (length(range.comp) == 1){
    A_mats = A_mats[[1]]
    S_mats = S_mats[[1]]
  }
  
  return(list("A" = A_mats, "S" = S_mats, "samples" = samples))
}
################################################################################

################################################################################
#'Find the best number of components
#'
#'@description
#' Analyze the ICA with different number of components to extract the best
#' related to one or more variables.
#'
#'@usage 
#'\code{Best_nc(icas_list, range.comp, metadata, metadata_id = 0, vars, is.categorical = FALSE)}
#'
#'@param icas_list output of `Range_ICA`
#'
#'@param range.comp a range in the format n:n to specify 
#'the range of number of components
#'
#'@param metadata a tibble containing the `vars` to associate and the
#' `metadata_id` column. Not restricted to only these.
#'
#'@param metadata_id a string indicating the name of the column 
#'containing the sample identifier used when running `Range_ICA`.
#'
#'@param vars a string/char vector indicating the variable/s for which you want
#' to maximize the association. For continuous variables more than one 
#' variable can be provided in form of a vector.
#' 
#'@param is.categorical boolean indicating whether the `var` provided is
#' categorical.
#'
#'@details
#' `Best_nc` takes the `Range_ICA` output and analyze it throughout
#'the different range of components  of each ICA analysis with respect to one or
#' more than one variable. The most associated component for each variable is 
#' selected as the representative of the number of components.
#'
#'If the target variable/s are continuous, spearman R is used to evaluate their
#'association.
#'
#'Else, the ratio between the each component Inter quartile range (IQR) and the
#' mean of each component splited by group IQRs is used for this purpose. This 
#'is needed given that sample weights (A matrix) don't usually fulfill normality
#'and homoscedasticity assumptions needed for most parametric tests. 
#'(CURRENTLY UNDER TESTING)
#'
#'@return
#' `Best_nc` returns a ggplot barplot ordered by the mean association of the 
#' number of components with the `vars` of your choice. ggplot syntax can be
#'used to directly alter or export the resulting plot.
#' 
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'@TBD 
#'- indpendence of range.ncomp
#'- output the ggplot and let them do the choice of plot
#'- option to order or not the results

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

  range.comp = names(icas_list$A)
  for (n.comp in range.comp) {
    ica.nc <- icas_list$A[[n.comp]] %>% as_tibble(rownames = "samples")
    test_vars <- dplyr::rename(metadata, samples = all_of(metadata_id)) %>%
      inner_join(ica.nc, by = "samples")
    
    if (is.categorical == FALSE){
      continuous_var <- c(names(ica.nc), vars)
      test_cont <- dplyr::select(test_vars, all_of(continuous_var)) %>%
        drop_na() %>%
        column_to_rownames("samples")
      
      res2 <- rcorr(x=as.matrix(test_cont), type = "spearman")
        
      res2$r <- res2$r[vars, names(ica.nc)[-1],drop=F]
      res2$P <- res2$P[vars, names(ica.nc)[-1],drop=F]
      
      for (cv in vars) { 
        correlations[[cv]]$Rs[[n.comp]] <- res2$r[cv,]
        correlations[[cv]]$Pvals[[n.comp]] <- res2$P[cv,]
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
        correlations[[dv]]$Rs[[n.comp]] <- as.numeric(disp_ratio) %>% set_names(.,names(disp_ratio))
        correlations[[dv]]$Pvals[[n.comp]] <- NULL
      }
    }

    # Check the best IC and save it as a representative of the nc
    for (cv in vars){
      best_c <- sort(abs(correlations[[cv]]$Rs[[n.comp]]),decreasing = T)[1] %>% names()
      to_plot %>% add_row(nc = n.comp,
                          var = cv,
                          IC = best_c,
                          R = abs(correlations[[cv]]$Rs[[n.comp]][[best_c]]),
                          Pval = correlations[[cv]]$Pvals[[n.comp]][[best_c]]) -> to_plot
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
    ggtitle(paste0("ICA Space correlation with variable/s: ", paste(vars, sep = ", "))) +
    rotate_x_text() +
    rremove("xlab")
}
################################################################################

################################################################################
#'Generate a bootstraped dataframe
#'
#'@description
#'bootstrap_expression bootstrap sample or genewise a given tibble containing 
#'expression data.
#'
#'@usage 
#'\code{}
#'
#'@param expression tibble with sample names in columns and genes/probes in
#' rows. Should be mean centered and fairly filtered.
#' 
#'@param expression_id name of the column containing the gene ids
#'
#'@param MARGIN integer either 1 or 2 determining whether rows (genes) or 
#'columns (samples) respectively should be bootstrapped.
#'
#'@param perc double beween 0 and 1 staying the proportion of the dataframe 
#'that should remain unaltered.
#'
#'@details
#' `bootstrap_expression` takes a tibble and bootstrap with replacement either
#'  genes (`MARGIN` = 1), or samples (`MARGIN` = 2) keeping a `perc` of the
#'  dataset unaltered
#'
#'@return
#' `bootstrap_expression` returns a ggplot barplot ordered by the mean association of the 
#' number of components with the `vars` of your choice. ggplot syntax can be
#'used to directly alter or export the resulting plot.
#' 
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'@TBD 
#'- do it the tidyverse way inside!!
bootstrap_expression <- function(expression,
                                 expression_id = "gene",
                                 MARGIN,
                                 perc = 0.9) {
  
  df <- expression %>% column_to_rownames(expression_id)
  
  if (MARGIN == 1) {
    # rows
    n <- round(nrow(df) * (1 - perc))
    
    out <- sample(rownames(df), n)
    df_boot <- df[!(rownames(df) %in% out), ]
    boot <- sample(rownames(df_boot), n)
    df_boot[out, ] <- df_boot[boot, ]
  }
  
  if (MARGIN == 2) {
    # cols
    n <- round(ncol(df) * (1 - perc))
    
    out <- sample(colnames(df), n)
    df_boot <- df[, !(colnames(df) %in% out)]
    boot <- sample(colnames(df_boot), n)
    df_boot[, out] <- df_boot[, boot]
  }
  return(df_boot %>% rownames_to_column(expression_id))
}

################################################################################

################################################################################
#'Find the most robust number of components (less outliers)
#'
#'@description
#'`Boot_ICA` compares the ICA results over a range of components with the same
#' results when bootrapping samples or genes to ensure the robustness of the
#' obtained components
#'@usage 
#'\code{Boot_ICA(expression, df_id, range.comp = 2:20, iterations = 1, seed = 0)}
#'
#'@param expression tibble with sample names in columns and genes/probes in
#' rows. Should be mean centered and fairly filtered.
#' 
#'@param df_id name of the column containing the gene ids
#' 
#'@param range.comp a range in the format n:n to specify the range of number of 
#'components to analyze. Not recomended to use more than 20 components. 
#'Specially with small datasets.
#'
#'@param iterations integer indicating how many bootstraps will be performed for
#' the comparison.
#' 
#'@param seed integer determining the random seed used for the bootstrrapping of 
#'genes/samples
#'
#'@param perc double vector containing the proportion of samples that will be 
#'kept unaltered in each bootstrap. Should be of length 2 and between 0 and 1.
#'
#'@details
#' `Boot_ICA` run  `Range_ICA` over the original ICA samples and for each
#'  iteration it generates a two bootstrapping datasets: one with bootstrapped 
#'  sample and other with bootstrapped genes. Bootstrapping its done with 
#'  replacement. The proportion of entries to replace its determined by the 
#'  `perc` parameter. When bootstrapping samples, the matrices containing gene
#'   weights (S) are compared. When bootstrapping genes, the matrices containing
#'   sample weights (A) are. For each number of components, the correlation between
#'   each original component and the absolute best of all the bootstrapping 
#'   components its saved. 
#'
#'@return
#' `Boot_ICA` returns a tibble with the absolute correlations of each ICA 
#' component with the best matching bootstrapping ICA component for each number
#' of components, iteration, and kind of bootstrapping (A samples, S for genes).
#' This can be a direct input for visualization via ggplot (see example).
#' 
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#' tcga_p <- read_rds("01_Input/tcga_icaready.RDS")
#' boot_res <- Boot_ICA(expression = tcga_p, df_id  =  "gene", range.comp = 2:20, iterations = 2, seed = 0)
#'
#'
#' boot_res %>% dplyr::filter(bootstrap == "A") %>% ggplot() + aes(x = nc, y = correlation) +
#' geom_boxplot(width = 0.5) +
#'   labs(y = "absolute pearson correlation", x = "number of components") +
#'   ggtitle("just sample bootstraps") +
#'   coord_cartesian(ylim = c(0.8, 1)) +
#'   theme_bw()
#' 
#' ggplot(boot_res) + aes(x = nc, y = correlation, fill = bootstrap) +
#'   geom_boxplot(width = 0.5) +
#'   labs(y = "absolute pearson correlation", x = "number of components") +
#'   ggtitle("All bootstraps") +
#'   coord_cartesian(ylim = c(0.8, 1)) +
#'   theme_bw()

#'@TBD 
#'- paralelization via library(foreach) https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
#'- Add translation of A and S to gene and sample

Boot_ICA <- function(expression,
                            df_id = "gene",
                            range.comp = 2:20,
                            iterations = 1,
                            seed = 0,
                            perc = c(0.9, 0.9)) {
  
  original_ICA <- Range_ICA(expression, df_id, range.comp)
  
  set.seed(seed)
  listof_correlations <- list()
  listof_correlations_id <- list()
  results <- list()
  
  for (i in 1:iterations) {
  bootgene_expression <- bootstrap_expression(expression, MARGIN = 1, perc = perc[1])
  bootsample_expression <- bootstrap_expression(expression, MARGIN = 2, perc = perc[2])
  bootstraps <- list(S = bootsample_expression, A = bootgene_expression)

  for (what_boot in c("A", "S")){
    boot_expression <- bootstraps[[what_boot]]

      # loop through iterations
      boot_ICA <- Range_ICA(boot_expression, df_id, range.comp)
      res <- boot_ICA[[what_boot]]
      base_res <- original_ICA[[what_boot]]

      correlations_id <- list()

      for (nc in names(base_res)) {
        # loop through ICA runs
        base_ic <- base_res[[nc]]
        res_ic <- res[[nc]]
        cor_list <- as.list(rep(0, ncol(base_ic)))
        cor_id_list <- as.list(rep(0, ncol(base_ic)))
        names(cor_list) <- colnames(base_ic)
        names(cor_id_list) <- colnames(base_ic)

        for (b_ic in colnames(base_ic)) {
          # loop through boot ICA components
          for (r_ic in colnames(res_ic)) {
            # loop through base ICA components
            pears_c <- cor(base_ic[, b_ic], res_ic[, r_ic])

            if (abs(pears_c) > abs(cor_list[[b_ic]])) {
              # save the highest correlation between each
              # base ICA component and every bootstrap ICA
              cor_list[[b_ic]] <- pears_c
              cor_id_list[[b_ic]] <- r_ic
            }
          }
        }
        # stopifnot(length(unique(unlist(cor_id_list))) == ncol(base_ic))
        listof_correlations[[nc]][[i]] <- cor_list
        # correlations_id[[nc]] <- cor_id_list
      }
      # listof_correlations_id[[i]] <- correlations_id
      results[[what_boot]] <- listof_correlations
    }
  }
  results_tibble <- melt(results) %>%
    as_tibble() %>%
    dplyr::rename(correlation = value, Component = L4, iteration = L3, nc = L2, bootstrap = L1) %>%
    mutate(correlation = abs(correlation),
           nc = factor(nc, level = names(results$A))) %>%
    dplyr::arrange()

  return(results_tibble)
  }

################################################################################

################################################################################
#'Explore sample weights distributions
#'
#'@description
#'`Sampleweights_distribution` explores the component space sample weight in detail to 
#'check the distribution of sample weights with respect to some variables.
#'This can give a clue on what a component is detecting
#'(binary classification, outliers ...).
#'
#'@usage 
#'\code{Sampleweights_distribution(ica, df, df_id, plotfile = NA)}
#'
#'@param ica output of JADE
#'
#'@param df tibble containing the `df_id` and the variables to be tested
#'
#'@param df_id string indicating the name of the `df` sample id column 
#'
#'@param plotfile string indicating the path and name to the pdf where the 
#'function should export the plots. If not provided it will use stdout.
#'
#'@details
#'`Sampleweights_distribution` will plot the sample weights with respect to the desired 
#'variables included in the `df`. For automatic selection of the best 
#'color scale possible ensure your columns are coded according to their content 
#'(continuous as numeric).
#'
#'@return
#'`Sampleweights_distribution`  will return a plot per column of the `df`.
#'
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'#'@TBD 
#'- availability of choosing the vars like the other functions
#'- independence of ggarrange

Sampleweights_distribution <- function(ica, df, df_id, plotfile = NA){
  
  A <- as_tibble(ica[["A"]], rownames = df_id)
  
  stopifnot("identifiers must mach between the df and the ica object" = all(A[[df_id]] == df[[df_id]]))

  if (!is.na(plotfile)){pdf(file=paste(plotfile, ".pdf", sep=""))}
  
  # Sample weights
  for (ann in dplyr::select(df, -all_of(df_id)) %>% names()){
    
    rug_aes <- df[[ann]]
    rug_name <- ann
    comps_plots <- lapply(colnames(A[-1]), function(ic){
      p <- 
        ggplot(A) +
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
    ggarrange(plotlist = comps_plots, common.legend = T,  legend = "bottom") %>% 
      annotate_figure(top = text_grob(ann, color = "black", face = "bold", size = 14)) %>% 
      print()
  }
  if (!is.na(plotfile)) {dev.off()}
}

################################################################################

################################################################################
#'Explore the components of an ICA
#'
#'@description
#'`ICA_explorator` allows the exploration of the component space as a whole, 
#'checking separately the continuous and discrete variables to evaluate the
#'components for associations to whatever metadata or analysis results provided.
#'
#'@usage 
#'\code{ICA_explorator(ica, df, df_cont = NA, df_disc = NA, df_id)}
#'
#'@param ica output of JADE
#'
#'@param df tibble containing the `df_id` and the variables to be tested
#'
#'@param df_cont string or character vector containing the variables you want to
#' test as continuous. Should be column names of `df`.
#' 
#'@param df_disc string or character vector containing the variables you want to
#' test as discrete Should be column names of `df`.
#' 
#'@param df_id string indicating the column name containing the sample ids of
#' `df`.
#'
#'@details
#'`ICA_explorator` will compute different plots depending in whether `df_cont` 
#'and/or `df_disc` are provided. When only one is it will avoid the calculation 
#'of the absent. For `df_cont` Spearman correlation coeficients are computed.
#'For `df_disc` the sample weights of each component are standarized and 
#'represented as the y axis
#'
#'@return
#'`ICA_explorator` will return a list containing one or two ggplot objects that
#' can be easily toggled afterwards using ggplot sintax to match the desired 
#' format. 
#' 
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'@TBD
#'- remove warning NANs when plotting R in correlogram
#'- clean undesired output to stdout
#'- fix ordering of the components TFs

ICA_explorator <- function(ica,
                             df,
                             df_cont = NA,
                             df_disc = NA,
                             df_id) {
  
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

#'Analyze a component sample weights indepth
#'
#'@description
#'`Sampleweights_indepth` allows the testing of a specific component with 
#'relation to a variable.
#'
#'@usage 
#'\code{Sampleweights_indepth(ica, interest_IC, df, df_id = 0, var, test = c(NULL, "pearson", "spearman", "ttest", "wilcox", "anova", "welch", "kruskal"))}
#'
#'@param ica output of JADE
#'
#'@param interest_IC string indicating the component to test. i.e. "IC.1"
#'
#'@param df tibble containing the `df_id` and the variables to be tested
#'
#'@param df_id string indicating the column name containing the sample ids of
#' `df`.
#'
#'@param var string indicating the variable you want to
#' test against. Should be in the column names of `df`.
#' 
#'@param test string indicating the statistical test you would want to perform 
#'to your data if any. Accepted options are (NULL, "pearson", "spearman",
#' "ttest", "wilcox", "anova", "welch", "kruskal"). 
#' 
#'@details
#'Details of the test might change depending in the fullfiled by the IC sample
#'weights. For pearson and spearman data must be in a double format. Normality 
#'is asserted using shapiro.test and a IC of 95%. Homoscedasticity is asserted
#'using fligner.test as it works better in cases where the normality
#'assumption is unfulfilled.
#'
#'@return
#'`Sampleweights_indepth` will return a ggplot object containing the analysis 
#' carried so the user can choose the final layout using ggplot syntax. 
#' For further details check the vignette (02_TCGA/Example_TCGA.Rmd)
#'
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'ic17_OS <- Sampleweights_indepth(ica = best_ica,
#'                      interest_IC = "IC.17",
#'                    df = metadata_pdac_death,
#'                    df_id = "bcr_patient_barcode",
#'                    var = "OS.time",
#'                    test = "spearman")
#'
#'ic17_OS + geom_point() + geom_smooth(method=lm) +
#'  scale_colour_gradientn(colours = terrain.colors(10)) +
#'  theme_bw()
#'  
#'@TBD

Sampleweights_indepth <- function(ica,
                                  interest_IC,
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
  if(vartype == "discrete"){
    v.pv <- fligner.test(get(interest_IC) ~ get(var), data = to_test)$p.value
    if (v.pv < 0.05) {
      eqvar <- FALSE
      } else {
        eqvar <- TRUE
      }
  } else {
    eqvar <- NaN
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
#'Explore Component Genes
#'
#'@description
#' `PlotGeneWeights` allow the exploration of the genes that have the most 
#' impact in the component.
#'
#'#'@usage 
#'\code{PlotGeneWeights(ica, interest_IC, expression, df_id, n_genes = 25, column_annotation = NA)}
#'
#'
#'@param ica output of JADE
#'
#'@param interest_IC string indicating the component to test. i.e. "IC.1"
#'
#'@param expression tibble containing the expression data used for ICA before 
#'standarization `df_id` containing the genes should be a column.
#'
#'@param df_id string indicating the column name containing the gene ids of
#' `df`. Should match those of the `ica`.
#'
#'@param n_genes integer indicating the number of up and down regulated genes 
#'you want to visualize
#'
#'@param column_annotation argument inherited by annotation_col inside pheatmap.
#'Must have the same sample ordering as the `df`.
#' 
#'@details
#'gene abreviations are recommended as df_ids. heatmap computed with pheatmap 
#'and ggplotified via as.ggplot()
#'
#'@return
#'a list of plots containing:
#'- "densplot" which is the density plot of the gene weights and can give useful
#'   information about how many genes are important and how.
#'   
#'- "heatmap" row standardize gene expression of the most positive and most 
#'negative `n_genes`
#'
#'@author 
#'Jacobo Solorzano
#'Network Biology for Immuno-oncology
#'Centre de Recherches en Cancerologie de Toulouse
#'
#'@examples
#'
#'gw <- PlotGeneWeights(ica = best_ica,
#' interest_IC = "IC.17",
#' expression = tcga_gn,
#' df_id = "gene",
#' n_genes  = 25,
#' column_annotation = NA)
#' 
#' ggarrange(gw[["densplot"]] + rremove("xylab"), gw[["heatmap"]], heights = c(1.5, 10), widths = c(0.2,1),
#'           labels = c("A", "B"),
#'           ncol = 2, nrow = 1)

PlotGeneWeights <- function(ica, interest_IC, expression, df_id, n_genes = 25, column_annotation = NA, annotation_colors = NA){
  return_plots <- list()
  
  S <- as_tibble(ica[["S"]], rownames = df_id)
  S_sort <- arrange(S, desc(get(interest_IC)))
  S_mneg <- tail(S_sort, n_genes)
  S_mpos <- head(S_sort, n_genes)
  
  return_plots[["densplot"]] <- ggplot(S) +
    aes_string(y = interest_IC) +
    geom_density(alpha=.5, fill="#AED3FA") +
      geom_hline(yintercept = 0, colour = "black") +
      scale_x_reverse() +
      annotate("rect",xmin = -Inf, xmax = Inf,   ymin = min(S_mneg[interest_IC]), ymax = max(S_mneg[interest_IC]),   fill = "blue", alpha = 0.5) +
      annotate("rect",xmin = -Inf, xmax = Inf,   ymin =  min(S_mpos[interest_IC]), ymax = max(S_mpos[interest_IC]),   fill = "red", alpha = 0.5) +
      theme_classic()

annot_row <- tibble(name = S_mpos[[df_id]], class = "most possitive") %>% add_row(name = S_mneg[[df_id]], class = "most negative") %>% column_to_rownames("name") 
    
    return_plots[["heatmap"]] <- dplyr::filter(expression, get(df_id) %in% rownames(annot_row)) %>%
      dplyr::arrange(match(get(df_id), rownames(annot_row))) %>%
      column_to_rownames(df_id) %>% 
      pheatmap(scale = "row", annotation_row = annot_row, cluster_rows = F, 
               annotation_col = column_annotation, annotation_colors = annotation_colors, show_rownames = F, show_colnames =  F) %>% as.ggplot()
    dev.off()
    
    return(return_plots)
  }
################################################################################