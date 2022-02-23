library(tidyverse)
library(AnnotationDbi)
library(hgu219.db)
library(gedepir) #remotes::install_github("bcm-uga/gedepir",build_vignettes = TRUE)
library(msigdbr)
library(GSVA)
library(fgsea)
library(hgu219.db)
setwd("/home/jacobo/Documents/05_DeconvolutionReview/01_Puleo/")
source(file = "../R/functions.R" )
################################################################################

puleo <- read_tsv("01_Input/ProcessedExpression.tsv") %>% column_to_rownames("rows")

# Filtering for the ICA
## Filter out of xy chr by probes
x <- hgu219CHR
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
to_keep <- xx[rownames(puleo)] %in% c(1:23, "NULL", "MT") # could change leaving autosomic XY
puleo <- puleo[to_keep, ]


## mean center gene_wise and selection of half mos variable genes
puleo_p <- apply(puleo, 1, function(x) x - mean(x))
puleo_p <- data.frame(t(puleo_p))

### ab_dv function test
# abdv_t <- matrix(c(seq(1,3,1),
#                     seq(1,6,2),
#                     seq(1,9,3)),
#                   nrow = 3, ncol = 3)
#
# apply(abdv_t,1, function(x) sum(abs(x-mean(x)))/length(x))
## out should be 0,2/3,4/3

abdv <- apply(puleo_p, 1, function(x) {
  sum(
    abs(
      x - mean(x)
    )
  ) / length(x)
})

median_abdv <- median(abdv) # median is the half
puleo_p <- puleo_p[abdv > median_abdv, ]
saveRDS(puleo_p, file = "01_Input/puleo_icaready.RDS")


# Annotation of dataset to have Ensembl IDs and Symbols
probes <- row.names(puleo)
probe_info <- AnnotationDbi::select(hgu219.db,
                                    probes,
                                    c("SYMBOL", "ENSEMBL", "GENENAME"))


ensembl_puleo <- Probes_to_whatever(puleo,
                                    probe_info, 
                                    "ENSEMBL")

sym_puleo <- Probes_to_whatever(puleo,
                                probe_info, 
                                "SYMBOL")


## Export both datasets for outside usage
write_tsv(ensembl_puleo %>% rownames_to_column(),
          "01_Input/ProcessedExpression_ensembl.tsv")

write_tsv(sym_puleo %>% rownames_to_column(),
          "01_Input/ProcessedExpression_symbol.tsv")

# Iterative ICA (checked)
icas_list <- Range_ICA(puleo_p, 2:20)
icas_list %>% write_rds("02_Output/puleo_bestfitICs.RDS")


# Iterative ICA Analysis
## Data
icas_list <- read_rds("02_Output/puleo_bestfitICs.RDS")

## Metadata
### Load
metadata <- read_tsv("01_Input/E-MTAB-6134.sdrf.txt")

### Processing
metadata <- as.data.frame(metadata)
metadata <- metadata[,c("Source Name",
                       "Characteristics[sex]",
                       "Characteristics[tumor grading]",
                       "Characteristics[TNM tumour grading]",
                       "Characteristics[resection margin]",
                       "Characteristics[clinical center]",
                       "Characteristics[ffpeblock age]",
                       "Characteristics[os.delay]",
                       "Characteristics[os.event]",
                       "Characteristics[dfs.delay]",
                       "Characteristics[dfs.event]",
                       "Characteristics[hightumcellclassif]",
                       "Characteristics[wholetumclassif]",
                       "Characteristics[average vaf]",
                       "Characteristics[krasmut]",
                       "Characteristics[tp53mut]",
                       "Characteristics[cdkn2amut]")]

colnames(metadata) <- c("sample_name",
                        "sex",
                        "tg",
                        "tnm_tg",
                        "resection_margin",
                        "center",
                        "block_age",
                        "os.delay",
                        "os.event",
                        "dfs.delay",
                        "dfs.event",
                        "hightumcellclassif",
                        "wholetumclassif",
                        "average_vaf",
                        "krasmut",
                        "tp53mut",
                        "cdkn2amut")

metadata$sex[metadata$sex %in% "male"] <- "M"
metadata$sex[metadata$sex %in% "female"] <- "F"

metadata$tg <- gsub("tumour grading G","", metadata$tg)
metadata$tg[metadata$tg %in% "not available"] <- NA

metadata$average_vaf[metadata$average_vaf %in% "not available"] <- NA

metadata$tg <- gsub("tumour grading G","", metadata$tg)
metadata$tg[metadata$tg %in% "not available"] <- NA

metadata$resection_margin <- gsub("resection margin R","",
                                  metadata$resection_margin)
metadata$resection_margin[metadata$resection_margin %in% "not available"] <- NA

metadata$krasmut[metadata$krasmut %in% "no mutation in KRas"] <- 0
metadata$krasmut[metadata$krasmut %in% "mutation in KRas"] <- 1

metadata$tp53mut[metadata$tp53mut %in% "no mutation in TP53"] <- 0
metadata$tp53mut[metadata$tp53mut %in% "mutation in TP53"] <- 1

metadata$cdkn2amut[metadata$cdkn2amut %in% "no mutation in CDKN2a"] <- 0
metadata$cdkn2amut[metadata$cdkn2amut %in% "mutation in CDKN2a"] <- 1

metadata$os.delay[metadata$os.delay %in% "not available"] <- NA
metadata$dfs.delay[metadata$dfs.delay %in% "not available"] <- NA

## get the component most associated to desired vars
Best_nc(icas_list = icas_list,
        range.comp = 2:20,
              metadata = metadata,
              metadata_id = "sample_name",
              cont_vars = c("os.delay","dfs.delay"))

# Best ICA analysis 
elected_ncomp <- 12
elected_ncomp_str <- paste("nc", elected_ncomp, sep = "")

best_ica <- list()
best_ica["S"] <- icas_list[["S"]][elected_ncomp_str]
best_ica["A"] <- icas_list[["A"]][elected_ncomp_str]
best_ica[["samples"]] <- icas_list[["samples"]]


## clinical data
metadata_cont <- c("os.delay", "dfs.delay", "average_vaf", "block_age")
metadata_disc <- colnames(metadata[, !colnames(metadata) %in% metadata_cont])
metadata_disc <- metadata_disc[-1] # exclude sample_name

ICA_explorator(
  ica = best_ica,
  df = metadata,
  df_cont = metadata_cont,
  df_disc = metadata_disc,
  df_id = "sample_name",
  plot_dir = "02_Output/ICA10_plots",
  analysis_name = "clinical_data",
  interest_IC = "IC.9"
)

## Puleo original classification + PAMG (even if data is ICGC or TCGA. Puleo is the latest)
prior_c <- read_tsv("01_Input/subDatasetMolCharac_puleo.tsv")
interest_class <- grep("^Puleo.*?", colnames(prior_c))


interest_class <- c(which(colnames(prior_c) %in% c("sample", "PAMG")),
                    interest_class)
puleo_c <- prior_c[, interest_class]


ICA_explorator(
  ica = best_ica,
  df = puleo_c,
  df_id = "sample",
  df_cont = colnames(puleo_c)[-1],
  df_disc = FALSE,
  plot_dir = "02_Output/ICA10_plots",
  analysis_name = "Classifications"
)

### Double check Supplementary data (ONLY PULEO)
gene_corrs <- read_tsv("01_Input/AllICA_GeneCor.tsv", col_names = TRUE,
                       col_types = "ccnnnnnnnnnn")
gene_corrs <- as.data.frame(gene_corrs)
rownames(gene_corrs) <- gene_corrs[,"probe"]
gene_corrs <- subset(gene_corrs, select = -c(GeneSymbol, probe))
gene_corrs <- gene_corrs[rownames(best_ica[["S"]]),]
res2 <- rcorr(x=as.matrix(best_ica[["S"]]), y=as.matrix(gene_corrs))

res2$r <- res2$r[colnames(best_ica[["S"]]), colnames(gene_corrs)]
res2$P <- res2$P[colnames(best_ica[["S"]]), colnames(gene_corrs)]

corrplot(res2$r, type = "full", order="original", method = "circle",
         p.mat = res2$P, sig.level = 0.05, insig = "blank", is.corr = T,
         mar=c(0,0,1,0)) # http://stackoverflow.com/a/14754408/54964

## Deconvolution
puleo_decon <- read_rds("01_Input/puleo_deconvolutions_m.RDS")
res <- Corr_by_ct(best_ica[["A"]], puleo_decon)



# Component characterization
## Translation of probes to gene IDs of S_mat
S_mat <- best_ica[["S"]]
S_ensembl <- Probes_to_whatever(S_mat,
                            probe_info,
                            "ENSEMBL")

S_sym <- Probes_to_whatever(S_mat,
                                probe_info,
                                "SYMBOL")



GW <- Explore_GW(S_sym = S_sym,
           S_ensembl = S_ensembl,
           elected_ncomp = elected_ncomp,
           interest_IC = "IC.9")

