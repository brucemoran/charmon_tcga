library(TCGAbiolinks)
library(tidyverse)
library(biomaRt)

##see https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html

##define a function to run GDCquery, GDCdownload, GDCprepare and return data
gdc_query_download_prepare <- function(PROJECT, DATA_CAT, DATA_TYPE, FILE_TYPE, LEGACY){

  ##make dir as hyphen to underscore
  proj <- gsub("-", "_", PROJECT)
  
  ##if OUT_DIR does not exist, create it
  dir.create(paste0(proj, "/data"), recursive = TRUE, showWarnings = FALSE)

  qry <- GDCquery(project = PROJECT,
                  data.category = DATA_CAT,
                  data.type = DATA_TYPE,
                  file.type = FILE_TYPE,
                  legacy = LEGACY)

  GDCdownload(qry,
              directory = paste0(proj, "/data"),
              method = "api",
              files.per.chunk = 50)

  out_data <- GDCprepare(qry, directory = paste0(proj, "/data"))
  
  ##tell function to send output as 'out_data'
  return(out_data)
}

gdc_qdp_ns <- function(list){
lapply(list, function(f){ 
  tcga_rna <- gdc_query_download_prepare(PROJECT = f,
                                         DATA_CAT = "Transcriptome Profiling",
                                         DATA_TYPE = "Gene Expression Quantification",
                                         FILE_TYPE = "FPKM.txt.gz",
                                         LEGACY = FALSE)
  
  proj <- gsub("-", "_", PROJECT)
  
  ##make tibbles
  ##fpkms
  tcga_rna_fpkm_tb <- tibble::as_tibble(assays(tcga_rna)$`HTSeq - FPKM`, rownames = "ensembl_gene_id") %>%
    dplyr::mutate(across(is.numeric, round, 3))
  
  ##biomart
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  biomart_tb <- tibble::as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                                "external_gene_name"),
                                                 mart = mart)) %>%
                dplyr::filter(ensembl_gene_id %in% tcga_rna_fpkm_tb$ensembl_gene_id)
                   
  tcga_rna_fpkm_bm_tb <- dplyr::left_join(biomart_tb, tcga_rna_fpkm_tb) %>%
                         dplyr::arrange(ensembl_gene_id)
  
  ##log2fpkms
  tcga_rna_log2fpkm_bm_tb <- tcga_rna_fpkm_bm_tb %>%
    dplyr::mutate(across(is.numeric, ~ .x + 0.0001)) %>%
    dplyr::mutate(across(is.numeric, log2)) %>%
    dplyr::mutate(across(is.numeric, round, 3))
  
  ##clinical
  tcga_rna_clin <- tibble::as_tibble(colData(tcga_rna))
  
  tcga_list <- list(tcga_rna, tcga_rna_fpkm_bm_tb, tcga_rna_log2fpkm_bm_tb, tcga_rna_clin)
  list_stubs <- c("tcga_rna", "tcga_rna_fpkm_bm_tb", "tcga_rna_log2fpkm_bm_tb", "tcga_rna_clin")
  list_names <- gsub("tcga", proj, list_stubs)
  names(tcga_list) <- list_names
  
  save_file <- paste0(paste0(proj, "/data/"), proj, ".list.RDS")
  saveRDS(tcga_list, file = save_file)
})
}

##run the function using different inputs
PROJECT_LIST <- list("TCGA-BRCA", "TCGA-COAD", "TCGA-LUAD", "TCGA-OV", "TCGA-UCEC")
gdc_qdp_ns(PROJECT_LIST)
