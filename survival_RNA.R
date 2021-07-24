library(tidyverse)
library(devtools)
devtools::install_github("https://github.com/DonaghEgan/rpartSurvivalClassifier", force = TRUE, ref = "dev")
#detach('package:rpartSurvivalClassifier', unload = TRUE)
library(rpartSurvivalClassifier)
library(readxl)

##download clinical survival data
url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6066282/bin/NIHMS978596-supplement-1.xlsx"
temps <- tempfile("./tmp.xlsx")
utils::download.file(url, temps)
tcga_surv_tb <- readxl::read_xlsx(temps, sheet = "TCGA-CDR") %>%
  dplyr::rename(patient = "bcr_patient_barcode")

PROJECT_LIST <- dir(pattern = "TCGA")

##run per project, i.e. TCGA disease type
project_output_list <- lapply(PROJECT_LIST, function(proj){
  
  print(paste0("Working on: ", proj))
  
  ##create output dir, load RDS data, create expr_df
  lapply(c("rpart", "median"), function(dii){
    lapply(c(paste0(proj, "/output/"), "plots"), function(di){
      dir.create(paste(di, dii, sep = "/"), recursive = TRUE, showWarnings = FALSE)
    })
  })
  lapply(c("rpart", "median"), function(di){
    dir.create(paste("plots", di, sep = "/"), recursive = TRUE, showWarnings = FALSE)
  })
  proj_list <- readRDS(file.path(proj, "data", paste0(proj, ".list.RDS")))
  expr_tb <- proj_list[[3]]
  clin_df <- proj_list[[4]]
  rm_cols <- colnames(tcga_surv_tb)[colnames(tcga_surv_tb) %in% colnames(clin_df)]
  rm_cols <- rm_cols[!rm_cols %in% "patient"]
  surv_clin_tb <- tcga_surv_tb %>% dplyr::filter(type %in% gsub("TCGA_", "", proj)) %>%
    dplyr::select(-rm_cols) %>%
    dplyr::left_join(., clin_df, by = "patient") %>%
    dplyr::select(-"...1")
  
  ##find ensembl_gene_id for known gene names and parse from expr_tb
  TRDVs <- unlist(lapply(c("TRDV1", "TRDV2", "TRDV3"), function(f){
    expr_tb[grep(f, unlist(expr_tb[,"external_gene_name"])),1]
  }))
  names(TRDVs) <- c("TRDV1", "TRDV2", "TRDV3")
  expr_df <- expr_tb %>% dplyr::select(-external_gene_name) %>%
    dplyr::filter(ensembl_gene_id %in% TRDVs) %>%
    t() %>%
    tibble::as_tibble(., rownames = "barcode")
  colnames(expr_df) <- expr_df[1,]
  expr_df <- expr_df[-1,]
  colnames(expr_df)[1] <- "barcode"
  expr_df[,2:dim(expr_df)[2]] <- lapply(expr_df[,2:dim(expr_df)[2]], as.numeric)
  
  ##survival using DSS, PFI
  surv_vec <- c("DSS", "PFI")
  
  ################
  ## run rpart ##
  ##############
  
  
  rpart_out_list <- lapply(surv_vec, function(surv){
    print(paste0("Running rpart on: ", surv))
    rpart_surv_tb <- rpartSurvivalClassifier::run_rpart(expr_df = expr_df, 
                                                        gene_ids = TRDVs, 
                                                        clin_df = surv_clin_tb, 
                                                        surv_event = surv, 
                                                        surv_time = paste0(surv, ".time"), 
                                                        join_el = "barcode",
                                                        title_text = proj)
    
    det <- dplyr::select(.data = rpart_surv_tb[[1]], -where(is.list))
    readr::write_csv(det, paste0(proj, "/output/", proj, ".", surv, ".rpart_surv_tb.csv"))
    
    if(!is.null(rpart_surv_tb)){
      print(paste0("Running rpart survival on: ", surv))
      rpart_lrt_list <- rpartSurvivalClassifier::run_surv_plot(clin_tb = rpart_surv_tb[[1]], 
                                                               gene_ids = TRDVs, 
                                                               group_name = "_rpart_group",
                                                               surv_event = surv, 
                                                               surv_time = paste0(surv, ".time"),
                                                               title_text = proj,
                                                               sub_text = "rpart stratification")
      
      #plot outputs from rpart splits
      lapply(rpart_surv_tb, function(p){
        
        ##remove NULL
        rpart_plots_list <- rpart_surv_tb[[2]][!sapply(rpart_surv_tb[[2]], is.null)]
        
        ##pdf outputs
        pdf(paste0(proj, "/output/rpart/rpart_", proj, "_", surv, ".pdf"))
        pdf(paste0("plots/rpart/rpart_", proj, "_", surv, ".pdf"))  
        
        ##apply over list of plot functions
        lapply(rpart_plots_list, function(pp){
          pp()
        })
        dev.off()
        dev.off()
      })
      
      #plot outputs from rpart survival
      rpart_surv_plots_list <- lapply(rpart_lrt_list, function(p){
        p[[1]]
      })
      
      ##again apply over plot list
      lapply(seq_along(rpart_surv_plots_list), function(pp){
        if(class(rpart_surv_plots_list[[pp]]) %in% "ggsurvplot"){
          rpart_surv_plots_list[[pp]]
          ggplot2::ggsave(filename = paste0(proj, "/output/rpart/rpart_", proj, "_", names(rpart_surv_plots_list)[pp], ".", surv, ".pdf"))
          ggplot2::ggsave(filename = paste0("plots/rpart/rpart_", proj, "_", names(rpart_surv_plots_list)[pp], ".", surv, ".pdf"))
          
        }
      })
      return(list(rpart_surv_tb = rpart_surv_tb, 
                  rpart_lrt_list = rpart_lrt_list, 
                  rpart_surv_plots_list = rpart_surv_plots_list))
    }
  })
  
  names(rpart_out_list) <- surv_vec
  
  #####################
  ## median survival ##
  ####################
  
  ##function
  med_func <- function(df){
    tibble::as_tibble(as.data.frame(lapply(seq_along(df), function(x){
      vec <- df[[x]]
      if(is.numeric(vec[1])){
        medv <- median(vec)
        ifelv <- ifelse(vec > medv, "High", "Low")
        tbo <- tibble::tibble(ifelv, vec)
        colnames(tbo) <- paste(names(df)[x], c("median_group", "log2tpm"), sep="_")
        return(tbo)
      } else {
        tbo <- tibble::tibble(vec)
        colnames(tbo) <- names(df)[x]
        return(tbo)
      }
    })))
  }
  surv_median_tb <- med_func(expr_df)
  
  ##join with clinical survival data
  median_surv_tb <- dplyr::left_join(surv_median_tb, surv_clin_tb)
  det <- dplyr::select(.data = median_surv_tb, -where(is.list))
  readr::write_csv(det, paste0(proj, "/output/", proj, ".", surv, ".median_surv_tb.csv"))
  
  median_out_list <- lapply(surv_vec, function(surv){
    print(paste0("Running median survival on: ", surv))
    
    ##run survival         
    median_lrt_list <- rpartSurvivalClassifier::run_surv_plot(clin_tb = median_surv_tb, 
                                                              gene_ids = TRDVs, 
                                                              group_name = "_median_group",
                                                              surv_event = surv, 
                                                              surv_time = paste0(surv, ".time"),
                                                              title_text = proj,
                                                              sub_text = "Median stratification")
    
    median_surv_plots_list <- lapply(median_lrt_list, function(p){
      p[[1]]
    })
    lapply(seq_along(median_surv_plots_list), function(pp){
      if(class(median_surv_plots_list[[pp]]) %in% "ggsurvplot"){
        median_surv_plots_list[[pp]]
        ggplot2::ggsave(filename = paste0(proj, "/output/median/median_", proj, "_", names(median_surv_plots_list)[pp], ".", surv, ".pdf"))
        ggplot2::ggsave(filename = paste0("plots/median/median_", proj, "_", names(median_surv_plots_list)[pp], ".", surv, ".pdf"))
      }
    })
    return(list(median_surv_tb = median_surv_tb, 
                median_lrt_list = median_lrt_list, 
                median_surv_plots_list = median_surv_plots_list))
  })
  
  names(median_out_list) <- surv_vec
  
  ##save generated data
  print("Saving and returning data...")
  save_file <- paste0(paste0(proj, "/output/"), proj, ".output_list.RDS")
  saveRDS(list(rpart_out_list = rpart_out_list, 
               median_out_list = median_out_list), 
          file = save_file)
  
  return(list(rpart_out_list = rpart_out_list, 
              median_out_list = median_out_list))
})

names(project_output_list) <- PROJECT_LIST
save_file <- paste0("survival_RNA.output_list.RDS")
saveRDS(project_output_list, file = save_file)

##write rmarkdown HTMLs
rmarkdown::render("./plots/rmds/TCGA_BRCA.Rmd", output = "plots/htmls/TCGA_BRCA.html")
rmarkdown::render("./plots/rmds/TCGA_COAD.Rmd", output = "plots/htmls/TCGA_COAD.html")
rmarkdown::render("./plots/rmds/TCGA_LUAD.Rmd", output = "plots/htmls/TCGA_LUAD.html")
rmarkdown::render("./plots/rmds/TCGA_OV.Rmd", output = "plots/htmls/TCGA_OV.html")
rmarkdown::render("./plots/rmds/TCGA_UCEC.Rmd", output = "plots/htmls/TCGA_UCEC.html")

