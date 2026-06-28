#' -------------------------------------------------------------------------
#' Author: Stephen Bradshaw
#' Contact: [stephen.bradshaw@dpird.wa.gov.au]
#' Date: Nov 2025
#' Title: Demonstration code for publication
#' Outline: Utilises a continuous time series of vessel activity (RCN) and periodic survey data from Phone Diary Surveys (PDS) to impute the latter.
#' 
#' Details
#'    - PDS is disparate TS data
#'    - RCN is a continuous TS dataset
#'    - Demonstrate imputation method using rSTAN
#'    - Demonstrate visualisation of results
#' 
#' Resources:
#'    - Demo data provided as an RDS file
#'     
#' Version:
#'  - v1.00 (20251111):
#'  - v1.10 (20260628): Altered model based on reviewer comments to include a logit transformation of the response variable (PDS) to ensure that the imputed values remain within the bounds of 0 and 1. This change was made to improve the model's performance and accuracy in predicting PDS values.
#' ------------------------------------------------------------------------- 

#### Package / Library Requirements ####
rm(list=ls())

#' Packages (not in environment) -------------------------------------------
list.of.packages <- c("magrittr", "tidyr", "dplyr", "stringr", "purrr"
                      , "rebus", "data.table", "lubridate"
                      , "foreach", "doParallel"
                      ,"httr","jsonlite"
                      , "ggplot2", "gridExtra"
                      , "RColorBrewer"
                      , "readxl"
                      , "data.table"
                      , "openxlsx"
                      # , "rstudioapi", "shinystan", "StanHeaders"
                      , "rstan")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#' Libraries ---------------------------------------------------------------
req_packages <- list.of.packages 
sapply(req_packages,require, character.only = TRUE, quietly=TRUE)

rm(list.of.packages, new.packages)
#####

#### Housekeeping ####
#--> Set TZ ####
Sys.setenv(TZ = "Australia/Perth")

#--> Assign numCores ####
numCores = detectCores()-2  # this will use all available cores minus 2
#--> Set Directories ####
tmp_openFile <- dirname(rstudioapi::getActiveDocumentContext()$path) %>% str_split("/") %>% unlist() %>% tail(1)
dirParent <- dirname(rstudioapi::getActiveDocumentContext()$path) %>% str_remove(tmp_openFile)
rm(tmp_openFile)
dirData <- paste0(dirParent, "01_data")
dirScripts <- paste0(dirParent, "02_scripts")
dirOutputs <- paste0(dirParent, "03_outputs")
setwd(dirParent)

#' Checks if directory present, creates if not
#' @param directory_name directory name as a string
#' @returns print statement of outcome
func_checkCreateDirectory <- function(directory_name) {
  if (!dir.exists(directory_name)) {
    dir.create(directory_name)
    print("Directory created successfully.")
  } else {
    print("Directory already exists.")
  }
}

#--> Create Directories ####
dirOutputs <- "03_outputs"
func_checkCreateDirectory(dirOutputs)

#--> Functions ####
source(paste0(dirParent, "/00_src/functions.R"))

#--> Misc Options ####
options(stringsAsFactors = FALSE)
options(scipen = 999)
#####

################################################################################
################################## DATA IMPORT #################################
################################################################################

#### LOAD DATA ####
#--> Read ####
df <- readRDS(dir(dirData, full.names = TRUE)[1])

#--> Review Dataset ####
glimpse(df)
#####

################################################################################
#################################### RUN STAN ##################################
################################################################################

#--> Set up process ####
## Specify data storage
store_imputations <- list()

##Specify specs for model
mod_model <- paste0(dirScripts, "/demo_rSTAN_logit.stan")

ramp_subset <- df$site %>% unique()

#--> Conditional run of stan ####
run_rSTAN <- TRUE
if (run_rSTAN){
  
  iter = 2e3
  
  model <- rstan::stan_model(mod_model, auto_write=FALSE)
  
  for (ramp in ramp_subset){
    
    #### TEST ####
    # ramp <- "DEMO BOAT RAMP"
    #####
    
    #### SET UP DATA #####
    site_list <- list()
    
    df_site <- df %>% 
      filter(site == ramp) %>%
      arrange(date_ym) %>%
      mutate(
        month = as.numeric(format(date_ym, "%m")),  # Extract month for seasonal effects
        year = as.numeric(format(date_ym, "%Y"))
      ) 
    
    year_levels <- df_site$year %>% unique() %>% sort()
    
    print("Original data structure (shows 3 rows):")
    print(head(df_site, 3))
    print(paste("Total observations:", nrow(df_site)))
    print(paste("Missing PDS values:", sum(is.na(df_site$retrievals_pds))))
    
    ## Create missing data indicators
    observed_idx <- which(!is.na(df_site$retrievals_pds))
    missing_idx  <- which(is.na(df_site$retrievals_pds))
    n_observed   <- length(observed_idx)
    n_missing    <- length(missing_idx)
    
    if (n_observed == 0) {
      stop(paste0("No observed PDS values available for site: ", ramp))
    }
    
    ## replace zero or negative SEs with a tiny positive epsilon
    eps <- 1e-6
    df_site <- df_site %>%
      mutate(
        se_retrievals_pds = ifelse(is.na(se_retrievals_pds) | se_retrievals_pds <= 0, eps, se_retrievals_pds),
        se_retrievals_rcn = ifelse(is.na(se_retrievals_rcn) | se_retrievals_rcn <= 0, eps, se_retrievals_rcn)
      )
    
    ## compute summaries
    pds_obs_vec <- df_site$retrievals_pds[observed_idx]
    
    N <- dim(df_site)[1]

    ## Prepare data for Stan
    sdata <- list(
      N = nrow(df_site),
      N_obs = n_observed,
      N_miss = n_missing,
      
      ## Observed data
      obs_idx = observed_idx,
      miss_idx = missing_idx,
      
      ## RCN (predictor, with measurement error)
      rcn_obs = df_site$retrievals_rcn,
      rcn_se_obs = df_site$se_retrievals_rcn,
      
      ## PDS (response)
      pds_obs = df_site$retrievals_pds[observed_idx],
      pds_se_obs = df_site$se_retrievals_pds[observed_idx],
      
      # Month info
      month = as.integer(df_site$month),
      N_months = 12,
      
      # Year info
      year_idx = match(df_site$year, year_levels),
      N_years = length(year_levels)
    )
    
    sdata %>% str()
    #####
    
    #### EXECUTION OF STAN MODEL #####
    #--> Running ####
    set.seed(1)
    fit_mcmc = rstan::sampling(
      object=model
      , data = sdata
      , chains = 4
      , cores = 1 #4 (warning in positron for URL. No effect on output eitherway)
      , iter = iter
      , warmup = iter/2
      , thin = 4
      , control = list(adapt_delta = 0.9, max_treedepth = 12)
    )
    # print(fit_mcmc)
    
    #--> Extracting and checking ####
    ## After fitting the model
    print(fit_mcmc)
    
    ## Extract posterior samples
    fit <- rstan::extract(fit_mcmc)
    
    ## Check structure
    cat("\n=== Posterior Structure ===\n")
    cat("pds_est_constrained dimensions:", dim(fit$pds_est_constrained), "\n")
    cat("  (iterations x N_all)\n")
    cat("pds_full_constrained dimensions:", dim(fit$pds_full_constrained), "\n")
    cat("  (iterations x N_all)\n")
    cat("month_effect dimensions:", dim(fit$month_effect), "\n")
    cat("  (iterations x 12 months)\n\n")
    
    ## [THINKING I SHOULD SAVE THESE?] Visualize key parameters
    check_Posteriors <- TRUE
    if ("code" == "checkPosteriors") {
      #--> [PLOT] MONTHLY BETAS FOR SITE ####
      ## Titles
      month_labels <- paste("Month effect -", month.abb)
      
      ## Set layout for 12 plots
      par(mfrow = c(3, 4),
          mar = c(2, 2, 4, 1),    # smaller margins for individual plots
          oma = c(5, 5, 2, 2))    # outer margins for shared labels
      
      ## Loop over each month
      for (i in 1:12) {
        hist(fit$month_effect[, i],
             main = month_labels[i],
             xlab = "",       # turn off per-plot x label
             ylab = "",       # turn off per-plot y label
             col = "darkgrey",
             border = "white",
             axes = TRUE)     # show axes ticks but not labels
      }
      
      ## Add shared axis labels
      mtext("Frequency", side = 2, outer = TRUE, line = 3, cex = 1.2)   # left y-axis
      mtext("Month effect on logit(q)", side = 1, outer = TRUE, line = 3, cex = 1.2)
      
      ## Reset plotting layout
      par(mfrow = c(1, 1))
      
      #--> [PLOT] Plot Distributions for mean site Beta and variation (SD) by month ####
      ## Reset margins to default
      par(mfrow = c(1, 2),
          mar = c(5, 4, 4, 2) + 0.1,   # default R margins: bottom, left, top, right
          oma = c(0, 2, 0, 0))          # no outer margin
      
      ## BETA CHOSEN FOR EACH TIME STEP --> MEAN
      hist(fit$alpha_0,
           main = "Overall mean (alpha_0)",
           xlab = "",
           ylab = ""
           )
      abline(v = mean(fit$alpha_0), col = "red", lwd = 2)
      
      ## RESIDUAL VARIATION ON LOG1P SCALE
      hist(fit$sigma_log,
           main = "Residual variation (sigma_log)",
           xlab = "",
           ylab = ""
      )
      abline(v = mean(fit$sigma_log), col = "red", lwd = 2)
      
      ## Add shared axis labels
      mtext("Frequency", side = 2, outer = TRUE, line = -0.5, cex = 1.2)   # left y-axis
      
      ## Reset plotting layout
      par(mfrow = c(1, 1))
      #####      
    } #end checking posteriors
    #####
    
    #### Compute Posterior Summaries ####
    #--> Get Posteriors ####
    cat("\n=== Extracting posterior draws ===\n")
  
    ## Extract posterior draws
    pds_mean_draws <- func_selectPDSDraws(fit = fit, site_name = ramp)

    ## Compute posterior mean
    pds_posterior_mean <- apply(pds_mean_draws, 2, mean)
    
    ## Or use median (more robust to outliers)
    pds_posterior_median <- apply(pds_mean_draws, 2, median)
    
    ## Compute posterior SE for each time point
    pds_se <- apply(pds_mean_draws, 2, sd)
    
    ## Compute empirical posterior intervals
    pds_posterior_ll <- apply(pds_mean_draws, 2, quantile, probs = 0.025)
    pds_posterior_ul <- apply(pds_mean_draws, 2, quantile, probs = 0.975)
    
    df_site_out <- df_site
    df_site_out$pds_full_mean <- pds_posterior_mean
    df_site_out$pds_full_se   <- pds_se  # Call it SE (estimation SE)
    
    ## For plotting: empirical 95% posterior interval
    df_site_out$pds_full_ll   <- pds_posterior_ll
    df_site_out$pds_full_ul   <- pds_posterior_ul
    
    ## Imputation flag
    df_site_out <- df_site_out %>%
      mutate(pds_imputed_flag = row_number() %in% sdata$miss_idx)
    
    #--> Store Results ####
    site_list <- list(
      sdata = sdata,
      fit_mcmc = fit_mcmc,
      df_site_out = df_site_out
    )
    store_imputations[[ramp]] <- site_list
    
  } #end loop for stan processing for all ramps (ONLY 1 IN DEMO)
  
  #--> SAVE OUTPUT ####
  saveRDS(store_imputations, file = paste0(dirOutputs, "/stanOuts_", Sys.Date() %>% format("%Y%m%d"), ".rds"))
  #####
  
}
#####

################################################################################
################################# VISUALISATION ################################
################################################################################

#--> Demo output ####
output_files <- dir(dirOutputs, full.names = TRUE)
output_files <- output_files[basename(output_files) %>% str_detect("^stanOuts_\\d{8}\\.rds$")]

if (length(output_files) == 0) {
  stop(paste0("No stan output files found in ", dirOutputs))
}

output_dates <- as.Date(
  sub(".*_(\\d{8})\\.rds$", "\\1", basename(output_files)),
  format = "%Y%m%d"
)

latest_output_file <- output_files[which(output_dates == max(output_dates)) %>% tail(1)]
store_imputations <- readRDS(latest_output_file)

store_imputations %>% lapply(glimpse)


#--> Collate Summary Statistics ####
tablist_summStats <- list()

for (n in 1:length(store_imputations)){
  # n<-1
  summary(store_imputations[[n]]$fit_mcmc)$summary %>% tail(10)
  summary(store_imputations[[n]]$fit_mcmc)$summary[,1] %>% names()
  summary(store_imputations[[n]]$fit_mcmc)$summary %>% colnames()
  
  sum_df <- summary(store_imputations[[n]]$fit_mcmc)$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter")
  
  tablist_summStats[[n]] <- sum_df %>%
    filter(grepl("^(alpha_0|sigma_month|sigma_year|sigma_log|month_effect\\[|year_effect\\[|mean_q|mean_month_effect|mean_year_effect)", parameter))  %>%
    mutate(site = store_imputations %>% names() %>% pluck(n)) %>% 
    select("site", "parameter","n_eff", "Rhat")
  
}

## Combine all sites into one table
tab_summStats <- bind_rows(tablist_summStats)


## Round to 3dp
tab_summStats <- tab_summStats %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

## Save to CSV
write.csv(tab_summStats, paste0(dirOutputs, "/tab_summStats.csv"), row.names = FALSE)


#--> Demo plot ####
ramp <- "DEMO BOAT RAMP"
cat("=== Generating Plot ===\n")
plot_obj <- func_plot_rcn_pds_quality(
  df = store_imputations[[ramp]]$df_site_out,
  site_name = store_imputations[[ramp]]$df_site_out$site[1]
)
print(plot_obj)

#--> Save plots for all sites ####
## A4 Landscape dimensions in inches (minus margins)
width_inches <- (297 - 40)/25.4  # 257mm to inches
height_inches <- (210 - 60)/25.4  # 170mm to inches

## Convert to pixels at 300 DPI
width_px <- width_inches * 300   # ≈ 3035 pixels
height_px <- height_inches * 300  # ≈ 2008 pixels

## Use in your code
jpeg(paste0(dirOutputs,"/plot_imputed_", ramp,".png"),
     width = width_px,    # ≈ 3035 pixels
     height = height_px,  # ≈ 2008 pixels
     units = "px",
     res = 300)

print(
  func_plot_rcn_pds_quality(
    df = store_imputations[[ramp]]$df_site_out,
    site_name = store_imputations[[ramp]]$df_site_out$site[1],
    
    impute_threshold_pct = 90,
    good_label = "<=90% imputed",
    bad_label = ">90% imputed",
    filter_date_from = "2011-09-01",
    filter_date_to   = "2022-06-01",
  )
)
dev.off()

#####


################################################################################
################################################################################
################################################################################
