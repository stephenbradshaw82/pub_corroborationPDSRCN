#' -------------------------------------------------------------------------
#' Project:     Camera Data Analysis Functions
#'
#' Created by:  Stephen Bradshaw
#' Contact: [stephen.bradshaw@dpird.wa.gov.au]
#' Modified:    11/11/2025
#' Version:     1.01 - Initial Release
#'              
#' Purpose:     These functions support PDS-RCN data analysis and plotting
#'              Scripts call functions listed below
#' -------------------------------------------------------------------------

###############################################################
########################## FUNCTIONS ##########################
###############################################################

#' Modification to the standard "in" function
#' Use: x %!in% vector etc
`%!in%` = Negate(`%in%`)


#' Selects posterior PDS draws from an extracted Stan fit object
#' @param fit extracted Stan fit list
#' @param site_name optional site name for clearer error messages
#' @returns matrix/array of posterior draws for the full PDS series
func_selectPDSDraws <- function(fit, site_name = NA_character_) {
  if ("pds_full_constrained" %in% names(fit)) {
    return(fit$pds_full_constrained)
  }
  
  if ("pds_est_constrained" %in% names(fit)) {
    return(fit$pds_est_constrained)
  }
  
  stop(
    paste0(
      "Could not find posterior PDS draws",
      ifelse(is.na(site_name), "", paste0(" for site: ", site_name)),
      ". Available draws: ", paste(names(fit), collapse = ", ")
    )
  )
}
#####


#' Plot function showing RCN and PDS records
#' @param df
#' @param site_name string for plot title (used for labelling)
#' @param impute_threshold_pct numeric threshold for proportion imputed to flag poor quality data (default 80)
#' @param good_label string label for good quality data (default "Good RCN data")
#' @param bad_label string label for poor quality data (default "Poor RCN data (>80% imputed)")
#' @param show_legacy_pds logical whether to show legacy PDS data (default TRUE)
#' @param legacy_pds_color string color for legacy PDS data (default "blue")
#' @param model_pds_color string color for PDS model data (default "darkorange2")
#' @param pds_label_legacy string label for legacy PDS data (default "PDS legacy (obs)")
#' @param pds_label_model string label for PDS model data (default "PDS model (Stan)")
#' @param filter_date_from string date to filter from (default "2011-09-01")
#' @param filter_date_to string date to filter to (default "2024-09-01")
#' @param date_col string name of date column in df (default "date_ym")
#' @param breaks_by string for x-axis breaks (default "3 months")
#' @returns ggplot object
func_plot_rcn_pds_quality <- function(
    df,
    site_name,
    impute_threshold_pct = 80,
    good_label = "Good RCN data",
    bad_label  = "Poor RCN data (>80% imputed)",
    show_legacy_pds = TRUE,
    legacy_pds_color = "blue",
    model_pds_color  = "darkorange2",
    pds_label_legacy = "PDS legacy (obs)",
    pds_label_model  = "PDS model (Stan)",
    filter_date_from = "2011-09-01",
    filter_date_to   = "2022-06-01",
    date_col = "date_ym",
    breaks_by = "3 months") {
  
  # ##### TEST ####
  # df <- df_site_out
  # site_name <- df_site$site[1]
  # impute_threshold_pct <- 70
  # good_label <- "Good RCN data"
  # bad_label  <- "Poor RCN data (>70% imputed)"
  # show_legacy_pds <- TRUE
  # legacy_pds_color <- "blue"
  # model_pds_color  <- "darkorange2"
  # pds_label_legacy <- "PDS legacy (obs)"
  # pds_label_model  <- "PDS model (Stan)"
  # filter_date_from <- "2011-09-01"
  # filter_date_to   <- "2022-06-01"
  # date_col <- "date_ym"
  # breaks_by <- "3 months"
  # ##############
  
  ## Filter date
  date_from <- as.Date(filter_date_from)
  date_to   <- as.Date(filter_date_to)
  df <- df %>% filter(.data[[date_col]] >= date_from & .data[[date_col]] < date_to)
  
  ## X breaks
  month_breaks <- seq(date_from, date_to, by = breaks_by)
  
  ## Threshold for RCN imputation
  max_val <- max(df$prop_imputed_rcn, na.rm = TRUE)
  threshold <- if (max_val <= 1) impute_threshold_pct/100 else impute_threshold_pct
  
  ## Shaded NA regions
  shade_ranges <- df %>%
    arrange(.data[[date_col]]) %>%
    mutate(na_group = cumsum(c(TRUE, diff(has_iSurvey) != FALSE))) %>%
    filter(has_iSurvey==FALSE) %>%
    group_by(na_group) %>%
    summarise(start = min(.data[[date_col]]),
              end   = max(.data[[date_col]]), .groups = "drop")
  
  ## Legacy PDS grouping to prevent joining across gaps
  if (show_legacy_pds) {
    df <- df %>%
      mutate(legacy_pds_flag = !pds_imputed_flag & !is.na(retrievals_pds)) %>%
      mutate(legacy_pds_group = cumsum(c(1, diff(as.numeric(legacy_pds_flag)) != 0)))
  }
  
  ## Build legend color lists dynamically
  color_values <- c(
    `RCN line` = "grey40",
    `PDS model` = model_pds_color,
    `FALSE` = "black",
    `TRUE` = "red"
  )
  color_labels <- c(
    `RCN line` = "RCN line",
    `PDS model` = pds_label_model,
    `FALSE` = good_label,
    `TRUE` = bad_label
  )
  if (show_legacy_pds) {
    color_values <- c(color_values, `PDS legacy` = legacy_pds_color)
    color_labels <- c(color_labels, `PDS legacy` = pds_label_legacy)
  }
  
  ## Plot
  ggplot(df, aes(x = .data[[date_col]], group = 1)) +
    
    ## Shaded NA region
    geom_rect(data = shade_ranges,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.3) +
    
    ## MODELLED PDS (always shown)
    geom_ribbon(aes(ymin = pds_full_ll, ymax = pds_full_ul),
                fill = model_pds_color, alpha = 0.15) +
    
    geom_line(aes(y = pds_full_mean, color = "PDS model"), linewidth = 0.9) +
    
    ## LEGACY PDS (optional)
    {if (show_legacy_pds)
      geom_line(data = df %>% filter(legacy_pds_flag),
                aes(y = retrievals_pds, color = "PDS legacy", group = legacy_pds_group),
                linewidth = 0.6, alpha = 0.8, na.rm = TRUE)
    } +
    {if (show_legacy_pds)
      geom_point(data = df %>% filter(legacy_pds_flag),
                 aes(y = retrievals_pds, color = "PDS legacy"),
                 size = 1.6, alpha = 0.7, na.rm = TRUE)
    } +
    {if (show_legacy_pds)
      geom_errorbar(data = df %>% filter(legacy_pds_flag),
                    aes(ymin = retrievals_pds - se_retrievals_pds,
                        ymax = retrievals_pds + se_retrievals_pds,
                        color = "PDS legacy"),
                    width = 5, alpha = 0.5, na.rm = TRUE)
    } +
    
    ##RCN (unchanged)
    geom_line(aes(y = retrievals_rcn, color = "RCN line"),
              linewidth = 0.4, alpha = 0.8) +
    geom_point(aes(y = retrievals_rcn, color = prop_imputed_rcn > threshold),
               shape = 16, size = 2) +
    geom_errorbar(aes(ymin = retrievals_rcn - se_retrievals_rcn,
                      ymax = retrievals_rcn + se_retrievals_rcn,
                      color = prop_imputed_rcn > threshold),
                  width = 5, alpha = 0.5) +
    
    ## Annotation (top-left inside panel)
    annotate(
      "text",
      x = date_from,#min(df[[date_col]], na.rm = TRUE), 
      y = max(c(df$retrievals_rcn + df$se_retrievals_rcn,
                df$retrievals_pds + df$se_retrievals_pds,
                df$pds_full_ll), na.rm = TRUE),
      
      label = "Error bars = SE; Shaded ribbon = 95% CI (model)",
      hjust = 0, vjust = -2,
      size = 3.5, color = "black"
    ) +
    
    ## Scales
    scale_color_manual(values = color_values, labels = color_labels) +
    scale_x_date(breaks = month_breaks,
                 labels = scales::date_format("%Y-%m"),
                 expand = c(0.01, 0),
                 limits = c(date_from, date_to)) +
    
    labs(title = paste0("RCN & PDS Counts Over Time at ", site_name),
         x = "Date", y = "Monthly Count", color = "Data Source / Quality") +
    
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      legend.position = "top",
      axis.ticks.x = element_line(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
}
#####


###############################################################
###############################################################
###############################################################
