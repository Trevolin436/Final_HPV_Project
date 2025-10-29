################################################################################
#                                                                              #
########################## *** Data visualizations *** #########################
#                                                                              #
#       Bioinformatics in Project (721):                                       #
#                                                                              #
#       /|\   ____________________________________________________   /|\       #
#       |*|   "Type-Specific HPV Prevalence and Infection Dynamics   |*|       #
#       |*|   in a Pre-Vaccine Cohort"                               |*|       #
#       \|/   ____________________________________________________   \|/       #
#                                                                              #
#       R script: fig3_objective_1_viz                                         #
#       Plot: (a) plot_baseline_hr_hpv_types_per_participant.png               #
#             (b) plot_baseline_lr_hpv_types_per_participant.png               #
#       Objective 1: What is the baseline prevalence of each HPV type?         #
#       Description: Creating Pie Charts for hr & lr types per participant     #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)     # read_csv() to load the CSV
library(dplyr)     # mutate(), filter(), if_any(), rowSums()
library(tidyr)     # complete() for filling missing categories
library(tibble)    # tibble() used to build indicator columns
library(lubridate) # year() for date filtering
library(ggplot2)   # pie charts (geom_col + coord_polar) and saving figures
library(scales)    # percent() formatting in legend labels
library(glue)      # glue() for dynamic titles/labels


################################################################################

# --- STEP 2: READ CSV FILES ---

# 01_merged_FRESH_HPV_data.csv
# df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/2_Data Cleaning/merged_tables_csv/01_merged_FRESH_HPV_data.csv")

# 04_HPV_genotyping_complete.csv
df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/2_Data Cleaning/data_cleaning_csv/04_HPV_genotyping_complete.csv")

# set the new working directory
setwd("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/3_Statistical Analysis/statistical_analysis_outputs")
getwd()


################################################################################

# --- STEP 3: FILTER THE "large_blood_date" (2013-2022)  ---

# Filter and sort by study ID and large blood date
df <- df %>%
  filter(year(large_blood_date) <= 2022) %>%
  arrange(study_id, large_blood_date)


################################################################################

# --- STEP 4: DATA VISUALISATION AND PLOT 2. ---

# ============================================================================ #
#   Objective 1: What is the baseline prevalence of each HPV type?             #
#   Figure 3a: "Baseline" LR-HPV: number of types per participant"             #
#   Figure 3b: "Baseline HR-HPV: number of types per participant"              #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Define HR/LR genotype column names
# ------------------------------------------------------------------------------
hr_type_cols <- c("x16","x18","x26","x31","x33","x35","x39","x45","x51",
                  "x52","x53","x56","x58","x59","x66","x68","x73","x82")

lr_type_cols <- c("x6","x11","x40","x42","x43","x54","x61",
                  "x67","x69","x70","x71","x72","x84")


# ------------------------------------------------------------------------------
# 2) Filter for baseline (visit_num == 1)
# ------------------------------------------------------------------------------
df_base <- df %>% 
  filter(visit_num == 1)


# ------------------------------------------------------------------------------
# 3) Only keep sections & numeric counts
# ------------------------------------------------------------------------------
unique(df_base$hr_hpv_sum)
unique(df_base$lr_hpv_sum)

hr_vars <- c("hr_hpv_types", paste0("hr_n_", 1:7))
lr_vars <- c("lr_hpv_types", paste0("lr_n_", 1:5))

tab1 <- df_base %>%
  mutate(
    # Overall HPV status (for column stratification)
    hpv_positive = factor(hpv_positive, levels = c(1, 0), labels = c("HPV+","HPV-")),
    
    # Any HR / any LR flags
    hr_hpv_types = if_any(all_of(hr_type_cols), ~ .x == 1),
    lr_hpv_types = if_any(all_of(lr_type_cols), ~ .x == 1),
    
    # Count how many distinct HR/LR types are present per person
    hr_count = rowSums(across(all_of(hr_type_cols), ~ as.integer(.x == 1)), na.rm = TRUE),
    lr_count = rowSums(across(all_of(lr_type_cols), ~ as.integer(.x == 1)), na.rm = TRUE)
  ) %>%
  { 
    # Build exact-count indicator columns using the current data (.)
    hr_ind <- tibble::tibble(
      !!!setNames(
        lapply(1:7,  function(k) .[["hr_count"]] == k),
        paste0("hr_n_", 1:7)
      )
    )
    lr_ind <- tibble::tibble(
      !!!setNames(
        lapply(1:5, function(k) .[["lr_count"]] == k),
        paste0("lr_n_", 1:5)
      )
    )
    dplyr::bind_cols(., hr_ind, lr_ind)
  }


# ------------------------------------------------------------------------------
# 4) Generic builder for HR and LR
# ------------------------------------------------------------------------------
pie_legend_only <- function(data, count_col, k_seq,
                            title_text, legend_title,
                            pal_start, pal_end,
                            zero_label = "No type"
) {
  # (1) Summarize including 0 and NA:
  df_raw <- data %>%
    dplyr::mutate(k = {{ count_col }}) %>%
    dplyr::mutate(
      cat = dplyr::case_when(
        is.na(k)       ~ "Missing",
        k == 0         ~ zero_label,
        k %in% k_seq   ~ paste0(k, " type", ifelse(k == 1, "", "s")),
        TRUE           ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(cat)) %>%
    dplyr::count(cat, name = "n")
  
  # (2) Legend order:
  num_labels  <- paste0(k_seq, " type", ifelse(k_seq == 1, "", "s"))
  level_order <- c(zero_label, num_labels, "Missing")
  
  df <- df_raw %>%
    tidyr::complete(cat = factor(level_order, levels = level_order),
                    fill = list(n = 0)) %>%
    dplyr::mutate(
      total = sum(n),
      p     = n / total
    )
  
  # (3) Colors:
  ramp <- grDevices::colorRampPalette(c(pal_start, pal_end))(length(k_seq))
  names(ramp) <- num_labels
  fill_values <- c(setNames("#FFFFFF", zero_label), 
                   ramp 
                   # "Missing" = missing_color
  )
  
  # (4) Legend labels: 
  nice_labels <- stats::setNames(
    glue::glue("{df$cat}: {df$n} ({scales::percent(df$p, accuracy = 0.1)})"),
    as.character(df$cat)
  )
  
  # (5) Plot: 
  ggplot2::ggplot(df, ggplot2::aes(x = "", y = n, fill = cat)) +
    ggplot2::geom_col(width = 1, color = "black", linewidth = 0.4) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = level_order,                        
      labels = nice_labels[level_order],
      name   = legend_title,
      drop   = FALSE,                                 
      guide  = ggplot2::guide_legend(
        override.aes = list(color = "black", linewidth = 0.4),
        title.theme  = ggplot2::element_text(size = 14, face = "bold"),
        label.theme  = ggplot2::element_text(size = 12)
      )
    ) +
    ggplot2::labs(
      title    = title_text,
      subtitle = glue::glue("Number of women = {sum(df$n)}"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
}


# ------------------------------------------------------------------------------
# 4) A. Pie chart: High-risk (0–7 types)
# ------------------------------------------------------------------------------

# Includes "No type" (0) and "Missing" for NA
p_hr <- pie_legend_only(
  data         = tab1,
  count_col    = hr_count,
  k_seq        = 1:7,
  title_text   = "Baseline HR-HPV: number of HPV types detected per participant",
  legend_title = "HR types per participant",
  pal_start    = "#FEE5D9",   
  pal_end      = "#A50F15",   
  # missing_color= "#7B1FA2",
  zero_label   = "No HR-type"
) 

print(p_hr)


# ------------------------------------------------------------------------------
# 5) B. Pie chart: Low-risk (0–5 types)
# ------------------------------------------------------------------------------
p_lr <- pie_legend_only(
  data         = tab1,
  count_col    = lr_count,
  k_seq        = 1:5,
  title_text   = "Baseline LR-HPV: number of HPV types detected per participant",
  legend_title = "LR types per participant",
  pal_start    = "#DEEBF7",   
  pal_end      = "#08519C",   
  # missing_color= "#7B1FA2",
  zero_label   = "No LR-type"
) 

print(p_lr)


# ------------------------------------------------------------------------------
# 6) Save PNG file
# ------------------------------------------------------------------------------
# ggplot2::ggsave("plot_baseline_hr_hpv_types_per_participant.png", plot = p_hr, width = 12, height = 7, dpi = 300)
# ggplot2::ggsave("plot_baseline_lr_hpv_types_per_participant.png", plot = p_lr, width = 12, height = 7, dpi = 300)

