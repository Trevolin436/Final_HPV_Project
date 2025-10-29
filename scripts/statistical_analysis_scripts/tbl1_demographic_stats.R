################################################################################
#                                                                              #
######################### *** Statistical Analysis *** #########################
#                                                                              #
#       Bioinformatics in Project (721):                                       #
#                                                                              #
#       /|\   ____________________________________________________   /|\       #
#       |*|   "Type-Specific HPV Prevalence and Infection Dynamics   |*|       #
#       |*|   in a Pre-Vaccine Cohort"                               |*|       #
#       \|/   ____________________________________________________   \|/       #
#                                                                              #
#       R script: tbl1_demographic_stats                                       #
#       HTML: Table 1. Baseline Characteristics                                #
#       Descriptive statistics: Baseline Characteristics (Demographic)         #
#       Description: Creating Baseline Characteristics Table                   #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(tidyverse)  # Includes dplyr, tidyr, readr, stringr, purrr — for data manipulation, reshaping, and reading CSVs
library(lubridate)  # For working with and filtering date variables
library(forcats)    # For ordering and relabeling categorical factors (HPV+, HPV−)
library(glue)       # For dynamic text and labels inside tables
library(gt)         # For creating the publication-ready Table 1 output


################################################################################

# --- STEP 2: READ CSV FILES ---

# 01_merged_FRESH_HPV_data.csv
df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/HPV_Project_Objectives/2_Data Cleaning/01_merged_FRESH_HPV_data.csv")

# 03_FRESH_HPV_data_cleaning.csv
# df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/HPV_Project_Objectives/2_Data Cleaning/03_FRESH_HPV_data_cleaning.csv")

# set the new working directory
setwd("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/HPV_Project_Objectives/5_Final Results/tables")
getwd()


################################################################################

# --- STEP 3: FILTER THE "large_blood_date" (2013-2022)  ---

# Filter and sort by study ID and large blood date
df <- df %>%
  filter(year(large_blood_date) <= 2022) %>%
  arrange(study_id, large_blood_date)


################################################################################

# --- STEP 4: FILTER THE "visit_num" (BASELINE) ---

# Filter and sort by visit_num
df <- df %>% 
  filter(visit_num == 1)


################################################################################

# --- STEP 5: PRE-VAX BASELINE CHARACTERISTICS TABLE 1. ---

# ============================================================================ #
#   Table 1: Baseline Characteristics                                          #
#   - Demographic Data                                                          #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Recode variables: Keep numeric vars numeric for mean/SD and median/IQR.
# ------------------------------------------------------------------------------
tab1 <- df %>%
  mutate(
    live_with_partner_baseline = factor(live_with_partner_baseline,
                                        levels = c(1, 2),
                                        labels = c("Yes", "No")),
    
    sti_symptoms_now = factor(sti_symptoms_now,
                              levels = c(0, 1),
                              labels = c("No", "Yes")),
    
    condoms_used_30_days = case_when(
      condoms_used_30_days == 1 ~ "Never",
      condoms_used_30_days == 2 ~ "Occasional",
      condoms_used_30_days == 3 ~ "Always",
      TRUE ~ NA_character_
    ) |> factor(levels = c("Never", "Occasional", "Always")),
    
    hpv_positive = factor(hpv_positive, levels = c(1, 0), labels = c("HPV+", "HPV-"))
  )

tab1 <- tab1 %>%
  mutate(hpv_positive = forcats::fct_relevel(hpv_positive, "HPV+", "HPV-"))

grp_levels <- c("HPV+", "HPV-")
g <- tab1$hpv_positive


# ------------------------------------------------------------------------------
# 2) Numeric formatters 
# ------------------------------------------------------------------------------
calc_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[!is.na(x)]
  if (!length(x)) return("—")
  sprintf("%.1f (%.1f)", mean(x), stats::sd(x))
}

calc_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[!is.na(x)]
  if (!length(x)) return("—")
  q <- stats::quantile(x, c(.25,.5,.75), na.rm = TRUE)
  sprintf("%.1f (%.1f–%.1f)", q[2], q[1], q[3])
}


# Grouped numeric summaries > character strings for "Total" / "HPV+" / "HPV-"
fmt_num_by_group <- function(x, g, fun) {
  out_total <- fun(x)
  parts <- lapply(grp_levels, function(lv) fun(x[g == lv]))
  c(out_total, setNames(parts, grp_levels))
}


# ------------------------------------------------------------------------------
# 3) Categorical block (subheading + child rows) for 3 columns
#    - Percent is out of non-missing within each column
#    - Missing row (count only)
# ------------------------------------------------------------------------------
pct_str <- function(n, d) if (d > 0) sprintf("%d (%.1f%%)", n, 100*n/d) else "0 (0.0%)"

expand_levels_block3 <- function(x, label, level_order, level_labels = NULL,
                                 indent = "    ", include_missing = TRUE,
                                 missing_label = "Missing") {
  x_chr <- as.character(x)
  f_all <- factor(x_chr, levels = level_order)
  
  # Totals:
  tab_all <- table(f_all, useNA = "no")
  d_all   <- sum(tab_all)
  if (is.null(level_labels)) level_labels <- level_order
  cnt_all <- as.numeric(tab_all)
  total_col <- if (d_all > 0) sprintf("%d (%.1f%%)", cnt_all, 100*cnt_all/d_all) else sprintf("%d (0.0%%)", cnt_all)
  
  
  # By group:
  build_group_col <- function(level_value) {
    f_g   <- factor(x_chr[g == level_value], levels = level_order)
    tab_g <- table(f_g, useNA = "no")
    d_g   <- sum(tab_g)
    cnt_g <- as.numeric(tab_g)
    if (d_g > 0) sprintf("%d (%.1f%%)", cnt_g, 100*cnt_g/d_g) else sprintf("%d (0.0%%)", cnt_g)
  }
  col_pos <- build_group_col("HPV+")
  col_neg <- build_group_col("HPV-")
  
  
  # Core rows (subheading & levels):
  core <- tibble::tibble(
    Characteristic = c(label, paste0(indent, level_labels)),
    Total = c("", total_col),
    `HPV+` = c("", col_pos),
    `HPV-` = c("", col_neg)
  )
  
  
  # missing row (counts by column, no %):
  if (include_missing) {
    miss_total <- sum(is.na(x))
    miss_pos   <- sum(is.na(x[g == "HPV+"]))
    miss_neg   <- sum(is.na(x[g == "HPV-"]))
    miss_row <- tibble::tibble(
      Characteristic = paste0(indent, missing_label),
      Total = if (miss_total > 0) sprintf("%d", miss_total) else "0",
      `HPV+` = if (miss_pos   > 0) sprintf("%d", miss_pos)   else "0",
      `HPV-` = if (miss_neg   > 0) sprintf("%d", miss_neg)   else "0"
    )
    dplyr::bind_rows(core, miss_row)
  } else {
    core
  }
}

# Convenience wrappers
expand_yes_no_3   <- function(x, label, ...) expand_levels_block3(x, label, c("No","Yes"), ...)
expand_never_3    <- function(x, label, ...) expand_levels_block3(x, label, c("Never","Occasional","Always"), ...)


# ------------------------------------------------------------------------------
# 4) Assemble all rows
# ------------------------------------------------------------------------------
sum_rows3 <- dplyr::bind_rows(
  # Numeric
  {
    vals <- fmt_num_by_group(tab1$age_baseline, g, calc_mean_sd)
    tibble::tibble(Characteristic = "Age [mean (SD)]",
                   Total = vals[[1]], `HPV+` = vals[[2]], `HPV-` = vals[[3]])
  },
  {
    vals <- fmt_num_by_group(tab1$age_first_sex_baseline_demo, g, calc_mean_sd)
    tibble::tibble(Characteristic = "Age at sexual debut [mean (SD)]",
                   Total = vals[[1]], `HPV+` = vals[[2]], `HPV-` = vals[[3]])
  },
  {
    vals <- fmt_num_by_group(tab1$number_sexual_partners, g, calc_median_iqr)
    tibble::tibble(Characteristic = "Number of sexual partners [median (IQR)]",
                   Total = vals[[1]], `HPV+` = vals[[2]], `HPV-` = vals[[3]])
  },
  {
    vals <- fmt_num_by_group(tab1$num_sex_episodes_7_day, g, calc_median_iqr)
    tibble::tibble(Characteristic = "Sex events in past 7 days [median (IQR)]",
                   Total = vals[[1]], `HPV+` = vals[[2]], `HPV-` = vals[[3]])
  },
  {
    vals <- fmt_num_by_group(tab1$num_sex_partners_30_day, g, calc_median_iqr)
    tibble::tibble(Characteristic = "Number of sex partners in past 30 days [median (IQR)]",
                   Total = vals[[1]], `HPV+` = vals[[2]], `HPV-` = vals[[3]])
  },
  
  # Categorical blocks (Missing counts)
  expand_never_3 (tab1$condoms_used_30_days,          "Condom use in past month [n (%)]", include_missing = FALSE),
  expand_yes_no_3(tab1$sti_symptoms_now,              "Current STI symptoms [n (%)]", include_missing = FALSE)
)


# Flags for styling
sum_rows3 <- sum_rows3 %>%
  dplyr::mutate(
    is_subheading = Total == "",
    is_child      = Total != "" & stringr::str_detect(Characteristic, "^\\s+")
  ) %>%
  dplyr::select(-is_subheading, -is_child)


# N for header (overall row-counts per group - denominators vary per row)
N_total <- nrow(tab1)
N_pos   <- sum(g == "HPV+")
N_neg   <- sum(g == "HPV-")


# ------------------------------------------------------------------------------
# 5) Render with gt
# ------------------------------------------------------------------------------
gt_baseline_1 <-
  sum_rows3 %>%
  gt::gt() %>%
  gt::tab_header(
    title    = gt::md("**Table 1. Baseline Characteristics**"),
  ) %>%
  gt::cols_label(
    Characteristic = gt::md("**Characteristic**"),
    Total          = gt::md(glue::glue("**Total<br>(N = {N_total})**")),
    `HPV+`         = gt::md(glue::glue("**HPV+<br>(N = {N_pos})**")),
    `HPV-`         = gt::md(glue::glue("**HPV−<br>(N = {N_neg})**"))
  ) %>%
  gt::tab_options(
    table.font.size = 20,
    heading.padding = gt::px(2),
    data_row.padding = gt::px(2),
    table_body.hlines.style = "solid",
    table_body.hlines.width = gt::px(0.5)
  ) %>%
  gt::tab_style(
    style = gt::cell_text(size = gt::px(20)),
    locations = gt::cells_title(groups = "title")
  )

print(gt_baseline_1)


# ------------------------------------------------------------------------------
# 6) Save HTML & PNG file
# ------------------------------------------------------------------------------
# gtsave(gt_baseline_1, "Table 1. Baseline Characteristics.html")
# gtsave(gt_baseline_1, "Table 1. Baseline Characteristics.png", vwidth = 1200)

