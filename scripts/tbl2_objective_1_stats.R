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
#       R script: tbl2_objective_1_stats                                       #
#       HTML: Table 2. Distribution of HPV genotypes at baseline               #
#       Objective 1: What is the baseline prevalence of each HPV type?         #
#       Descriptive statistics: HPV Genotype Distribution                      #
#       Description: Creating HPV Genotype Distribution Table                  #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)      # read_csv() to load the CSV
library(dplyr)      # filter(), mutate(), if_any(), across(), arrange()
library(stringr)    # str_detect() to find genotype columns
library(lubridate)  # year() for date filtering
library(glue)       # glue() for dynamic headers/labels
library(gtsummary)  # tbl_summary(), add_overall(), modify_*(), as_gt()
library(gt)         # Styling/formatting the final table (tab_header, tab_style, etc.)
library(gtExtras)   # Adds gtsave() and gt preview helpers


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

# --- STEP 4: FILTER THE "visit_num" (BASELINE) ---

# Filter and sort by visit_num
df <- df %>% 
  filter(visit_num == 1)


################################################################################

# --- STEP 5: PRE-VAX HPV GENOTYPE DISTRIBUTION TABLE 2. ---

# ============================================================================ #
#   Objective 1: What is the baseline prevalence of each HPV type?             #
#   Table 2: Distribution of HPV genotypes at baseline                         #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Define the hr and lr genotypes and vaccine sets
# ------------------------------------------------------------------------------
hr_vars <- c("x16","x18","x26","x31","x33","x35","x39","x45","x51",
             "x52","x53","x56","x58","x59","x66","x68","x73","x82")

lr_vars <- c("x6","x11","x40","x42","x43","x54","x61",
             "x67","x69","x70","x71","x72","x84")

# vaccines (presence of ANY type in that vaccine set)
bivalent   <- c("x16","x18")
quadriv    <- c("x6","x11","x16","x18")
nonavalent <- c("x6","x11","x16","x18","x31","x33","x45","x52","x58")

# Genotype columns
geno_cols <- names(df)[str_detect(names(df), "^x\\d+$")]


# ---- Prevalence-based ordering for HR/LR types ----
prev_hr <- vapply(hr_vars, function(v) mean(df[[v]] == 1, na.rm = TRUE), numeric(1))
prev_lr <- vapply(lr_vars, function(v) mean(df[[v]] == 1, na.rm = TRUE), numeric(1))

hr_vars_ord <- names(sort(prev_hr, decreasing = TRUE))
lr_vars_ord <- names(sort(prev_lr, decreasing = TRUE))


# ------------------------------------------------------------------------------
# 2) Build tab1 with all flags needed (Yes/No factors)
# ------------------------------------------------------------------------------
tab1 <- df %>%
  mutate(
    # Overall HPV (any genotype detected at baseline)
    hpv_sum     = rowSums(across(all_of(geno_cols)), na.rm = TRUE),
    any_hpv = factor(if_else(hpv_sum > 0, "Yes", "No"), levels = c("Yes","No")),
    
    # HR/LR (any)
    hr_hpv_types = factor(if_any(all_of(hr_vars), ~ .x == 1), levels = c(TRUE,FALSE),
                          labels = c("Yes","No")),
    lr_hpv_types = factor(if_any(all_of(lr_vars), ~ .x == 1), levels = c(TRUE,FALSE),
                          labels = c("Yes","No")),
    
    # Vaccines (any type in set)
    bivalent_any   = factor(if_any(all_of(bivalent),   ~ .x == 1), levels = c(TRUE,FALSE),
                            labels = c("Yes","No")),
    quadriv_any    = factor(if_any(all_of(quadriv),    ~ .x == 1), levels = c(TRUE,FALSE),
                            labels = c("Yes","No")),
    nonaval_any    = factor(if_any(all_of(nonavalent), ~ .x == 1), levels = c(TRUE,FALSE),
                            labels = c("Yes","No")),
    
    # Single vs Multiple (based on hpv_sum)
    single_hpv    = factor(hpv_sum == 1, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    multiple_hpv  = factor(hpv_sum >= 2, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    
    # Exact number of types (2–8)
    type_2 = factor(hpv_sum == 2, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_3 = factor(hpv_sum == 3, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_4 = factor(hpv_sum == 4, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_5 = factor(hpv_sum == 5, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_6 = factor(hpv_sum == 6, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_7 = factor(hpv_sum == 7, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    type_8 = factor(hpv_sum == 8, levels = c(TRUE,FALSE), labels = c("Yes","No")),
    
    # Each individual genotype is a "Yes/No" factor
    across(all_of(c(hr_vars, lr_vars)),
           ~ factor(if_else(.x == 1, "Yes", "No"), levels = c("Yes","No")))
  ) %>%
  mutate(
    hpv_positive = factor(hpv_positive, levels = c(1,0), labels = c("HPV+","HPV-"))
  )


# ------------------------------------------------------------------------------
# 3) Labels and row order
# ------------------------------------------------------------------------------
labels_list <- list(
  # General HPV block
  any_hpv      ~ "Any HPV",
  hr_hpv_types ~ "HR-HPV types (any)",
  lr_hpv_types ~ "LR-HPV types (any)",

  # HPV Vaccines block
  bivalent_any   ~ "Bivalent types (HPV16, HPV18)",
  quadriv_any    ~ "Quadrivalent types (HPV06, HPV11, HPV16, HPV18)",
  nonaval_any    ~ "Nonavalent types (HPV06, HPV11, HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58)",
  
  # Single/Multiple block
  single_hpv   ~ "Single HPV infection",
  multiple_hpv ~ "Multiple HPV infections",
  
  # Exact counts 2–8
  type_2 ~ "2 types detected",
  type_3 ~ "3 types detected",
  type_4 ~ "4 types detected",
  type_5 ~ "5 types detected",
  type_6 ~ "6 types detected",
  type_7 ~ "7 types detected",
  type_8 ~ "8 types detected"
)

row_order <- c(
  # General HPV block
  "any_hpv","hr_hpv_types","lr_hpv_types",
  
  # HPV Vaccines block
  "bivalent_any","quadriv_any","nonaval_any",
  
  # Single/Multiple block
  "single_hpv","multiple_hpv",
  
  # Exact counts 2–8
  "type_2","type_3","type_4","type_5","type_6","type_7","type_8"
)

N_baseline <- sum(!is.na(tab1$hpv_positive))


# ------------------------------------------------------------------------------
# 4) Build gtsummary with sections
# ------------------------------------------------------------------------------
# build the right-side "n (%)"
pad_1  <- 120
pad_2 <- 92
pad_npc   <- function(lbl, n) paste0(lbl, strrep("\u00A0", n), "n (%)")


# --- build the gtsummary table and create padded 'section' labels ---
tbl_genotypes_fixed <-
  tab1 %>%
  dplyr::select(all_of(c(row_order, "hpv_positive"))) %>%
  gtsummary::tbl_summary(
    by        = hpv_positive,
    type      = list(all_dichotomous() ~ "dichotomous"),
    value     = list(all_dichotomous() ~ "Yes"),
    statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
    missing   = "no",
    label     = labels_list
  ) %>%
  gtsummary::add_overall(last = TRUE) %>%
  gtsummary::modify_header(all_stat_cols() ~ "**{level}**") %>%
  gtsummary::modify_header(stat_0 ~ glue("**N = {scales::comma(N_baseline)}**")) %>%
  gtsummary::modify_table_body(~ .x %>%
                                 dplyr::mutate(
                                   section_base = dplyr::case_when(
                                     variable %in% c("any_hpv","hr_hpv_types","lr_hpv_types") ~ "General HPV types",
                                     variable %in% c("bivalent_any","quadriv_any","nonaval_any") ~ "HPV vaccine types",
                                     variable %in% c("single_hpv","multiple_hpv",
                                                     "type_2","type_3","type_4","type_5","type_6","type_7","type_8") ~
                                       "Single and multiple HPV infections",
                                     TRUE ~ NA_character_
                                   ),
                                   # apply *different* right-padding per group
                                   section = dplyr::case_when(
                                     section_base == "General HPV types" ~ pad_npc("General HPV types", pad_1),
                                     section_base == "HPV vaccine types" ~ pad_npc("HPV vaccine types", pad_1),
                                     section_base == "Single and multiple HPV infections" ~
                                       pad_npc("Single and multiple HPV infections", pad_2),
                                     TRUE ~ section_base
                                   )
                                 )
  ) %>%
  gtsummary::modify_header(label ~ "") %>%
  gtsummary::modify_table_styling(columns = stat_1, hide = TRUE) %>%
  gtsummary::modify_table_styling(columns = stat_2, hide = TRUE) %>%
  gtsummary::modify_footnote(all_stat_cols() ~ NA_character_)


# --- pull the row-group labels that exist after padding ---
grp_final <- tbl_genotypes_fixed$table_body %>%
  dplyr::filter(!is.na(section)) %>%
  dplyr::pull(section) %>%
  unique()


# ------------------------------------------------------------------------------
# 5) Convert to gt and style
# ------------------------------------------------------------------------------
gt_genotypes <-
  tbl_genotypes_fixed %>%
  gtsummary::as_gt(groupname_col = "section") %>%
  gt::tab_style(
    style = gt::cell_text(style = "italic"),
    locations = gt::cells_row_groups(groups = grp_final)
  ) %>%
  gt::tab_style(
    style = gt::cell_text(size = gt::px(20)),
    locations = gt::cells_title(groups = "title")
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_row_groups(groups = grep("n \\(%\\)", grp_final, value = TRUE))
  ) %>%
  gt::tab_header(
    title = gt::md("**Table 2. Distribution of HPV genotypes at baseline**")
  ) %>%
  gt::tab_options(
    heading.padding  = gt::px(2),
    table.font.size  = 20,
    data_row.padding = gt::px(2)
  )

print(gt_genotypes)


# ------------------------------------------------------------------------------
# 6) Save HTML & PNG file
# ------------------------------------------------------------------------------
# gtsave(gt_genotypes, "Table 2. Distribution of HPV genotypes at baseline.html")
# gtsave(gt_genotypes, "Table 2. Distribution of HPV genotypes at baseline.png", vwidth = 1200)

