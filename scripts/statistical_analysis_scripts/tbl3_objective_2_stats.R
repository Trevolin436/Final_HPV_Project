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
#       R script: tbl3_objective_2_stats                                       #
#       HTML: Table 3. HPV Outcomes at Follow-up                               #
#       Objective 2: The percentage of women who clear, acquire or persist     #
#                    with individual and grouped (high-risk vs low-risk)       #
#                    HPV types.                                                #
#       Descriptive statistics: HPV Outcomes at Follow-up                      #                         
#       Description: Creating HPV Outcomes at Follow-up Table                  #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)      # read_csv()
library(dplyr)      # mutate(), filter(), select(), arrange(), group_by(), summarise(), rowwise()
library(tidyr)      # pivot_longer(), pivot_wider()
library(stringr)    # str_detect()
library(lubridate)  # year()
library(gt)         # gt table formatting & gtsave()


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

# --- STEP 4: PRE-VAX HPV OUTCOMES AT FOLLOW-UP TABLE 3. ---

# ============================================================================ #
#   Objective 2: The percentage of women who clear, acquire or persist with    #
#                individual and grouped (high-risk vs low-risk) HPV types.     #
#   Table 3: "HPV Outcomes at Follow-up"                                       #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Define genotype sets
# ------------------------------------------------------------------------------
ID    <- "study_id"
VISIT <- "visit_num"
stopifnot(ID %in% names(df), VISIT %in% names(df))

df <- df %>%
  mutate(visit_num = suppressWarnings(as.numeric(.data[[VISIT]]))) %>%
  filter(visit_num %in% c(1, 2, 3))


# ------------------------------------------------------------------------------
# 2) Genotype columns
# ------------------------------------------------------------------------------
geno_cols <- names(df)[str_detect(names(df), "^x\\d+$")]

# force 0/1
df <- df %>%
  mutate(across(all_of(geno_cols),
                ~ as.integer(pmin(pmax(suppressWarnings(as.numeric(.x)), 0), 1)),
                .names = "{.col}"))

# hpv_sum if missing
if (!"hpv_sum" %in% names(df)) {
  df <- df %>% mutate(hpv_sum = rowSums(across(all_of(geno_cols))))
}


# ------------------------------------------------------------------------------
# 3) Define genotype sets
# ------------------------------------------------------------------------------
hr_full <- c("x16","x18","x26","x31","x33","x35","x39","x45","x51",
             "x52","x53","x56","x58","x59","x66","x68","x73","x82")

lr_full <- c("x6","x11","x40","x42","x43","x54","x61",
             "x67","x69","x70","x71","x72","x84")

hr_types <- intersect(hr_full, geno_cols)
lr_types <- intersect(lr_full, geno_cols)
bivalent <- intersect(c("x16","x18"), geno_cols)
quadriv  <- intersect(c("x16","x18","x6","x11"), geno_cols)
nonaval  <- intersect(c("x6","x11","x16","x18","x31","x33","x45","x52","x58"), geno_cols)


# ------------------------------------------------------------------------------
# 4) Row-wise flags
# ------------------------------------------------------------------------------
df <- df %>%
  mutate(
    hpv_sum   = rowSums(across(all_of(geno_cols)), na.rm = TRUE),
    any_hpv   = as.integer(hpv_sum > 0),
    any_hr    = if (length(hr_types))  as.integer(rowSums(across(all_of(hr_types)),  na.rm = TRUE) > 0) else 0L,
    any_lr    = if (length(lr_types))  as.integer(rowSums(across(all_of(lr_types)),  na.rm = TRUE) > 0) else 0L,
    any_biv   = if (length(bivalent))  as.integer(rowSums(across(all_of(bivalent)),  na.rm = TRUE) > 0) else 0L,
    any_quad  = if (length(quadriv))   as.integer(rowSums(across(all_of(quadriv)),   na.rm = TRUE) > 0) else 0L,
    any_nona  = if (length(nonaval))   as.integer(rowSums(across(all_of(nonaval)),   na.rm = TRUE) > 0) else 0L,
    type_single = as.integer(hpv_sum == 1),
    type_multi  = as.integer(hpv_sum >= 2)
  )


# ------------------------------------------------------------------------------
# 5) Participant & visit wide for a flag
# ------------------------------------------------------------------------------
pivot_flag <- function(flag) {
  df %>%
    select(all_of(c(ID, VISIT, flag))) %>%
    group_by(.data[[ID]], .data[[VISIT]]) %>%
    summarise(
      value = if (all(is.na(.data[[flag]]))) NA_integer_ else max(.data[[flag]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = all_of(VISIT), values_from = value) %>%
    rename(v1 = `1`, v2 = `2`, v3 = `3`)
}

# Compute denominators (using ANY HPV as the cohort anchor)
W_any       <- pivot_flag("any_hpv")
baseline_N  <- sum(!is.na(W_any$v1))                                        # for prevalence
inc_clr_N   <- sum(!is.na(W_any$v1) & !is.na(W_any$v3))                     # for incidence & clearance
pers_N      <- sum(!is.na(W_any$v1) & !is.na(W_any$v2) & !is.na(W_any$v3))  # for persistence
follow_N    <- min(inc_clr_N, pers_N, na.rm = TRUE)                         # for follow-up (all inc/clr/pers)


# ------------------------------------------------------------------------------
# 5.a) Define the complete-follow-up cohort
# ------------------------------------------------------------------------------
W_any       <- pivot_flag("any_hpv")
cohort_ids  <- W_any %>%
  dplyr::filter(!is.na(v1), !is.na(v2), !is.na(v3)) %>%
  dplyr::pull(!!rlang::sym(ID)) %>%
  unique()

follow_N <- length(cohort_ids)


# ------------------------------------------------------------------------------
# 5.b) Participant & visit wide for a flag (restricted to cohort)
# ------------------------------------------------------------------------------
pivot_flag_cohort <- function(flag) {
  pivot_flag(flag) %>%
    dplyr::filter(.data[[ID]] %in% cohort_ids)
}


# ------------------------------------------------------------------------------
# 6) Type-level matrices & baseline denominators (restricted to cohort)
# ------------------------------------------------------------------------------
# Build matrices on the complete-follow-up cohort
make_type_matrix_13 <- function(type_cols) {
  df %>%
    dplyr::filter(.data[[ID]] %in% cohort_ids, .data[[VISIT]] %in% c(1, 3)) %>%
    dplyr::select(all_of(c(ID, VISIT)), all_of(type_cols)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(type_cols),
                        names_to = "type", values_to = "val") %>%
    dplyr::group_by(.data[[ID]], .data[[VISIT]], type) %>%
    dplyr::summarise(
      val = if (all(is.na(val))) NA_integer_ else max(val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(visit_lab = paste0("v", .data[[VISIT]])) %>%
    dplyr::select(-all_of(VISIT)) %>%
    tidyr::pivot_wider(names_from = visit_lab, values_from = val)
}


make_type_matrix_123 <- function(type_cols) {
  df %>%
    dplyr::filter(.data[[ID]] %in% cohort_ids, .data[[VISIT]] %in% c(1, 2, 3)) %>%
    dplyr::select(all_of(c(ID, VISIT)), all_of(type_cols)) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(type_cols),
                        names_to = "type", values_to = "val") %>%
    dplyr::group_by(.data[[ID]], .data[[VISIT]], type) %>%
    dplyr::summarise(
      val = if (all(is.na(val))) NA_integer_ else max(val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(visit_lab = paste0("v", .data[[VISIT]])) %>%
    dplyr::select(-all_of(VISIT)) %>%
    tidyr::pivot_wider(names_from = visit_lab, values_from = val)
}

type_mat_13  <- make_type_matrix_13(geno_cols)    # V1 vs V3
type_mat_123 <- make_type_matrix_123(geno_cols)   # V1, V2, V3

# format
fmt <- function(n, N) sprintf("%d/%d (%.1f%%)", n, N, if (N > 0) 100*n/N else NA_real_)


# ------------------------------------------------------------------------------
# 7) Persistence
# ------------------------------------------------------------------------------
# All three visits
persistence_rule <- function(v1, v2, v3) { (v1 == 1) & (v2 == 1) & (v3 == 1) }


# Baseline presence/absence filters per group (restricted to cohort)
baseline_status_by_group <- function(type_set) {
  if (!length(type_set)) {
    tibble::tibble(!!ID := character(0), present_v1 = logical(0))
  } else {
    v1_mat <- type_mat_123 %>%
      dplyr::filter(type %in% type_set) %>%
      dplyr::group_by(.data[[ID]]) %>%
      dplyr::summarise(any_v1 = any(v1 == 1, na.rm = TRUE), .groups = "drop")
    v1_mat %>% dplyr::rename(!!ID := .data[[ID]], present_v1 = any_v1)
  }
}


# ------------------------------------------------------------------------------
# 8) Core calculators with denominators
# ------------------------------------------------------------------------------
# For Any HPV & HR HPV: denominator = full follow-up cohort (follow_N) for all outcomes
# For LR & Vaccine groups:
#   - Incidence denom = cohort_ids who were baseline ABSENT for that group
#   - Clearance & Persistence denom = cohort_ids who were baseline PRESENT for that group


# calculation INCIDENCE group:
calc_group_incidence <- function(type_set, flag) {
  sub <- type_mat_13 %>% dplyr::filter(type %in% type_set)
  
  if (flag %in% c("any_hpv", "any_hr")) {
    denom_ids <- cohort_ids
  } else {
    base <- baseline_status_by_group(type_set)
    denom_ids <- base %>% dplyr::filter(!present_v1) %>% dplyr::pull(!!rlang::sym(ID))
  }
  
  if (!length(denom_ids) || !nrow(sub)) return(fmt(0, length(denom_ids)))
  
  n_inc <- sub %>%
    dplyr::filter(.data[[ID]] %in% denom_ids) %>%
    dplyr::group_by(.data[[ID]]) %>%
    dplyr::summarise(any_inc = any(v1 == 0 & v3 == 1, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(n = sum(any_inc, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(n)
  
  fmt(n_inc, length(denom_ids))
}


# calculation CLEARANCE group:
calc_group_clearance <- function(type_set, flag) {
  sub <- type_mat_13 %>% dplyr::filter(type %in% type_set)
  
  if (flag %in% c("any_hpv", "any_hr")) {
    denom_ids <- cohort_ids
  } else {
    base <- baseline_status_by_group(type_set)
    denom_ids <- base %>% dplyr::filter(present_v1) %>% dplyr::pull(!!rlang::sym(ID))
  }
  
  if (!length(denom_ids) || !nrow(sub)) return(fmt(0, length(denom_ids)))
  
  n_clr <- sub %>%
    dplyr::filter(.data[[ID]] %in% denom_ids) %>%
    dplyr::group_by(.data[[ID]]) %>%
    dplyr::summarise(any_clr = any(v1 == 1 & v3 == 0, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(n = sum(any_clr, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(n)
  
  fmt(n_clr, length(denom_ids))
}


# calculation PERSISTENCE group:
calc_group_persistence <- function(type_set, flag) {
  sub <- type_mat_123 %>% dplyr::filter(type %in% type_set)
  
  if (flag %in% c("any_hpv", "any_hr")) {
    denom_ids <- cohort_ids
  } else {
    base <- baseline_status_by_group(type_set)
    denom_ids <- base %>% dplyr::filter(present_v1) %>% dplyr::pull(!!rlang::sym(ID))
  }
  
  if (!length(denom_ids) || !nrow(sub)) return(fmt(0, length(denom_ids)))
  
  n_pers <- sub %>%
    dplyr::filter(.data[[ID]] %in% denom_ids) %>%
    dplyr::group_by(.data[[ID]]) %>%
    dplyr::summarise(any_pers = any(persistence_rule(v1, v2, v3), na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(n = sum(any_pers, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(n)
  
  fmt(n_pers, length(denom_ids))
}


# Mapping flags -> type sets 
is_type_set_flag <- function(flag)
  flag %in% c("any_hpv","any_hr","any_lr","any_biv","any_quad","any_nona")

flag_to_types <- function(flag) {
  switch(flag,
         any_hpv = geno_cols,
         any_hr  = hr_types,
         any_lr  = lr_types,
         any_biv = bivalent,
         any_quad= quadriv,
         any_nona= nonaval,
         character(0))
}


# ------------------------------------------------------------------------------
# 9) Row labels & W matrices 
# ------------------------------------------------------------------------------
rows <- dplyr::bind_rows(
  tibble::tibble(
    ROW  = c("Any HPV","HR HPV","LR HPV",
             "Bivalent types (HPV16, HPV18)",
             "Quadrivalent types (HPV16, HPV18, HPV06, HPV11)",
             "Nonavalent types (HPV06, HPV11, HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58)"),
    FLAG = c("any_hpv","any_hr","any_lr","any_biv","any_quad","any_nona")
  )
)


# ------------------------------------------------------------------------------
# 10) Compute table with denominators per the study design
# ------------------------------------------------------------------------------
out <- rows %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    follow_incidence   = calc_group_incidence(flag_to_types(FLAG), FLAG),
    follow_clearance   = calc_group_clearance(flag_to_types(FLAG), FLAG),
    follow_persistence = calc_group_persistence(flag_to_types(FLAG), FLAG)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(row_label = ROW, follow_incidence, follow_clearance, follow_persistence)


# ------------------------------------------------------------------------------
# 11) HTML table with gt
# ------------------------------------------------------------------------------
gt_tbl_3 <-
  out %>%
  gt::gt(rowname_col = "row_label") %>%
  gt::cols_label(
    follow_incidence   = gt::md("**Incidence**"),
    follow_clearance   = gt::md("**Clearance**"),
    follow_persistence = gt::md("**Persistence**")
  ) %>%
  gt::cols_align("center", columns = c(follow_incidence, follow_clearance, follow_persistence)) %>%
  gt::tab_header(
    title = gt::md("**Table 3. HPV Outcomes at Follow-up\\***")
  ) %>%
  gt::tab_footnote(
    footnote = gt::md("**The study followed up all participants with any HR-HPV at baseline, N=324. 
                      Denominators for Any HPV and HR-HPV are 324 as per study design. 
                      Denominators for LR-HPV and vaccine groups are informed by baseline detection (clearance/persistence) or absence (incidence).
                      Incidence: HPV types absent at baseline but present at Visit 3; 
                      Clearance: HPV types present at baseline but absent at Visit 3; 
                      Persistence: HPV types present at all three visits.*"),
  ) %>%
  gt::tab_style(
    style = gt::cell_text(size = gt::px(20)),
    locations = gt::cells_title(groups = "title")
  ) %>%
  gt::tab_options(
    heading.padding  = gt::px(2),
    table.font.size  = 20,
    data_row.padding = gt::px(2)
  )

print(gt_tbl_3)


# ------------------------------------------------------------------------------
# 12) Save HTML & PNG file
# ------------------------------------------------------------------------------
# gtsave(gt_tbl_3, "Table 3. HPV Outcomes at Follow-up.html") 
# gtsave(gt_tbl_3, "Table 3. HPV Outcomes at Follow-up.png", vwidth = 1200)

