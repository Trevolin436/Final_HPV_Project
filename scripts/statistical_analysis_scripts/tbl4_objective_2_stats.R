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
#       R script: tbl4_objective_2_stats                                       #
#       HTML: Table 4. HPV incidence, clearance and persistence across visits  #
#       Descriptive statistics: HPV ICP Across Visits                          #                         
#       Description: Creating HPV ICP Across Visits Table                      #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)      # read_csv()
library(dplyr)      # mutate(), filter(), select(), across(), if_any(), group_by(), summarise(), arrange()
library(tidyr)      # pivot_longer(), pivot_wider()
library(stringr)    # str_detect() to find genotype columns (x##)
library(lubridate)  # year() for date filtering
library(purrr)      # map(), pmap_dfr() for row-wise/list ops
library(rlang)      # sym() for tidy-eval in functions
library(gt)         # build and export the final GT table (gtsave)


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

# --- STEP 4: PRE-VAX HPV PERSISTENCE VS CLEARANCE TABLE 4. ---

# ============================================================================ #
#   Objective 2: The percentage of women who clear, acquire or persist with    #
#                individual and grouped (high-risk vs low-risk) HPV types.     #
#   Table 3: "HPV incidence, clearance and persistence across visits"          #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Define genotype sets
# ------------------------------------------------------------------------------
hr_cols <- c("x16","x18","x26","x31","x33","x35","x39","x45","x51",
             "x52","x53","x56","x58","x59","x66","x68","x73","x82")
lr_cols <- c("x6","x11","x40","x42","x43","x54","x61",
             "x67","x69","x70","x71","x72","x84")
type_cols <- c(hr_cols, lr_cols)

# keep only columns that exists
type_cols <- intersect(type_cols, names(df))
hr_cols   <- intersect(hr_cols, names(df))
lr_cols   <- intersect(lr_cols, names(df))


# ------------------------------------------------------------------------------
# 2) Visit date
# ------------------------------------------------------------------------------
if ("large_blood_date" %in% names(df)) {
  df <- df %>% mutate(visit_date = as.Date(large_blood_date))
} else if ("visit_date" %in% names(df)) {
  df <- df %>% mutate(visit_date = as.Date(visit_date))
} else {
  stop("Please provide a visit date column (e.g., large_blood_date or visit_date).")
}


# ------------------------------------------------------------------------------
# 2) Collapse within-visit replicates (OR multiple rows same date)
#     - OR the same participant has multiple rows per date; "any=1" wins
# ------------------------------------------------------------------------------
df_vis <- df %>%
  dplyr::select(study_id, visit_date, dplyr::all_of(type_cols)) %>%
  dplyr::group_by(study_id, visit_date) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(type_cols),
                  ~ as.integer(any(. == 1, na.rm = TRUE))),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    hpv_sum = rowSums(dplyr::across(dplyr::all_of(type_cols)), na.rm = TRUE)
  ) %>%
  dplyr::arrange(study_id, ) %>%
  dplyr::mutate(
    hr_any = if (length(hr_cols)) dplyr::if_any(dplyr::all_of(hr_cols), ~ .x == 1) else NA,
    lr_any = if (length(lr_cols)) dplyr::if_any(dplyr::all_of(lr_cols), ~ .x == 1) else NA
  )


# ------------------------------------------------------------------------------
# 3) VISIT ORDER & WIDE FLAGS
# ------------------------------------------------------------------------------
df_vis_ord <- df %>%
  dplyr::arrange(study_id, visit_date) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(v_idx = dplyr::row_number()) %>%
  dplyr::ungroup()

pull_v123 <- function(data, var_name) {
  v <- rlang::sym(var_name)
  data %>%
    dplyr::select(study_id, v_idx, flag = !!v) %>%
    tidyr::pivot_wider(
      names_from  = v_idx,
      values_from = flag,
      names_prefix = "V"
    )
  # NOTE: Keep NA to avoid inflating denominators (do not replace NA with 0).
}


# ------------------------------------------------------------------------------
# 4) Build row keys/labels/sections (only for columns that exist)
# ------------------------------------------------------------------------------
has_lr_any <- "lr_any" %in% names(df_vis_ord)

# Build row keys and labels (remove hr_any to avoid duplicate HR row)
row_keys <- c(
  intersect(hr_cols, names(df_vis_ord)),
  if (has_lr_any) "lr_any",
  intersect(lr_cols, names(df_vis_ord))
)

row_labels <- c(
  paste0("HPV", sub("^x", "", intersect(hr_cols, names(df_vis_ord)))),
  if (has_lr_any) "LR-HPV types",
  paste0("HPV", sub("^x", "", intersect(lr_cols, names(df_vis_ord))))
)

row_section <- c(
  rep("High-risk types", length(intersect(hr_cols, names(df_vis_ord)))),
  if (has_lr_any) "Low-risk types",
  rep("Low-risk types", length(intersect(lr_cols, names(df_vis_ord))))
)

# Safety: if nothing to compute, create an empty tibble
if (length(row_keys) == 0) {
  v123_list <- tibble::tibble(study_id=character(), V1=integer(), V2=integer(), V3=integer(), type=character())
} else {
  v123_list <- purrr::map(row_keys, ~{
    pull_v123(df_vis_ord, .x) %>% dplyr::mutate(type = .x)
  }) %>% dplyr::bind_rows()
}


# Convenience sets
all_types <- type_cols
hr_set    <- hr_cols
lr_set    <- lr_cols


# Pre-compute group-wide “any-type per visit” (keep NA as NA)
any_any <- df_vis_ord %>%
  dplyr::group_by(study_id, v_idx) %>%
  dplyr::summarise(
    any_type = if (length(all_types)) as.integer(rowSums(dplyr::across(dplyr::all_of(all_types)), na.rm = TRUE) > 0) else NA_integer_,
    .groups="drop"
  ) %>%
  tidyr::pivot_wider(names_from = v_idx, values_from = any_type, names_prefix = "V") %>%
  dplyr::mutate(type = "any_any")

hr_any_v <- df_vis_ord %>%
  dplyr::group_by(study_id, v_idx) %>%
  dplyr::summarise(
    hr_any = if (length(hr_set)) as.integer(rowSums(dplyr::across(dplyr::all_of(hr_set)), na.rm = TRUE) > 0) else NA_integer_,
    .groups="drop"
  ) %>%
  tidyr::pivot_wider(names_from = v_idx, values_from = hr_any, names_prefix = "V") %>%
  dplyr::mutate(type = "HR_GROUP")

lr_any_v <- df_vis_ord %>%
  dplyr::group_by(study_id, v_idx) %>%
  dplyr::summarise(
    lr_any = if (length(lr_set)) as.integer(rowSums(dplyr::across(dplyr::all_of(lr_set)), na.rm = TRUE) > 0) else NA_integer_,
    .groups="drop"
  ) %>%
  tidyr::pivot_wider(names_from = v_idx, values_from = lr_any, names_prefix = "V") %>%
  dplyr::mutate(type = "LR_GROUP")

# Vaccine groups per visit
grp_flags_v <- function(cols, lab){
  if (!length(cols)) return(tibble(study_id = character(), V1=integer(), V2=integer(), V3=integer(), type = character()))
  df_vis_ord %>%
    dplyr::group_by(study_id, v_idx) %>%
    dplyr::summarise(flag = as.integer(rowSums(dplyr::across(dplyr::all_of(cols)), na.rm = TRUE) > 0), .groups="drop") %>%
    tidyr::pivot_wider(names_from = v_idx, values_from = flag, names_prefix = "V") %>%
    dplyr::mutate(type = lab)
}

bivalent  <- intersect(c("x16","x18"), names(df_vis_ord))
quadriv   <- intersect(c("x6","x11","x16","x18"), names(df_vis_ord))
nonaval   <- intersect(c("x6","x11","x16","x18","x31","x33","x45","x52","x58"), names(df_vis_ord))

g_biv  <- grp_flags_v(bivalent,  "VAX_BIVALENT")
g_quad <- grp_flags_v(quadriv,   "VAX_QUADRIV")
g_nona <- grp_flags_v(nonaval,   "VAX_NONAVAL")


# ------------------------------------------------------------------------------
# 5) Denominator rules (require appropriate visit observed)
# ------------------------------------------------------------------------------
.valid_v3 <- function(tab) !is.na(tab$V3)

# Incidence denom = baseline negative AND has V3 observed
den_inc_base0 <- function(tab) sum(tab$V1 == 0 & .valid_v3(tab), na.rm = TRUE)

# Clearance denom = baseline positive AND has V3 observed
den_pos_base1 <- function(tab) sum(tab$V1 == 1 & .valid_v3(tab), na.rm = TRUE)

# Persistence denom = baseline positive AND both V2 & V3 observed
den_pers_base1_all3obs <- function(tab) {
  sum(tab$V1 == 1 & !is.na(tab$V2) & !is.na(tab$V3), na.rm = TRUE)
}


# ------------------------------------------------------------------------------
# 6) Event counters (missingness - NA excluded)
# ------------------------------------------------------------------------------
# Incidence: negative at V1 -> positive at V3
calc_inc  <- function(tab) sum(tab$V1 == 0 & tab$V3 == 1, na.rm = TRUE)

# Clearance: positive at V1 -> negative at V3
calc_clr  <- function(tab) sum(tab$V1 == 1 & tab$V3 == 0, na.rm = TRUE)

# Persistence (type/group): present at all three visits
calc_pers_all3 <- function(tab) sum(tab$V1 == 1 & tab$V2 == 1 & tab$V3 == 1, na.rm = TRUE)

fmt_nN <- function(n, d) {
  if (is.na(d) || d <= 0) return("0/0 (NA%)")
  sprintf("%d/%d (%.1f%%)", n, d, 100*n/d)
}


# ------------------------------------------------------------------------------
# 7) Build per-row results using correct denominators
# ------------------------------------------------------------------------------
row_compute <- function(type_key, label, section){
  # source table
  src <- if (type_key == "any_any") {
    any_any
  } else if (type_key == "HR_GROUP") {
    hr_any_v
  } else if (type_key == "LR_GROUP") {
    lr_any_v
  } else if (type_key == "VAX_BIVALENT") {
    g_biv
  } else if (type_key == "VAX_QUADRIV") {
    g_quad
  } else if (type_key == "VAX_NONAVAL") {
    g_nona
  } else {
    dplyr::filter(v123_list, type == type_key)
  }
  
  # Numerators
  inc_n <- calc_inc(src)
  clr_n <- calc_clr(src)
  prs_n <- calc_pers_all3(src)
  
  # Denominators
  inc_d <- den_inc_base0(src)
  clr_d <- den_pos_base1(src)
  prs_d <- den_pers_base1_all3obs(src)
  
  tibble::tibble(
    section = section,
    hpv_type = label,
    Incidence   = fmt_nN(inc_n, inc_d),
    Clearance   = fmt_nN(clr_n, clr_d),
    Persistence = fmt_nN(prs_n, prs_d)
  )
}


# ------------------------------------------------------------------------------
# 8) Display keys/labels
# ------------------------------------------------------------------------------
display_keys   <- c("any_any", "HR_GROUP", row_keys, "VAX_BIVALENT", "VAX_QUADRIV", "VAX_NONAVAL")
display_labels <- c("Any HPV", "HR-HPV types", row_labels,
                    "Bivalent types", "Quadrivalent types", "Nonavalent types")
display_secs   <- c("High-risk types", "High-risk types", row_section,
                    "Vaccine groups", "Vaccine groups", "Vaccine groups")

table_dat <- purrr::pmap_dfr(
  list(display_keys, display_labels, display_secs),
  row_compute
)


# ------------------------------------------------------------------------------
# 9) Render GT table
# ------------------------------------------------------------------------------
table_dat_f <- table_dat %>%
  dplyr::filter(!hpv_type %in% c(
    "Any HPV",
    "HR-HPV types",
    "LR-HPV types",
    "Bivalent types",
    "Quadrivalent types",
    "Nonavalent types"
  )) %>%
  # Replace HPV6 with HPV06 in the display labels
  dplyr::mutate(hpv_type = dplyr::if_else(hpv_type == "HPV6", "HPV06", hpv_type))

present_groups <- intersect(
  c("High-risk types", "Low-risk types", "Vaccine groups"),
  unique(table_dat_f$section)
)

tbl_4 <- table_dat_f %>%
  gt::gt(groupname_col = "section", rowname_col = "hpv_type") %>%
  {
    if (length(present_groups)) {
      gt::tab_style(
        .,
        style = gt::cell_text(style = "italic"),
        locations = gt::cells_row_groups(groups = present_groups)
      )
    } else {
      .
    }
  } %>%
  gt::tab_header(title = gt::md("**Table 4. HPV incidence, clearance, and persistence across visits\\***")) %>%
  gt::cols_label(
    Incidence   = gt::md("Incidence"),
    Clearance   = gt::md("Clearance"),
    Persistence = gt::md("Persistence")
  ) %>%
  gt::fmt_markdown(columns = dplyr::everything()) %>%
  gt::tab_stubhead(label = gt::md("**HPV type**")) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels(dplyr::everything())
  ) %>%
  gt::tab_style(
    style = gt::cell_text(size = gt::px(15)),
    locations = gt::cells_title(groups = "title")
  ) %>%
  gt::tab_options(table.font.size = 15, data_row.padding = gt::px(2)) %>%
  gt::tab_footnote(
    footnote = gt::md("**Incidence: HPV types absent at baseline but present at Visit 3.  
                       Clearance: HPV types present at baseline but absent at Visit 3.  
                       Persistence: HPV types present at all three visits.*"),
  )

print(tbl_4)


# ------------------------------------------------------------------------------
# 10) Quick check on V3 cohort size
# ------------------------------------------------------------------------------
# Number of participants with any V3 row (after collapsing)
n_with_v3 <- df_vis_ord %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(has_v3 = any(v_idx == 3), .groups="drop") %>%
  dplyr::summarise(N = sum(has_v3)) %>% dplyr::pull(N)
print(n_with_v3)


# ------------------------------------------------------------------------------
# 11) Save HTML & PNG file
# ------------------------------------------------------------------------------
# gtsave(tbl_4, "Table 4. HPV incidence, clearance, and persistence across visits.html")
# gtsave(tbl_4, "Table 4. HPV incidence, clearance, and persistence across visits.png", vwidth = 1200)

