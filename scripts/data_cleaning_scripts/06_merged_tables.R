################################################################################
#                                                                              #
############################ *** Data Cleaning *** #############################
#                                                                              #
#       Bioinformatics in Project (721):                                       #
#                                                                              #
#       /|\   ____________________________________________________   /|\       #
#       |*|   "Type-Specific HPV Prevalence and Infection Dynamics   |*|       #
#       |*|   in a Pre-Vaccine Cohort"                               |*|       #
#       \|/   ____________________________________________________   \|/       #
#                                                                              #
#       R Script: 06_merged_tables                                             #
#       CSV: 01_merged_FRESH_HPV_data.csv                                      # 
#       Description: Merge the 03. and 04. into one entire table and clean     #
#                    the data.                                                 #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# LIST OF CSV/TSV FILES:
# 1. 01_caprisa_data_cleaning.csv
# 2. 02_CVL_Pellets_data_cleaning.csv
# -> 3. 03_FRESH_HPV_data_cleaning.csv
# -> 4. 04_HPV_genotyping_data_cleaning.csv
# 5. 05_HLA_results_for_FRESH_updated.csv


################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(tidyverse)  # Collection of R packages for data manipulation
library(tidyr)
library(dplyr)      # For data manipulation
library(janitor)    # For cleaning column names and checking data
library(readr )     # For reading CSV file
library(readxl)     # For reading Excel files
library(openxlsx)   # For opening excel file


################################################################################

# --- STEP 2: READ CSV FILES ---

# 03_FRESH_HPV_data_cleaning.csv
df_03 <- read_csv("03_FRESH_HPV_data_cleaning.csv")
colSums(is.na(df_03))
table(duplicated(df_03))
# View(df_03)

# 04_HPV_genotyping_data_cleaning.csv
df_04 <- read_csv("04_HPV_genotyping_complete.csv")
colSums(is.na(df_04))
table(duplicated(df_04)) 
# View(df_04)


################################################################################

# --- STEP 3: FILTER THE "vst_dt" (2013-2022)  ---

# Filter and sort by study ID and vst_dt
df_03 <- df_03 %>%
  filter(year(vst_dt) <= 2022) %>%
  arrange(study_id, vst_dt)

table(df_03$vst_dt)

################################################################################

# --- STEP 4: FILTER THE "visit_num" (BASELINE) ---

# Filter and sort by visit_num
df_04 <- df_04 %>%
  filter(visit_num == 1)


################################################################################

# --- STEP 5: MERGE BY STUDY ID ---


# Merge datasets 03 and 04 by "study_id"
# Alternative (kept for reference): full_join
# merged_df_01 <- dplyr::full_join(df_03, df_04,
#                                  by = "study_id",
#                                  suffix = c(".x", ".y"))

# Inner join used for analysis
merged_df_01 <- dplyr::inner_join(
  df_04, df_03,
  by = "study_id",
  suffix = c(".x", ".y")
)

# Check concordance of visit_code columns from both sources
identical(merged_df_01$visit_code.x, merged_df_01$visit_code.y)

# NA counts for visit_code in each source
sum(is.na(merged_df_01$visit_code.x))  # real NAs = 46
sum(is.na(merged_df_01$visit_code.y))  # real NAs = 0


# ------------------------------------------------------------------------------
# VISIT CODE RECONCILIATION
# ------------------------------------------------------------------------------

# - If both .x and .y are present and different, combine as "x / y"
# - If .x is NA, use .y
# - Otherwise, keep .x
merged_df <- merged_df_01 %>%
  dplyr::mutate(
    visit_code = dplyr::case_when(
      !is.na(visit_code.x) & !is.na(visit_code.y) & visit_code.x != visit_code.y ~
        paste(visit_code.x, visit_code.y, sep = " / "),
      is.na(visit_code.x) ~ visit_code.y,
      TRUE ~ visit_code.x
    )
  ) %>%
  dplyr::select(-visit_code.x, -visit_code.y)

view(merged_df$visit_code)
view(merged_df)

# Quick diagnostics
sapply(merged_df, class)
summary(merged_df)
glimpse(merged_df)
colSums(is.na(merged_df))
table(duplicated(merged_df))
unique(merged_df)


################################################################################

# --- STEP 6: REMOVE DUPLICATES --- 

# Duplicate checks
duplicated(merged_df)
sum(duplicated(merged_df))  # number of row duplicates

# Drop duplicate rows while preserving the first occurrence
remove_duplicate_merged_df_01 <- merged_df %>% dplyr::distinct()

table(duplicated(remove_duplicate_merged_df_01))  # FALSE = 1439


################################################################################

# --- STEP 7: CHECKING AND CLEANING THE MERGED TABLE ---

# NA diagnostics across the full table
# Check and see if there are any real NAs in all columns
is.na(remove_duplicate_merged_df_01) 
# Get row numbers where run_no is NA AFTER coerce
which(is.na(remove_duplicate_merged_df_01))
# Count the number of NAs AFTER coerce
sum(is.na(remove_duplicate_merged_df_01)) 


# --------------------------------------------------------------------------------
# FILTERING AND SORTING
# --------------------------------------------------------------------------------

updated_merged_df_01 <- remove_duplicate_merged_df_01

# Filter by year and then sort by participant ("study_id") and date ("large_blood_date")
updated_merged_df_01 <- updated_merged_df_01 %>%
  dplyr::filter(lubridate::year(large_blood_date) <= 2022) %>%
  dplyr::arrange(study_id, large_blood_date)

table(updated_merged_df_01$large_blood_date)

# Post-filter diagnostics:
# sapply to check the data types
sapply(updated_merged_df_01, class)
# glimpse
glimpse(updated_merged_df_01)
# str
str(updated_merged_df_01)
# check if there are any NAs available in the columns
colSums(is.na(updated_merged_df_01))
# view the merged table
view(updated_merged_df_01)


################################################################################

# --- STEP 8: SAVE THE CSV/TSV ---

# library(readr) # For CSV export

setwd("~/Dropbox/Honours_CBCB/Hons2025/Trevi/")
# getwd()
# 
# # Save merged dataset (by visit_code)
write.csv(updated_merged_df_01, "01_merged_FRESH_HPV_data.csv", row.names = FALSE)

