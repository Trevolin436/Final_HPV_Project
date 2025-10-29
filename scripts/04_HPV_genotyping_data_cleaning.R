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
#       R script: 04_HPV_genotyping_data_cleaning                              #
#       CSV: (1) 04_HPV_genotyping_final_data_cleaning.csv                     # 
#            (2) 04_HPV_genotyping_complete.csv                                #
#       Description: Cleaning the entire table                                 #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(tidyverse)  # Collection of R packages for data science and manipulation
library(dplyr)      # For data manipulation
library(janitor)    # For cleaning column names and checking data
library(tidyr)      # For reshaping and handling missing values
library(readxl)     # For reading Excel files
library(openxlsx)   # For only excels


################################################################################

# --- STEP 2: READ & LOAD THE EXCEL " "04. HPV_genotyping_data" ---

# HPV_genotyping_data_04 <- read.xlsx("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/1_HPV Data/HPV_raw_data/04. HPV_genotyping_data.xlsx")
# excel_sheets("04.HPV_genotyping_data.xlsx") # List all sheets
HPV_genotyping_data <- read_excel("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/1_HPV Data/HPV_raw_data/04. HPV_genotyping_data.xlsx", sheet = 2)  # Or use sheet name
# view(HPV_genotyping_data)

# set the new working directory
setwd("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/HPV_Project_Objectives/2_Data Cleaning/data_cleaning_csv")
getwd()


################################################################################

# --- STEP 3: FILTER THE "large_blood_date" (2013-2022)  ---

# Filter and sort by study ID and large blood date
HPV_genotyping_data <- HPV_genotyping_data %>%
  filter(year(large_blood_date) <= 2022) %>%
  arrange(study_id, large_blood_date)
# view(HPV_genotyping_data)
table(HPV_genotyping_data$large_blood_date)


################################################################################

# --- STEP 3: COMPREHENSIVE DATA PREPROCESSING AND QUALITY ASSURANCE --- 


# Initialize data subset with essential variables for analysis pipeline
HPV_genotyping_data <- HPV_genotyping_data %>%
  select(
    study_id,
    visit_code,
    large_blood_date,
    nugent_gardnerella,
    nugent_mobiluncus,
    nugent_lactobacillus,
    nugent_total_score,
    nugent_interpretation,
    sti_mycoplasma_genitalium_pcr,
    sti_trichomonas_vaginalis_pcr,
    sti_neisseria_gonorrhoeae_pcr,
    sti_chlamydia_trachomatis_pcr,
    `16`, `18`, `26`, `31`, `33`, `35`, `39`, `45`, `51`, `52`, `53`,
    `56`, `58`, `59`, `66`, `68`, `73`, `82`, `6`, `11`, `40`, `42`,
    `43`, `54`, `61`, `67`, `69`, `70`, `71`, `72`, `84`,
    REPORT,
    COMMENT)


# ------------------------------------------------------------------------------
# DATA INTEGRITY ASSESSMENT: IDENTIFYING AND REMOVING NULL RECORDS
# ------------------------------------------------------------------------------

# Comprehensive identification of completely empty observations
null_observations <- HPV_genotyping_data %>%
  dplyr::filter(dplyr::if_all(dplyr::everything(), is.na))

# Display count of null observations for quality control
cat("Quality Control Report: Found", nrow(null_observations), "completely empty records")

# Remove null observations from primary dataset
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::filter(!dplyr::if_all(dplyr::everything(), is.na))


# ------------------------------------------------------------------------------
# DUPLICATE RECORD MANAGEMENT: DETECTION AND ELIMINATION
# ------------------------------------------------------------------------------

# Systematic identification of duplicate entries
duplicate_entries <- HPV_genotyping_data %>%
  dplyr::filter(duplicated(.))

# Quality control reporting for duplicate records
cat("Quality Control Report: Identified", nrow(duplicate_entries), "duplicate records")

# Remove all duplicate entries while preserving first occurrence
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::distinct()

# Verification of duplicate removal
duplicate_verification <- if (any(duplicated(HPV_genotyping_data))) {
  cat("Warning: Duplicate records still present in dataset\n")
} else {
  cat("Success: All duplicate records have been successfully removed\n")
}


# ------------------------------------------------------------------------------
# COLUMN NAMING STANDARDIZATION
# ------------------------------------------------------------------------------

# Standardize column names for consistency across analysis pipeline
HPV_genotyping_data <-HPV_genotyping_data %>% 
  clean_names()


# ------------------------------------------------------------------------------
# TEMPORAL DATA ORGANIZATION
# ------------------------------------------------------------------------------

# Sort dataset chronologically by participant and visit date for longitudinal analysis
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::arrange(study_id, large_blood_date)


# ------------------------------------------------------------------------------
# VISIT FREQUENCY FILTERING: EXCLUDING EXCESSIVE VISITORS
# ------------------------------------------------------------------------------

# Identify participants with excessive visit frequency (>3 visits)
excessive_visitors <- HPV_genotyping_data %>%
  dplyr::distinct(study_id, large_blood_date) %>%
  dplyr::count(study_id, name = "visit_frequency") %>%
  dplyr::filter(visit_frequency > 3)

# Remove participants with excessive visits from analysis dataset
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::filter(!study_id %in% excessive_visitors$study_id)


# ------------------------------------------------------------------------------
# COMPREHENSIVE DATA EXPLORATION: UNIQUE VALUE ANALYSIS
# ------------------------------------------------------------------------------

# Define columns for detailed value inspection (excluding ID and date columns)
inspection_columns <- setdiff(names(HPV_genotyping_data),
                              c("study_id", "large_blood_date", "visit_code"))

# Systematic examination of unique values in each column
for (column_name in inspection_columns) {
  cat("\n=== UNIQUE VALUES IN", toupper(column_name), "===\n")
  print(unique(HPV_genotyping_data[[column_name]]))
}


# ------------------------------------------------------------------------------
# TEXT DATA CORRECTION: SPELLING AND FORMATTING STANDARDIZATION
# ------------------------------------------------------------------------------

# Following are changes based on unique (from above)
# Change "Blood contaminatated" to "Blood contaminated":

# Correct spelling error in comment field
HPV_genotyping_data$comment <- stringr::str_replace_all(
  HPV_genotyping_data$comment,
  "Blood contaminatated",
  "Blood contaminated"
)

# Verification of spelling correction
spelling_verification <- HPV_genotyping_data %>%
  dplyr::filter(stringr::str_detect(comment, "Blood contaminatated")) %>%
  dplyr::select(comment)

if (nrow(spelling_verification) == 0) {
  cat("Spelling correction successful: No instances of 'Blood contaminatated' remain")
} else {
  cat("Warning: Instances of 'Blood contaminatated' still present:")
  print(spelling_verification)
}


# ------------------------------------------------------------------------------
# MISSING DATA STANDARDIZATION: CONVERTING STRING NA TO PROPER NULL VALUES
# ------------------------------------------------------------------------------

# Convert string representations of NA to proper missing values
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(., "NA"))) %>%
  dplyr::mutate(dplyr::across(where(is.factor), ~ dplyr::na_if(., "NA")))

# Display data structure for verification
str(HPV_genotyping_data)


# ------------------------------------------------------------------------------
# DATA TYPE CONVERSION: ENSURING PROPER VARIABLE CLASSES
# ------------------------------------------------------------------------------

# Convert variables to appropriate data types for analysis
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(
    large_blood_date = as.Date(large_blood_date),
    dplyr::across(c(nugent_gardnerella, nugent_mobiluncus, nugent_lactobacillus, 
                    nugent_total_score), as.numeric))


################################################################################

# --- STEP 5: CHECK ALL THE COLUMNS ---


# ------------------------------------------------------------------------------
# COMPREHENSIVE DATA VALIDATION: SYSTEMATIC QUALITY CHECKS
# ------------------------------------------------------------------------------

# Validation 1: Study ID Pattern Verification ("study_id")
invalid_study_ids <- HPV_genotyping_data %>%
  dplyr::filter(!stringr::str_detect(study_id, "^\\d{3}-\\d{2}-\\d{4}-\\d{3,4}$"))

print("Invalid Study IDs:")
print(invalid_study_ids$study_id)

# Validation 2: Visit Code Pattern Verification ("visit_code")
invalid_visit_codes <- HPV_genotyping_data %>%
  dplyr::filter(!stringr::str_detect(visit_code, "^\\d{1}-\\d{3}-\\d{1}$"))

print("Invalid Visit Codes:")
print(invalid_visit_codes$visit_code)


# ------------------------------------------------------------------------------
# NUGENT SCORE VALIDATION AND CORRECTION
# ------------------------------------------------------------------------------

# Calculate and correct Nugent total score discrepancies
nugent_discrepancies <- sum(!is.na(HPV_genotyping_data$nugent_total_score) &
                              HPV_genotyping_data$nugent_total_score != 
                              (HPV_genotyping_data$nugent_gardnerella +
                                 HPV_genotyping_data$nugent_mobiluncus +
                                 HPV_genotyping_data$nugent_lactobacillus))

# Apply corrections to Nugent total scores
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(nugent_total_score = dplyr::if_else(
    !is.na(nugent_total_score) &
      nugent_total_score != (nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus),
    nugent_gardnerella + nugent_mobiluncus + nugent_lactobacillus,
    nugent_total_score))

cat("Nugent total score corrections applied to", nugent_discrepancies, "records")


# ------------------------------------------------------------------------------
# NUGENT INTERPRETATION VALIDATION AND CORRECTION
# ------------------------------------------------------------------------------

# Count interpretation mismatches based on total score
interpretation_mismatches <- HPV_genotyping_data %>%
  dplyr::filter(!is.na(nugent_total_score) & !is.na(nugent_interpretation) &
                  (
                    (nugent_total_score < 3 & nugent_interpretation != "NO BV") |
                      (nugent_total_score >= 4 & nugent_total_score <= 6 & nugent_interpretation != "INTERMEDIATE") |
                      (nugent_total_score >= 7 & nugent_interpretation != "BV")
                  )
  ) %>% nrow()

# Correct Nugent interpretations based on total scores
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(nugent_interpretation = dplyr::case_when(
    nugent_total_score < 3 ~ "NO BV",
    nugent_total_score >= 4 & nugent_total_score <= 6 ~ "INTERMEDIATE",
    nugent_total_score >= 7 ~ "BV",
    TRUE ~ nugent_interpretation  
  ))

cat("Nugent interpretation corrections applied to", interpretation_mismatches, "records")


# ------------------------------------------------------------------------------
# STI DATA BINARY CONVERSION: STANDARDIZING DETECTION RESULTS
# ------------------------------------------------------------------------------

# Convert STI detection results to binary format for consistent analysis
# 0 = NOT DETECTED (no infection)
# 1 = DETECTED (infection present)

sti_columns <- grep("^sti_", names(HPV_genotyping_data), value = TRUE)

HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(dplyr::across(all_of(sti_columns),
                              ~ dplyr::case_when(. == "DETECTED" ~ 1,
                                                 . == "NOT DETECTED" ~ 0,
                                                 TRUE ~ NA_real_)))


################################################################################

# --- STEP 6: REPLICATE SAMPLE MANAGEMENT: CONCORDANCE AND DISCORDANCE ANALYSIS ---


# Create replicate dataset for participants with multiple samples
replicate_dataset <- HPV_genotyping_data %>%
  dplyr::group_by(study_id, large_blood_date) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup()

# Identify HPV type columns for concordance analysis
hpv_genotype_columns <- grep("^x\\d+$", names(HPV_genotyping_data), value = TRUE)

# Classify replicate visits as concordant or discordant
replicate_concordance <- replicate_dataset %>%
  dplyr::group_by(study_id, large_blood_date) %>%
  dplyr::summarise(is_concordant = n_distinct(dplyr::across(all_of(hpv_genotype_columns))) == 1, 
                   .groups = "drop")

# Separate concordant and discordant replicate visits
concordant_replicates <- replicate_dataset %>%
  dplyr::semi_join(replicate_concordance %>%
                     dplyr::filter(is_concordant),
                   by = c("study_id", "large_blood_date"))

discordant_replicates <- replicate_dataset %>%
  dplyr::semi_join(replicate_concordance %>%
                     dplyr::filter(!is_concordant),
                   by = c("study_id", "large_blood_date"))

# Define columns for replicate processing
nugent_columns <- grep("nugent_\\D", names(replicate_dataset), value = TRUE)
comprehensive_columns <- c(nugent_columns, sti_columns, hpv_genotype_columns)

# Process concordant replicates: retain row with minimum missing values
processed_concordant_replicates <- concordant_replicates %>%
  dplyr::group_by(study_id, large_blood_date) %>%
  dplyr::mutate(missing_count = rowSums(is.na(dplyr::across(all_of(comprehensive_columns))))) %>%
  dplyr::slice_min(order_by = missing_count, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-missing_count)

# Process discordant replicates: preserve positive HPV results and first non-NA values
processed_discordant_replicates <- discordant_replicates %>%
  dplyr::group_by(study_id, visit_code, large_blood_date) %>%
  dplyr::summarise(
    
    # Process Nugent and STI columns (retain first non-NA value)
    dplyr::across(all_of(c(nugent_columns, sti_columns)), ~ {
      valid_values <- .[!is.na(.)]
      if (length(valid_values) > 0) valid_values[1] else NA
    }, .names = "{.col}"),
    
    # Process HPV columns (retain maximum value, handling NA safely)
    dplyr::across(all_of(hpv_genotype_columns), ~ {
      if (all(is.na(.))) NA_integer_ else max(., na.rm = TRUE)
    }, .names = "{.col}"),
    
    # Process report field (retain "RESULT" if present)
    report = {
      result_values <- report[report == "RESULT"]
      if (length(result_values) > 0) result_values[1] else NA_character_
    },
    
    # Process comment field (retain only NA values)
    comment = {
      if (all(is.na(comment))) NA_character_ else NA_character_
    },
    .groups = "drop"
  )

# Reconstruct dataset by removing replicates and adding processed versions
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::anti_join(replicate_dataset) %>%
  dplyr::bind_rows(processed_concordant_replicates) %>%
  dplyr::bind_rows(processed_discordant_replicates) %>%
  dplyr::arrange(study_id, large_blood_date)


# ------------------------------------------------------------------------------
# CATEGORICAL DATA ENCODING: NUGENT INTERPRETATION STANDARDIZATION
# ------------------------------------------------------------------------------

# Encode Nugent interpretation as numerical categories
# 0 = "NO BV"
# 1 = "BV" 
# 2 = "INTERMEDIATE"

HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(nugent_interpretation = dplyr::case_when(
    is.na(nugent_interpretation) ~ NA_real_,
    nugent_interpretation == "NO BV" ~ 0,
    nugent_interpretation == "INTERMEDIATE" ~ 2,
    nugent_interpretation == "BV" ~ 1))


################################################################################

# --- STEP 7: DERIVED VARIABLE CONSTRUCTION: COMPREHENSIVE HPV ANALYSIS METRICS ---


# 1. Total HPV Infection Count
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(hpv_sum = rowSums(dplyr::select(., starts_with("x")), na.rm = FALSE))

# 2. High-Risk HPV Infection Count
high_risk_hpv_types <- c("x16", "x18", "x26", "x31", "x33", "x35", "x39", "x45", "x51",
                         "x52", "x53", "x56", "x58", "x59", "x66", "x68", "x73", "x82")

HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(hr_hpv_sum = rowSums(dplyr::select(., all_of(high_risk_hpv_types)), na.rm = FALSE))

# 3. Low-Risk HPV Infection Count
low_risk_hpv_types <- c("x6", "x11", "x40", "x42", "x43", "x54", "x61",
                        "x67", "x69", "x70", "x71", "x72", "x84")

HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(lr_hpv_sum = rowSums(dplyr::select(., all_of(low_risk_hpv_types)), na.rm = FALSE))

# 4. HPV Infection Category Classification
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(hpv_group = dplyr::case_when(
    is.na(hpv_sum) ~ NA_character_,
    hpv_sum == 0 ~ "No infection",
    hpv_sum == 1 ~ "Single infection",
    hpv_sum >= 2 ~ "Multiple infection"))

# 5. HPV Positivity Status
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(hpv_positive = dplyr::case_when(
    rowSums(dplyr::select(., all_of(hpv_genotype_columns)) == 1, na.rm = TRUE) > 0 ~ 1,
    rowSums(!is.na(dplyr::select(., all_of(hpv_genotype_columns)))) == 0 ~ NA_real_,
    TRUE ~ 0))

# 6. High-Risk and Low-Risk HPV Positivity Status
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::mutate(
    hr_hpv_positive = dplyr::case_when(
      rowSums(!is.na(dplyr::select(., all_of(high_risk_hpv_types)))) == 0 ~ NA_real_,
      rowSums(dplyr::select(., all_of(high_risk_hpv_types)) == 1, na.rm = TRUE) > 0 ~ 1,
      TRUE ~ 0
    ),
    
    lr_hpv_positive = dplyr::case_when(
      rowSums(!is.na(dplyr::select(., all_of(low_risk_hpv_types)))) == 0 ~ NA_real_,
      rowSums(dplyr::select(., all_of(low_risk_hpv_types)) == 1, na.rm = TRUE) > 0 ~ 1,
      TRUE ~ 0
    ))


# ------------------------------------------------------------------------------
# FINAL DATASET PREPARATION: ANALYSIS-READY DATA STRUCTURE
# ------------------------------------------------------------------------------

# str(HPV_genotyping_data$)

# Create final analysis dataset with selected variables
HPV_genotyping_data <- HPV_genotyping_data %>%
  dplyr::select(
    study_id,
    visit_code,
    large_blood_date,
    nugent_gardnerella,
    nugent_mobiluncus,
    nugent_lactobacillus,
    nugent_total_score,
    nugent_interpretation,
    sti_mycoplasma_genitalium_pcr,
    sti_trichomonas_vaginalis_pcr,
    sti_neisseria_gonorrhoeae_pcr,
    sti_chlamydia_trachomatis_pcr,
    hpv_positive,
    x16, x18, x26, x31, x33, x35, x39, x45, x51, 
    x52, x53, x56, x58, x59, x66, x68, x73, x82,
    hr_hpv_positive,
    hr_hpv_sum,
    x6, x11, x40, x42, x43, x54, x61,
    x67, x69, x70, x71, x72, x84,
    lr_hpv_positive,
    lr_hpv_sum,
    hpv_sum,
    hpv_group,
    report)


# Create complete genotyping dataset (participants with full HPV data)
HPV_genotyping_complete <- HPV_genotyping_data %>%
  dplyr::filter(dplyr::if_all(any_of(hpv_genotype_columns), ~ !is.na(.))) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(visit_num = dplyr::row_number()) %>%
  dplyr::ungroup()

# Create comprehensive dataset (includes all participants regardless of HPV data completeness)
# Note: This dataset should not be used for HPV-specific analyses
HPV_genotyping_data_final_clean <- HPV_genotyping_data %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(visit_num = dplyr::row_number()) %>%
  dplyr::ungroup()


################################################################################

# --- STEP 8: SAVE CSV/TSV FILES ---

# This dataset includes all participants with missing HPV test results and those with complete results  
# write_csv(HPV_genotyping_data, "04_HPV_genotyping_final_data_cleaning.csv")

# This dataset serves as the primary analytical dataset; excludes missing HPV values and represents the actual tested population 
# write_csv(HPV_genotyping_complete, "04_HPV_genotyping_complete.csv")

