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
#       R script: 03_FRESH_HPV_sexual_behavior_data_cleaning                   #
#       CSV: 03_FRESH_HPV_data_cleaning.csv                                    # 
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

# --- STEP 2: READ & LOAD THE EXCEL "03. FRESH_HPV_sexual_behavior_metadata_04.06.2025" --- 

# excel_sheets("03. FRESH_HPV_sexual_behavior_metadata_04.06.2025.xlsx") # List all sheets
FRESH_HPV_data_03 <- read_excel("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/1_HPV Data/HPV_raw_data/03. FRESH_HPV_sexual_behavior_metadata_04.06.2025.xlsx", sheet = 2)  # Or use sheet name
view(FRESH_HPV_data_03)


################################################################################

# --- STEP 3: COLUMN NAME STANDARDIZATION AND DATA CLEANING --- 

# Replace the old name with new name with underscore and lowercase
FRESH_HPV_data_03 <- FRESH_HPV_data_03 %>%
  rename(
    study_id = subjid,
    visit_id = visitid,
    vst_dt = vstdt,
    num_sex_episodes_7_day = numsexepisodes7day,
    num_sex_partners_7_days = numsexpartners7days,
    num_sex_episodes_30_day = numsexepisodes30day,
    condoms_used_30_days = condomsused30days,
    num_sex_partners_30_day = numsexpartners30day,
    family_planning_type = familyplanningtype,
    sti_symptoms_now = stisymptomsnow
  )

# View the renamed dataset
view(FRESH_HPV_data_03)

# Clean the names:
clean_FRESH_HPV_data <- FRESH_HPV_data_03 %>% clean_names()

# View the cleaned dataset
view(clean_FRESH_HPV_data)


################################################################################

# --- STEP 4: DATA STRUCTURE EXPLORATION AND ASSESSMENT --- 

# Create working copy for structure analysis
check_FRESH_HPV_data_1 <- clean_FRESH_HPV_data

# a) Data overview using glimpse (cleaner alternative to str())
glimpse(check_FRESH_HPV_data_1)
# Dataset dimensions: 902 rows × 21 columns

# b) Comprehensive data summary
summary(check_FRESH_HPV_data_1)

# Data structure summary:
# All variables are currently stored as character type
# Variables include: study_id, visit_id, vst_dt, visit_code, 
# visit_date_hiv_risk_3month, sexual behavior variables, 
# demographic variables, and partner information


################################################################################

# --- STEP 5: DUPLICATE RECORD IDENTIFICATION AND REMOVAL --- 

# Check for duplicates
duplicate_check <- duplicated(check_FRESH_HPV_data_1)
duplicate_count <- sum(duplicate_check)

cat("Duplicate records found:", duplicate_count, "\n")

# Remove duplicate records while preserving first occurrence
remove_duplicate_FRESH_HPV_data <- check_FRESH_HPV_data_1 %>% 
  distinct()


################################################################################

# --- STEP 6: MISSING DATA ASSESSMENT AND HANDLING --- 

# a) Comprehensive missing value analysis
missing_values_summary <- colSums(is.na(remove_duplicate_FRESH_HPV_data))

cat("Missing values per column:\n")
print(missing_values_summary)

# All columns show 0 missing values, indicating no true NAs present

# Create working copy for missing data processing
df_FRESH <- remove_duplicate_FRESH_HPV_data


# col: "sti_symptoms_now"
# b) Handle string representations of missing values
# Define common string representations of missing data
unique(df_FRESH$num_sex_partners_7_days)
class(df_FRESH$num_sex_partners_7_days)

fake_na_values <- c("", "NA")

# Convert string NAs to proper R NA values across specified columns
coerced_FRESH_HPV_data <- remove_duplicate_FRESH_HPV_data %>%
  dplyr::mutate(dplyr::across(
    c(num_sex_episodes_7_day, num_sex_partners_7_days, 
      num_sex_episodes_30_day, condoms_used_30_days, num_sex_partners_30_day,
      oral_sex_last_7_days_yn, penile_vaginal_last_7_days_yn, anal_last_7_days_yn,
      family_planning_type, sti_symptoms_now,
      age_baseline, age_first_sex_baseline_demo,
      number_sexual_partners, live_with_partner_baseline),
    ~ Reduce(function(x, val) dplyr::na_if(x, val), fake_na_values, init = .)
  ))

# Verification of NA conversion
cat("NA conversion verification:\n")
cat("Unique values in num_sex_partners_7_days after conversion:\n")
print(unique(coerced_FRESH_HPV_data$num_sex_partners_7_days))

# View the processed dataset
View(coerced_FRESH_HPV_data)

# Check missing values after string-to-NA conversion
missing_after_conversion <- colSums(is.na(coerced_FRESH_HPV_data))
cat("Missing values after string-to-NA conversion:\n")
print(missing_after_conversion)


################################################################################

# --- STEP 7: DATA TYPE CONVERSION AND STANDARDIZATION ---

# Load required library for date handling
library(lubridate)

# Check current data types
current_data_types <- sapply(coerced_FRESH_HPV_data, class)
cat("Current data types:\n")
print(current_data_types)

# Convert data types to appropriate formats
convert_dt_FRESH_HPV <- coerced_FRESH_HPV_data %>%
  dplyr::mutate(
    # Date conversions (format: "24/07/2025")
    vst_dt = as.Date(vst_dt),
    visit_date_hiv_risk_3month = as.Date(visit_date_hiv_risk_3month),
    date_lower = as.Date(date_lower),
    date_upper = as.Date(date_upper),
    
    # Numeric conversions for sexual behavior variables
    num_sex_episodes_7_day = as.numeric(num_sex_episodes_7_day),
    num_sex_partners_7_days = as.numeric(num_sex_partners_7_days),
    num_sex_episodes_30_day = as.numeric(num_sex_episodes_30_day),
    condoms_used_30_days = as.numeric(condoms_used_30_days),
    num_sex_partners_30_day = as.numeric(num_sex_partners_30_day),
    
    # Numeric conversions for sexual activity indicators
    oral_sex_last_7_days_yn = as.numeric(oral_sex_last_7_days_yn),
    penile_vaginal_last_7_days_yn = as.numeric(penile_vaginal_last_7_days_yn),
    anal_last_7_days_yn = as.numeric(anal_last_7_days_yn),
    
    # Numeric conversions for family planning and symptoms
    family_planning_type = as.numeric(family_planning_type),
    sti_symptoms_now = as.numeric(sti_symptoms_now),
    
    # Numeric conversions for demographic variables
    age_baseline = as.numeric(age_baseline),
    age_first_sex_baseline_demo = as.numeric(age_first_sex_baseline_demo),
    number_sexual_partners = as.numeric(number_sexual_partners),
    live_with_partner_baseline = as.numeric(live_with_partner_baseline)
  )


################################################################################

# --- STEP 8: FINAL DATA VALIDATION AND QUALITY ASSURANCE ---

# Create final validation dataset
check_FRESH_HPV_data_2 <- convert_dt_FRESH_HPV

# a) Final data structure overview
glimpse(check_FRESH_HPV_data_2)
# Final dataset dimensions: 902 rows × 21 columns

# b) Comprehensive summary of cleaned data
summary(check_FRESH_HPV_data_2)

# c) Final quality checks
# Check for missing values after type conversion
final_missing_values <- colSums(is.na(check_FRESH_HPV_data_2))
cat("Final missing values summary:\n")
print(final_missing_values)

# Verify no duplicates remain
final_duplicate_check <- duplicated(check_FRESH_HPV_data_2)
final_duplicate_count <- sum(final_duplicate_check)
duplicate_table <- table(final_duplicate_check)

cat("Final duplicate check:\n")
cat("Duplicates found:", final_duplicate_count, "\n")
print(duplicate_table)  # FALSE = 902, TRUE = 0


# ------------------------------------------------------------------------------
# Final output: Cleaned dataset
# ------------------------------------------------------------------------------

# Create final cleaned dataset
final_FRESH_HPV_data <- check_FRESH_HPV_data_2

# View the final cleaned dataset
View(final_FRESH_HPV_data)

# Data cleaning summary
cat("\n=== DATA CLEANING SUMMARY ===\n")
cat("Original dataset: 902 rows × 21 columns\n")
cat("Final dataset: 902 rows × 21 columns\n")
cat("Duplicates removed: 0\n")
cat("Missing values handled: String NAs converted to proper R NAs\n")
cat("Data types standardized: Dates and numeric variables properly formatted\n")
cat("Column names standardized: Lowercase with underscores\n")
cat("Dataset ready for analysis: TRUE\n")


################################################################################

# --- STEP 9: SAVE CSV/TSV FILES ---
library(writexl) # For writing excel file

# filepath <- "C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/2_Data Cleaning/data_cleaning_csv/03_FRESH_HPV_data_cleaning.csv"
write.csv(final_FRESH_HPV_data, filepath, row.names = FALSE)

