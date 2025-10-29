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
#       R script: fig2_objective_1_viz                                         #
#       Objective 1: What is the baseline prevalence of each HPV type?         #   
#       Plots: plot_baseline_genotype_distribution.png                         #
#       Description: Creating Bar chart (Horizontal) for each HPV type         #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)      # read_csv() to load the CSV
library(dplyr)      # filter(), mutate(), summarise(), across(), if_else(), case_when()
library(tidyr)      # pivot_longer() to reshape for plotting
library(forcats)    # fct_inorder() to lock bar order by prevalence
library(lubridate)  # year() for date filtering
library(ggplot2)    # plotting (geom_col, scales, theme, ggsave)
library(ggpattern)  # patterned bars & pattern legends used in geom_col_pattern()


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

# --- STEP 4: DATA VISUALISATION AND PLOT 1. ---

# ============================================================================ #
#   Objective 1: What is the baseline prevalence of each HPV type?             #
#   Figure 2: "HPV Genotype Distribution at Baseline"                          #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Define genotype sets 
# ------------------------------------------------------------------------------
hr_set <- c("x16","x18","x26","x31","x33","x35","x39","x45","x51",
            "x52","x53","x56","x58","x59","x66","x68","x73","x82")

lr_set <- c("x6","x11","x40","x42","x43","x54","x61",
            "x67","x69","x70","x71","x72","x84")

# Vaccine-type sets per your instruction
hr_vax <- c("x16","x18","x31","x45","x52","x58")
lr_vax <- c("x6","x11")

all_types <- c(hr_set, lr_set)


# ------------------------------------------------------------------------------
# 2) Filter for baseline (visit_num == 1)
# ------------------------------------------------------------------------------
df_base <- df %>% 
  filter(visit_num == 1)


# ------------------------------------------------------------------------------
# 2) Compute % positive for each genotype (0/1 columns)
# ------------------------------------------------------------------------------
# Columns x16, x18, ... are coded 0/1 (or NA).
sum_tbl <- df_base |>
  summarise(across(all_of(all_types),
                   ~ mean(.x == 1, na.rm = TRUE) * 100)) |>
  pivot_longer(everything(),
               names_to = "type_var", values_to = "percent")


# Add labels (HPV16, HPV18, ...) and grouping (HR/LR)
sum_tbl <- sum_tbl |>
  mutate(
    group = if_else(type_var %in% hr_set, "HR-HPV", "LR-HPV"),
    vaccine_class = case_when(
      type_var %in% hr_vax ~ "HR-HPV vaccine type",
      type_var %in% lr_vax ~ "LR-HPV vaccine type",
      TRUE                 ~ "Non-vaccine type"
    ),
    type_label = gsub("^x", "HPV", type_var),
    type_label = dplyr::if_else(type_label == "HPV6", "HPV06", type_label),
    # Replace NaN (all-NA columns) with 0
    percent = ifelse(is.nan(percent), 0, percent)
  ) 

# Order types by percent (descending), keeping HR then LR together
sum_tbl <- sum_tbl |>
  arrange(desc(percent), group) |>
  mutate(type_label = fct_inorder(type_label))


# ------------------------------------------------------------------------------
# 3) Pattern & color settings
# ------------------------------------------------------------------------------
hr_red    <- "#D62728"  # bar fill
lr_blue   <- "#1F77B4"  # bar fill
hr_red_dk <- "#D62728"  # pattern stroke (HR)
lr_blue_dk<- "#1F77B4"  # pattern stroke (LR)

# targeted vaccine genotypes get patterns
patterned_types <- c("x6","x11","x16","x18","x31","x33","x45","x52","x58")

sum_tbl <- sum_tbl |>
  mutate(
    vaccine_class = dplyr::case_when(
      type_var %in% patterned_types & group == "HR-HPV" ~ "HR-HPV vaccine type",
      type_var %in% patterned_types & group == "LR-HPV" ~ "LR-HPV vaccine type",
      TRUE ~ "Non-vaccine type"
    ),
    pattern_key = dplyr::case_when(
      vaccine_class == "HR-HPV vaccine type" ~ "stripe",
      vaccine_class == "LR-HPV vaccine type" ~ "circle",
      TRUE                                   ~ "none"
    ),
    # pattern stroke color (no white/NA)
    pat_col = dplyr::case_when(
      group == "HR-HPV" ~ hr_red_dk,
      TRUE              ~ lr_blue_dk
    )
  )


# ------------------------------------------------------------------------------
# 4) Bar chart (Horizontal, baseline only)
# ------------------------------------------------------------------------------
hpv_plot <- ggplot(sum_tbl, aes(x = percent, y = type_label, fill = group)) +
  ggpattern::geom_col_pattern(
    aes(
      pattern = pattern_key,       # "stripe","circle","none"
      pattern_colour = pat_col     # darker shade of the same color
    ),
    width  = 0.8,
    colour = "black",
    linewidth = 0.3,
    pattern_alpha   = 1,
    pattern_density = 0.6,
    pattern_spacing = 0.035,
    pattern_size    = 0.28,
    pattern_key_scale_factor = 0.9
  ) +
  # HR/LR color legend
  scale_fill_manual(
    values = c("HR-HPV" = hr_red, "LR-HPV" = lr_blue),
    breaks = c("HR-HPV", "LR-HPV"),
    labels = c("HR-HPV", "LR-HPV"),
    name   = NULL
  ) +
  # Vaccine-type legend
  ggpattern::scale_pattern_manual(
    values = c("none"="none", "stripe"="stripe", "circle"="circle"),
    breaks = c("stripe","circle"),
    labels = c("HR-HPV vaccine type","LR-HPV vaccine type"),
    name   = NULL
  ) +
  ggpattern::scale_pattern_colour_identity(guide = "none") +
  ggpattern::scale_pattern_fill_identity(guide = "none") +
  
  geom_text(aes(label = sprintf("%.1f%%", percent)), hjust = -0.1, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "HPV Genotype Distribution at Baseline",
    x = "Prevalence (%)", 
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.direction= "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank()
  ) +
  guides(
    # Color legend: solid red/blue
    fill = guide_legend(
      order = 1, nrow = 1, byrow = TRUE,
      override.aes = list(
        fill           = c("#D62728", "#1F77B4"),  
        pattern        = "none",                   
        pattern_colour = NA,                       
        pattern_fill   = NA,
        colour         = "black"
      )
    ),
    # Pattern legend: show red stripes + blue dots
    pattern = guide_legend(
      order = 2, nrow = 1, byrow = TRUE,
      override.aes = list(
        fill           = c("#D62728", "#1F77B4"),   # pure red/blue boxes
        pattern_colour = c("#D62728", "#1F77B4")),  # colored stripes/dots
      colour         = "black"
    )
  )

print(hpv_plot)


# ------------------------------------------------------------------------------
# 5) Save PNG file
# ------------------------------------------------------------------------------
# ggplot2::ggsave("HPV Genotype Distribution at Baseline.png", plot = hpv_plot, width = 12, height = 7, dpi = 300)

