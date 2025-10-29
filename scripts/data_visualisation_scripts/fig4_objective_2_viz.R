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
#       R script: fig4_objective_2_viz                                         #
#       Plot: plot_hpv_status_by_study_end.png                                 #
#       Objective 2: The percentage of women who clear, acquire or persist     #
#                    with individual and grouped (high-risk vs low-risk)       #
#                    HPV types.                                                #
#       Description: Creating Grouped Bar Chart (Vertical)                     #
#                                                                              #  
#       Code Author: "Trevolin Pillay" (29717051)                              #
#       Edited/Review by Dr.TJ Sanko                                           #
#                                                                              #
################################################################################

# --- STEP 1: LOAD NECESSARY LIBRARIES ---

library(readr)     # read_csv(), write_csv() for loading/saving data
library(dplyr)     # mutate(), filter(), arrange(), group_by(), across()
library(tidyr)     # pivot_longer(), complete() if needed
library(lubridate) # year() for date filtering
library(stringr)   # str_extract() to pull percentages from "n/N (%)"
library(xml2)      # read_html() to parse the HTML file
library(rvest)     # html_table() to extract the GT table into a data frame
library(ggplot2)   # plotting (geom_col, coord, scales, theme, ggsave)
library(ggtext)    # Enables rich text formatting in ggplot2
library(tibble)    # Creates modern, tidy data frames


################################################################################

# --- STEP 2: READ CSV FILES ---

# 01_merged_FRESH_HPV_data.csv
# df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/2_Data Cleaning/merged_tables_csv/01_merged_FRESH_HPV_data.csv")

# 04_HPV_genotyping_complete.csv
df <- read_csv("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/2_Data Cleaning/data_cleaning_csv/04_HPV_genotyping_complete.csv")

# set the new working directory
setwd("C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/4_Visualization and Interpretation/data_visualisation_plots")
getwd()


################################################################################

# --- STEP 3: FILTER THE "large_blood_date" (2013-2022)  ---

# Filter and sort by study ID and large blood date
df <- df %>%
  filter(year(large_blood_date) <= 2022) %>%
  arrange(study_id, large_blood_date)


################################################################################

# --- STEP 4: DATA VISUALISATION AND PLOT 3. ---

# ============================================================================ #
#   Objective 2: The percentage of women who clear, acquire or persist with    #
#                individual and grouped (high-risk vs low-risk) HPV types.     #
#   Based on "Table 3. HPV Outcomes at Follow-up"                              #
#   Figure 4: "HPV status by study end"                                        #
# ============================================================================ #


# ------------------------------------------------------------------------------
# 1) Read the HTML table
# ------------------------------------------------------------------------------
html_path <- "C:/Users/Trevolin Pillay/Desktop/BScHons Research Project/Project (721)/Final Results/HPV_Project_Objectives/5_Final Results/tables/Table 3. HPV Outcomes at Follow-up.html"

# Parse HTML tables
tables <- rvest::html_table(xml2::read_html(html_path), fill = TRUE)
df_raw <- tables[[1]]


# ------------------------------------------------------------------------------
# 2) Clean structure: header rows and data rows
# ------------------------------------------------------------------------------

# Row 2 contains column headers & rows 3+ contain data
colnames(df_raw) <- c("row", "Incidence", "Clearance", "Persistence")

# Remove empty or NA rows
df_raw <- df_raw %>%
  filter(!is.na(row) & row != "")


# ------------------------------------------------------------------------------
# 3) Extract numeric percentages
# ------------------------------------------------------------------------------
extract_pct <- function(x) {
  x <- stringr::str_extract(x, "[0-9.]+(?=%)")
  as.numeric(x)
}

df_clean <- df_raw %>%
  mutate(across(c(Incidence, Clearance, Persistence), extract_pct))


# ------------------------------------------------------------------------------
# 4) Reshape to long format 
# ------------------------------------------------------------------------------
df_long <- df_clean %>%
  pivot_longer(
    cols = c(Incidence, Clearance, Persistence),
    names_to = "Metric", values_to = "Percent"
  ) %>%
  mutate(
    # Keep the order as in the HTML table
    row    = factor(row, levels = df_clean$row),
    Metric = factor(Metric, levels = c("Incidence","Clearance","Persistence"))
  ) %>%
  filter(!is.na(Percent))


df_plot <- df_long %>%
  group_by(row) %>%
  mutate(
    # Order metrics by descending percent *within each row*
    metric_order = rank(-Percent, ties.method = "first")  # numeric rank
  ) %>%
  ungroup()


# ------------------------------------------------------------------------------
# 5) Plot grouped bar chart
# ------------------------------------------------------------------------------
df_plot <- df_plot %>%
  mutate(Metric = factor(Metric, levels = c("Clearance", "Persistence", "Incidence")))

label_map <- c(
  "Bivalent types (HPV16, HPV18)" =
    "Bivalent&nbsp; types",
  "Quadrivalent types (HPV16, HPV18, HPV06, HPV11)" =
    "Quadrivalent&nbsp; types",
  "Nonavalent types (HPV06, HPV11, HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58)" =
    "Nonavalent&nbsp; types"
)

# Data used to create the vaccine legend
vax_legend_df <- tibble(
  vaccine = factor(
    c(
      "Bivalent types = (HPV16, HPV18)",
      "Quadrivalent types = (HPV16, HPV18, HPV06, HPV11)",
      "Nonavalent types = (HPV06, HPV11, HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58)"
    ),
    levels = c(
      "Bivalent types = (HPV16, HPV18)",
      "Quadrivalent types = (HPV16, HPV18, HPV06, HPV11)",
      "Nonavalent types = (HPV06, HPV11, HPV16, HPV18, HPV31, HPV33, HPV45, HPV52, HPV58)"
    )
  )
)


ggplot(
  df_plot,
  aes(x = row, y = Percent, fill = Metric, group = interaction(row, metric_order))
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percent)),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 2.6
  ) +
  geom_point(
    data = vax_legend_df,
    aes(x = Inf, y = Inf, shape = vaccine),
    inherit.aes = FALSE,
    alpha = 0, size = 3, show.legend = TRUE
  ) +
  labs(
    title = "HPV status by study end",
    x = NULL, y = "Percentage (%)",
    fill = NULL,   
    shape = NULL   
  ) +
  scale_x_discrete(labels = label_map) +
  scale_fill_manual(
    values = c("Clearance"="#009E73","Persistence"="#CC79A7","Incidence"="#E69F00"),
    breaks = c("Clearance","Persistence","Incidence"),
    drop = FALSE
  ) +
  guides(
    fill  = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)  
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .08))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = ggtext::element_markdown(angle = 35, hjust = 1),
    legend.position = "top",
    legend.box = "vertical",          
    legend.direction = "horizontal",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.margin = margin(t = -5, b = 10)
  )


# ------------------------------------------------------------------------------
# 6) Save PNG & CSV file
# ------------------------------------------------------------------------------
# ggsave("HPV status by study end.png", width = 12, height = 6, dpi = 300)

# Save the clean data
# write_csv(df_long, "Metrics_parsed_percentages.csv")

