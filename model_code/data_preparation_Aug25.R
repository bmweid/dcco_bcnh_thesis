# Data Preparation Script

# Load required library
library(tidyverse)

# Function to safely read CSV files
safe_read_csv <- function(file_path) {
  tryCatch(
    read_csv(file_path, col_types = cols(.default = "c")),
    error = function(e) {
      warning(paste("Error reading file:", file_path, "\nError:", e$message))
      return(NULL)
    }
  )
}

# Import data
dataset_paths <- c(
  nest_density = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/pivoted_morans_i_results.csv",
  colony_growth_index = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/pivoted_growth_index.csv",
  road_proximity = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/summed_road_proximity_bcnh.csv",
  bcnh_nest_success = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/bcnh_nest_success.csv",
  dcco_usurpation = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/dcco_usurpation.csv",
  dcco_management = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/dcco_management.csv"
)

datasets <- lapply(dataset_paths, safe_read_csv)

# Remove NULL entries (failed imports)
datasets <- datasets[!sapply(datasets, is.null)]

# Print information about imported datasets
cat("Imported datasets:\n")
for (name in names(datasets)) {
  df <- datasets[[name]]
  cat(sprintf("%s: %d rows, %d columns\n", name, nrow(df), ncol(df)))
  cat("Columns:", paste(names(df), collapse = ", "), "\n\n")
}

# Combine data & select specific columns
combined_data <- datasets$nest_density %>%
  select(year, bcnh_nest_density = morans_i_BCNH, dcco_nest_density = morans_i_DCCO) %>%
  mutate(across(c(bcnh_nest_density, dcco_nest_density), as.numeric),
         total_nest_density = bcnh_nest_density + dcco_nest_density) %>%
  left_join(
    datasets$colony_growth_index %>%
      select(year, growth_index_BCNH, growth_index_DCCO) %>%
      mutate(across(c(growth_index_BCNH, growth_index_DCCO), as.numeric)),
    by = "year"
  ) %>%
  left_join(
    datasets$road_proximity %>%
      select(year, bcnh_road_proximity) %>%
      mutate(bcnh_road_proximity = as.numeric(bcnh_road_proximity)),
    by = "year"
  ) %>%
  left_join(
    datasets$bcnh_nest_success %>%
      select(year, bcnh_nest_success) %>%
      mutate(bcnh_nest_success = as.numeric(bcnh_nest_success)),
    by = "year"
  ) %>%
  left_join(
    datasets$dcco_usurpation %>%
      select(year, dcco_usurpation) %>%
      mutate(dcco_usurpation = as.numeric(dcco_usurpation)),
    by = "year"
  ) %>%
  left_join(
    datasets$dcco_management %>%
      select(year, winter_nestremoval, `deterrence_ activenestremoval`, deterrence_night) %>%
      rename(deterrence_activenestremoval = `deterrence_ activenestremoval`) %>%
      mutate(across(c(winter_nestremoval, deterrence_activenestremoval, deterrence_night), as.numeric)),
    by = "year"
  )

# Check the structure of the combined dataset
str(combined_data)
summary(combined_data)

# Print the first few rows of the combined dataset
print(head(combined_data))

# Save the combined dataset as a CSV file
output_path <- "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/current working files/model datasets/combined_dataset.csv"
write_csv(combined_data, output_path)
cat(sprintf("Combined dataset saved to: %s\n", output_path))