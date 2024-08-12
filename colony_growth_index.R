library (tidyverse)
library (dplyr)
library (snakecase)
library (tidyr)

#import & clean data
#dcco
data_dcco <- read.csv("C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/data/post-clean/merged_filtered_dcco_all1.csv")


# long pivots data so that each year per tag will have it's own row with all data connected to the tag available.
#data_dcco_pivot <- data_dcco |>
#  pivot_longer(
#    cols = x_1992:x_2023,
#    names_to = "year",
#    values_to = "count"
#  )

# Save cleaned data as CSV to post-clean data folder
#write.csv(data_dcco_pivot, file = "C:\\Users\\baill\\OneDrive\\r-projects\\cormorant-colony-analysis\\data\\post-clean\\data_dcco_cleaned.csv", 
         # row.names = FALSE)

#prep dcco data
data_dcco_sum <- data_dcco |>
  group_by(year, peninsula)|>
  summarize(
    total_count = sum(nest_count, na.rm= TRUE),
  )
  

#import & clean data
#bcnh
#data_bcnh <- read.csv("C:\\Users\\baill\\OneDrive\\r-projects\\cormorant-colony-analysis\\data\\raw\\data_set_raw_bcnh.csv")
#janitor::clean_names(data_bcnh, )
#colnames(data_bcnh)[2] <- "coord_x"
#colnames(data_bcnh) <- to_snake_case(colnames(data_dcco))

# long pivots data so that each year per tag will have it's own row with all data connected to the tag available.
#data_bcnh_pivot <- data_bcnh |>
#  pivot_longer(
#   cols = x_1992:x_2023,
#    names_to = "year",
#    values_to = "count"
# )

# Save cleaned data as CSV to post-clean data folder
#write.csv(data_bcnh_pivot, file = "C:\\Users\\baill\\OneDrive\\r-projects\\cormorant-colony-analysis\\data\\post-clean\\data_bcnh_cleaned.csv", 
#          row.names = FALSE)

data_bcnh_sum <- data_bcnh_cleaned |>
  group_by(year) |>
  summarize(
    total_count = sum(count, na.rm= TRUE)
  )

#colony growth index
#>gi <- population colony i year t+1- population year t / max (larger of two values) (population colony i year t+1- population year t)
#needs to be calculated for each year; add as column in table? maybe doesn't need to be column...
#using total counts each year

colony_growth_index_dcco <- data_dcco_sum |>
  mutate(growth_index = (lag(total_count) + total_count) / pmax(lag(total_count), 
                                                                total_count, na.rm = TRUE))

colony_growth_index_bcnh <- data_bcnh_sum |>
  mutate(growth_index = (lag(total_count) + total_count) / pmax(lag(total_count), 
                                                                total_count, na.rm = TRUE))

#nest success
janitor::clean_names(bcnh_nest_success_cleaned, )
colnames(bcnh_nest_success_cleaned) <- to_snake_case(colnames(bcnh_nest_success_cleaned))
bcnh_nest_success <- bcnh_nest_success_cleaned$nest_productivity


#prep site data- divide by peninsula?

#prep survey data- 