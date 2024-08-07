library(dplyr)

merged_df <- read.csv ("cormorant-colony-analysis/data/raw/exportmergedfromarcgis.csv")
class(merged_df)
help("group_by")
# Merging the data
merged_df <- merged_df |>
  dplyr::group_by(treeid, year) |>
  dplyr::summarise(
    xcoord = first(xcoord),
    ycoord = first(ycoord),
    tag = first(tag),
    peninsula = first(peninsula),
    northing = first(northing),
    easting = first(easting),
    bcnh_count = sum(bcnh_count, na.rm = TRUE),
    dcco_count = sum(dcco_count, na.rm = TRUE)
  )

# Display the merged data
print(merged_df)

write.csv(merged_df, file = "C:/Users/baill/OneDrive/r-projects/cormorant-colony-analysis/data/post-clean/arcgismergedcleanedbcnhdcco.csv", row.names = FALSE)