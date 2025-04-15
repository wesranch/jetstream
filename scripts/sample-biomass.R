#install.packages("loon.data")
library(loon.data)
library(dplyr)
library(ggplot2)


data("alaska_forest")
head(alaska_forest)

#group
alaska_forest2 <- alaska_forest %>%
  #group_by(year, common) %>%
  #mutate(observations = n(), .groups = "drop")%>%
  #ungroup()%>%
  filter(common %in% c("Alaskan birch", "Black spruce", 
                       "Trembling aspen", "White spruce"))%>%
  filter(!year %in% 1994:1999)


split_dfs <- split(alaska_forest2 , f = alaska_forest2$year )

for (df in split_dfs){
  year <- as.character(unique(df$year))
  df <- df %>%
    select(site, mass, year, common, longitude, latitude)
  write.csv(df, paste0("../thesis/data/input/biomass_", year, ".csv"), row.names = FALSE)
}

