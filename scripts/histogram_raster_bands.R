#histogram of RS inputs

library(terra)
library(tidyverse)
library(ggridges)
library(dplyr)
library(RColorBrewer)

# path to data
df <- read.csv("/Users/wancher/Documents/thesis/data/output/pixel-vals-indices-topo.csv")
df <- df %>% select(-common, -site, -mass)
# df base r
par(mfrow = c(6, 6), mar = c(3, 3, 2, 1))

# loop over each var
for(i in 1:ncol(df)) {
  plot(1:nrow(df), df[[i]], 
       main = colnames(df)[i], 
       xlab = "Index", 
       ylab = colnames(df)[i], 
       col = "red", 
       pch = 16, 
       cex = 0.7)
}

par(mfrow = c(1, 1)) 

#example file
r <- rast(path)
band_names <- names(r)
normalize_layer <- function(layer) {
  min_val <- min(layer[], na.rm = TRUE)
  max_val <- max(layer[], na.rm = TRUE)
  normalized <- (layer - min_val) / (max_val - min_val)
  return(normalized)
}

normalized_rasters <- lapply(1:nlyr(r), function(i) {
  band <- r[[i]]
  normalized_band <- normalize_layer(band)
  return(normalized_band)
})

# Convert the list of normalized and reduced layers into a SpatRaster object
normalized_r <- rast(normalized_rasters)

#convert to df
df <- as.data.frame(normalized_r, xy = FALSE, na.rm = TRUE)
df_long <- df %>%
  pivot_longer(cols = everything(), 
               names_to = "Band", 
               values_to = "Value")
df_long_sampled <- df_long %>%
  sample_n(500)

#ridge plot
display.brewer.all()
colors <- colorRampPalette(brewer.pal(11, "Paired"))(57)  # 57 colors

plt <- ggplot(df_long_sampled, aes(x = Value, y = Band, fill = Band)) +
  geom_density_ridges_gradient(scale = 10, rel_min_height = 0.002, linewidth = 0.35) +
  scale_fill_manual(values = colors)+
  #scale_fill_viridis_d(option = "plasma") +  
  theme(
    axis.title = element_blank(),  
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.ticks = element_blank()   
  )
plt #this one for rs
ggsave(
  filename = "C:/Users/wesra/Downloads/Histogram_RS_Inputs_2015.jpeg",
  plot = plt,
  width = 10,
  height =10,
  dpi = 600
)



colors <- colorRampPalette(brewer.pal(8, "Accent"))(17) #landis

plt2 <- ggplot(df_long_sampled, aes(x = Value, y = Band, fill = Band)) +
  geom_density_ridges_gradient(scale = 10, rel_min_height = 0.002, linewidth = 0.35) +
  scale_fill_manual(values = colors)+
  #scale_fill_viridis_d(option = "plasma") +  
  theme(
    axis.title = element_blank(),  
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.ticks = element_blank()   
  )
plt2
ggsave(
  filename = "C:/Users/wesra/Downloads/Histogram_Landis_Input_Maps.png",
  plot = plt2,
  width = 8,
  height = 6,
  dpi = 600
)
