# Setup
library(tidyverse)
library(FactoMineR)
library(vegan)
library(maps)
setwd("C:/Users/David/Desktop/STAT 3280/Project")

# Import data and perform PCA
# Code from: https://github.com/DEST-bio/data-paper/blob/main/Figures7_and_9/code/PCA_InbreedvsPools_addtoFolder.R
# Data adapted from: https://github.com/DEST-bio/data-paper/blob/main/Figures7_and_9/data/DEST_DGN_AllSNPs_Metadata.Rdata
load("DEST_DGN_AllSNPs_Metadata.Rdata")
DEST_DGN_metadata$country = gsub("USA","United States", DEST_DGN_metadata$country)

dat_filt_maf_LD500_naimp[which(DEST_DGN_metadata$continent == "Europe" & 
                                 DEST_DGN_metadata$Continental_clusters == "Europe_W" &
                                 !DEST_DGN_metadata$country %in% c("Egypt", "Cyprus", "Turkey", "Finland", "United Kingdom", "Sweden", "Denmark")),] %>%  
  PCA(scale.unit = F, graph = F, ncp = 50) -> LD500_naimp_PCA_50PCs_object

LD500_naimp_PCA_50PCs_object$ind$coord %>% 
  as.data.frame() %>% 
  mutate(sampleId = rownames(.)) %>% 
  left_join(., DEST_DGN_metadata) -> PCA_coords_metadata

# Filter data for European countries and collapse data into one point per country
PC_EU <- PCA_coords_metadata %>%
  filter(continent == "Europe" & Continental_clusters != "North_America") %>%
  group_by(country) %>%
  summarize(lat = mean(lat),
            long = mean(long),
            Dim.1 = mean(Dim.1),
            Dim.2 = mean(Dim.2))

# Perform Procrustes analysis to compare PC coordinates with geographic coordinates
pro <- procrustes(X = PC_EU %>% select(long, lat),
                  Y = PC_EU %>% select(Dim.2, Dim.1))
PC_EU <- cbind(PC_EU, pro$Yrot) %>%
  select(country, 
         sample_lat=lat, sample_long=long,
         proc_lat="1", proc_long="2") %>%
  mutate(proc_lat = proc_lat + mean(sample_lat),
         proc_long = proc_long + mean(sample_long)) %>%
  mutate(sample_long_text_adj = case_when(
           country == "Germany" ~ -0.5,
           country == "France" ~ -1.2,
           country == "Switzerland" ~ 0.2,
           TRUE ~ 0),
         sample_lat_text_adj = case_when(
           country == "Spain" ~ -0.5,
           country == "Germany" ~ -0.7,
           country == "Portugal" ~ 0.8,
           country == "Italy" ~ -0.8,
           country == "Austria" ~ -0.8,
           country == "Netherlands" ~ 0.8,
           country == "Switzerland" ~ -0.8,
           country == "France" ~ -0.2,
           TRUE ~ 0),
         pc_long_text_adj = case_when(
           country == "Germany" ~ -2.7,
           country == "Switzerland" ~ -2.7,
           country == "Austria" ~ 2.2,
           TRUE ~ 0),
         pc_lat_text_adj = case_when(
           country == "Spain" ~ 0.8,
           country == "Germany" ~ -0.4,
           country == "Portugal" ~ 0.8,
           country == "Italy" ~ -0.8,
           country == "France" ~ -0.8,
           country == "Austria" ~ -0.6,
           country == "Netherlands" ~ -0.5,
           country == "Switzerland" ~ 0.4,
           TRUE ~ 0),)

map_eu <- map_data("world") %>%
  mutate(hide = region %in% c("UK", "Ireland", "Algeria", "Tunisia", "Morocco", "Isle of Man"))

pc_col <- "#C21010"
land_col <- "#CFE8A9"
sea_col <- "#FFFDE3"

# Plot Procrustes analysis results with map in background
ggplot(PC_EU) + 
  geom_polygon(map_eu,
               mapping=aes(x=long, y=lat, group=group, 
                           fill=hide, color=hide)) +
  scale_fill_manual(values=c(land_col, sea_col)) +
  scale_color_manual(values=c("black", sea_col)) +
  geom_point(aes(x=sample_long, y=sample_lat), size=4) +
  geom_text(aes(x=sample_long + sample_long_text_adj, 
                y=sample_lat + sample_lat_text_adj,
                label=country),
            size=6, fontface="bold") +
  geom_point(aes(x=proc_long, y=proc_lat), size=4, color=pc_col) +
  geom_text(aes(x=proc_long + pc_long_text_adj, 
                y=proc_lat + pc_lat_text_adj, 
                label=country),
            color=pc_col, size=6, fontface="bold") +
  geom_curve(aes(x=sample_long, y=sample_lat,
                 xend=proc_long, yend=proc_lat),
             arrow=arrow(length=unit(0.3, "cm")),
             curvature=0.2, size=1) +
  geom_text(aes(label="Geography of samples", x=-9.5, y=54), 
            size=7, hjust=0) +
  geom_text(aes(label="Geography reconstructed\nfrom genetics", x=-9.5, y=51.9), 
            color=pc_col, size=7, hjust=0) +
  theme_classic() +
  theme(panel.background = element_rect(fill=sea_col),
        legend.position = "none") +
  coord_fixed(xlim=c(-9, 15), ylim=c(36, 54)) + 
  xlab("Longitude / PC 2") +
  ylab("Latitude / PC 1") +
  ggtitle(expression(paste("Reconstructing",
                           italic(" D. melanogaster "),
                           "geography from genetics via PCA and Procrustes analysis")))

