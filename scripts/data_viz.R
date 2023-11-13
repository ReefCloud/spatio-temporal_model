# Codes associated with visualizations of model inputs 
# Author: Julie Vercelloni
# Contact: j.vercelloni@aims.gov.au
# Last update: November 2023 

# Clean the workspace
rm(list=ls())

# Load R package 
source("R/packages.R")
source("R/functions.R")

## Import benthic data and filter by focus tier4

TIER4 <<- "1808"
X <- read.csv("data/data_ready_new.csv") %>%
  filter(tier4 == TIER4)

# Import associated tier4 
tier4 <- geojson_sf("data/tier-4.json")

# Import predictive layer

HexPred_sf_raw <- st_read("data/predictive_layer_v3.shp", quiet = TRUE) %>%
  dplyr::select(Tier5:L_5____2,reef_ar,reefid,geometry)%>%
  dplyr::rename(MaxDHW= DHW_t5, 
                LagMaxDHW.1 = LDHW_5_1,
                LagMaxDHW.2 = LDHW_5_2,
                `Wave_hours(weighted)` = st_5___,
                `LagWave_hours(weighted).1` = L_5____1,
                `LagWave_hours(weighted).2` = L_5____2) 

X$Tier5 <- as.character(X$tier5)
X$fYEAR <- as.character(X$fYEAR)

HexPred_sf <- HexPred_sf_raw

# Making reefid 
HexPred_sf$reefid <- as.factor(HexPred_sf$reefid)

HexPred_reefid <- HexPred_sf %>%
  filter(fYEAR == "2004") %>%
  group_by(Tier5) %>% 
  mutate(reefid_merged = as.factor(paste0(reefid, collapse = "_"))) %>%
  dplyr::select(Tier5, reefid_merged) %>%
  distinct() %>%
  st_drop_geometry()

HexPred_reefid2 <-inner_join(HexPred_sf %>% data.frame() , HexPred_reefid) %>%
  group_by(Tier5, fYEAR) %>% 
  filter(row_number()==1) %>%
  replace(is.na(.), 0) %>%
  st_as_sf(sf_column_name = "geometry")%>%
  dplyr::select(., - reefid) %>%
  rename(reefid = reefid_merged) 

tal_reefid <- HexPred_sf %>%
  filter(fYEAR == "2004") %>%
  group_by(Tier5) %>% 
  count() %>%
  ungroup() %>%
  dplyr::rename(reef_n = n) %>%
  group_by(reef_n) %>%
  count() %>%
  mutate(prop = (n / length(unique(HexPred_reefid2$Tier5))*100)) %>%
  arrange(desc(prop)) %>%
  st_drop_geometry()

### Coral data 
# For the visualization, keep only the locations with more than 2 observations 

reef_tal <- X %>%
  group_by(REEF, TRANSECT_NO) %>% count() %>%
  filter(n>2)

X_vis <- X %>%
  filter(REEF %in% reef_tal$REEF)

p_vis_data <- ggplot(X_vis) + geom_line(aes(x = fYEAR, y = (COUNT/TOTAL)*100, colour=  as.factor(TRANSECT_NO), group = interaction(as.factor(TRANSECT_NO), REEF)), 
                                        show.legend = FALSE) + 
  facet_wrap(~REEF, ncol=4) + theme_bw() +
  labs(x = "Year", y = "Coral cover", subtitle = paste("Tier",TIER4, sep="")) +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 8, margin = margin()))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

ggsave(plot =  p_vis_data, width=8, height=10, file = "extra/viz_coral_data.png")

### Spatial scales of aggregation
#Tier4 

map_tier4 <- ggplot() + 
  geom_sf(data = tier4, fill="transparent") +
  geom_sf(data = tier4 %>% filter(tier_id == TIER4), fill="#0072B2FF", alpha=.3) +
  geom_label_repel(data = tier4 %>% filter(tier_id == TIER4), aes(label = tier_id, geometry = geometry),
                   stat = "sf_coordinates",
                   min.segment.length = 0) + theme_bw() +
  xlab("Longitude") + ylab("Latitude")

ggsave(plot =  map_tier4, width=4, height=4, file = "extra/viz_tier4_loc.png")

# Tier5
map_tier5 <- ggplot() + 
  geom_sf(data = tier4 %>% filter(tier_id == TIER4), fill="transparent")+
  geom_sf(data = HexPred_sf_raw %>% filter(fYEAR == "2004"), col="gray90", size=.8) +
  geom_sf(data = HexPred_sf_raw %>% filter(fYEAR == "2004" & Tier5 %in% unique(X$tier5)), fill = "#0072B2FF", col= "black", size = 4.8)+ theme_bw() + 
  ggtitle(TIER4)


ggsave(plot =  map_tier5, width=6, height=6, file = "extra/viz_coral_data_loc.png")

### Disturbance data

#### Cyclone exposure

yy <- unique(HexPred_reefid2$fYEAR)
plot.list <- list()

for ( i in 1: length(unique(HexPred_sf_raw$fYEAR))){
  plot.list[[i]] <- tm_shape(HexPred_sf_raw %>% filter(fYEAR == yy[i])) +
    tm_polygons("Wave_hours(weighted)",title=paste0("Wave intensity in ", yy[i]),legend.hist=TRUE,
                breaks = seq(0,plyr::round_any(max(HexPred_sf_raw$`Wave_hours(weighted)`), 1, f = ceiling), by = 4), palette = "-magma") +
    tm_layout(legend.outside = TRUE, legend.outside.position = "right", legend.hist.size = 1, legend.hist.height = 0.4, 
              legend.hist.width = 0.7)
  tmap_save(plot.list[[i]], file = sprintf('extra/Cycl%02d.png', i) )
} 


png.files_cycl <- sprintf("extra/Cycl%02d.png", 1:length(unique(HexPred_reefid2$fYEAR))) 
GIF.convert(png.files_cycl, output = "extra/animation_cyclone.gif")

# Remove png files from the extra folder 
del_files <- list.files(path = "extra/", pattern='Cycl', all.files= T, full.names= T)
unlink(del_files)

#### Heat stress

plot.list <- list()

for ( i in 1: length(unique(HexPred_sf_raw$fYEAR))){
  plot.list[[i]] <- tm_shape(HexPred_sf_raw %>% filter(fYEAR == yy[i])) +
    tm_polygons("MaxDHW",title=paste0("DHW in ", yy[i]),legend.hist=TRUE,
                breaks = seq(0,plyr::round_any(max(HexPred_sf_raw$MaxDHW), 1, f = ceiling), by = 4), palette = "-viridis") +
    tm_layout(legend.outside = TRUE, legend.outside.position = "right", legend.hist.size = 1, legend.hist.height = 0.4, 
              legend.hist.width = 0.7)
  tmap_save(plot.list[[i]], file = sprintf('extra/DHW%02d.png', i) )
} 

png.files_dhw <- sprintf("extra/DHW%02d.png", 1:length(unique(HexPred_reefid2$fYEAR))) 
GIF.convert(png.files_dhw, output = "extra/animation_dhw.gif")

# Remove png files from the extra folder 
del_files <- list.files(path = "extra/", pattern='DHW', all.files= T, full.names= T)
unlink(del_files)


#### QA/QC

# Applying control quality - finding extreme events values for tier5 without observations 
out_cycl <- quantile(HexPred_sf_raw$`Wave_hours(weighted)`, probs = 0.975)
out_bleach <- quantile(HexPred_sf_raw$MaxDHW, probs = 0.975)

HexPred_sf <- HexPred_sf_raw %>%  
  mutate(As.Data = ifelse(Tier5 %in% X$Tier5, "Yes", "No"))%>%
  mutate(across(c(`Wave_hours(weighted)`,`LagWave_hours(weighted).1`,`LagWave_hours(weighted).2`), ~ifelse(.x >= out_cycl & As.Data == "No", NA, .x))) %>%
  mutate(across(c(MaxDHW,LagMaxDHW.1,LagMaxDHW.2), ~ifelse( .x >= out_bleach & As.Data == "No", NA, ifelse(.x < out_bleach, .x, .x)))) 

na_cyclone <- HexPred_sf %>% filter(is.na(`Wave_hours(weighted)`))

# For viz 
tier5_all_centr <- st_centroid(HexPred_sf %>% filter(fYEAR == "2004"))

map_na_cycl <- ggplot() + 
  geom_sf(data = tier4 %>% filter(tier_id == TIER4), fill="transparent")+
  geom_sf(data = tier5_all_centr %>% dplyr::select(., -fYEAR), col="gray90", size=.8) +
  geom_sf(data = st_as_sf(na_cyclone), fill = "#CCBA72", col= "black", size = 4.3) +
  theme_bw() + facet_wrap(~ fYEAR, ncol = 2)
ggsave(plot =  map_na_cycl , width=8, height=6, file = "extra/viz_cycl_control.png")

na_dhw <-  HexPred_sf %>% filter(is.na(`MaxDHW`)) 

map_na_dhw <- ggplot() + 
  geom_sf(data = tier4 %>% filter(tier_id == TIER4), fill="transparent")+
  geom_sf(data = tier5_all_centr %>% dplyr::select(., -fYEAR), col="gray90", size=.8) +
  geom_sf(data = st_as_sf(na_dhw), fill = "#79402E", col= "black", size = 4.3) +
  theme_bw() + facet_wrap(~ fYEAR, ncol = 2)
ggsave(plot =  map_na_dhw , width=8, height=6, file = "extra/viz_dhw_control.png")



