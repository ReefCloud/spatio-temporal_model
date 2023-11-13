# Codes associated with the spatio-temporal model 
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

X$Tier5 <- as.character(X$tier5)
X$fYEAR <- as.character(X$fYEAR)

# Import predictive layer 

HexPred_sf_raw <- st_read("data/predictive_layer_v3.shp", quiet = TRUE) %>%
  dplyr::select(Tier5:L_5____2,reef_ar,reefid,geometry)%>%
  dplyr::rename(MaxDHW = DHW_t5, # rename to align with covariate.hexpred table
                LagMaxDHW.1 = LDHW_5_1,
                LagMaxDHW.2 = LDHW_5_2,
                `Wave_hours(weighted)` = st_5___,
                `LagWave_hours(weighted).1` = L_5____1,
                `LagWave_hours(weighted).2` = L_5____2) %>%
  filter(fYEAR %in% X$fYEAR)

# Applying control quality - finding extreme events values for tier5 without observations 
out_cycl <- quantile(HexPred_sf_raw$`Wave_hours(weighted)`, probs = 0.975)
out_bleach <- quantile(HexPred_sf_raw$MaxDHW, probs = 0.975)

HexPred_sf <- HexPred_sf_raw %>%  
  mutate(As.Data = ifelse(Tier5 %in% X$Tier5, "Yes", "No"))%>%
  mutate(across(c(`Wave_hours(weighted)`,`LagWave_hours(weighted).1`,`LagWave_hours(weighted).2`), ~ifelse(.x >= out_cycl & As.Data == "No", NA, .x))) %>%
  mutate(across(c(MaxDHW,LagMaxDHW.1,LagMaxDHW.2), ~ifelse( .x >= out_bleach & As.Data == "No", NA, ifelse(.x < out_bleach, .x, .x)))) %>%
  mutate(across(c(MaxDHW:`LagWave_hours(weighted).2`), ~ as.numeric(scale(.)))) 

# Making reefid 

HexPred_sf$reefid <- as.factor(HexPred_sf$reefid)

HexPred_reefid <- HexPred_sf %>%
  filter(fYEAR == min(HexPred_sf$fYEAR)) %>%
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


X$Year <- as.Date(paste0(as.character(X$fYEAR),"-01-01"))  # needs to be a Date object
X$k_Z <- X$TOTAL                                           # this is ntrials
lon_idx <- which(names(X) == "LONGITUDE")                  
lat_idx <- which(names(X) == "LATITUDE")
STObj <- stConstruct(x = X,                               
                     space = c(lon_idx, lat_idx), 
                     time = "Year",                      
                     interval = TRUE)      # time reflects an interval

# Making BAUs 

HexPred_sp <- as_Spatial(HexPred_reefid2)                                    # convert to sp
nHEX <- nrow(subset(HexPred_sp, fYEAR == min(HexPred_sf$fYEAR)))       # no. of hexagons
nYEAR <- length(unique(HexPred_sp@data$fYEAR))        # no. of years     

HexPred_sp@data$n_spat <- rep(1:nHEX, each = nYEAR)   # index for each spatial BAU 
BAUs_spat <- subset(HexPred_sp, fYEAR == min(HexPred_sf$fYEAR))        # extract spatial grid (first year)
coordnames(BAUs_spat) <- c("LONGITUDE", "LATITUDE")

# Construct spatio-temporal BAUs (will not contain covariate information for now)

ST_BAUs <- auto_BAUs(manifold = STplane(),
                     data = STObj,
                     spatial_BAUs = BAUs_spat,
                     tunit = "years")

ST_BAUs <- ST_BAUs[, 1:nYEAR, 1:2]                 # remove last year (automatically inserted by FRK)
ST_BAUs$fYEAR <- as.character(ST_BAUs$t + (min(as.numeric(HexPred_sf$fYEAR))-1))    # create fYEAR variable 
ST_BAUs$n_spat <- rep(1:nHEX, nYEAR)               # create (spatial) index for each BAU 

# Update BAUs with covariate information
ST_BAUs@data <- left_join(ST_BAUs@data, HexPred_sp@data , by = c("fYEAR","n_spat")) 
ST_BAUs$fs <- 1                   # scalar fine-scale variance matrix
ST_BAUs@sp@proj4string  <- CRS()  # set CRS to NA

# Update BAUs with random effects and change class 

ST_BAUs@data$reefid <- as.character(ST_BAUs@data$reefid) 
ST_BAUs@data$reefid <- as.factor(ST_BAUs@data$reefid)
ST_BAUs@data$yearid <- as.factor(ST_BAUs@data$fYEAR)

# Covariates must only be in BAUs, so remove covariates associated with data

overlapping_fields <- intersect(names(STObj@data), names(ST_BAUs@data))
STObj@data[,overlapping_fields] <- NULL

# Create basis functions

basis <- auto_basis(STplane(),
                    ST_BAUs,
                    tunit = "years",
                   #  nres = 2L, # for development (model runs in ~2min)
                     nres = 3L, # for final run (model runs in ~30min)
                    regular = TRUE)

# Fit FRK model 
M <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
           Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2 + (1 | reefid), 
         data = list(STObj), 
         BAUs = ST_BAUs, 
         basis = basis, 
         response = "binomial", 
         link = "logit", 
         K_type = "precision", 
         method = "TMB", 
         est_error = FALSE)

pred <- predict(M, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                      names_to = "draw", 
                      values_to = "pred"
  )

# Summary predictions at tier5 

pred_sum_sf <- post_dist_df %>% group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% group_by(Tier5) %>% slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower) %>%
  mutate(Tier5_fYEAR = paste0(Tier5,fYEAR))

## Tier5 

### Spatio-temporal predictions and uncertainty

yy <- unique(pred_sum_sf$fYEAR)
plot.list <- list()

pred_sum_sf <- pred_sum_sf %>%
  mutate(pred_perc = pred *100) %>%
  mutate(Unc_perc = Unc * 100)

# Coral cover 
for ( i in 1: length(unique(pred_sum_sf$fYEAR))){
  plot.list[[i]] <- tm_shape(pred_sum_sf %>% filter(fYEAR == yy[i])) +
    tm_polygons("pred_perc",title=paste0("Coral cover in ", yy[i]),legend.hist=TRUE,
                breaks = seq(0,plyr::round_any(max(pred_sum_sf$pred_perc), 10, f = ceiling), by = 10)) +
    tm_layout(legend.outside = TRUE, legend.outside.position = "right", legend.hist.size = 1, legend.hist.height = 0.4, 
              legend.hist.width = 0.7)
  tmap_save(plot.list[[i]], file = sprintf('extra/Pred%02d.png', i) )
} 

# Associated uncertainty 
for ( i in 1: length(unique(pred_sum_sf$fYEAR))){
  plot.list[[i]] <- tm_shape(pred_sum_sf %>% filter(fYEAR == yy[i])) +
    tm_polygons("Unc_perc",title=paste0("Unc. in ", yy[i]),
                breaks = seq(0, plyr::round_any(max(pred_sum_sf$Unc_perc), 10, f = ceiling), by = 10), palette = "BrBG") +
    tm_layout(legend.outside = TRUE, legend.outside.position = "right", legend.hist.size = 1, legend.hist.height = 0.4, 
              legend.hist.width = 0.7)
  tmap_save(plot.list[[i]], file = sprintf('extra/Pred_unc%02d.png', i) )
} 

# Create the GIF

png.files_cover <- sprintf("extra/Pred%02d.png", 1:length(unique(pred_sum_sf$fYEAR))) 
png.files_unc <- sprintf("extra/Pred_unc%02d.png", 1:length(unique(pred_sum_sf$fYEAR))) 

GIF.convert(png.files_cover, output = "extra/animation_cover.gif")
GIF.convert(png.files_unc, output = "extra/animation_cover_unc.gif")

# Remove png files from the extra folder 
del_files <- list.files(path = "extra/", pattern='Pred', all.files= T, full.names= T)
unlink(del_files)

### Coral cover trajectories of tier5 with data 

pred_FRK_data <- pred_sum_sf %>% filter(Tier5 %in% unique(X$tier5))

X$Tier5 <- as.character(X$tier5)
X$fYEAR <- as.character(X$fYEAR)

# Keep tier5 with more than three observations only 

X_tal <- X %>% group_by(tier5) %>% tally()%>%
  filter(n<3)

pred_with_data <- pred_FRK_data %>% data.frame() %>% full_join(X) %>%
  filter(!Tier5 %in% X_tal$tier5) %>% droplevels()

p_data <- ggplot(pred_with_data) + 
  geom_line(aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF)), 
            show.legend = FALSE, linewidth=.1, col="grey30") + 
  geom_ribbon(aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.5, fill ="#0072B2FF") +
  geom_line(aes(x=fYEAR, y=pred*100, group=1),size=.6) +
  facet_wrap(~Tier5, ncol=3) +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

ggsave(plot =  p_data , width=8, height=10, file = "extra/viz_pred_tier5_data.png")

### Coral cover trajectories of tier5 with and without data 

# Pick 9 random tiers for visualization purpose 

random_tier5 <- sample(unique(pred_with_data$Tier5), 9, replace = FALSE)

#p_data_cut <- ggplot(pred_with_data %>% 
#                       filter(Tier5 %in% random_tier5)) + 
#  geom_line(aes(x = fYEAR, y = COUNT/TOTAL, group = interaction(as.factor(TRANSECT_NO), REEF)), 
#            show.legend = FALSE, linewidth=.1, col="grey30") + 
#  geom_ribbon(aes(x=fYEAR,ymin=.lower, ymax=.upper, group=1),alpha=.5, fill ="#0072B2FF") +
#  geom_line(aes(x=fYEAR, y=pred, group=1),size=.6) +
#  facet_wrap(~Tier5, ncol=3) + ylim(0,1) +
#  labs(x = "Year", y = "Coral cover", subtitle = "with data") +
#  theme_bw()+
#  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
#        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
#        axis.title.x=element_text(size=11),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        strip.background = element_rect(fill = 'white'))+
#  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

# Centroid locations

tier5_nodata_centr <- st_centroid(pred_sum_sf %>% filter(fYEAR == "2004" & !Tier5 %in% unique(X$tier5) ))
tier5_data_centr <- st_centroid(pred_sum_sf%>% filter(fYEAR == "2004" & Tier5 %in% unique(X$tier5) ))
tier5_all_centr <- st_centroid(pred_sum_sf %>% filter(fYEAR == "2004"))

# Get distances
# find the nearest neighbor to the centroid of the each polygon with data

dist_hex <- t(st_distance(tier5_data_centr, tier5_all_centr)) %>% data.frame()
colnames(dist_hex) <- tier5_data_centr$Tier5
rownames(dist_hex) <- tier5_all_centr$Tier5

dist_hex_long <- dist_hex %>%
  tibble::rownames_to_column()%>%
  rename("Tier5_all" = rowname) %>%
  tidyr::pivot_longer(!Tier5_all, names_to = "Tier5_data", values_to = "distance")

dist_hex_long$distance <- as.numeric(dist_hex_long$distance) 

# Get quantiles to select close and far tier5 cells from data 

quant_dist <- quantile(dist_hex_long$distance)

# Tier 5 with no data and further from the data locations 
# Using distance smaller than 25th quantile and pick 9 random  

dist_min <- dist_hex_long %>% group_by(Tier5_data) %>%
  filter(!distance == 0) %>%
  filter(distance < quant_dist[2]) %>%
  data.frame() %>%
  sample_n(9)

closest_tier5 <- pred_sum_sf %>% filter(Tier5 %in% dist_min$Tier5_all)

p_close <- ggplot(closest_tier5) + 
  geom_ribbon(aes(x=fYEAR,ymin=.lower, ymax=.upper, group=1),alpha=.9, fill= "#D68A8A") +
  geom_line(aes(x=fYEAR, y=pred, group=1), size=.6) +
  facet_wrap(~Tier5) + ylim(0,1) +
  labs(x = "", y = "Coral cover", subtitle = "close cells") + theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

# Tier 5 with no data and further from the data locations 
# Using distance greater than the 75% quantile and pick 9 random  

dist_max <- dist_hex_long %>%
  filter(distance > quant_dist[4]) %>%
  data.frame() %>%
  sample_n(9)

far_tier5 <- pred_sum_sf %>% filter(Tier5 %in% dist_max$Tier5_all)

p_far <- ggplot(far_tier5) + 
  geom_ribbon(aes(x=fYEAR,ymin=.lower, ymax=.upper, group=1),alpha=.3, fill= "#83BFA9") +
  geom_line(aes(x=fYEAR, y=pred, group=1), size=.6) +
  facet_wrap(~Tier5) + ylim(0,1) +
  labs(x = "Year", y = "Coral cover", subtitle = "far cells") + theme_bw() +
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

# Common viz 

tier4 <- geojson_sf("data/tier-4.json")

map_tier5 <- ggplot() + 
  geom_sf(data = tier4 %>% filter(tier_id == TIER4), fill="transparent")+
  geom_sf(data = tier5_all_centr, col="gray90", size=.8) +
  geom_sf(data = st_as_sf(pred_with_data), fill = "#0072B2FF", col= "black", size = 4.8) +
  geom_sf(data = closest_tier5, fill = "#D68A8A", col = "black", size = 4.8) +
  geom_sf(data = far_tier5, fill = "#83BFA9", col = "black", size = 4.8) + theme_bw()

p_traj <- p_close + p_far + plot_layout(nrow = 2)

ggsave(plot =  map_tier5 , width=6, height=6, file = "extra/viz_map_pred_alltier5.png")
ggsave(plot =  p_traj , width=8, height=10, file = "extra/viz_pred_alltier5.png")

### Tier4

# Weight by reef areas

Reef_Areas <- HexPred_sf %>% distinct(Tier5, .keep_all = TRUE) %>%
  mutate(reef_ar = reef_ar / 1000000) %>% # areas in square kilometer
  dplyr::select(Tier5, reef_ar) %>% 
  st_drop_geometry()

Sim_Tier4 <- left_join(post_dist_df,Reef_Areas)

tot_area <- sum(Reef_Areas$reef_ar) 

Sim_Tier4_coverage <- Sim_Tier4 %>% 
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(fYEAR, draw) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE)) 

pred_tier4 <-  Sim_Tier4_coverage %>% group_by(fYEAR) %>% 
  ggdist::median_hdci(cover)%>%
  dplyr::select(fYEAR:.upper)%>%
  mutate(tier4 = TIER4) %>%
  data.frame() 

colnames(pred_tier4) <- c("Year", "Mean" ,"Lower", "Upper", "tier4")

p_tier4 <- ggplot(pred_tier4 %>% data.frame()) +
  geom_ribbon(aes(x = Year, ymin=Lower, ymax=Upper, group=1), alpha=.4, fill="orange")+
  geom_line(aes(x=Year, y=Mean,group=1), col="black", linewidth=1.1)+
  xlab("Year") +ylab("Coverage (sq km)")+theme_bw()+
  theme(axis.text.x = element_text(size=13),legend.position = "none",
        legend.text = element_text(colour = "black", size = 9), 
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15))+
  scale_x_discrete(breaks=c(2004,2006,2008,2010,2012,2014,2016,2018,2020,2022))

ggsave(plot =  p_tier4 , width=6, height=4, file = "extra/viz_pred_tier4.png")

### Effect of disturbances

# Full table 
coef_table_all <- coef_uncertainty(M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE)%>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value")%>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value)

# Fixed effects only 

coef_table_fixed <- coef_table_all %>%
  filter(Type == "fixed")  %>% 
  mutate(param = case_when(param == "MaxDHW" ~ "Heat stress",
                                                    param == "LagMaxDHW.1" ~ "Heat stress (lag1)",
                                                    param == "LagMaxDHW.2" ~ "Heat stress (lag2)",
                                                    param == "Wave_hours.weighted." ~ "Cyclone exposure",
                                                    param == "LagWave_hours.weighted..1" ~ "Cyclone exposure (lag1)",
                                                    param == "LagWave_hours.weighted..2" ~ "Cyclone exposure (lag2)"))


p_coef <- ggplot(coef_table_fixed[-1,], aes(y=param, x=`50%`))+ geom_point() +
  geom_errorbar(aes(y = param, xmin = `2.5%`, xmax = `97.5%`), width=.1) + 
  geom_vline(xintercept = 0, linetype = "dashed") +theme_bw() +
  xlab("Effect size") + ylab("")

ggsave(plot =  p_coef,width=6, height=4, file = "extra/fixed_effects_model.png")

### Model validation
# Create DHARMA residuals 
res <- make_dharma_res(M)

# Visual diagnostics 
p_res <- gg_dharma(res, integerResponse = FALSE)
ggsave(plot =  p_res,width=6, height=7, file = "extra/dharma_model.png")



