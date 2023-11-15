# Codes associated with the basis function exploratory analysis
# Author: Julie Vercelloni
# Contact: j.vercelloni@aims.gov.au
# Last update: November 2023 

# Clean the workspace
rm(list=ls())

# Load R package and functions
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

## Basis function

# Model full
basis <- auto_basis(STplane(),
                    ST_BAUs,
                    tunit = "years",
                    #nres = 2L, # for development (model runs in ~2min each)
                    nres = 3L, # for final run (model runs in ~10-15min each)
                    regular = TRUE)

p2_1 <- show_basis(basis@Basis2)

# Find locations of temporal knots 
orign_knots <- basis@Basis2@df     
 
# fit
start_time <- Sys.time()
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
end_time <- Sys.time()
time <- end_time - start_time

pred_10K <- predict(M, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_10K <- as.data.frame(pred_10K$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model full") %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

# Model 5K
# Pick 5 knots

# Spatial remains similar 
G_spatial <- auto_basis(manifold = plane(), 
                        data=ST_BAUs,
                        #nres = 2L, # for development (model runs in ~2min)
                        nres = 3L, # for final run (model runs in ~10-15min)
                        prune= 0)

# Pick 5 random knots
knots <- sample(orign_knots$loc1,5)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)

p2_2 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_5K <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
              Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2 + (1 | reefid), 
            data = list(STObj), 
            BAUs = ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_5K <- end_time - start_time

pred_5K <- predict(M_5K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_5K <- as.data.frame(pred_5K$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 5K") %>%
tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

# Model 3K
# Pick 3 knots
knots <- sample(orign_knots$loc1,3)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)

p2_3 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_3K <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
              Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2 + (1 | reefid), 
            data = list(STObj), 
            BAUs = ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_3K <- end_time - start_time

pred_3K <- predict(M_3K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_3K <- as.data.frame(pred_3K$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 3K") %>%
tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

# Model 1K
# Pick 1 knot
knots <- sample(orign_knots$loc1,1)

G_temporal <- local_basis(manifold=real_line(),      # functions on real line
                          loc = matrix(knots),   # location parameter
                          scale = rep(2,length(knots)),          # scale parameter
                          type = "Gaussian")
basis <- TensorP(G_spatial,G_temporal)

p2_4 <- show_basis(basis@Basis2)

start_time <- Sys.time()
M_1K <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
              Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2 + (1 | reefid), 
            data = list(STObj), 
            BAUs = ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_1K <- end_time - start_time

pred_1K <- predict(M_1K, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_1K <- as.data.frame(pred_1K$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model 1K") %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

# Model noK
G_spatial <- auto_basis(manifold = plane(), 
                        data=ST_BAUs,
                        nres = 2L, 
                        prune= 0)

start_time <- Sys.time()

M_noK <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
               Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2 + (1 | reefid), 
            data = list(STObj), 
            BAUs = ST_BAUs, 
            basis = G_spatial, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)
end_time <- Sys.time()
time_noK <- end_time - start_time

pred_noK <- predict(M_noK, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_noK <- as.data.frame(pred_noK$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  mutate(knots = "Model noK") %>%
tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc, knots),
                    names_to = "draw", 
                    values_to = "pred"
)

### Summary predictions at tier5 

post_dist_df <- rbind(post_dist_df_10K, post_dist_df_5K, post_dist_df_3K, post_dist_df_1K, post_dist_df_noK)

pred_sum_sf <- post_dist_df %>% group_by(fYEAR,Tier5, knots) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% group_by(Tier5) %>% slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower) %>%
  mutate(Tier5_fYEAR = paste0(Tier5,fYEAR))

pred_FRK_data <- pred_sum_sf %>% filter(Tier5 %in% unique(X$tier5))

X$Tier5 <- as.character(X$tier5)
X$fYEAR <- as.character(X$fYEAR)

X_tal <- X %>% group_by(tier5) %>% tally()%>%
  filter(n<3)

pred_with_data <- pred_FRK_data %>% data.frame() %>% full_join(X) %>%
  filter(!Tier5 %in% X_tal$tier5) %>% droplevels()

pred_with_data$knots <- as.factor(pred_with_data$knots)
pred_with_data$knots <- factor(pred_with_data$knots, levels = c("Model full", "Model 5K",
                                                                "Model 3K", "Model 1K", "Model noK"))
  
# Data fit 
p_data <- ggplot(pred_with_data %>% filter(Tier5 == "10064")) + 
  geom_line(aes(x = fYEAR, y = (COUNT/TOTAL)*100, group = interaction(as.factor(TRANSECT_NO), REEF)), 
            show.legend = FALSE, linewidth=.1, col="grey30") + 
  geom_ribbon(aes(x=fYEAR,ymin=.lower*100, ymax=.upper*100, group=1),alpha=.5, fill ="#0072B2FF") +
  geom_line(aes(x=fYEAR, y=pred*100, group=1),size=.6) +
  facet_wrap(~ knots, ncol=3) +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

ggsave(plot =  p_data , width=8, height=6, file = "extra/viz_pred_tier5_data_basisfunction.png")

### Summary predictions at tier4

Reef_Areas <- HexPred_sf %>% distinct(Tier5, .keep_all = TRUE) %>%
  mutate(reef_ar = reef_ar / 1000000) %>% # areas in square kilometer
  dplyr::select(Tier5, reef_ar) %>% 
  st_drop_geometry()

Sim_Tier4 <- inner_join(post_dist_df, Reef_Areas)
tot_area <- sum(Reef_Areas$reef_ar) 

Sim_Tier4_coverage <- Sim_Tier4 %>% 
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(fYEAR, draw, knots) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE),
            cover_prop = cover / tot_area) 

pred_tier4 <-  Sim_Tier4_coverage %>% group_by(knots, fYEAR) %>% 
  ggdist::median_hdci(cover_prop)%>%
  dplyr::select(knots:.upper)%>%
  mutate(tier4 = TIER4) %>%
  data.frame() 

colnames(pred_tier4) <- c("knots","Year", "Mean" ,"Lower", "Upper", "tier4")

pred_tier4$knots <- as.factor(pred_tier4$knots)
pred_tier4$knots <- factor(pred_tier4$knots, levels = c("Model full", "Model 5K",
                                                                "Model 3K", "Model 1K", "Model noK"))

p_tier4 <- ggplot(pred_tier4 %>% data.frame()) +
  geom_ribbon(aes(x = Year, ymin=Lower, ymax=Upper, group=1), alpha=.2, fill="#0072B2FF")+
  geom_line(aes(x=Year, y=Mean,group=1), col="black", linewidth=1.1)+
  geom_point(aes(x=Year, y=Mean), col="black", size=2.1)+
  xlab("Year") +ylab("Coral cover")+theme_bw()+ facet_wrap(~knots)+
  theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=8),axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'))+
  scale_x_discrete(breaks=c(2004,2008,2012,2016,2020))

ggsave(plot =  p_tier4 , width=8, height=6, file = "extra/viz_pred_tier4_basisfunction.png")

# Viz basis function locations 

p_basis <- p2_1 + ggtitle("Model full") + p2_2 + ggtitle("Model 5K") + 
           p2_3 + ggtitle("Model 3K") + p2_4 + ggtitle("Model 1K")

ggsave(plot =  p_basis , width=6, height=5, file = "extra/viz_basisfunction.png")

#### AIC table
model_name <- c("Model full", "Model 5K",
                "Model 3K", "Model 1K", "Model noK")
AIC_values <- c(AIC(M), AIC(M_5K), AIC(M_3K), AIC(M_1K), AIC(M_noK))

AIC_table <- data.frame(cbind(model_name, round(AIC_values,2))) 
colnames(AIC_table) <- c("Model", "AIC")
AIC_table <- write.csv(AIC_table, file = "extra/table_aic_basisfunction.csv", row.names = F)

#### Effect of disturbances
coef_table_all <- coef_uncertainty(M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[1])

coef_table_5K <- coef_uncertainty(M_5K, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[2])


coef_table_3K <- coef_uncertainty(M_3K, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[3])

coef_table_1K<- coef_uncertainty(M_1K, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[4])

coef_table_noK<- coef_uncertainty(M_noK, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[5])

# All table together 

coef_table <- rbind(coef_table_all[-1,], coef_table_5K[-1,],
                    coef_table_3K[-1,], coef_table_1K[-1,],
                    coef_table_noK[-1,]) %>%
  mutate(param = case_when(param == "MaxDHW" ~ "Heat stress",
                           param == "LagMaxDHW.1" ~ "Heat stress (lag1)",
                           param == "LagMaxDHW.2" ~ "Heat stress (lag2)",
                           param == "Wave_hours.weighted." ~ "Cyclone exposure",
                           param == "LagWave_hours.weighted..1" ~ "Cyclone exposure (lag1)",
                           param == "LagWave_hours.weighted..2" ~ "Cyclone exposure (lag2)"))

coef_table$Name <- as.factor(coef_table$Name)
coef_table$Name<- factor(coef_table$Name, levels = c("Model full", "Model 5K",
                                                        "Model 3K", "Model 1K", "Model noK"))
# Make viz
p_coef <- ggplot()+ geom_point(data = coef_table, aes(y=param, x=`50%`)) +
  geom_errorbar(data = coef_table, aes(y = param, xmin = `2.5%`, xmax = `97.5%`), width=.1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~Name) +
  theme_bw() + labs(x = "Effect size", y = "")

ggsave(plot =  p_coef,width=6, height=4, file = "extra/fixed_effects_basisfunction.png")

#### Computing time 

time_values <- c(time/60, time_5K/60, time_3K/60, time_1K/60, time_noK/60)
time_table <- data.frame(cbind(model_name, round(time_values,2)))
colnames(time_table) <- c("Model", "time")
write.csv(time_table, file = "extra/table_time_basisfunction.csv", row.names = F)
