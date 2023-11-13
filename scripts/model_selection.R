## Appendix A: Codes associated with the model selection
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

# Import predictive layer and make reefid

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

# Applying control quality - removing values associated with extreme disturbance events for tier5 without observations 
out_cycl <- quantile(HexPred_sf_raw$`Wave_hours(weighted)`, probs = 0.975)
out_bleach <- quantile(HexPred_sf_raw$MaxDHW, probs = 0.975)

HexPred_sf <- HexPred_sf_raw %>%  
  mutate(As.Data = ifelse(Tier5 %in% X$Tier5, "Yes", "No"))%>%
  mutate(across(c(`Wave_hours(weighted)`,`LagWave_hours(weighted).1`,`LagWave_hours(weighted).2`), 
                ~ifelse(.x >= out_cycl & As.Data == "No", NA, .x))) %>%
  mutate(across(c(MaxDHW,LagMaxDHW.1,LagMaxDHW.2), ~ifelse( .x >= out_bleach & As.Data == "No", NA, ifelse(.x < out_bleach, .x, .x)))) %>%
  mutate(across(c(MaxDHW:`LagWave_hours(weighted).2`), ~ as.numeric(scale(.))))

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

####################################################
#################################################### MODEL SELECTION 
####################################################

# Making objects for the models - same for all models   

## Construct STIDF object from benthic data

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
nHEX <- nrow(subset(HexPred_sp, fYEAR == 2004))       # no. of hexagons
nYEAR <- length(unique(HexPred_sp@data$fYEAR))        # no. of years     

HexPred_sp@data$n_spat <- rep(1:nHEX, each = nYEAR)   # index for each spatial BAU 
BAUs_spat <- subset(HexPred_sp, fYEAR == 2004)        # extract spatial grid (first year)
coordnames(BAUs_spat) <- c("LONGITUDE", "LATITUDE")

# Construct spatio-temporal BAUs (will not contain covariate information for now)

ST_BAUs <- auto_BAUs(manifold = STplane(),
                     data = STObj,
                     spatial_BAUs = BAUs_spat,
                     tunit = "years")

ST_BAUs <- ST_BAUs[, 1:nYEAR, 1:2]                 # remove last year (automatically inserted by FRK)
ST_BAUs$fYEAR <- as.character(ST_BAUs$t + 2003)    # create fYEAR variable 
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
                    # nres = 2L, # for development (model runs in ~2min)
                    nres = 3L, # for final run (model runs in ~30min per model)
                    regular = TRUE)

###############################
############################### Run model full
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

# Extract model predictions 

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

pred_sum_sf <- post_dist_df %>% 
  group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% 
               group_by(Tier5) %>% 
               slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower)

###############################
############################### Run model random only
start_time <- Sys.time()

M_no <- FRK(f = COUNT ~ 1 + (1 | reefid), 
            data = list(STObj), 
            BAUs = ST_BAUs, 
            basis = basis, 
            response = "binomial", 
            link = "logit", 
            K_type = "precision", 
            method = "TMB", 
            est_error = FALSE)

end_time <- Sys.time()
time_no <- end_time - start_time

# Extract model predictions and focus on data

pred_no <- predict(M_no, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_no <- as.data.frame(pred_no$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                      names_to = "draw", 
                      values_to = "pred"
  )

# Summary predictions at tier5 

pred_sum_sf_no <- post_dist_df_no %>% 
  group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% 
               group_by(Tier5) %>% 
               slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower)

###############################
############################### Run model covariates only
start_time <- Sys.time()

M_no_random <- FRK(f = COUNT ~ 1 + MaxDHW + LagMaxDHW.1 + LagMaxDHW.2 + 
                     Wave_hours.weighted. + LagWave_hours.weighted..1 + LagWave_hours.weighted..2, 
                   data = list(STObj), 
                   BAUs = ST_BAUs, 
                   basis = basis, 
                   response = "binomial", 
                   link = "logit", 
                   K_type = "precision", 
                   method = "TMB", 
                   est_error = FALSE)

end_time <- Sys.time()
time_no_random  <- end_time - start_time

# Extract model predictions and focus on data

pred_no_random <- predict(M_no_random, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_no_random <- as.data.frame(pred_no_random$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                      names_to = "draw", 
                      values_to = "pred"
  )

# Summary predictions at tier5 

pred_sum_sf_no_random <- post_dist_df_no_random %>% 
  group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% 
               group_by(Tier5) %>% 
               slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower)

###############################
############################### Run model covariates (no lags) + random
start_time <- Sys.time()

M_no_lag <- FRK(f = COUNT ~ 1 + MaxDHW + Wave_hours.weighted. + (1 | reefid), 
                data = list(STObj), 
                BAUs = ST_BAUs, 
                basis = basis, 
                response = "binomial", 
                link = "logit", 
                K_type = "precision", 
                method = "TMB", 
                est_error = FALSE)

end_time <- Sys.time()
time_no_lag <- end_time - start_time

# Extract model predictions and focus on data

pred_no_lag <- predict(M_no_lag, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_no_lag <- as.data.frame(pred_no_lag$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                      names_to = "draw", 
                      values_to = "pred"
  )

# Summary predictions at tier5 

pred_sum_sf_no_lag <- post_dist_df_no_lag %>% 
  group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% 
               group_by(Tier5) %>% 
               slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower)

###############################
############################### Run model covariates only (no lags)
start_time <- Sys.time()

M_no_lag_no_random <- FRK(f = COUNT ~ 1 + MaxDHW +  Wave_hours.weighted., 
                          data = list(STObj), 
                          BAUs = ST_BAUs, 
                          basis = basis, 
                          response = "binomial", 
                          link = "logit", 
                          K_type = "precision", 
                          method = "TMB", 
                          est_error = FALSE)

end_time <- Sys.time()
time_no_lag_no_random <- end_time - start_time

# Extract model predictions and focus on data

pred_no_lag_no_random <- predict(M_no_lag_no_random, type = c("mean"))

# Extracting posterior distributions of predictive locations 

post_dist_df_no_lag_no_random <- as.data.frame(pred_no_lag_no_random$MC$mu_samples) %>% 
  mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
  mutate(Tier5 = ST_BAUs@data$Tier5) %>%
  mutate(id_loc = row_number()) %>%
  tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                      names_to = "draw", 
                      values_to = "pred"
  )

# Summary predictions at tier5 

pred_sum_sf_no_lag_no_random <- post_dist_df_no_lag_no_random %>% 
  group_by(fYEAR,Tier5) %>% 
  ggdist::median_hdci(pred)%>%
  inner_join(HexPred_reefid2 %>% 
               group_by(Tier5) %>% 
               slice(1) %>% dplyr::select(geometry,Tier5)) %>% 
  st_as_sf(sf_column_name = "geometry") %>%
  mutate(Unc = .upper - .lower)

######## Results 

#### AIC table
model_name <- c("Model full", "Model random only", "Model covariates only", "Model covariates (no lags) + random", "Model covariates only (no lags)")
AIC_values <- c(AIC(M), AIC(M_no), AIC(M_no_random), AIC(M_no_lag), AIC(M_no_lag_no_random))

AIC_table <- data.frame(cbind(model_name, round(AIC_values,2))) 
colnames(AIC_table) <- c("Model", "AIC")
AIC_table <- write.csv(AIC_table, file = "extra/table_aic_modelchoice.csv", row.names = F)


#### Effect of disturbances
coef_table_all <- coef_uncertainty(M, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[1])

coef_table_all_no_random <- coef_uncertainty(M_no_random, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[3])


coef_table_all_no_lag <- coef_uncertainty(M_no_lag, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = TRUE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[4])

coef_table_all_no_lag_no_random <- coef_uncertainty(M_no_lag_no_random, percentiles = c(2.5, 50, 97.5), nsim = 400, random_effects = FALSE) %>%
  data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols =  !rowname, names_to = "param", values_to = "value") %>%
  mutate(Type = ifelse(str_starts(param, "X.Intercept"), "random", "fixed")) %>%
  tidyr::pivot_wider(names_from = rowname, values_from = value) %>%
  filter(Type == "fixed") %>%
  mutate(Name = model_name[5])

# All table together 

coef_table <- rbind(coef_table_all[-1,], coef_table_all_no_random[-1,],
                    coef_table_all_no_lag[-1,], coef_table_all_no_lag_no_random[-1,]) %>%
  mutate(param = case_when(param == "MaxDHW" ~ "Heat stress",
                           param == "LagMaxDHW.1" ~ "Heat stress (lag1)",
                           param == "LagMaxDHW.2" ~ "Heat stress (lag2)",
                           param == "Wave_hours.weighted." ~ "Cyclone exposure",
                           param == "LagWave_hours.weighted..1" ~ "Cyclone exposure (lag1)",
                           param == "LagWave_hours.weighted..2" ~ "Cyclone exposure (lag2)"))

coef_table$Name <- as.factor(coef_table$Name)
coef_table$Name <- factor(coef_table$Name, levels = c("Model full", "Model covariates only", "Model covariates (no lags) + random", "Model covariates only (no lags)"))

# Make viz
p_coef <- ggplot()+ geom_point(data = coef_table, aes(y=param, x=`50%`)) +
  geom_errorbar(data = coef_table, aes(y = param, xmin = `2.5%`, xmax = `97.5%`), width=.1) + 
  geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~Name) +
  theme_bw() + labs(x = "Effect size", y = "")

ggsave(plot =  p_coef,width=6, height=4, file = "extra/fixed_effects_modelchoice.png")

#### Predictions range and uncertainty 

pred_sum_all_pred <- rbind(pred_sum_sf %>% dplyr::select(pred) %>% mutate(`Model name` = model_name[1]),
                           pred_sum_sf_no %>% dplyr::select(pred) %>% mutate(`Model name` = model_name[2]),
                           pred_sum_sf_no_random %>% dplyr::select(pred) %>% mutate(`Model name` = model_name[3]),
                           pred_sum_sf_no_lag %>% dplyr::select(pred) %>% mutate(`Model name` = model_name[4]),
                           pred_sum_sf_no_lag_no_random %>% dplyr::select(pred) %>% mutate(`Model name` = model_name[5]))

p.range <- ggplot(pred_sum_all_pred) + aes(y = factor(`Model name`, level = model_name), x= pred) +
  stat_halfeye(alpha = .6, .width = .95) + 
  theme_bw() +
  ylab("") + xlab("Prediction range") 

ggsave(plot =  p.range, width=5, height=5, file = "extra/range_modelchoice.png")

pred_sum_all <- rbind(pred_sum_sf %>% dplyr::select(Unc) %>% mutate(`Model name` = model_name[1]),
                      pred_sum_sf_no %>% dplyr::select(Unc) %>% mutate(`Model name` = model_name[2]),
                      pred_sum_sf_no_random %>% dplyr::select(Unc) %>% mutate(`Model name` = model_name[3]),
                      pred_sum_sf_no_lag %>% dplyr::select(Unc) %>% mutate(`Model name` = model_name[4]),
                      pred_sum_sf_no_lag_no_random %>% dplyr::select(Unc) %>% mutate(`Model name` = model_name[5]))

p.unc <- ggplot(pred_sum_all) + aes(y = factor(`Model name`, level = model_name), x= Unc) +
  stat_halfeye(alpha = .6, .width = .95) + 
  theme_bw() +
  ylab("") + xlab("Uncertainty range") 

ggsave(plot =  p.unc, width=5, height=5, file = "extra/uncertainty_modelchoice.png")

#### Tier4 trajectories 

# Weight by reef areas
Reef_Areas <- HexPred_sf %>% distinct(Tier5, .keep_all = TRUE) %>%
  mutate(reef_ar = reef_ar / 1000000) %>% # areas in square kilometer
  dplyr::select(Tier5, reef_ar) %>% 
  st_drop_geometry()

Sim_Tier4_full <- left_join(post_dist_df,Reef_Areas) %>%
  mutate(Name = model_name[1])

Sim_Tier4_no <- left_join(post_dist_df_no,Reef_Areas) %>%
  mutate(Name = model_name[2])

Sim_Tier4_no_random <- left_join(post_dist_df_no_random,Reef_Areas) %>%
  mutate(Name = model_name[3])

Sim_Tier4_no_lag <- left_join(post_dist_df_no_lag,Reef_Areas) %>%
  mutate(Name = model_name[4])

Sim_Tier4_no_lag_no_random <- left_join(post_dist_df_no_lag_no_random,Reef_Areas) %>%
  mutate(Name = model_name[5])

Sim_Tier4 <- rbind(Sim_Tier4_full, Sim_Tier4_no, Sim_Tier4_no_random, Sim_Tier4_no_lag, Sim_Tier4_no_lag_no_random)

Sim_Tier4_coverage <- Sim_Tier4 %>% 
  mutate(weighted_pred = pred * reef_ar) %>%
  group_by(Name, fYEAR, draw) %>%
  summarize(cover = sum(weighted_pred, na.rm = TRUE)) 

pred_tier4 <-  Sim_Tier4_coverage %>% group_by(Name, fYEAR) %>% 
  ggdist::median_hdci(cover)%>%
  dplyr::select(Name:.upper)%>%
  mutate(tier4 = TIER4) %>%
  data.frame() 

colnames(pred_tier4) <- c("Name","Year", "Mean" ,"Lower", "Upper", "tier4")

pred_tier4$Name <- as.factor(pred_tier4$Name)
pred_tier4$Name <- factor(pred_tier4$Name, levels = c("Model full", "Model random only" ,"Model covariates only", "Model covariates (no lags) + random", "Model covariates only (no lags)"))


p_tier4 <- ggplot(pred_tier4 %>% data.frame()) +
  geom_ribbon(aes(x = Year, ymin=Lower, ymax=Upper, group=1), alpha=.4, fill="orange")+
  geom_line(aes(x=Year,y=Mean,group=1),col="black",size=1.1)+ facet_wrap(~Name)+
  xlab("Year") + ylab("Coverage (sq km)")+ theme_bw()+
  theme(axis.text.x = element_text(size=11, angle = 90, hjust = 1),legend.position = "none",
        legend.text = element_text(colour = "black", size = 9), 
        axis.text.y = element_text(size=11),axis.title.y=element_text(size=13),
        axis.title.x=element_text(size=13))+
  scale_x_discrete(breaks=c(2004,2006,2008,2010,2012,2014,2016,2018,2020,2022))

ggsave(plot =  p_tier4, width=8, height=5, file = "extra/tier4_modelchoice.png")

#### Computing time 

time_values <- c(time, time_no, time_no_random, time_no_lag, time_no_lag_no_random)
time_table <- data.frame(cbind(model_name, round(time_values,2)))
colnames(time_table) <- c("Model", "time")
time_table <- write.csv(time_table, file = "extra/table_time_modelchoice.csv", row.names = F)

