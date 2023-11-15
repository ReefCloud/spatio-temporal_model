## Appendix A: Codes associated with the leave-out data analysis
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
#################################################### LEAVE-OUT DATA ANALYSIS 
####################################################
X <- X %>%
  mutate(id = row_number()) %>%
  mutate(Transect = paste(REEF,TRANSECT_NO, sep="_"),
         Reef = str_replace(sub(" Site.*", "", REEF), " ", "_"))


## Preparing the lists 

list.indicators <- list()
pred_sum_list <- list()
test.validation <- c("1 - rm(20% obs)", "2 - rm(20% reef)", "3 - rm(20% site)", "4 - rm(20% transect)", "5 - rm(4YRS)")
computing_time <- list()

# Make training and testing datasets for each validation test 

for (i in 1:length(test.validation)){
  test.name <- test.validation[i]
  if (test.name == test.validation[1]) {
    train <- X %>% sample_frac(.80)
    test  <- anti_join(X, train, by = 'id')
  } else if (test.name == test.validation[2]){
    reef_pick <-  sample(unique(X$Reef), ceiling(length(unique(X$Reef))*.8), replace = F) 
    train <-  X %>% filter(Reef %in% reef_pick)
    test <- anti_join(X, train, by = 'id')
  } else if (test.name == test.validation[3]){
    site_pick <-  sample(unique(X$REEF), ceiling(length(unique(X$REEF))*.8)) 
    train <-  X %>% filter(REEF %in% site_pick)
    test <-  anti_join(X, train, by = 'id')
  } else if(test.name == test.validation[4]){
    transect_pick <-  sample(unique(X$Transect), ceiling(length(unique(X$Transect))*.8)) 
    train <-  X %>% filter(Transect %in% transect_pick)
    test <-  anti_join(X, train, by = 'id')  
  }else{
    tal_transect <- X %>% group_by(Transect) %>% tally() %>% mutate(year_keep = n - 4) %>% filter(!year_keep < 1) %>% arrange(Transect)
    
    train <- X %>% inner_join(tal_transect) %>%   arrange(Transect) %>%
      tidyr::nest(.by = Transect)  %>% 
      mutate(n = tal_transect$year_keep) %>%
      mutate(samp = map2(data, n, sample_n)) %>%
      dplyr::select(-c(data,n)) %>%
      tidyr::unnest(samp) %>%
      data.frame()
    
    test <-  anti_join(X, train, by = 'id')   
  }
  
  ## Construct STIDF object from data
  train$Year <- as.Date(paste0(as.character(train$fYEAR),"-01-01"))  # needs to be a Date object
  train$k_Z <- train$TOTAL                                           # this is ntrials
  lon_idx <- which(names(train) == "LONGITUDE")                  
  lat_idx <- which(names(train) == "LATITUDE")
  STObj <- stConstruct(x = train,                               
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
                      nres = 3L, # for final run (model runs in ~10-15min per model)
                      regular = TRUE)
  
  ## Fit FRK model
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
  computing_time[i] <- end_time - start_time
  
  ## Predict on test data
  pred <- predict(M, type = c("mean"))
  
  ## Extracting posterior distributions of predictive locations 
  
  post_dist_df <- as.data.frame(pred$MC$mu_samples) %>% 
    mutate(fYEAR = ST_BAUs@data$fYEAR) %>%
    mutate(Tier5 = ST_BAUs@data$Tier5) %>%
    mutate(id_loc = row_number()) %>%
    tidyr::pivot_longer(!c(fYEAR,Tier5,id_loc),
                        names_to = "draw", 
                        values_to = "pred"
    )
  
  test$Tier5 <- as.character(test$tier5)
  test$fYEAR <- as.character(test$fYEAR)
  
  test_pred <- inner_join(test, post_dist_df, by = c("Tier5","fYEAR"))
  
  test_pred_sum <- test_pred %>% group_by(fYEAR,Tier5) %>% 
    ggdist::median_hdci(pred)%>%
    inner_join(test_pred %>% group_by(Tier5, fYEAR) %>% slice(1) %>%
                 mutate(COVER = COUNT / TOTAL) %>%
                 dplyr::select(Tier5,fYEAR,COVER)) %>%
    mutate(p_se = .upper - .lower) %>%
    mutate(Test = test.name)
  
  pred_sum_list[[i]] <- test_pred_sum 
  
  ## Measures
  
  list.indicators[[i]]<- c(test.name, round(coverage95(test_pred_sum$COVER, test_pred_sum$.lower, test_pred_sum$.upper),2),
                           round(IS95(test_pred_sum$COVER, test_pred_sum$.lower, test_pred_sum$.upper),2),
                           round(RMSPE(test_pred_sum$COVER, test_pred_sum$pred),2),
                           round(crps(test_pred_sum$COVER, data.frame(test_pred_sum$pred, test_pred_sum$p_se))$CRPS,2), round(AIC(M),2))
}

indicators_table <- data.frame(do.call(rbind,list.indicators))
colnames(indicators_table) <- c("Test","95% coverage", "95% interval score", "Root-mean-squared prediction error", "Continuous ranked probability score", "AIC")

indicator_table_leaveout <- write.csv(indicators_table, file = "extra/indicator_table_leaveout_modelchoice.csv", row.names = F)

time_values <- data.frame(do.call(rbind,computing_time))
time_table <- data.frame(cbind(test.validation, round(time_values,2)))
colnames(time_table) <- c("Test", "time")

time_table_leaveout <- write.csv(time_table, file = "extra/table_time_leaveout_modelchoice.csv", row.names = F)


