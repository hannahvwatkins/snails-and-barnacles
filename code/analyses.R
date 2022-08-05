library(tidyverse) #to wrangle data and plot
library(glmmTMB) #to run models
library(ggeffects) #to generate model predictions
library(labelled) #to remove attributes from scaled columns
library(performance) #model checks
library(DHARMa) #model checks
library(rsample) #for generating CIs for the hurdle models

#load in data-------------------------------------------------------------------
snail_transects_full <- read_csv("data/Snail_Data_Transects.csv") %>%
  #create a unique quad id
  unite("quadrat", c(Transect,distance_from_water), remove = FALSE) 

snail_quad_level <- snail_transects_full %>%
  #only include one line per quadrat, where the total counts are included 
  #instead of the individual-level data
  filter(!is.na(total_num)) %>%
  #scale continuous variables
  mutate(scale_distance_from_water = scale(distance_from_water),
         scale_length = scale(length)) %>% 
  #remove the attributes from the scaled columns so they're just normal cols
  remove_attributes(c("scaled:center","scaled:scale")) 

snail_quad_nooutliers <- snail_transects_full %>%
  #remove the three outliers to compare the models with and without them
  filter(total_num < 20) %>%
  #only include one line per quadrat, where the total counts are included 
  #instead of the individual-level data
  filter(!is.na(total_num)) %>%
  #scale continuous variables
  mutate(scale_distance_from_water = scale(distance_from_water),
         scale_length = scale(length)) %>% 
  #remove the attributes from the scaled columns so they're just normal cols
  remove_attributes(c("scaled:center","scaled:scale"))

snail_individuals <- snail_transects_full %>% 
  #only take observations where there is a length recording (i.e., when there
  #was actually a snail observed)
  filter(!is.na(length))  %>%
  #scale continuous variables
  mutate(scale_distance_from_water = scale(distance_from_water),
         scale_length = scale(length),
  #turn "barnacled" into a factor and a binary response
         barnacled = as.factor(barnacled), 
         barnacled_binary = case_when(barnacled == "no" ~ 0, 
                                      barnacled == "yes" ~ 1)) %>% 
  #remove the attributes from the scaled columns so they're just normal cols
  remove_attributes(c("scaled:center","scaled:scale")) 

snail_movement <- read_csv("data/Snail Data - Behavioural_observations.csv") %>% 
  #remove the three observations with trematodes, as these very clearly impact
  #snail size but are too correlated with other variables to go into the mdoel
  filter(trematodes == "no") %>% 
  mutate(wgt_ratio = wet_wgt_barnacles/wet_wgt_snail,
         scale_length = scale(length),
         wet_wgt_barnacles = case_when(barnacled == "no" ~ 0,
                                       TRUE ~ wet_wgt_barnacles),
         scale_wgt_barnacles = scale(wet_wgt_barnacles),
         wgt_ratio = case_when(barnacled == "no" ~ 0,
                               TRUE ~ wgt_ratio),
         scale_wgt_ratio = scale(wgt_ratio),
         #create a new variables of zeros and ones for speed, where 0 = did not 
         #move and 1 = moved
         distance_logistic = case_when(distance == 0 ~ 0,
                                              TRUE ~ 1)) %>% 
  remove_attributes(c("scaled:center","scaled:scale")) 

snail_movement_gamma <- snail_movement %>% 
  #get rid of all the zeros for the analysis of the non-zero values
  filter(distance != 0)

snail_movement_nooutliers <- read_csv("data/Snail Data - Behavioural_observations.csv") %>% 
  #remove the three observations with trematodes, as these very clearly impact
  #snail size but are too correlated with other variables to go into the mdoel
  filter(trematodes == "no") %>% 
  filter(length > 20) %>% 
  mutate(wgt_ratio = wet_wgt_barnacles/wet_wgt_snail,
         scale_length = scale(length),
         wet_wgt_barnacles = case_when(num_barnacles == 0 ~ 0,
                                       TRUE ~ wet_wgt_barnacles),
         scale_wgt_barnacles = scale(wet_wgt_barnacles),
         wgt_ratio = case_when(num_barnacles == 0 ~ 0,
                               TRUE ~ wgt_ratio),
         scale_wgt_ratio = scale(wgt_ratio),
         #create a new variables of zeros and ones for speed, where 0 = did not 
         #move and 1 = moved
         distance_logistic = case_when(distance == 0 ~ 0,
                                       TRUE ~ 1)) %>% 
  remove_attributes(c("scaled:center","scaled:scale")) 

snail_recapture <- read_csv("data/snail_recapture.csv") %>% 
  #relevel the months
  mutate(date = factor(date, levels = c("June", "July"))) %>% 
  #switch to df to get rid of a weird tibble issue in the for loop
  as.data.frame()
#we can also rearrange this dataset so each individual row is a snail, rather 
#than a group of snails - this'll make it easier to plot and analyze
snail_recap_long <- NULL
df <- NULL
for (i in 1:nrow(snail_recapture)){
  df <- tibble(observer = rep(snail_recapture[i,1], times = 10),
               date = rep(snail_recapture[i,2], times = 10),
               pond = rep(snail_recapture[i,3], times = 10),
               barnacled = rep(snail_recapture[i,4], times = 10),
               recap = c(rep("yes", times = snail_recapture[i,5]),
                         rep("no", times = (10 - snail_recapture[i,5]))),
               recap_binary = case_when(recap == "yes" ~ 1,
                                        recap == "no" ~ 0))
  snail_recap_long <- rbind(snail_recap_long, df)
  }

#-------------------------------------------------------------------------------
#Q1: Is there a relationship between snail abundance per quadrat and ---- 
#distance from the water?
#since we're looking at counts, we'll try a Poisson model first
total_num <- glmmTMB(total_num ~ scale_distance_from_water + (1|Transect), 
                         data = snail_quad_level, family = poisson) 
#Check for overdispersion, because that's common with Poisson
check_overdispersion(total_num)

# Overdispersion detected so use negative binomial distribution instead 
total_num_nb <- glmmTMB(total_num ~ scale_distance_from_water + 
                          I(scale_distance_from_water^2) +
                           (1|Transect), 
                         family = nbinom2,
                         data = snail_quad_level)
summary(total_num_nb)

#and we'll do a quick model check on that with DHARMa
total_num_check <- simulateResiduals(total_num_nb)
plot(total_num_check)

#and compare to the data without the outliers
total_num_nb_out <- glmmTMB(total_num ~ scale_distance_from_water + 
                              I(scale_distance_from_water^2) +
                              (1|Transect), 
                            family = nbinom2,
                            data = snail_quad_nooutliers) 
summary(total_num_nb_out)
total_num_out_check <- simulateResiduals(total_num_nb_out)
plot(total_num_out_check)

#and we can make the same two models with a linear effect of disatnce
total_num_nb_linear <- glmmTMB(total_num ~ scale_distance_from_water + 
                          (1|Transect), 
                        family = nbinom2,
                        data = snail_quad_level)
summary(total_num_nb_linear)
#and we'll do a quick model check on that with DHARMa
total_num_check_linear <- simulateResiduals(total_num_nb_linear)
plot(total_num_check_linear)

#and compare to the data without the outliers
total_num_nb_out_linear <- glmmTMB(total_num ~ scale_distance_from_water + 
                              (1|Transect), 
                            family = nbinom2,
                            data = snail_quad_nooutliers) 
summary(total_num_nb_out_linear)
total_num_out_check_linear <- simulateResiduals(total_num_nb_out_linear)
plot(total_num_out_check_linear)

#Q2: Is there a relationship between snail length and distance from the --------
#water?
snail_length <- glmmTMB(length ~ scale_distance_from_water + 
                          I(scale_distance_from_water^2) +
                          (1|Transect/quadrat), 
                         data = snail_individuals)
summary(snail_length)
#check model fit
snail_length_check <- simulateResiduals(snail_length)
plot(snail_length_check)

#linear option instead
snail_length_linear <- glmmTMB(length ~ scale_distance_from_water + 
                          (1|Transect/quadrat), 
                        data = snail_individuals)
summary(snail_length_linear)
#check model fit
snail_length_check_linear <- simulateResiduals(snail_length_linear)
plot(snail_length_check_linear)

#Q3: Does the probability of being barnacled change with distance from the------ 
#water and with snail length?
mod_prob_barn <- glmmTMB(barnacled_binary ~ scale_distance_from_water * 
                           scale_length +
                           (1|Transect/quadrat), 
                         data = snail_individuals, 
                         family = binomial)
summary(mod_prob_barn)

#check model fit
mod_prob_barn_check <- simulateResiduals(mod_prob_barn)
plot(mod_prob_barn_check)

#Q4.1: Do barnacles affect snail speed?-----------------------------------------
#Hurdle models
# This is two submodels:
# one a logistic regression that deals with snails that moved vs didn't move and 
# the other a Gamma  model that just considers snails that moved. 
# Use with three barnacle metrics: presence/absence, weight of barnacles, ratio 
# of barnacle to snail weight

hurdle_mod1 <- glmmTMB(distance_logistic ~ barnacled + scale_length + (1|pond), 
                     family = binomial(link = "logit"),
                     data = snail_movement)
summary(hurdle_mod1) 
hurdle_mod1_check <- simulateResiduals(hurdle_mod1)
plot(hurdle_mod1_check)
#not perfect, but pretty good

hurdle_mod2 <- glmmTMB(distance ~ barnacled + scale_length + (1|pond),
                     family = Gamma(link = "log"),
                     data = snail_movement_gamma)
summary(hurdle_mod2) 
hurdle_mod2_check <- simulateResiduals(hurdle_mod2)
plot(hurdle_mod2_check)
#looks excellent

#Q4.2: Does the number of barnacles impact snail speed--------------------------
hurdle_mod3 <- glmmTMB(distance_logistic ~ scale_wgt_barnacles + scale_length + 
                         (1|pond), 
                       family = binomial(link = "logit"),
                       data = snail_movement)
summary(hurdle_mod3) 
hurdle_mod3_check <- simulateResiduals(hurdle_mod3)
plot(hurdle_mod3_check)
#looks fine

hurdle_mod4 <- glmmTMB(distance ~ scale_wgt_barnacles + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)
summary(hurdle_mod4) 
hurdle_mod4_check <- simulateResiduals(hurdle_mod4)
plot(hurdle_mod4_check)
#looks excellent

#Q4.3: Does the number of barnacles impact snail speed--------------------------
hurdle_mod5 <- glmmTMB(distance_logistic ~ scale_wgt_ratio + scale_length + (1|pond), 
                       family = binomial(link = "logit"),
                       data = snail_movement)
summary(hurdle_mod5) 
hurdle_mod5_check <- simulateResiduals(hurdle_mod5)
plot(hurdle_mod5_check)
#looks fine

hurdle_mod6 <- glmmTMB(distance ~ scale_wgt_ratio + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)
summary(hurdle_mod6) 
hurdle_mod6_check <- simulateResiduals(hurdle_mod6)
plot(hurdle_mod6_check)
#looks excellent

#Q5: Do barnacles impact snail allometry---------------------------------------
mod_allometry <- glmmTMB(wet_wgt_snail ~ barnacled * scale_length +
                           (1|pond),
                         data = snail_movement)
summary(mod_allometry)
allometry_check <- simulateResiduals(mod_allometry)
plot(allometry_check)
#looks great

#with no outlier
mod_allometry_nooutlier <- glmmTMB(wet_wgt_snail ~ barnacled * scale_length +
                                     (1|pond),
                                   data = snail_movement_nooutliers)
summary(mod_allometry_nooutlier)
allometry_nooutlier_check <- simulateResiduals(mod_allometry_nooutlier)
plot(allometry_nooutlier_check)
#looks great
#Q6: Do recapture rates vary depending on whether barnacles are present?-----
mod_recap <- glmmTMB(recap_binary ~ barnacled * date + 
                       (1|pond), 
                     data = snail_recap_long, 
                     family = binomial)
summary(mod_recap)
recap_check <- simulateResiduals(mod_recap)
plot(recap_check)
#nice