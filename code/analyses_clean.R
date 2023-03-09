library(tidyverse) #to wrangle data and plot
library(glmmTMB) #to run models
library(ggeffects) #to generate model predictions
library(performance) #model checks
library(DHARMa) #model checks
library(rsample) #for generating CIs for the hurdle models
library(car) #for AIC checks
#each section number refers to a corresponding figure in the manuscript

#load in data-------------------------------------------------------------------
snail_transects_full <- read_csv("data/Snail_Data_Transects.csv") %>%
  #create a unique quad id
  unite("quadrat", c(Transect,distance_from_marsh), remove = FALSE) 

snail_quad_level <- snail_transects_full %>%
  #only include one line per quadrat, where the total counts are included 
  #instead of the individual-level data
  filter(!is.na(total_num)) %>%
  #scale continuous variables - use c() to remove attributes
  mutate(scale_distance_from_water = c(scale(distance_from_water))) %>% 
  #remove individual-specific columns to avoid confusion
  select(-c(length, barnacled))

snail_quad_nooutliers <- snail_transects_full %>%
  #remove the three outliers to compare the models with and without them
  filter(total_num < 20) %>%
  #only include one line per quadrat, where the total counts are included 
  #instead of the individual-level data
  filter(!is.na(total_num)) %>%
  #scale continuous variables
  mutate(scale_distance_from_water = c(scale(distance_from_water)))  %>% 
  #remove individual-specific columns to avoid confusion
  select(-c(length, barnacled))

snail_individuals <- snail_transects_full %>% 
  #only take observations where there is a length recording (i.e., when there
  #was actually a snail observed)
  filter(!is.na(length))  %>%
  fill(total_num, .direction = "down") %>% 
  #scale continuous variables
  mutate(scale_distance_from_water = c(scale(distance_from_water)),
         scale_length = c(scale(length)),
         #turn "barnacled" into a factor and a binary response
         barnacled = as.factor(barnacled), 
         barnacled_binary = case_when(barnacled == "no" ~ 0, 
                                      barnacled == "yes" ~ 1)) 

snail_movement <- read_csv("data/Snail Data - Behavioural_observations.csv") %>% 
  #remove the three observations with trematodes, as these very clearly impact
  #snail size but are too correlated with other variables to go into the model
  filter(trematodes == "no") %>% 
  mutate(wgt_ratio = wet_wgt_barnacles/wet_wgt_snail,
         scale_length = c(scale(length)),
         #make sure there are zeros instead of NAs for the no barnacles
         wet_wgt_barnacles = case_when(barnacled == "no" ~ 0,
                                       TRUE ~ wet_wgt_barnacles),
         wgt_ratio = case_when(barnacled == "no" ~ 0,
                               TRUE ~ wgt_ratio),
         scale_wgt_barnacles = c(scale(wet_wgt_barnacles)),
         scale_wgt_ratio = c(scale(wgt_ratio)),
         #also log and scale the two continuos predictors because there are some
         #pretty strong outliers
         #use half of the lowest non-zero value for the constant to be able to 
         #transform the zeros
         log_wgt_barnacles = log(wet_wgt_barnacles + 0.015),
         log_wgt_ratio = log(wgt_ratio + 0.01239669),
         scale_log_wgt_barnacles = c(scale(log_wgt_barnacles)),
         scale_log_wgt_ratio = c(scale(log_wgt_ratio)),
         #create a new variables of zeros and ones for speed, where 0 = did not 
         #move and 1 = moved
         distance_logistic = case_when(distance == 0 ~ 0,
                                       TRUE ~ 1))

snail_movement_gamma <- snail_movement %>% 
  #get rid of all the zeros for the analysis of the non-zero values
  filter(distance != 0)

snail_movement_nooutliers <- read_csv("data/Snail Data - Behavioural_observations.csv") %>% 
  #remove the three observations with trematodes, as these very clearly impact
  #snail size but are too correlated with other variables to go into the mdoel
  filter(trematodes == "no") %>% 
  filter(length > 20) %>% 
  mutate(wgt_ratio = wet_wgt_barnacles/wet_wgt_snail,
         scale_length =c(scale(length)),
         wet_wgt_barnacles = case_when(num_barnacles == 0 ~ 0,
                                       TRUE ~ wet_wgt_barnacles),
         scale_wgt_barnacles = c(scale(wet_wgt_barnacles)),
         wgt_ratio = case_when(num_barnacles == 0 ~ 0,
                               TRUE ~ wgt_ratio),
         scale_wgt_ratio = c(scale(wgt_ratio)),
         #create a new variables of zeros and ones for speed, where 0 = did not 
         #move and 1 = moved
         distance_logistic = case_when(distance == 0 ~ 0,
                                       TRUE ~ 1))

snail_recapture <- read_csv("data/snail_recapture.csv") %>% 
  #relevel the months
  mutate(date = factor(date, levels = c("June", "July"))) %>% 
  #switch to df to get rid of a weird tibble issue in the for loop
  as.data.frame()
#we can also rearrange this dataset so each individual row is a snail, rather 
#than a group of snails - this'll make it easier to plot and analyze
#10 snails in each group were originally caputured, so the number of snails not
#recaptured is 10-num_recap
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

#Summary stats------------------------------------------------------------------
#total quadrats surveyed in population analysis
quads <- snail_transects_full %>% 
  distinct(quadrat) %>% 
  nrow()

#total snails seen across all quads
total_snails <- snail_quad_level %>% 
  summarize(n = sum(total_num))

#snails with individual measurements made
snails_detail <- snail_individuals %>% 
  nrow()

#number of snails per quadrat
range(snail_quad_level$total_num)
mean(snail_quad_level$total_num)
sd(snail_quad_level$total_num)

#snail sizes
range(snail_individuals$length)
mean(snail_individuals$length)
sd(snail_individuals$length)

#snails with barnacles - quadrat totals
snail_quad_level %>% 
  summarize(num_barnacled = sum(num_barnacled, na.rm = TRUE),
            num_total = sum(total_num)) %>% 
  mutate(prop_barnacled = num_barnacled/num_total)

#snails with barnacles - measured
snail_individuals %>% 
  group_by(barnacled) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = barnacled, values_from = n) %>% 
  mutate(prop_barnacled = yes/(yes+no))

#proportion of snails recaptured
snail_recap_long %>% 
  group_by(recap, date) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = recap, values_from = n) %>% 
  mutate(prop_recap = yes/(yes+no))

#Q1A: Is there a relationship between snail abundance per quadrat and ---- 
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

#and we can make the same two models with a linear effect of distance
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

#the models with a linear term are definitely a worse fit based on the residuals
#and from a biological sense the quadratic models are much more sound
#but we can also run a quick AIC comparison to confirm that they are
#statistically better
AIC(total_num_nb_out_linear, total_num_nb_out)
#quadratic is better!

#Q1B: Is there a relationship between snail length and distance from the --------
#water?
snail_length <- glmmTMB(length ~ scale_distance_from_water + 
                          I(scale_distance_from_water^2) +
                          (1|Transect/quadrat), 
                        data = snail_individuals)
summary(snail_length)
#check model fit
snail_length_check <- simulateResiduals(snail_length)
plot(snail_length_check)
#the residuals aren't perfect but pretty reasonable given the data

#linear option instead
snail_length_linear <- glmmTMB(length ~ scale_distance_from_water + 
                                 (1|Transect/quadrat), 
                               data = snail_individuals)
summary(snail_length_linear)
#check model fit
snail_length_check_linear <- simulateResiduals(snail_length_linear)
plot(snail_length_check_linear)
#much worse - quadratic is a lot better

#Q2: Does the probability of being barnacled change with distance from the------ 
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
#residuals look good!

#Q3: Does the weight of barnacles impact snail speed--------------------------
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

#Q4: Do recapture rates vary depending on whether barnacles are present?-----
mod_recap <- glmmTMB(recap_binary ~ barnacled * date + 
                       (1|pond), 
                     data = snail_recap_long, 
                     family = binomial)
summary(mod_recap)
recap_check <- simulateResiduals(mod_recap)
plot(recap_check)
#nice

#QS1: Do barnacles affect snail speed?-----------------------------------------
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
#QS2: Does the number of barnacles impact snail speed--------------------------
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

#QS3: Do barnacles impact snail allometry---------------------------------------
#with no outlier
mod_allometry_nooutlier <- glmmTMB(wet_wgt_snail ~ barnacled * scale_length +
                                     (1|pond),
                                   data = snail_movement_nooutliers)
summary(mod_allometry_nooutlier)
allometry_nooutlier_check <- simulateResiduals(mod_allometry_nooutlier)
plot(allometry_nooutlier_check)
#looks great