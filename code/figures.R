library(tidyverse) #for plotting and general data wrangling
library(ggeffects) #for plotting model output
library(rsample) #for generating CIs on hurdle model predictions
library(labelled) #for removing attributes from data
library(patchwork) #for combining panels
library(glmmTMB) #for running final models
library(brms) #for plotting hurdle model output
source("code/ggplot_paper_theme.R")
options(mc.cores = parallel::detectCores())

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
#run the final model from the analyses.R script
total_num_nb <- glmmTMB(total_num ~ scale_distance_from_water + 
                          I(scale_distance_from_water^2) +
                          (1|Transect), 
                        family = nbinom2,
                        data = snail_quad_level)
#predict 
total_num_predict <- ggpredict(total_num_nb, 
                               terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_quad_level$distance_from_water) +
           mean(snail_quad_level$distance_from_water))
#plot predictions on raw data
ggplot(total_num_predict, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_quad_level, 
             aes(distance_from_water, total_num),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Number of snails")

#ggsave("figures/Figure1A.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

total_num_nb_out <- glmmTMB(total_num ~ scale_distance_from_water + 
                              I(scale_distance_from_water^2) +
                              (1|Transect), 
                            family = nbinom2,
                            data = snail_quad_nooutliers) 
total_num_predict_out <- ggpredict(total_num_nb_out, 
                               terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_quad_nooutliers$distance_from_water) +
         mean(snail_quad_nooutliers$distance_from_water))

#plot predictions on raw data
ggplot(total_num_predict_out, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_quad_nooutliers, 
             aes(distance_from_water, total_num),
             size = 0.5) +
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Number of snails") 

#ggsave("figures/Figure1B.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#linear versions of both
total_num_nb_linear <- glmmTMB(total_num ~ scale_distance_from_water + 
                                 (1|Transect), 
                               family = nbinom2,
                               data = snail_quad_level)

#predict 
total_num_predict_linear <- ggpredict(total_num_nb_linear, 
                               terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_quad_level$distance_from_water) +
           mean(snail_quad_level$distance_from_water))
#plot predictions on raw data
ggplot(total_num_predict_linear, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_quad_level, 
             aes(distance_from_water, total_num),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Number of snails")

#ggsave("figures/Figure1C.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#linear no outliers
total_num_nb_out_linear <- glmmTMB(total_num ~ scale_distance_from_water + 
                                     (1|Transect), 
                                   family = nbinom2,
                                   data = snail_quad_nooutliers)
#predict 
total_num_predict_linear_out <- ggpredict(total_num_nb_out_linear, 
                               terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_quad_level$distance_from_water) +
           mean(snail_quad_level$distance_from_water))
#plot predictions on raw data
ggplot(total_num_predict_linear_out, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_quad_nooutliers, 
             aes(distance_from_water, total_num),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Number of snails")

#ggsave("figures/Figure1D.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q2: Is there a relationship between snail length and distance from the --------
#water?
snail_length <- glmmTMB(length ~ scale_distance_from_water + 
                          I(scale_distance_from_water^2) +
                          (1|Transect/quadrat), 
                        data = snail_individuals)
#predict 
snail_length_predict <- ggpredict(snail_length, 
                               terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_individuals$distance_from_water) +
           mean(snail_individuals$distance_from_water))
#plot predictions on raw data
ggplot(snail_length_predict, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_individuals, 
             aes(distance_from_water, length),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Snail length (mm)")

#ggsave("figures/Figure2A.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

snail_length_linear <- glmmTMB(length ~ scale_distance_from_water + 
                          (1|Transect/quadrat), 
                        data = snail_individuals)
#predict 
snail_length_predict_linear <- ggpredict(snail_length_linear, 
                                  terms="scale_distance_from_water [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(distance_from_water = x*sd(snail_individuals$distance_from_water) +
           mean(snail_individuals$distance_from_water))
#plot predictions on raw data
ggplot(snail_length_predict_linear, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_individuals, 
             aes(distance_from_water, length),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Snail length (mm)")

#ggsave("figures/Figure2B.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q3: Does the probability of being barnacled change with distance from the------ 
#water and with snail length?
mod_prob_barn <- glmmTMB(barnacled_binary ~ scale_distance_from_water * 
                           scale_length +
                           (1|Transect/quadrat), 
                         data = snail_individuals, 
                         family = binomial)
#predict 
mod_prob_predict <- ggpredict(mod_prob_barn, 
                                  terms="scale_length [all]") %>% 
  #convert scaled distances back to unscaled
  mutate(length = x*sd(snail_individuals$length) +
           mean(snail_individuals$length))
#plot predictions on raw data
ggplot(mod_prob_predict, aes(length, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_individuals, 
             aes(length, barnacled_binary),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Snail length (mm)",
       y = "Probability of carrying a barnacle")

#ggsave("figures/Figure3.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

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

hurdle_mod2 <- glmmTMB(distance ~ barnacled + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)

#predict 
hurdle1_predict <- ggpredict(hurdle_mod1, 
                              terms=c("barnacled")) %>% 
  mutate(barnacled = x)
#plot predictions on raw data
fig4A.1 <- ggplot(hurdle1_predict, aes(barnacled, predicted)) + 
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low), 
                width = 0.1) +
  geom_point(data = snail_movement, 
             aes(barnacled, distance_logistic),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  theme_paper() +
  labs(x = "Barnacle(s) present",
       y = "Probability of snail movement",
       title = "A") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))
#predict 
hurdle2_predict <- ggpredict(hurdle_mod2, 
                             terms=c("barnacled")) %>% 
  mutate(barnacled = x)
#plot predictions on raw data
fig4A.2 <-  ggplot(hurdle2_predict, aes(barnacled, predicted)) + 
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low), 
                width = 0.1) +
  geom_point(data = snail_movement_gamma, 
             aes(barnacled, distance),
             size = 0.5, position = position_jitter(width = 0.1)) + 
  theme_paper() +
  labs(x = "Barnacle(s) present",
       y = "Distance moved (cm)",
       title = "B") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

fig4A.1 + fig4A.2

#ggsave("figures/Figure4A_separate.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#having trouble with bootstrapping the combined model predictions of the two 
#submodels, so we'll rerun it in brms with flat priors to get an approximation
hurdle_mod12 <- brm(bf(distance ~ barnacled + scale_length + (1|pond),
                       hu ~ barnacled + scale_length + (1|pond)), 
                       family = hurdle_gamma(),
                       data = snail_movement, 
                    seed = 123,
                    iter = 5000,
                    control = list(adapt_delta = 0.9))
summary(hurdle_mod12)

hurdle12_predict <- ggpredict(hurdle_mod12, 
                             terms=c("barnacled")) %>% 
  mutate(barnacled = x)
#ignore uncertainty warning because we only want the global effect

ggplot(hurdle12_predict, aes(barnacled, predicted)) + 
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low), 
                width = 0.1) +
  geom_point(data = snail_movement, 
             aes(barnacled, distance),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  theme_paper() +
  labs(x = "Barnacle(s) present",
       y = "Distance moved (cm)")

#ggsave("figures/Figure4A_combined.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q4.2: Does the number of barnacles impact snail speed--------------------------
hurdle_mod3 <- glmmTMB(distance_logistic ~ scale_wgt_barnacles + scale_length + 
                         (1|pond), 
                       family = binomial(link = "logit"),
                       data = snail_movement)


hurdle_mod4 <- glmmTMB(distance ~ scale_wgt_barnacles + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)

#predict 
hurdle3_predict <- ggpredict(hurdle_mod3, 
                             terms=c("scale_wgt_barnacles [n=100]")) %>% 
  mutate(scale_wgt_barnacles = x,
         wet_wgt_barnacles = x*sd(snail_movement$wet_wgt_barnacles) +
           mean(snail_movement$wet_wgt_barnacles))
#plot predictions on raw data
fig4B.1 <- ggplot(hurdle3_predict, aes(wet_wgt_barnacles, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement, 
             aes(wet_wgt_barnacles, distance_logistic),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  theme_paper() +
  labs(x = "Wet weight of barnacles (g)",
       y = "Probability of snail movement",
       title = "A") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))
#predict 
hurdle4_predict <- ggpredict(hurdle_mod4, 
                             terms=c("scale_wgt_barnacles [n=100]")) %>% 
  mutate(scale_wgt_barnacles = x,
         wet_wgt_barnacles = x*sd(snail_movement$wet_wgt_barnacles) +
           mean(snail_movement$wet_wgt_barnacles))

#plot predictions on raw data
fig4B.2 <- ggplot(hurdle4_predict, aes(wet_wgt_barnacles, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement_gamma, 
             aes(wet_wgt_barnacles, distance),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Wet weight of barnacles (g)",
       y = "Distance moved (cm)",
       title = "B") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

fig4B.1 + fig4B.2

#ggsave("figures/Figure4B_separate.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#having trouble with bootstrapping the combined model predictions of the two 
#submodels, so we'll rerun it in brms with flat priors to get an approximation
hurdle_mod34 <- brm(bf(distance ~ scale_wgt_barnacles + scale_length + (1|pond),
                       hu ~ scale_wgt_barnacles + scale_length + (1|pond)), 
                    family = hurdle_gamma(),
                    data = snail_movement, 
                    seed = 345,
                    iter = 5000,
                    control = list(adapt_delta = 0.93))
summary(hurdle_mod34)

hurdle34_predict <- ggpredict(hurdle_mod34, 
                              terms=c("scale_wgt_barnacles [n=100]")) %>% 
  mutate(scale_wgt_barnacles = x,
         wet_wgt_barnacles = x*sd(snail_movement$wet_wgt_barnacles) +
           mean(snail_movement$wet_wgt_barnacles))
#ignore uncertainty warning because we only want the global effect

ggplot(hurdle34_predict, aes(wet_wgt_barnacles, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement, 
             aes(wet_wgt_barnacles, distance),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Wet weight of barnacles (g)",
       y = "Distance moved (cm)") 

#ggsave("figures/Figure4B_combined.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q4.3: Does the number of barnacles impact snail speed--------------------------
hurdle_mod5 <- glmmTMB(distance_logistic ~ scale_wgt_ratio + scale_length + (1|pond), 
                       family = binomial(link = "logit"),
                       data = snail_movement)

hurdle_mod6 <- glmmTMB(distance ~ scale_wgt_ratio + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)

#predict 
hurdle4_predict <- ggpredict(hurdle_mod5, 
                             terms=c("scale_wgt_ratio [n=100]")) %>% 
  mutate(scale_wgt_ratio = x,
         wgt_ratio = x*sd(snail_movement$wgt_ratio) +
           mean(snail_movement$wgt_ratio))
#plot predictions on raw data
fig4C.1 <- ggplot(hurdle4_predict, aes(wgt_ratio, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement, 
             aes(wgt_ratio, distance_logistic),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  theme_paper() +
  labs(x = "Ratio of barnacle weight to snail weight",
       y = "Probability of snail movement",
       title = "A") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))
#predict 
hurdle6_predict <- ggpredict(hurdle_mod6, 
                             terms=c("scale_wgt_ratio [n=100]")) %>% 
  mutate(scale_wgt_ratio = x,
         wgt_ratio = x*sd(snail_movement$wgt_ratio) +
           mean(snail_movement$wgt_ratio))

#plot predictions on raw data
fig4C.2 <- ggplot(hurdle6_predict, aes(wgt_ratio, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement_gamma, 
             aes(wgt_ratio, distance),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Ratio of barnacle weight to snail weight",
       y = "Distance moved (cm)",
       title = "B") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

fig4C.1 + fig4C.2

#ggsave("figures/Figure4C_separate.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#having trouble with bootstrapping the combined model predictions of the two 
#submodels, so we'll rerun it in brms with flat priors to get an approximation
hurdle_mod56 <- brm(bf(distance ~ scale_wgt_ratio + scale_length + (1|pond),
                       hu ~ scale_wgt_ratio + scale_length + (1|pond)), 
                    family = hurdle_gamma(),
                    data = snail_movement, 
                    seed = 5678,
                    iter = 5000,
                    control = list(adapt_delta = 0.93))
summary(hurdle_mod56)

hurdle56_predict <- ggpredict(hurdle_mod56, 
                              terms=c("scale_wgt_ratio [n=100]")) %>% 
  mutate(scale_wgt_ratio = x,
         wgt_ratio= x*sd(snail_movement$wgt_ratio) +
           mean(snail_movement$wgt_ratio))
#ignore uncertainty warning because we only want the global effect

ggplot(hurdle56_predict, aes(wgt_ratio, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_movement, 
             aes(wgt_ratio, distance),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Ratio of barnacle weight to snail weight",
       y = "Distance moved (cm)") 

#ggsave("figures/Figure4C_combined.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q5: Do barnacles impact snail allometry---------------------------------------
mod_allometry <- glmmTMB(wet_wgt_snail ~ barnacled * scale_length +
                           (1|pond),
                         data = snail_movement)
#predict 
mod_allo_predict <- ggpredict(mod_allometry, 
                              terms=c("scale_length", "barnacled")) %>% 
  #convert scaled distances back to unscaled
  mutate(length = x*sd(snail_movement$length) +
           mean(snail_movement$length),
         barnacled = group)
#plot predictions on raw data
ggplot(mod_allo_predict, aes(length, predicted, 
                             group = barnacled, fill = barnacled)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(colour = barnacled)) +
  geom_point(data = snail_movement, 
             aes(length, wet_wgt_snail, colour = barnacled),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Snail length (mm)",
       y = "Snail wet weight (g)",
       colour = "Barnacle(s) present?",
       fill = "Barnacle(s) present?")

#ggsave("figures/Figure5A.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#with no outlier
mod_allometry_nooutlier <- glmmTMB(wet_wgt_snail ~ barnacled * scale_length +
                           (1|pond),
                         data = snail_movement_nooutliers)
#predict 
mod_allo_predict_nooutlier <- ggpredict(mod_allometry_nooutlier, 
                              terms=c("scale_length", "barnacled")) %>% 
  #convert scaled distances back to unscaled
  mutate(length = x*sd(snail_movement_nooutliers$length) +
           mean(snail_movement_nooutliers$length),
         barnacled = group)
#plot predictions on raw data
ggplot(mod_allo_predict_nooutlier, aes(length, predicted, 
                             group = barnacled, fill = barnacled)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(colour = barnacled)) +
  geom_point(data = snail_movement_nooutliers, 
             aes(length, wet_wgt_snail, colour = barnacled),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Snail length (mm)",
       y = "Snail wet weight (g)",
       colour = "Barnacle(s) present?",
       fill = "Barnacle(s) present?")

#ggsave("figures/Figure5B.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)

#Q6: Do recapture rates vary depending on whether barnacles are present?-----
mod_recap <- glmmTMB(recap_binary ~ barnacled * date + 
                       (1|pond), 
                     data = snail_recap_long, 
                     family = binomial)

mod_recap_predict <- ggpredict(mod_recap, 
                               terms=c("date", "barnacled")) %>% 
  rename(date = x,
         barnacled = group)
#plot predictions on raw data
ggplot(mod_recap_predict, aes(date, predicted,  colour = barnacled)) + 
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymax = conf.high, ymin = conf.low), 
                position = position_dodge(width = 0.8),
                width = 0.1) +
  geom_point(data = snail_recap_long, 
             aes(date, recap_binary, group = fct_rev(barnacled)),
             size = 0.5, position = position_jitterdodge(dodge.width = 0.8,
                                                         jitter.height = 0.05)) + 
  theme_paper() +
  labs(x = "Month",
       y = "Probability of recapture",
       colour = "Barnacle(s) present?",
       fill = "Barnacle(s) present?")

#ggsave("figures/Figure6.png", device = "png",
#       height = 120, width = 180, units = "mm", dpi = 600)
