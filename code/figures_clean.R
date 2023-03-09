library(tidyverse) #for plotting and general data wrangling
library(ggeffects) #for plotting model output
library(patchwork) #for combining panels
library(glmmTMB) #for running final models
source("code/ggplot_paper_theme.R")
pal <- c("#71a89b", "#405e57")

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



#Q1A: Is there a relationship between snail abundance per quadrat and ---- 
#distance from the water?
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
fig1a <- ggplot(total_num_predict_out, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_quad_nooutliers, 
             aes(distance_from_water, total_num),
             size = 0.5) +
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Number of snails",
       title = "A") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 1, size = 12))
fig1a

#Q1B: Is there a relationship between snail length and distance from the --------
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
fig1b <- ggplot(snail_length_predict, aes(distance_from_water, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line() +
  geom_point(data = snail_individuals, 
             aes(distance_from_water, length),
             size = 0.5) + 
  theme_paper() +
  labs(x = "Distance from water (m)",
       y = "Snail length (mm)",
       title = "B") +
  theme(plot.title = element_text(hjust = -0.05, vjust = 1, size = 12))

fig1a/fig1b
#ggsave("figures/Figure1.png", device = "png",
#       height = 180, width = 135, units = "mm", dpi = 600)
#ggsave("figures/Figure1.pdf",
#       height = 180, width = 135, units = "mm", dpi = 600)

#Q2: Does the probability of being barnacled change with distance from the------ 
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

#ggsave("figures/Figure2.png", device = "png",
#       height = 90, width = 135, units = "mm", dpi = 600)
#ggsave("figures/Figure2.pdf",
#       height = 90, width = 135, units = "mm", dpi = 600)

#Q3: Does the weight of barnacles impact snail speed--------------------------
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
fig3A <- ggplot(hurdle3_predict, aes(wet_wgt_barnacles, predicted)) + 
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
fig3B <- ggplot(hurdle4_predict, aes(wet_wgt_barnacles, predicted)) + 
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

fig3A + fig3B

#ggsave("figures/Figure3.png", device = "png",
#       height = 90, width = 150, units = "mm", dpi = 600)
#ggsave("figures/Figure3.pdf",
#       height = 90, width = 150, units = "mm", dpi = 600)

#to present the results as percent declines for each unit increase of wet 
#weight:
#vals <- NULL
#for (i in 0:4){
#  scaled <- (i - mean(snail_movement$wet_wgt_barnacles))/
#    sd(snail_movement$wet_wgt_barnacles)
#  vals[[i+1]] <- scaled
#}
#
#hurdle3_predict_subsample <- ggpredict(hurdle_mod3, 
#                             terms=c("scale_wgt_barnacles [vals]")) %>% 
#  mutate(scale_wgt_barnacles = x,
#         wet_wgt_barnacles = x*sd(snail_movement$wet_wgt_barnacles) +
#           mean(snail_movement$wet_wgt_barnacles))
#percent_changes <- NULL
#for(i in 1:4){
#  percent_changes[i] <- (hurdle3_predict_subsample$predicted[i+1]- 
#    hurdle3_predict_subsample$predicted[i])/ 
#    hurdle3_predict_subsample$predicted[i]
#}


#Q4: Do recapture rates vary depending on whether barnacles are present?-----
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
  scale_colour_manual(values = rev(pal)) +
  scale_fill_manual(values = rev(pal)) +
  theme_paper() +
  labs(x = "Month",
       y = "Probability of recapture",
       colour = "Barnacle(s)\npresent?",
       fill = "Barnacle(s)\npresent?")

#ggsave("figures/Figure4.png", device = "png",
#       height = 90, width = 135, units = "mm", dpi = 600)
#ggsave("figures/Figure4.pdf", 
#       height = 90, width = 135, units = "mm", dpi = 600)

#QS1: Dos the presence of barnacles affect snail speed?-------------------------
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
figS1.A <- ggplot(hurdle1_predict, aes(barnacled, predicted)) + 
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
figS1.B <-  ggplot(hurdle2_predict, aes(barnacled, predicted)) + 
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

figS1.A + figS1.B

#ggsave("figures/FigureS1.png", device = "png",
#       height = 90, width = 150, units = "mm", dpi = 600)
#ggsave("figures/FigureS1.pdf",
#       height = 90, width = 150, units = "mm", dpi = 600)

#QS2: Does the number of barnacles impact snail speed--------------------------
hurdle_mod5 <- glmmTMB(distance_logistic ~ scale_wgt_ratio + scale_length + (1|pond), 
                       family = binomial(link = "logit"),
                       data = snail_movement)

hurdle_mod6 <- glmmTMB(distance ~ scale_wgt_ratio + scale_length + (1|pond),
                       family = Gamma(link = "log"),
                       data = snail_movement_gamma)

#predict 
hurdle5_predict <- ggpredict(hurdle_mod5, 
                             terms=c("scale_wgt_ratio [n=100]")) %>% 
  mutate(scale_wgt_ratio = x,
         wgt_ratio = x*sd(snail_movement$wgt_ratio) +
           mean(snail_movement$wgt_ratio))
#plot predictions on raw data
figS2A <- ggplot(hurdle5_predict, aes(wgt_ratio, predicted)) + 
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
figS2B <- ggplot(hurdle6_predict, aes(wgt_ratio, predicted)) + 
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

figS2A + figS2B

#ggsave("figures/FigureS2.png", device = "png",
#       height = 90, width = 150, units = "mm", dpi = 600)
#ggsave("figures/FigureS2.pdf", 
#       height = 90, width = 150, units = "mm", dpi = 600)
#QS3: Do barnacles impact snail allometry---------------------------------------
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
  scale_colour_manual(values = rev(pal)) +
  scale_fill_manual(values = rev(pal)) +
  theme_paper() +
  labs(x = "Snail length (mm)",
       y = "Snail wet weight (g)",
       colour = "Barnacle(s)\npresent?",
       fill = "Barnacle(s)\npresent?")

#ggsave("figures/FigureS3.png", device = "png",
#       height = 90, width = 135, units = "mm", dpi = 600)
#ggsave("figures/FigureS3.pdf", 
#       height = 90, width = 135, units = "mm", dpi = 600)
