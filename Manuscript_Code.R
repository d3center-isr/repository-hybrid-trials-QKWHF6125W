##################################
# Proximal and Distal Model Code #
##################################

# This file performs the proximal and distal model analysis for a hybrid SMART-MRT. It uses the following inputs:
  # finaldata_w_resp: data for the proximal model including responders and nonresponders
  # distal_nonresp: data for the distal model including only nonresponders
  # abar_plot: prediction data plotting distal model outcomes for values of abar at the mean, one standard deviation below the mean, and one standard deviation above the mean
# Tables and figures are numbered to match those in the manuscript
# All output tables are formatted using Kable

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "rootSolve", "kableExtra", "knitr", "geepack", "data.table")

# Set environments
path_estimators = Sys.getenv("path_estimators")
path_data = Sys.getenv("path_data")

# Import clean data
finaldata_w_resp = read.csv(file.path(path_data, "finaldata_w_resp.csv")) # proximal model data: responders and nonresponders
distal_nonresp = read.csv(file.path(path_data, "data_distal_nonresp.csv")) # distal model 2 data: nonresponders only
abar_plot = read.csv(file.path(path_data, "abar_plot.csv")) # prediction data plotting distal model outcomes for various abar values

# Source for estimators
source(file.path(path_estimators, "Hanna_Estimator_resp.R"))

#####
# Proximal Model
#####

# Make time.classified = 0 for responders (this avoids issues of colinearity)
finaldata_w_resp1 = within(finaldata_w_resp, time.classified[nonresponder == 0] <- 0)
finaldata_w_resp1$owner = factor(finaldata_w_resp1$owner)

# Control variable regression
resp_model_control = geeglm(outcome ~ a1 + a2 + time.classified + sex + a2days + bmi.bl, data = finaldata_w_resp1, id = owner)

# Control variable output table
resp_model_control_table = data.table(tidy(resp_model_control))
cilow = data.table(resp_model_control_table$estimate - 1.96 * resp_model_control_table$std.error)
cihigh = data.table(resp_model_control_table$estimate + 1.96 * resp_model_control_table$std.error)
resp_model_control_table = data.frame(cbind(resp_model_control_table[,2:3], cilow, cihigh, resp_model_control_table[,4:5]))
names(resp_model_control_table) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "Wald", "Pr>|W|")
rownames(resp_model_control_table) = c("Intercept", "$Z_{i1}$", "$Z_{i2}$", "Week Rerandomized", "Sex", "Days Since Rerandomization", "Baseline BMI")
resp_model_control_table = round(resp_model_control_table, digits = 4)
resp_model_control_table %>% kbl(caption = "Table 8.1: Control Variables") %>% kable_classic(full_width = F, html_font = "Times New Roman") %>% footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit")

# Add variable for interaction of a1 and a2
finaldata_w_resp1 = finaldata_w_resp %>%
  mutate(a1_a2 = a1*a2)

# Proximal model
proximal = binary_outcome_moderated_effect_resp(
  dta = finaldata_w_resp1,
  control_var = c("a1", "a2", "time.classified", "sex", "a2days", "bmi.bl"),
  moderator = c("a1", "a2", "a1_a2"),
  id_var = "owner",
  day_var = "studyday.0start",
  trt_var = "intervention",
  outcome_var = "outcome",
  avail_var = NULL,
  prob_treatment = "prob",
  significance_level = 0.05)

# Proximal model output table
proximal_model = as.data.frame(proximal[1])
names(proximal_model) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "T", "Pr>|T|")
rownames(proximal_model) = c("Intercept", "$Z_{i1}$", "$Z_{i2}$", "$Z_{i1} Z_{i2}$")
proximal_round = round(proximal_model, digits = 4)
proximal_round %>% kbl(caption = "Table 8.2: Results (Parameter Estimates) for Proximal Model") %>% kable_classic(full_width = F, html_font = "Times New Roman") %>% footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit")

#####
# Distal Model
#####

##### Distal Model : nonresponders only #####

# Distal model 
model_distal = geeglm(weight_bl_lbs - weight_6mo_lbs ~ a1 * a2 * abar_cnr + sex_cnr + bmi_bl_cnr, data = distal_nonresp, id = id, subset = rerand == 1)

# Distal model output table
model_distal_table = data.table(tidy(model_distal))
cilow = data.table(model_distal_table$estimate - 1.96 * model_distal_table$std.error)
cihigh = data.table(model_distal_table$estimate + 1.96 * model_distal_table$std.error)
model_distal_table = data.frame(cbind(model_distal_table[,2:3], cilow, cihigh, model_distal_table[,4:5]))
names(model_distal_table) = c("Estimate", "Robust SE", "95% CI LL", "95% CI UL", "Wald", "Pr>|W|")
rownames(model_distal_table) = c("Intercept", "$Z_{i1}$", "$Z_{i2}$", "$\\bar{A_i}^{(2)}$", "Sex", "Baseline BMI", "$Z_{i1} Z_{i2}$", "$Z_{i1} \\bar{A_i}^{(2)}$", "$Z_{i2}\\bar{A_i}^{(2)}$", "$Z_{i1}Z_{i2}\\bar{A_i}^{(2)}$")
model_distal_table = round(model_distal_table, digits = 4)
model_distal_table %>% kbl(caption = "Table 9: Results (Parameter Estimates) for Distal Model") %>% kable_classic(full_width = F, html_font = "Times New Roman") %>% footnote(general = "CI: confidence interval; LL: lower limit; UL: upper limit; all covariates centered")

##### Plot #####

ggplot(abar_plot, aes(x = factor(group))) +
  geom_point(aes(y = mean, color = "red", shape = '15'), size = 2) +
  geom_point(aes(y = sd_above, color = "green", shape = '16'), size = 2) +
  geom_point(aes(y = sd_below, color = "blue", shape = '17'), size = 2) +
  scale_colour_manual(name = 'Rate of message delivery', values =c('red'='red','green'='green', 'blue'='blue'), labels = c('Mean (0.6435)','1 SD Above Mean (0.7437)', '1 SD Below Mean (0.5432)')) +
  scale_shape_manual(name = 'Rate of message delivery', values =c('15'=15, '16'=16, '17'=17), labels = c('Mean (0.6435)','1 SD Above Mean (0.7437)', '1 SD Below Mean (0.5432)')) +
  # guides(colour = guide_legend(override.aes = list(shape = c(15,16,17)))) +
  labs(title = str_wrap("Figure 4: Estimated weight loss by the initial options, subsequent options and rate of message delivery", 60), x = "Sequence of initial and subsequent options", y = "Estimated weight loss") + 
  ylim(0, 10) 