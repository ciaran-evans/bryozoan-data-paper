# Installing packages
# uncomment if these are not installed
# install.packages("tidyverse")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("RLRsim")
# install.packages("latex2exp")
# install.packages("patchwork")
# install.packages("xtable")

library(tidyverse)
library(lme4)
library(lmerTest)
library(RLRsim)
library(latex2exp)
library(patchwork)
library(xtable)

# importing and organizing the raw data
# Here we do initial data tidying

metab <- read_csv("bryozoan_raw.csv")
metab_b <- metab[,1:4]
metab_w <- metab[,5:8]
colnames(metab_b) <- c("Stage", "Run", "Mass", "Metabolic")
colnames(metab_w) <- c("Stage", "Run", "Mass", "Metabolic")
metab_b <- metab_b %>%
  mutate(Species = "bugula")
metab_w <- metab_w %>%
  mutate(Species = "watersipora")

bryozoan <- metab_b %>%
  rbind(metab_w) %>%
  drop_na() %>%
  mutate(Mass = 10^Mass,
         Metabolic = 10^Metabolic,
         Stage = tolower(Stage)) %>%
  dplyr::select(Species, Run, Stage, Mass, Metabolic) %>%
  arrange(Species, Run, Stage) %>%
  mutate(Run = cumsum(c(1, diff(Run) != 0)))

### uncomment to save the data
# write_csv(bryozoan, file="bryozoan_data.csv")


# Section 3

### Number of runs
bryozoan %>%
  pull(Run) %>%
  unique() %>%
  length()

### Number of rows in the data
nrow(bryozoan)

# Section 4.1

### Number of individuals
bryozoan %>%
  filter(Stage != "late") %>%
  nrow()

### check the same number of individuals is 
### measured in early and late stage for each run
### (result is empty, so the numbers are always the same)
bryozoan %>%
  filter(Stage != "larvae") %>%
  group_by(Species, Run) %>%
  summarize(num_early = sum(Stage == "early"),
            num_late = sum(Stage == "late")) %>%
  ungroup() %>%
  filter(num_early != num_late)


# Section 4.2

## check whether mass and metabolic rate were measured once
## or twice

### For the Watersipora: there are no NAs after the full join,
### so each mass is measured once
bryozoan %>%
  filter(Species == "watersipora",
         Stage == "early") %>%
  full_join(bryozoan %>%
              filter(Species == "watersipora",
                     Stage == "late"),
            by = c("Species", "Run", "Mass")) %>%
  is.na() %>%
  sum()

### For the bugula:
### Lots of NAs
bryozoan %>%
  filter(Species == "bugula",
         Stage == "early") %>%
  full_join(bryozoan %>%
              filter(Species == "bugula",
                     Stage == "late"),
            by = c("Species", "Run", "Mass")) %>%
  is.na() %>%
  sum()

### NAs are 4*number of early bugula. So what we 
### would expect if each mass was measured twice
bryozoan %>%
  filter(Species == "bugula",
         Stage == "early") %>%
  nrow()

## Figure 1
### distributions of mass and metabolic rate
### (before fixing the error)
p1 <- bryozoan %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early", "late")) %>%
  ggplot(aes(x = Stage, y = Mass)) +
  geom_boxplot(lwd=0.7) +
  facet_wrap(~Species) +
  theme_bw() +
  labs(y = "Mass (micrograms)")

p2 <- bryozoan %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early", "late")) %>%
  ggplot(aes(x = Stage, y = Metabolic)) +
  geom_boxplot(lwd=0.7) +
  facet_wrap(~Species) +
  theme_bw() +
  labs(y = "Metabolic rate (mJ/hour)")

### Uncomment to save the plot
# pdf(file = "bryozoan_eda_1.pdf", width=10, height=4)
p1 + p2
# dev.off()


## fixing errors
bugula_early_mass <- bryozoan %>%
  filter(Species == "bugula",
         Stage == "early") %>%
  pull(Mass)

bryozoan$Mass[bryozoan$Mass < 1] <- bugula_early_mass

## Uncomment to save the corrected data
# write_csv(bryozoan, file="bryozoan_data_fixed.csv")



# Section 4.3

## Figure 2
### relationship between mass and metabolic rate

### Uncomment to save the plot
# pdf(file = "bryozoan_eda_2.pdf", width=9, height=4)
bryozoan %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early", "late")) %>%
  ggplot(aes(x = Mass, y = Metabolic, color = Stage, 
             shape = Stage)) +
  geom_point(alpha = 0.8, size=1.5) +
  geom_smooth(se=F, method="lm") +
  facet_wrap(~Species) +
  theme_bw() +
  labs(x = "Mass (micrograms)",
       y = "Metabolic rate (mJ/hour)") +
  guides(color = guide_legend(override.aes = list(linetype = 0)))
# dev.off()



## Figure 3
### differences by run

### Uncomment to save the plot
# pdf(file = "bryozoan_eda_3.pdf", width=10, height=4)
bryozoan %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early", "late")) %>%
  ggplot(aes(x = as.factor(Run), y = Metabolic)) +
  geom_boxplot() +
  facet_grid(Species~Stage) +
  theme_bw() +
  labs(x = "Run", y = "Metabolic rate (mJ/hour)")
# dev.off()


# Section 5.1

## subsetting the data to focus only on early-stage bugula

bugula_early <- bryozoan %>%
  filter(Species == "bugula",
         Stage == "early")


# Section 5.2

## Fitting a simple linear regression with early-stage bugula
## (untransformed data)

be_lm <- lm(Metabolic ~ Mass, data = bugula_early)

## Figure 4
### Diagnostic plots for be_lm

p1 <- bugula_early %>%
  mutate(residuals = residuals(be_lm),
         predicted = predict(be_lm)) %>%
  ggplot(aes(x = predicted, y = residuals)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0, 
              color="blue", lwd=1.2) +
  labs(x = "Predicted metabolic rate (mJ/hour)",
       y = "Residuals") +
  theme_bw()

p2 <- bugula_early %>%
  mutate(residuals = residuals(be_lm)) %>%
  ggplot(aes(sample = residuals)) +
  geom_qq() +
  geom_qq_line() +
  labs(x = "Theoretical normal quantiles",
       y = "Observed residual quantiles") +
  theme_bw()

### uncomment to save the plot
# pdf(file = "bryozoan_diagnostics_1.pdf", width=9, height=4)
p1 + p2
# dev.off()

### summary of fitted model
summary(be_lm)

### p-value for hypothesis test
### H0: beta1 = 0   HA: beta1 > 0
pt(4.946, df = 195, lower.tail=F)

### 95% confidence interval for the slope

0.006427 - qt(0.025, 195, lower.tail=F)*0.0013
0.006427 + qt(0.025, 195, lower.tail=F)*0.0013




# Section 6: log transformations

## Figure 5
### motivation for log-transforming both predictor and response
p1 <- data.frame(x = 1:100, y = 1:100) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  theme_bw() +
  labs(x = "Mass", 
       y = "Metabolism",
       title = TeX("Efficiency is the same ($\\beta = 1$)")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p2 <- data.frame(x = 1:100, y = (1:100)^2) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  theme_bw() +
  labs(x = "Mass", 
       y = "Metabolism",
       title = TeX("Larger individuals are less efficient ($\\beta > 1$)")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p3 <- data.frame(x = 1:100, y = (1:100)^0.25) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  theme_bw() +
  labs(x = "Mass", 
       y = "Metabolism",
       title = TeX("Larger individuals are more efficient ($0 < \\beta < 1$)")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

### uncomment to save the plot
# pdf(file="transformation_motivation.pdf", width = 12, height=4)
p1 + p2 + p3
# dev.off()


## fitting transformed model

### creating log-transformed variables
bugula_early <- bugula_early %>%
  mutate(log_mass = log(Mass),
         log_metabolic = log(Metabolic))


be_lm <- lm(log_metabolic ~ log_mass,
            data = bugula_early)

summary(be_lm)

## confidence interval for log-log slope

0.5991 - qt(0.025, 195, lower.tail=F)*0.1162
0.5991 + qt(0.025, 195, lower.tail=F)*0.1162

## p-value
### H0: beta = 1   HA: beta < 1
(0.5991 - 1)/0.1162
pt(-3.45, 195)



# Section 7: Multiple regression

## Creating log-transformed variables for the full data
bryozoan <- bryozoan %>%
  mutate(log_mass = log(Mass),
         log_metabolic = log(Metabolic))

## excluding late-stage measurements to avoid dependence issues
bryozoan_larvae_early <- bryozoan %>%
  filter(Stage != "late")

## Figure 6
### log metabolic rate vs. log mass for larval and early stages

### Uncomment to save plot
# pdf(file = "bryozoan_eda_4.pdf", width=9, height=4)
bryozoan_larvae_early %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early")) %>%
  ggplot(aes(x = log_mass, 
             y = log_metabolic,
             color = Stage,
             shape = Stage)) +
  geom_point(alpha = 0.8, size=1.5) +
  geom_smooth(se=F, method="lm") +
  facet_wrap(~Species) +
  theme_bw() +
  labs(x = "log(Mass)",
       y = "log(Metabolic rate)") +
  guides(color = guide_legend(override.aes = list(linetype = 0)))
# dev.off()

## Re-level the Stage variable so larval is baseline
## This is a little more intuitive, since larval comes
## chronologically before early stage
bryozoan_larvae_early <- bryozoan_larvae_early %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early"))

## Fitting model in Equation (5)
ble_lm <- lm(log_metabolic ~ Stage*Species + log_mass, 
             data = bryozoan_larvae_early)

summary(ble_lm)

## Check that residual plot for the multiple regression model looks reasonable
## (not shown in paper)
bryozoan_larvae_early %>%
  mutate(pred = predict(ble_lm),
         resid = residuals(ble_lm)) %>%
  ggplot(aes(x = pred, y = resid,
             color = Species)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0,
              color = "blue", lwd = 1.2) +
  labs(x = "Predicted log(metabolic rate)",
       y = "Residuals") +
  theme_bw()


## confidence interval

0.577 - qt(0.025, 563, lower.tail=F)*0.073
0.577 + qt(0.025, 563, lower.tail=F)*0.073

## test
### H0: beta = 1   HA: beta < 1

pt((0.577 - 1)/0.073, 563)


## does the slope vary? Answer appears to be no. Compare against a
## full model with all interaction terms

ble_lm_full <- lm(log_metabolic ~ Stage*Species*log_mass, 
                  data = bryozoan_larvae_early)

anova(ble_lm, ble_lm_full)




# Section 8: Mixed effects models

## simple random intercept model (Equation (7)) with REML
ble_lme_simple <- lmer(log_metabolic ~ (1|Run), 
                       data = bryozoan_larvae_early)

## Do we need to include Run in the model? This is a comparison
## between the random intercepts model, and an intercept-only linear model

### intercept-only linear model
ble_lm_null <- lm(log_metabolic ~ 1, 
                       data = bryozoan_larvae_early)

### Using likelihood ratio test with the anova function
### Possible issue: conservative because variance = 0 is on the edge
### of the parameter space
anova(ble_lme_simple, ble_lm_null)

### Using a simulation-based LRT from RLRsim
exactRLRT(ble_lme_simple)



## Incorporating Species, Stage, and mass (Equation (8))

ble_lme <- lmer(log_metabolic ~ (1|Run) + Stage*Species + 
                  log_mass, data = bryozoan_larvae_early)

summary(ble_lme)

### Testing whether the random effect variance is 0
exactRLRT(ble_lme)

### Figure 7
### checking assumptions for the model in Equation (8)

p1 <- bryozoan_larvae_early %>%
  mutate(pred = predict(ble_lme),
         resid = residuals(ble_lme)) %>%
  ggplot(aes(x = pred, y = resid,
             color = Species, shape = Species)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_abline(slope = 0, intercept = 0,
              color = "blue", lwd = 1.2) +
  labs(x = "Predicted log(metabolic rate)",
       y = "Residuals") +
  theme_bw()

p2 <- bryozoan_larvae_early %>%
  mutate(resid = residuals(ble_lme)) %>%
  ggplot(aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  labs(x = "Theoretical normal quantiles",
       y = "Observed residual quantiles") +
  theme_bw()

p3 <- data.frame(sample = coef(ble_lme)$Run[,1]) %>% 
  ggplot(aes(sample = sample)) + 
  geom_qq() + 
  geom_qq_line() +
  labs(x = "Theoretical normal quantiles",
       y = "Estimated random effect quantiles") +
  theme_bw()

### Uncomment to save the plot
# pdf(file="mixed_model_diagnostics.pdf", width=12, height=3)
p1 + p2 + p3
# dev.off()

### t confidence interval
0.67612 - qt(0.025, df=560.37, lower.tail = F)*0.07411
0.67612 + qt(0.025, df=560.37, lower.tail = F)*0.07411


### parametric bootstrap CI (simple percentile interval)
boot_fun <- function(fitted_mod){
  return(coef(fitted_mod)$Run$log_mass[1])
}

param_boot <- bootMer(ble_lme, boot_fun, nsim = 10000,
                      seed = 3, type = "parametric")

quantile(param_boot$t, 0.025)
quantile(param_boot$t, 0.975)


### Hypothesis testing
### H0: beta = 1   HA: beta < 1
(0.67612 - 1)/0.07411



## adding random slopes (Equation (9))

ble_lme_2 <- lmer(log_metabolic ~ Stage*Species + 
                  log_mass + (log_mass|Run), 
                  data = bryozoan_larvae_early)

summary(ble_lme_2)

### Making a nice table of coefficients (for Table 1)
summary(ble_lme_2)$coefficients[,c(1,2,4)] %>%
  xtable()

### Parametric bootstrap confidence interval for the slope
### (simple percentile interval)
boot_fun <- function(fitted_mod){
  return(mean(coef(fitted_mod)$Run$log_mass))
}

param_boot <- bootMer(ble_lme_2, boot_fun, nsim = 10000,
                      seed = 3, type = "parametric")

quantile(param_boot$t, 0.025)
quantile(param_boot$t, 0.975)


# Not included in paper: incorporating late-stage bugula
## Now what if we wanted to include the late-stage 
## bugula in our model too? 
## One way to do that could be to include a random effect 
## term for each individual

## first we need to make a term for the individual id

bryozoan_id <- bryozoan %>%
  filter(Stage != "larvae") %>%
  arrange(Species, Run, Mass, Stage) %>%
  mutate(id = cumsum(c(1, diff(Mass) != 0)) + 
           nrow(bryozoan %>% filter(Stage == "larvae"))) %>%
  rbind(bryozoan %>%
          filter(Stage == "larvae") %>%
          mutate(id = 1:n())) %>%
  arrange(Species, Run, Stage)

max(bryozoan_id$id)
max(table(bryozoan_id$id))


bryozoan_lme <- lmer(log_metabolic ~ (1|Run/id) + 
                       Stage*Species + log_mass, 
                     data = bryozoan_id)

bryozoan_lme_full <- lmer(log_metabolic ~ (1|Run/id) + 
                       Stage*Species*log_mass, 
                     data = bryozoan_id)

anova(bryozoan_lme, bryozoan_lme_full)



bryozoan %>%
  mutate(Stage = fct_relevel(Stage, "larvae", "early", "late")) %>%
  ggplot(aes(x = log_mass, 
             y = log_metabolic,
             color = Stage,
             shape = Stage)) +
  geom_point() +
  geom_smooth(se=F, method="lm") +
  facet_wrap(~Species) +
  theme_bw() +
  labs(x = "log(Mass)",
       y = "log(Metabolic rate)",
       color = "Stage", shape = "Stage")
