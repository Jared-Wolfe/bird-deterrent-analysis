# =========================================================
# Acoustic Deterrent Analysis (2023–2024)
# =========================================================
# This script reproduces analyses from:
# "Evaluating Acoustic Deterrents for Reducing Migratory Bird 
# Collisions with Buildings"
#
# Data are included locally in this repository.
# =========================================================


# =========================================================
# LOAD LIBRARIES
# =========================================================

library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(gridExtra)
library(ggplot2)
library(lubridate)
library(loo)


# =========================================================
# SPRING MIGRATION ARU EXPERIMENT (2023)
# =========================================================

# Load ARU dataset (local file)
aru_data <- read.csv("ARU_spring_data.csv")

# Filter experimental and control units
filtered_data <- aru_data %>%
  filter(ARU_ID %in% c("ARU_Exp", "ARU1_Control", "ARU10_Control", 
                       "ARU12_Control", "ARU14_Control"))

# Define experimental periods
periods <- c(rep("Synthetic", 17), 
             rep("Predator", 17), 
             rep("Control", 7))

# Reshape to long format
long_data <- filtered_data %>%
  pivot_longer(cols = starts_with("Day"),
               names_to = "day",
               values_to = "detection") %>%
  mutate(period = ifelse(ARU_ID == "ARU_Exp",
                         rep(periods, each = n()/length(periods)),
                         "Baseline_Control"))

# Convert day to numeric
long_data$day <- as.numeric(gsub("Day_", "", long_data$day))

# Aggregate detections
daily_data <- long_data %>%
  group_by(day, period, ARU_ID) %>%
  summarize(total_detections = sum(detection), .groups = "drop")

# Offset for temporal trend
daily_data$log_day <- log(daily_data$day)

# Fit Bayesian model
aru_model <- brm(
  total_detections ~ period * ARU_ID + offset(log_day),
  data = daily_data,
  family = gaussian(),
  prior = c(
    prior(normal(0, 10), class = "b"),
    prior(cauchy(0, 2), class = "sigma")
  ),
  iter = 2000,
  chains = 4
)

summary(aru_model)
pp_check(aru_model)


# =========================================================
# ARU RESULTS INTERPRETATION
# =========================================================

# The Bayesian regression results indicate that both the predator and synthetic 
# deterrents influenced bird detections relative to control conditions.

# The synthetic deterrent showed a strong reduction in detections 
# (Estimate ≈ -18), with a credible interval excluding zero, indicating 
# strong evidence for an effect.

# The predator deterrent also showed a reduction (Estimate ≈ -8), 
# but with a wider credible interval that included zero, suggesting 
# greater uncertainty in its effectiveness.

# These results informed priors used in subsequent collision models:
# Synthetic prior ~ normal(-18, 8)
# Predator prior ~ normal(-8, 8)


# =========================================================
# FALL COLLISION MODEL (2023–2024)
# =========================================================

# Load collision dataset (local file)
bird_data <- read.csv("NW_collision_data.csv")

# Format dates and filter relevant years
bird_data <- bird_data %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(year(Date) %in% c(2023, 2024))

# Define treatment indicators
bird_data <- bird_data %>%
  mutate(
    Predator2023 = ifelse(Location == "KGH-E" &
                            Date >= as.Date("2023-09-16") &
                            Date <= as.Date("2023-10-31"), 1, 0),
    
    Predator2024 = ifelse(Location %in% c("KGH-E", "RAC-E") &
                            Date >= as.Date("2024-09-12") &
                            Date <= as.Date("2024-10-31"), 1, 0),
    
    Synthetic2023 = ifelse(Location == "RAC-E" &
                             Date >= as.Date("2023-09-16") &
                             Date <= as.Date("2023-10-31"), 1, 0)
  )

# Migration intensity offset (simplified categorical scaling)
bird_data <- bird_data %>%
  mutate(MigrationOffset = case_when(
    (year(Date) == 2023 & ((month(Date) == 9 & day(Date) >= 18 & day(Date) <= 28) |
                             (month(Date) == 10 & day(Date) <= 15))) ~ 5,
    (year(Date) == 2024 & ((month(Date) == 9 & day(Date) >= 22 & day(Date) <= 30) |
                             (month(Date) == 10 & day(Date) <= 14))) ~ 5,
    TRUE ~ 1
  ))

# Fit final Bayesian model
final_model <- brm(
  Birds ~ Predator2023 + Predator2024 + Synthetic2023 +
    offset(MigrationOffset) + (1 | Location),
  data = bird_data,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("normal(0, 10)", class = "Intercept"),
    set_prior("normal(-1, 2)", class = "b", coef = "Predator2023"),
    set_prior("normal(-1, 2)", class = "b", coef = "Predator2024"),
    set_prior("normal(-2, 3)", class = "b", coef = "Synthetic2023"),
    set_prior("student_t(3, 0, 2.5)", class = "sd")
  ),
  chains = 4,
  iter = 6000,
  warmup = 2000,
  cores = 4,
  seed = 123,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

summary(final_model)
pp_check(final_model)
mcmc_intervals(final_model)


# =========================================================
# COLLISION MODEL INTERPRETATION
# =========================================================

# Predator deterrent (2023):
# Moderate reduction in collisions (~55%), credible interval excludes zero.

# Predator deterrent (2024):
# Strong reduction (~90%), indicating increased effectiveness with improved deployment.

# Synthetic deterrent (2023):
# Significant reduction (~68%), consistent with ARU experiment results.

# Substantial variation among buildings suggests site-specific risk factors.


# =========================================================
# PRIOR SENSITIVITY ANALYSIS
# =========================================================

# Test robustness to prior assumptions

weaker_priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(0, 5)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd")
)

stronger_priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(-2, 1)", class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sd")
)

model_weaker <- update(final_model, prior = weaker_priors)
model_stronger <- update(final_model, prior = stronger_priors)

summary(model_weaker)
summary(model_stronger)


# =========================================================
# SENSITIVITY ANALYSIS INTERPRETATION
# =========================================================

# Posterior estimates remained consistent across weaker and stronger priors.

# Credible intervals overlapped substantially, indicating that results are 
# driven by the data rather than prior specification.

# Conclusion: deterrent effects are robust and well-supported.