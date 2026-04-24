#################################################################################################
#################################################################################################
#################################################################################################
############################  SPRING MIGRATION ARU EXPERIMENT IN 2023 ###########################
#################################################################################################
#################################################################################################
#################################################################################################

# Loading data
aru_data <- read.csv("https://drive.google.com/uc?export=download&id=1CJQx_Oacr0KA2XccK7sZhO5HNNLdp_o0", header = TRUE, skip = 0)

# Load necessary libraries
# dplyr and tidyr for data manipulation, brms for Bayesian modeling
# bayesplot and gridExtra for visualization, ggplot2 for plotting

library(dplyr)
library(tidyr)
library(brms)
library(bayesplot)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(loo)    # For model comparison using LOO-CV

# Filter data to include only the ARU_Exp and control ARUs
filtered_data <- aru_data %>%
  filter(ARU_ID %in% c("ARU_Exp", "ARU1_Control", "ARU10_Control", "ARU12_Control", "ARU14_Control"))

# Define periods for ARU_Exp
dates <- names(filtered_data)[5:ncol(filtered_data)]  # Adjusted for the new dataset structure
periods <- c(rep("Synthetic", 17), 
             rep("Predator", 17), 
             rep("Control", 7))

# Reshape the data to long format and add period information
long_data <- filtered_data %>%
  gather(day, detection, -Score, -Label, -Migrant_songbird, -ARU_ID) %>%
  mutate(period = ifelse(ARU_ID == "ARU_Exp", rep(periods, each = n()/length(periods)), "Baseline_Control"))

# Ensure 'day' is numeric
# Remove the 'Day_' prefix if it exists, and convert the day column to numeric
long_data$day <- as.numeric(gsub("Day_", "", long_data$day))

# Aggregate data by day, period, and ARU_ID
daily_data <- long_data %>%
  group_by(day, period, ARU_ID) %>%
  summarize(total_detections = sum(detection),
            .groups = 'drop')  # This drops the grouping after summarizing

# Add log of day as an offset
# Now that 'day' is numeric, apply the log transformation
daily_data$log_day <- log(daily_data$day)

# Define the model with ARU_ID as an additional factor
model_formula <- bf(total_detections ~ period * ARU_ID + offset(log_day))

# Fit the Bayesian model using brm() from the brms package
fit <- brm(model_formula, data = daily_data, family = gaussian(), 
           prior = c(prior(normal(0, 10), class = "b"),
                     prior(cauchy(0, 2), class = "sigma")),
           iter = 2000, chains = 4, cores = 4)

# Summarize the model
summary(fit)

# Posterior predictive checks
pp_check(fit)

# Extract the posterior samples for visualization
posterior_samples <- as_draws_df(fit)

# Rename the parameters in the posterior samples for better labeling in plots
posterior_samples <- posterior_samples %>%
  rename(`Synthetic signal` = b_periodSynthetic, `Predator signal` = b_periodPredator)

# Posterior Density Plots
density_plot <- mcmc_dens_overlay(posterior_samples, pars = c("Predator signal", "Synthetic signal")) +
  ggtitle("Posterior Density Plots") +
  xlab("Effect Size") +
  ylab("Density") +
  theme_minimal()

# Posterior Interval Plots
interval_plot <- mcmc_intervals(posterior_samples, pars = c("Predator signal", "Synthetic signal")) +
  ggtitle("Posterior Interval Plots") +
  xlab("Effect Size") +
  theme_minimal()

# Comparison Plots
comparison_plot <- mcmc_areas(posterior_samples, pars = c("Predator signal", "Synthetic signal")) +
  ggtitle("Comparison of Posteriors") +
  xlab("Effect Size") +
  theme_minimal()

# Combine the plots using grid.arrange
combined_plot <- grid.arrange(density_plot, interval_plot, comparison_plot, ncol = 1)

# Print the combined plot to display it
print(combined_plot)

# The Bayesian regression results indicate that both the "predator" and "synthetic" 
# deterrents influenced bird detections compared to the control period and other 
# control ARUs. Specifically, the "synthetic" deterrent showed a significant 
# reduction in detections (Estimate = -18.13, 95% CI: [-35.02, -1.25]), 
# suggesting it was particularly effective at deterring birds compared to the 
# control. The "predator" deterrent also showed a reduction in detections 
# (Estimate = -8.54, 95% CI: [-25.35, 8.28]), but the credible interval 
# includes zero, indicating uncertainty about its effectiveness relative to the 
# control period. When examining the control ARUs, variability in detections 
# was observed, with some control ARUs like ARU10_Control and ARU12_Control 
# showing increases in detections compared to the baseline ARU_Exp, while 
# ARU1_Control exhibited a significant decrease. Over time, as indicated by 
# the offset (log_day), detections generally declined, aligning with 
# expectations of temporal trends in bird activity. These results were 
# instrumental in updating the priors for future analyses, with the "synthetic" 
# deterrent's prior centered on a stronger negative effect (normal(-18.13, 8.68)) 
# and the "predator" deterrent's prior reflecting a moderate effect with greater 
# uncertainty (normal(-8.54, 8.62)). These informed priors will enhance the 
# precision and relevance of subsequent models by incorporating empirical 
# evidence from the current analysis.



############################################################################################
############################################################################################
############################################################################################
##########################  UPDATE BAYESIAN MODEL ##########################################
############################################################################################
############################################################################################
############################################################################################

# Load required libraries explicitly for Bayesian modeling and data manipulation.
library(brms)       # Bayesian regression modeling (brms interface to Stan).
library(dplyr)      # Data manipulation and transformation.
library(lubridate)  # Robust and explicit handling of date operations.

# Load the bird collision dataset explicitly from your provided Google Drive URL for reproducibility.
bird_data <- read.csv("https://drive.google.com/uc?export=download&id=1ZmaL_Fh4UOmASRqZg0u_4gj1aP5yY0IM")

# Convert the raw 'Date' column explicitly into R date format to facilitate accurate date-based operations.
bird_data <- bird_data %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

# Filter data to include ONLY observations from the years 2023 and 2024.
# This ensures consistency by removing years with additional confounding deterrent treatments (e.g., window films).
bird_data <- bird_data %>%
  filter(year(Date) %in% c(2023, 2024))

# Create binary indicator variable ('Predator2023') set to 1 if Predator deterrent was active at "KGH-E" in Fall 2023; otherwise, explicitly set to 0.
bird_data <- bird_data %>%
  mutate(Predator2023 = ifelse(Location == "KGH-E" &
                                 Date >= as.Date("2023-09-16") &
                                 Date <= as.Date("2023-10-31"), 1, 0))

# Create binary indicator variable ('Predator2024') set to 1 if Predator deterrent was active at "KGH-E" or "RAC-E" in Fall 2024; otherwise, explicitly set to 0.
bird_data <- bird_data %>%
  mutate(Predator2024 = ifelse(Location %in% c("KGH-E", "RAC-E") &
                                 Date >= as.Date("2024-09-12") &
                                 Date <= as.Date("2024-10-31"), 1, 0))

# Create binary indicator variable ('Synthetic2023') set to 1 if Synthetic deterrent was active at "RAC-E" in Fall 2023; otherwise, explicitly set to 0.
bird_data <- bird_data %>%
  mutate(Synthetic2023 = ifelse(Location == "RAC-E" &
                                  Date >= as.Date("2023-09-16") &
                                  Date <= as.Date("2023-10-31"), 1, 0))

# Create a numeric 'MigrationOffset' variable controlling explicitly for daily variations in migration intensity (1 = low, 5 = highest).
bird_data <- bird_data %>%
  mutate(MigrationOffset = case_when(
    (year(Date) == 2023 & ((month(Date) == 9 & day(Date) >= 18 & day(Date) <= 28) |
                             (month(Date) == 10 & day(Date) <= 15))) ~ 5,
    (year(Date) == 2024 & ((month(Date) == 9 & day(Date) >= 22 & day(Date) <= 30) |
                             (month(Date) == 10 & day(Date) <= 14))) ~ 5,
    (month(Date) == 9 & day(Date) >= 15 & day(Date) <= 24) |
      (month(Date) == 10 & day(Date) >= 11 & day(Date) <= 15) ~ 4,
    (month(Date) == 9 & day(Date) >= 10 & day(Date) < 15) |
      (month(Date) == 10 & day(Date) >= 16 & day(Date) <= 20) ~ 3,
    (month(Date) == 9 & day(Date) >= 5 & day(Date) < 10) |
      (month(Date) == 10 & day(Date) >= 21 & day(Date) <= 25) ~ 2,
    TRUE ~ 1))

# Define and run the Bayesian zero-inflated Poisson model WITHOUT YEAR VARIABLE ('Year_f').
# Model explicitly estimates only the deterrent effects and building-level variability ('Location').
final_model <- brm(
  
  # Model formula includes ONLY deterrent effects explicitly separated by year/treatment,
  # plus an explicit offset for migration intensity and explicit random intercept for each building ('Location').
  formula = Birds ~ Predator2023 + Predator2024 + Synthetic2023 +
    offset(MigrationOffset) + (1 | Location),
  
  # Specify the prepared dataset containing only observations from 2023 and 2024.
  data = bird_data,
  
  # Specify a zero-inflated Poisson family because the bird collision data are counts with many zeroes.
  family = zero_inflated_poisson(),
  
  # Sensible middle-ground priors, justified explicitly by earlier analyses, to balance conservative and empirical expectations:
  prior = c(
    set_prior("normal(0, 10)", class = "Intercept"),                # Explicit weak prior on intercept allows flexibility.
    set_prior("normal(-1, 2)", class = "b", coef = "Predator2023"), # Explicit moderately informed prior (Predator deterrent 2023).
    set_prior("normal(-1, 2)", class = "b", coef = "Predator2024"), # Explicit moderately informed prior (Predator deterrent 2024).
    set_prior("normal(-2, 3)", class = "b", coef = "Synthetic2023"),# Explicit moderate prior (Synthetic deterrent 2023).
    set_prior("student_t(3, 0, 2.5)", class = "sd")),               # Explicit weakly informative prior for building variability.
  
  # Configuration of robust Bayesian sampling:
  chains = 4,               # Run four independent Markov chains for robust inference.
  iter = 6000,              # Specify 6000 total iterations per chain for precision.
  warmup = 2000,            # Warm-up phase of 2000 iterations for sampler optimization.
  cores = 4,                # Parallelize sampling across four CPU cores.
  seed = 123,               # Set random seed ensuring reproducibility.
  
  # Stringent sampler control parameters to ensure model convergence:
  control = list(adapt_delta = 0.999, max_treedepth = 15))

# Summarize and output detailed posterior estimates, credible intervals, and diagnostic metrics.
summary(final_model)

# Define the file path
save_path <- "C:\\Users\\jdwolfe\\Dropbox\\Analyses\\Deterrent analysis"

# Save models to .rds files in the specified directory
saveRDS(final_model, file = file.path(save_path, "final_model.rds"))

# Visualize posterior parameter credible intervals to facilitate clear interpretation.
library(bayesplot)
mcmc_intervals(final_model, prob = 0.95)

# Perform posterior predictive checks to ensure model predictions match observed data distributions.
pp_check(final_model)

# ===============================================================================
# DETAILED BAYESIAN MODEL SUMMARY AND DIAGNOSTICS
# (Zero-inflated Poisson model explicitly assessing bird deterrent effectiveness)
# ===============================================================================

# FAMILY AND FORMULA:
# Family: zero_inflated_poisson (appropriate explicitly for count data with many zeros)
# Formula explicitly used:
# Birds ~ Predator2023 + Predator2024 + Synthetic2023 + offset(MigrationOffset) + (1 | Location)

# DATA DETAILS:
# Number of explicit observations used: 1980
# Data explicitly restricted to years 2023 and 2024 only.

# -------------------------------------------------------------------------------
# MULTILEVEL HYPERPARAMETERS (BUILDING-LEVEL RANDOM EFFECT: LOCATION):
# -------------------------------------------------------------------------------
# sd(Intercept): Estimate = 1.62; Est.Error = 0.56; 95% CI = [0.83, 3.01]
# EXPLICIT INTERPRETATION:
# Clear and substantial variability among buildings. Some locations explicitly 
# experienced inherently higher or lower bird collision rates due to unmodeled 
# site-specific characteristics (e.g., building structure, vegetation, location).

# -------------------------------------------------------------------------------
# REGRESSION COEFFICIENTS:
# -------------------------------------------------------------------------------

# INTERCEPT (BASELINE RATE without deterrents):
# Estimate = -3.68; Est.Error = 0.55; 95% CI = [-4.86, -2.65]
# EXPLICIT INTERPRETATION:
# Represents baseline bird collision log-count (no deterrent applied). Serves explicitly 
# as reference to evaluate explicit deterrent effects.

# PREDATOR DETERRENT EFFECT (2023):
# Estimate = -0.81; Est.Error = 0.35; 95% CI = [-1.47, -0.12]
# EXPLICIT INTERPRETATION:
# Credible interval explicitly excludes zero. Predator deterrent explicitly reduced 
# bird collisions significantly (~55% reduction, exp(-0.81) ≈ 0.45) during Fall 2023.

# PREDATOR DETERRENT EFFECT (2024):
# Estimate = -2.26; Est.Error = 0.39; 95% CI = [-3.02, -1.51]
# EXPLICIT INTERPRETATION:
# Credible interval explicitly excludes zero. Strong evidence predator deterrent 
# substantially reduced collisions (~90% reduction, exp(-2.26) ≈ 0.10) during Fall 2024.
# The predator deterrent explicitly showed increased effectiveness compared to 2023.

# SYNTHETIC DETERRENT EFFECT (2023):
# Estimate = -1.14; Est.Error = 0.46; 95% CI = [-2.02, -0.23]
# EXPLICIT INTERPRETATION:
# Credible interval explicitly excludes zero. Clear evidence synthetic deterrent 
# significantly reduced collisions (~68% reduction, exp(-1.14) ≈ 0.32) during Fall 2023.

# -------------------------------------------------------------------------------
# ZERO-INFLATION PARAMETER ('zi'):
# Estimate = 0.84; Est.Error = 0.01; 95% CI = [0.82, 0.86]
# EXPLICIT INTERPRETATION:
# Model explicitly indicates approximately 84% probability of excess zeros in collision counts,
# explicitly validating choice of zero-inflated Poisson model.

# -------------------------------------------------------------------------------
# MODEL DIAGNOSTICS AND SAMPLING QUALITY METRICS:
# -------------------------------------------------------------------------------

# RHAT (POTENTIAL SCALE REDUCTION FACTOR):
# Explicitly Rhat = 1.00 for all parameters, explicitly indicating excellent chain convergence.

# EFFECTIVE SAMPLE SIZE (ESS):
# Bulk and Tail ESS explicitly exceed 3000 for most parameters, explicitly ensuring 
# posterior estimates are robust, reliable, and stable.

# POSTERIOR PREDICTIVE CHECKS ('pp_check'):
# Explicit posterior predictive checks indicate the model adequately reproduces the 
# observed collision data distribution. Visual inspection explicitly confirms good fit.

# -------------------------------------------------------------------------------
# EXPLICIT SUMMARY OF KEY FINDINGS:
# -------------------------------------------------------------------------------
# 1. Predator deterrent explicitly effective in both 2023 (~55% reduction) and 2024 (~90% reduction).
# 2. Synthetic deterrent explicitly effective in 2023 (~68% reduction).
# 3. Significant and explicit building-level variability indicates the need for targeted mitigation strategies.

# =============================================================================== 
# END OF DETAILED MODEL SUMMARY AND DIAGNOSTICS
# ===============================================================================



############################################################################################
############################################################################################
############################################################################################
##########################  SENSITIVITY ANALYSIS ON PRIORS #################################
############################################################################################
############################################################################################
############################################################################################

# Load necessary library for Bayesian modeling if not already loaded.
library(brms)  # brms is the package that runs Bayesian models via Stan.

# -------------------------------------------------------------------------------
# GOAL OF THIS SECTION:
# Test how sensitive your results are to the choice of prior distributions.
# Specifically, fit new models with (1) weaker priors and (2) stronger priors,
# then compare parameter estimates to the original model.
# If posterior estimates are similar, your results are robust to prior choice.
# -------------------------------------------------------------------------------

# Define ALTERNATIVE PRIOR SETS for Sensitivity Analysis:

# 1. Define a set of LESS INFORMATIVE PRIORS (i.e., wider spread, less constraining)
weaker_priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),           # Intercept prior stays weak
  set_prior("normal(0, 5)", class = "b", coef = "Predator2023"), # Prior centered at 0 (neutral), wider SD (5)
  set_prior("normal(0, 5)", class = "b", coef = "Predator2024"),
  set_prior("normal(0, 5)", class = "b", coef = "Synthetic2023"),
  set_prior("student_t(3, 0, 2.5)", class = "sd")             # Random effect prior remains weak
)

# 2. Define a set of STRONGER INFORMATIVE PRIORS (i.e., narrower spread, stronger beliefs)
stronger_priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),            # Intercept prior stays weak
  set_prior("normal(-2, 1)", class = "b", coef = "Predator2023"), # Assume stronger negative effect (mean = -2, SD = 1)
  set_prior("normal(-2, 1)", class = "b", coef = "Predator2024"),
  set_prior("normal(-2, 1)", class = "b", coef = "Synthetic2023"),
  set_prior("student_t(3, 0, 2.5)", class = "sd")             # Random effect prior remains weak
)

# -------------------------------------------------------------------------------
# FIT THE MODELS USING THE DIFFERENT PRIORS:
# -------------------------------------------------------------------------------

# Fit a model using the weaker priors
model_weaker <- brm(
  formula = Birds ~ Predator2023 + Predator2024 + Synthetic2023 + offset(MigrationOffset) + (1 | Location),
  data = bird_data,
  family = zero_inflated_poisson(),
  prior = weaker_priors,
  chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 123,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

# Fit a model using the stronger priors
model_stronger <- brm(
  formula = Birds ~ Predator2023 + Predator2024 + Synthetic2023 + offset(MigrationOffset) + (1 | Location),
  data = bird_data,
  family = zero_inflated_poisson(),
  prior = stronger_priors,
  chains = 4, iter = 6000, warmup = 2000, cores = 4, seed = 123,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

# Define the file path again (if needed)
save_path <- "C:\\Users\\jdwolfe\\Dropbox\\Analyses\\Deterrent analysis"

# Load the saved model from the .rds file
final_model <- readRDS(file = file.path(save_path, "final_model.rds"))

# -------------------------------------------------------------------------------
# SUMMARIZE RESULTS:
# -------------------------------------------------------------------------------

# Summarize and inspect posterior means and credible intervals
summary(final_model)      # Original model
summary(model_weaker)     # Model with weaker priors
summary(model_stronger)   # Model with stronger priors

# -------------------------------------------------------------------------------
# OPTIONAL: CREATE VISUAL COMPARISON OF POSTERIOR DISTRIBUTIONS:
# -------------------------------------------------------------------------------

# Extract posterior samples for each model
posterior_original <- as_draws_df(final_model) %>%
  rename(`Synthetic signal` = b_Synthetic2023,
         `Predator signal (2023)` = b_Predator2023,
         `Predator signal (2024)` = b_Predator2024)

posterior_weaker <- as_draws_df(model_weaker) %>%
  rename(`Synthetic signal` = b_Synthetic2023,
         `Predator signal (2023)` = b_Predator2023,
         `Predator signal (2024)` = b_Predator2024)

posterior_stronger <- as_draws_df(model_stronger) %>%
  rename(`Synthetic signal` = b_Synthetic2023,
         `Predator signal (2023)` = b_Predator2023,
         `Predator signal (2024)` = b_Predator2024)

# Create a combined dataset for plotting
combined_posterior <- bind_rows(
  posterior_original %>% mutate(Prior = "Original"),
  posterior_weaker %>% mutate(Prior = "Weaker"),
  posterior_stronger %>% mutate(Prior = "Stronger")
) %>%
  pivot_longer(cols = c(`Synthetic signal`, `Predator signal (2023)`, `Predator signal (2024)`),
               names_to = "Parameter", values_to = "Estimate")

# Plot density comparisons
library(ggplot2)

ggplot(combined_posterior, aes(x = Estimate, fill = Prior)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior Distributions under Different Priors",
       x = "Effect Size (log scale)",
       y = "Density") +
  scale_fill_manual(values = c("black", "blue", "red"))

# -------------------------------------------------------------------------------
# RESULTS SUMMARY: Prior Sensitivity Analysis
# -------------------------------------------------------------------------------

# This block summarizes the results of a formal prior sensitivity analysis
# conducted to evaluate the robustness of posterior estimates to different
# specifications of prior distributions.

# Three versions of the Bayesian zero-inflated Poisson model were compared:
# (1) the original model using moderately informed priors,
# (2) a model with weaker priors centered at zero with wider standard deviations,
# and (3) a model with stronger priors centered on more negative effects
# with narrower standard deviations.

# Posterior estimates for all three deterrent coefficients remained highly
# consistent across prior conditions, with overlapping credible intervals
# and minimal shifts in posterior means.

# Specifically, the effect of the predator deterrent in 2023 was estimated as:
# -0.81 [-1.47, -0.12] in the original model,
# -0.82 [-1.48, -0.11] under weaker priors, and
# -0.96 [-1.59, -0.32] under stronger priors.

# The effect of the predator deterrent in 2024 was:
# -2.26 [-3.02, -1.51] in the original model,
# -2.28 [-3.07, -1.50] with weaker priors, and
# -2.38 [-3.09, -1.68] with stronger priors.

# The synthetic deterrent effect (2023) was estimated at:
# -1.14 [-2.02, -0.23] in the original model,
# -1.12 [-2.01, -0.18] under weaker priors, and
# -1.31 [-2.12, -0.48] under stronger priors.

# All three models yielded parameter estimates that were credibly different from zero,
# regardless of prior specification. Visual comparison of the posterior densities
# confirmed that distributions overlapped nearly completely, with only minor
# sharpening under stronger priors.

# All chains converged successfully (Rhat = 1.00 for all parameters), with effective
# sample sizes (Bulk_ESS and Tail_ESS) exceeding 3,000 in all cases, indicating that
# posterior inferences were statistically stable.

# Conclusion: the deterrent effects inferred by the model are robust to prior assumptions,
# and the magnitude and direction of effects are supported by the data rather than driven
# by prior specification.


