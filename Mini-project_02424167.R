setwd("C:/Users/xiao/Desktop")

#Install the necessary packages
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
library(readr)
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")
library(emmeans)
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
library(lme4)
if (!requireNamespace("lmerTest", quietly = TRUE)) install.packages("lmerTest")
library(lmerTest)
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
library(Matrix)
if (!require(dplyr)) install.packages("dplyr")
library(dplyr)
library(ggplot2)
# Read data
file_path <- "TAB_Ukaleqarteq2017_2022F.csv"
tab_data <- read_csv(file_path)

# Data preprocessing
# Logarithmic conversion of SST, SIC and Hg
tab_data <- tab_data %>%
  mutate(SST_log = log(SST + 1),
         SIC_log = log(SIC + 1),
         Hg_log = log(Hg_blood + 1),
         year = as.factor(year)) # Ensure that the year is a factor type

# Analysis of variance (ANOVA) and Tukey post hoc comparisons
# ANOVA analysis of SSTs
anova_sst <- aov(SST_log ~ year, data = tab_data)
print(summary(anova_sst))

# Perform Tukey ex post facto comparisons
tukey_sst <- TukeyHSD(anova_sst)
print(tukey_sst)

# Linear mixed effects modeling (LMM)
# Use of LMM to assess the relationship between mercury concentrations with changes in SST/SIC
lmm_model <- lmer(Hg_log ~ SST_log * SIC_log + (1 | year), data = tab_data)
print(summary(lmm_model))


#####Warm ocean conditions affect TAB and foraging patterns

# Read data
tab_data <- read_csv("TAB_Ukaleqarteq2017_2022F.csv")

# Grouping by year and calculating SST correlation statistics
sst_summary <- tab_data %>%
  group_by(year) %>%
  summarise(
    mean = mean(SST, na.rm = TRUE),
    SE = sd(SST, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, n() - 1) * SE,
    CI_upper = mean + qt(0.975, n() - 1) * SE,
    range_lower = min(SST, na.rm = TRUE),
    range_upper = max(SST, na.rm = TRUE),
    N = n()
  ) %>%
  ungroup() 

# Calculate the overall mean, standard error and confidence interval
overall_summary <- summarise(tab_data,
                             mean = mean(SST, na.rm = TRUE),
                             SE = sd(SST, na.rm = TRUE) / sqrt(n()),
                             CI_lower = mean - qt(0.975, df = n() - 1) * SE,
                             CI_upper = mean + qt(0.975, df = n() - 1) * SE,
                             range_lower = min(SST, na.rm = TRUE),
                             range_upper = max(SST, na.rm = TRUE),
                             N = n()
)

# Combined year statistics and aggregate statistics
final_summary <- bind_rows(sst_summary, overall_summary)

# Print results
print(final_summary)


# # Grouped by year and calculated SST correlation statistics
sic_summary <- tab_data %>%
  group_by(year) %>%
  summarise(
    mean = mean(SIC, na.rm = TRUE),
    SE = sd(SIC, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, n() - 1) * SE,
    CI_upper = mean + qt(0.975, n() - 1) * SE,
    range_lower = min(SIC, na.rm = TRUE),
    range_upper = max(SIC, na.rm = TRUE),
    N = n()
  ) %>%
  ungroup() 

# Calculate the overall mean, standard error, and confidence intervals
overall_summary1 <- summarise(tab_data,
                             mean = mean(SIC, na.rm = TRUE),
                             SE = sd(SIC, na.rm = TRUE) / sqrt(n()),
                             CI_lower = mean - qt(0.975, df = n() - 1) * SE,
                             CI_upper = mean + qt(0.975, df = n() - 1) * SE,
                             range_lower = min(SIC, na.rm = TRUE),
                             range_upper = max(SIC, na.rm = TRUE),
                             N = n()
)

# Combined year statistics and aggregate statistics
final_summary1 <- bind_rows(sic_summary, overall_summary1)

# Print results
print(final_summary1)


# Read in the data from the CSV file
data <- read_csv("TAB_Ukaleqarteq2017_2022F.csv")

# Calculate the mean, SE, and CI by year for blood Hg
hg_stats_by_year <- data %>%
  group_by(year) %>%
  summarize(
    mean = mean(Hg_blood, na.rm = TRUE),
    SE = sd(Hg_blood, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE,
    N = n(),
    .groups = 'drop'  # Drop grouping for further operations
  )

# Ensure year is a character before combining
hg_stats_by_year$year <- as.character(hg_stats_by_year$year)

# Calculate overall statistics for blood Hg
overall_hg_stats <- data %>%
  summarize(
    mean = mean(Hg_blood, na.rm = TRUE),
    SE = sd(Hg_blood, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE,
    N = n()
  )

# Create a one-row data frame for the overall stats to match the structure
overall_hg_stats_df <- data.frame(
  year = "Overall",
  overall_hg_stats
)

# Bind the yearly stats and the overall stats
hg_stats <- bind_rows(hg_stats_by_year, overall_hg_stats_df)

# Print the results
print(hg_stats)

# Install and load the necessary packages
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(readr)) install.packages("readr")
library(ggplot2)
library(readr)

# Plotting bar graphs with error lines
plot_with_error_bars <- function(data, y_var, title, y_label) {
  ggplot(data, aes(x = as.factor(year), y = !!sym(y_var), group = year)) +
    geom_bar(stat = "summary", fun = "mean", fill = "blue", color = "black", width = 0.7) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    labs(title = title, x = "Year", y = y_label)
}

# SST bar graph
sst_barplot <- plot_with_error_bars(data, "SST", "Yearly Average Sea Surface Temperature", "Temperature (°C)")

# SIC Bar Chart
sic_barplot <- plot_with_error_bars(data, "SIC", "Yearly Average Sea Ice Coverage", "Coverage (%)")

print(sst_barplot)
print(sic_barplot)

# Installation and loading of necessary packages
if (!require(lme4)) install.packages("lme4")
if (!require(ggplot2)) install.packages("ggplot2")
library(lme4)
library(ggplot2)


# Convert year data into factors to facilitate subsequent analysis
data$year <- as.factor(data$year)

# Categorize the data according to the SST conditions, 'Warm' for 2018 and 2021, and 'Cold' for the others
data$temperature_condition <- ifelse(data$year %in% c("2018", "2021"), "Warm", "Cold")

# Generate summary statistics for 'prop.flightT'
flight_summary <- data %>%
  group_by(year, temperature_condition) %>%
  summarise(
    mean = mean(prop.flightT, na.rm = TRUE),
    SE = sd(prop.flightT, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE
  ) %>%
  ungroup()
print(flight_summary)


# Plot containing actual observations, predicted means and 95% confidence intervals
plot_flight <- ggplot() +
  geom_point(data = data, aes(x = year, y = prop.flightT, color = temperature_condition), alpha = 0.5) +
  geom_point(data = flight_summary, aes(x = year, y = mean, color = temperature_condition), shape = 17, size = 3) +
  geom_errorbar(data = flight_summary, aes(x = year, ymin = CI_lower, ymax = CI_upper, color = temperature_condition), width = 0.1) +
  scale_color_manual(values = c("Cold" = "blue", "Warm" = "red")) +
  labs(title = "Proportion of Time Spent in Flight",
       x = "Year",
       y = "Proportion of Time",
       color = "SST Condition") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    )+
  ylim(0,0.6)


# Print Graphics
print(plot_flight)


# Generate summary statistics for 'prop.colony'
colony_summary <- data %>%
  group_by(year, temperature_condition) %>%
  summarise(
    mean = mean(prop.colony, na.rm = TRUE),
    SE = sd(prop.colony, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE
  ) %>%
  ungroup()
print(colony_summary)

# Plot containing actual observations, predicted means and 95% confidence intervals
plot_colony <- ggplot() +
  geom_point(data = data, aes(x = year, y = prop.colony, color = temperature_condition), alpha = 0.5) +
  geom_point(data = colony_summary, aes(x = year, y = mean, color = temperature_condition), shape = 17, size = 3) +
  geom_errorbar(data = colony_summary, aes(x = year, ymin = CI_lower, ymax = CI_upper, color = temperature_condition), width = 0.1) +
  scale_color_manual(values = c("Cold" = "blue", "Warm" = "red")) +
  labs(title = "Proportion of Time Spent in Colony",
       x = "Year",
       y = "Proportion of Time",
       color = "SST Condition") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,0.6)

# Print Graphics
print(plot_colony)


# Generate summary statistics for 'prop.giving'
dive_summary <- data %>%
  group_by(year, temperature_condition) %>%
  summarise(
    mean = mean(prop.totaldive, na.rm = TRUE),
    SE = sd(prop.totaldive, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE
  ) %>%
  ungroup()
print(dive_summary)

plot_dive <- ggplot() +
  geom_point(data = data, aes(x = year, y = prop.totaldive, color = temperature_condition), alpha = 0.5) +
  geom_point(data = dive_summary, aes(x = year, y = mean, color = temperature_condition), shape = 17, size = 3) +
  geom_errorbar(data = dive_summary, aes(x = year, ymin = CI_lower, ymax = CI_upper, color = temperature_condition), width = 0.1) +
  scale_color_manual(values = c("Cold" = "blue", "Warm" = "red")) +
  labs(title = "Proportion of Time Spent in Diving",
       x = "Year",
       y = "Proportion of Time",
       color = "SST Condition") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,0.6)

# Print Graphics
print(plot_dive)


# Generate summary statistics for 'prop.ice'
ice_summary <- data %>%
  group_by(year, temperature_condition) %>%
  summarise(
    mean = mean(prop.ice, na.rm = TRUE),
    SE = sd(prop.ice, na.rm = TRUE) / sqrt(n()),
    CI_lower = mean - qt(0.975, df = n() - 1) * SE,
    CI_upper = mean + qt(0.975, df = n() - 1) * SE
  ) %>%
  ungroup()

print(ice_summary)

plot_ice <- ggplot() +
  geom_point(data = data, aes(x = year, y = prop.ice, color = temperature_condition), alpha = 0.5) +
  geom_point(data = ice_summary, aes(x = year, y = mean, color = temperature_condition), shape = 17, size = 3) +
  geom_errorbar(data = ice_summary, aes(x = year, ymin = CI_lower, ymax = CI_upper, color = temperature_condition), width = 0.1) +
  scale_color_manual(values = c("Cold" = "blue", "Warm" = "red")) +
  labs(title = "Proportion of Time Spent in Ice",
       x = "Year",
       y = "Proportion of Time",
       color = "SST Condition") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,0.6)

# Print Graphics
print(plot_ice)


#####Effects of mercury on the flight of small puffins Modeling
library(lme4)
library(lmerTest) 
library(readr)

# Read data
data <- read_csv("TAB_Ukaleqarteq2017_2022F.csv") 

# Preprocess the data to make sure all the variables are of the right type
data$year <- as.factor(data$year)

# Interaction model
model_with_interaction <- lmer(prop.flightT ~ Hg_blood * SST * SIC + (1|year), data = data)
print(model_with_interaction)
plot(residuals(model_with_interaction))
drop1(model_with_interaction, test="Chisq")
AIC(model_with_interaction)

# Build the model and output a summary
# Impact on flight time (prop.flightT)
model_flightT1 <- lmer(prop.flightT ~ Hg_blood + SST + SIC  + (1|year), data = data)
summary(model_flightT1)
anova(model_with_interaction, model_flightT1)
model_flightT2 <- lmer(prop.flightT ~ SST + SIC  + (1|year), data = data)
summary(model_flightT2)
AIC(model_with_interaction, model_flightT1)
AIC(model_flightT2, model_flightT1)

# Mean SST values for high SST-low SIC years
mean_sst_high_sst_low_sic <- mean(c(3.62, 5.89))

# Mean SST values for low SST-high SIC years
mean_sst_low_sst_high_sic <- mean(c(2.47, 1.42, -0.055))

# Using the SST coefficients from the model
sst_coef <- 0.023096

# Calculation of expected time-of-flight changes under high SST and low SIC conditions
predicted_change_high_sst_low_sic <- sst_coef * (mean_sst_high_sst_low_sic - mean_sst_low_sst_high_sic)

# Converting changes to percentages
predicted_percentage_change_high_sst_low_sic <- predicted_change_high_sst_low_sic * 100

# output result
cat("In years with high SST and low SIC (2018 and 2021), compared to years with low SST and high SIC (2017, 2019, and 2020),",
    "Proportion of flight time expected to increase for little auks", predicted_percentage_change_high_sst_low_sic, "%。\n")


install.packages("car")
library(car)  # to calculate the VIF
# Build a general linear model instead of a mixed effects model
model_for_vif <- lm(prop.flightT ~ SST + SIC, data = data)

# Calculate VIF
vif_values <- vif(model_for_vif)
print(vif_values)

# Checking the VIF value
if(any(vif_values > 5)) {
  cat("Some variables have a VIF greater than 5, indicating potential multicollinearity.\n")
} else if(any(vif_values > 10)) {
  cat("Some variables have a VIF greater than 10, indicating serious multicollinearity.\n")
} else {
  cat("No significant multicollinearity detected among the variables.\n")
}

# Creating residual plots
plot(residuals(model_flightT1) ~ fitted(model_flightT1))
abline(h = 0, col = "red")  
plot(residuals(model_flightT2) ~ fitted(model_flightT2))
abline(h = 0, col = "red") 
# Creating a Q-Q diagram
qqnorm(residuals(model_flightT1))
qqline(residuals(model_flightT1), col = "red")  
qqnorm(residuals(model_flightT2))
qqline(residuals(model_flightT2), col = "red") 

# # Impacts on colony time (prop. colony)
model_colony1 <- lmer(prop.colony ~ Hg_blood + SST + SIC + (1|year), data = data)
summary(model_colony1)
model_colony2 <- lmer(prop.colony ~ SST + SIC + (1|year), data = data)
summary(model_colony2)
AIC(model_colony2, model_colony1)

# Creating residual plots
plot(residuals(model_colony1) ~ fitted(model_colony1))
abline(h = 0, col = "red")  
plot(residuals(model_colony2) ~ fitted(model_colony2))
abline(h = 0, col = "red") 
# Creating a Q-Q diagram
qqnorm(residuals(model_colony1))
qqline(residuals(model_colony1), col = "red")  
qqnorm(residuals(model_colony2))
qqline(residuals(model_colony2), col = "red")


# Impact on dive time (prop.totaldive)
model_totaldive1 <- lmer(prop.totaldive ~ Hg_blood + SST + SIC + (1|year), data = data)
summary(model_totaldive1)
model_totaldive2 <- lmer(prop.totaldive ~ SST + SIC + (1|year), data = data)
summary(model_totaldive2)
AIC(model_totaldive1, model_totaldive2)
# Creating residual plots
plot(residuals(model_totaldive1) ~ fitted(model_totaldive1))
abline(h = 0, col = "red")  
plot(residuals(model_totaldive2) ~ fitted(model_totaldive2))
abline(h = 0, col = "red") 
# Creating a Q-Q diagram
qqnorm(residuals(model_totaldive1))
qqline(residuals(model_totaldive1), col = "red")  
qqnorm(residuals(model_totaldive2))
qqline(residuals(model_totaldive2), col = "red")

# Calculation of averages based on data 
mean_sic_high_sst <- mean(c(0.09, 0.04))  # Mean SIC for years with high SST and low SIC
mean_sic_other_years <- mean(c(1.83, 2.99, 10.4))  # Average SIC for other years

# Calculate the average change in SIC
sic_change <- mean_sic_high_sst - mean_sic_other_years

# Calculate the percentage change in dive time using the model's estimated SIC coefficients
sic_effect_on_dive <- abs(sic_change) * 0.503561

# Percentage change in output forecast
cat("In years with high SST and low SIC (2018 and 2021), the proportion of dive time 
    for little auks is expected to increase compared to other years", sic_effect_on_dive, "%。\n")



# Impact on ice time (prop.ice)
model_ice1 <- lmer(prop.ice ~ Hg_blood + SST + SIC + (1|year), data = data)
summary(model_ice1)
model_ice2 <- lmer(prop.ice ~ SST + SIC + (1|year), data = data)
summary(model_ice2)
model_ice3 <- lmer(prop.ice ~  SIC + (1|year), data = data)
summary(model_ice3)
AIC(model_ice1, model_ice2)
# Creating residual plots
plot(residuals(model_ice1) ~ fitted(model_ice1))
abline(h = 0, col = "red")  
plot(residuals(model_ice2) ~ fitted(model_ice2))
abline(h = 0, col = "red") 
plot(residuals(model_ice3) ~ fitted(model_ice3))
abline(h = 0, col = "red")
# Creating a Q-Q diagram
qqnorm(residuals(model_ice1))
qqline(residuals(model_ice1), col = "red")  
qqnorm(residuals(model_ice2))
qqline(residuals(model_ice2), col = "red")
qqnorm(residuals(model_ice3))
qqline(residuals(model_ice3), col = "red")

# Define the value of the SIC
sic_high_sst <- c(0.09, 0.04)  # SIC values for years with high SST and low SIC
sic_low_sst <- c(1.83, 2.99, 10.4)  # SIC values for low SST high SIC years

# Calculate the average SIC value
mean_sic_high_sst <- mean(sic_high_sst)
mean_sic_low_sst <- mean(sic_low_sst)

# SIC effects estimated using the model
sic_effect_estimate <- 1.62221

# Calculating the predicted effect of a change in SIC on prop.ice
sic_effect_on_ice <- sic_effect_estimate * (mean_sic_high_sst - mean_sic_low_sst)

# output result
cat("In years with high SST and low SIC (2018 and 2021), compared to years with low SST and high SIC (2017, 2019 and 2020),",
    "Little auks expected to spend more time on the ice", sic_effect_on_ice, "%。\n")
