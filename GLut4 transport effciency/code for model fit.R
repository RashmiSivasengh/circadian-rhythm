library(dplyr)
library(ggplot2)
library(lme4)

# Import your data first (assuming you saved it as "glut4_data.csv")
df <- read.csv("glut4_data.csv")

# Calculate the GLUT4 ratio (Surface/Total)
df <- df %>%
  mutate(GLUT4_Ratio = GLUT4_Surface / Glut4.Total)

# Convert TimePoint to numeric if needed
df$TimePoint <- as.numeric(df$TimePoint)


ggplot(df, aes(x = TimePoint, y = GLUT4_Ratio)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "GLUT4 Translocation over Time",
       x = "Time (min)", y = "Surface/Total GLUT4 Ratio") +
  theme_minimal()

# Mixed-effects model with Time as a fixed effect and Glut4 Total as a covariate
model <- lmer(GLUT4_Ratio ~ TimePoint + Glut4.Total + (1 | Sample), data = df)
summary(model)

# Residual plot
plot(model)

# Predicted vs Actual
df$Predicted <- predict(model)

ggplot(df, aes(x = Predicted, y = GLUT4_Ratio)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title = "Predicted vs Actual GLUT4 Ratios",
       x = "Predicted", y = "Actual") +
  theme_minimal()
