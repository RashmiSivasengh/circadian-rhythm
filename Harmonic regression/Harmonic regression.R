# Load necessary libraries
library(ggplot2)
library(minpack.lm)
library(dplyr)

# Sample Time Data
Time <- c(0, 6, 18, 24)

# PER3 Data
Per3_Group1 <- c(79.3388, 306.8659, 150.8181, 129.7069)
Per3_Group2 <- c(86.3280, 225.8811, 173.1701, 178.6724)
Per3_Group3 <- c(96.1552, 58.6567, 219.1544, 193.6247)

# Combine PER3 Data into a Single Data Frame
data <- data.frame(
  Time = rep(Time, 3),
  Expression = c(Per3_Group1, Per3_Group2, Per3_Group3),
  Condition = 'PER3'
)

# Define the Harmonic Function
harmonic_func <- function(t, A, omega, phi, C) {
  return(A * sin(omega * t + phi) + C)
}

# Attempt to Fit the Harmonic Regression Model for PER3
fit_per3 <- tryCatch(
  nlsLM(
    Expression ~ harmonic_func(Time, A, omega, phi, C),
    data = data,
    start = list(A = 100, omega = 2 * pi / 24, phi = 0, C = mean(data$Expression)),
    lower = c(A = 0, omega = 0, phi = -pi, C = 0),
    upper = c(A = 1000, omega = 2 * pi, phi = pi, C = 1000)
  ),
  error = function(e) {
    message("PER3 fit failed, using mean line instead.")
    return(NULL)
  }
)

# Generate Fitted Curve Data
time_fine <- seq(min(Time), max(Time), length.out = 1000)

if (!is.null(fit_per3)) {
  # If the model fit succeeds
  fit_per3_values <- predict(fit_per3, newdata = list(Time = time_fine))
} else {
  # If the model fit fails, use a straight mean line
  fit_per3_values <- rep(mean(data$Expression), length(time_fine))
}

# Plot PER3 Data and Fitted Curve
ggplot() +
  geom_point(data = data, aes(x = Time, y = Expression), color = 'black', size = 2) +
  geom_line(aes(x = time_fine, y = fit_per3_values), color = 'black', size = 1) +
  labs(
    title = 'PER3 Harmonic Regression Fit',
    x = 'Time (Hours)',
    y = 'Expression (logCPM Z-Scores)'
  ) +
  theme_minimal()
