library(tidyverse)
library(minpack.lm)
library(ggpubr)
library(broom)

# Reshape and annotate
long_df <- df %>%
  pivot_longer(-Time, names_to = "Condition", values_to = "Luminescence") %>%
  mutate(
    Insulin = factor(str_extract(Condition, "\\d+nM"),
                     levels = c("5nM", "10nM", "20nM", "30nM", "40nM", "50nM"))
  ) %>%
  # Remove 48-hour time point for 5nM only
  filter(!(Insulin == "5nM" & Time == 48))

# Constrained harmonic regression + bootstrap CI function
fit_harmonic_ci <- function(df) {
  model <- nlsLM(
    Luminescence ~ a + b * sin(2 * pi * Time / period + c),
    data = df,
    start = list(a = mean(df$Luminescence), b = sd(df$Luminescence), period = 24, c = 0),
    lower = c(-Inf, -Inf, 20, -Inf),
    upper = c(Inf, Inf, 26, Inf),
    control = nls.lm.control(maxiter = 500)
  )
  
  pred_time <- seq(min(df$Time), max(df$Time), length.out = 100)
  preds <- predict(model, newdata = data.frame(Time = pred_time))
  
  boot_preds <- replicate(500, {
    boot_sample <- df[sample(nrow(df), replace = TRUE), ]
    boot_model <- tryCatch(
      nlsLM(
        Luminescence ~ a + b * sin(2 * pi * Time / period + c),
        data = boot_sample,
        start = coef(model),
        lower = c(-Inf, -Inf, 20, -Inf),
        upper = c(Inf, Inf, 26, Inf),
        control = nls.lm.control(maxiter = 200)
      ),
      error = function(e) NULL
    )
    
    if (inherits(boot_model, "nls")) {
      predict(boot_model, data.frame(Time = pred_time))
    } else {
      rep(NA, length(pred_time))
    }
  })
  
  ci_lower <- apply(boot_preds, 1, quantile, 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_preds, 1, quantile, 0.975, na.rm = TRUE)
  
  tibble(Time = pred_time, Fit = preds, Lower = ci_lower, Upper = ci_upper)
}

# Fit models per insulin condition
harmonic_results <- long_df %>%
  group_by(Insulin) %>%
  nest() %>%
  mutate(preds_ci = map(data, possibly(fit_harmonic_ci, NULL))) %>%
  filter(!map_lgl(preds_ci, is.null)) %>%
  unnest(preds_ci)

# Plot with 12-hour axis
ggplot() +
  geom_point(data = long_df, aes(x = Time, y = Luminescence), alpha = 0.5, size = 1.5, color = "grey20") +
  geom_ribbon(data = harmonic_results, aes(x = Time, ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.3) +
  geom_line(data = harmonic_results, aes(x = Time, y = Fit), color = "lightblue", linewidth = 1) +
  facet_wrap(~Insulin, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +
  theme_pubclean(base_size = 14) +
  labs(
    x = "Time (Hours)",
    y = "% Change in Luminescence",
    title = "Harmonic Regression (Period constrained 20–26 hrs)",
    subtitle = "Solid lines: Harmonic fit; Shaded area: 95% CI; Points: Observations"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey70", fill = NA)
  )

ggplot() +
  geom_point(data = long_df, aes(x = Time, y = Luminescence), alpha = 0.5, size = 1.5, color = "grey20") +
  geom_ribbon(data = harmonic_results, aes(x = Time, ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.3) +
  geom_line(data = harmonic_results, aes(x = Time, y = Fit), color = "lightblue", linewidth = 1) +
  facet_wrap(~Insulin) +  # All panels now share the same y-axis scale
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +
  theme_pubclean(base_size = 14) +
  labs(
    x = "Time (Hours)",
    y = "% Change in Luminescence",
    title = "Harmonic Regression (Period constrained 20–26 hrs)",
    subtitle = "Solid lines: Harmonic fit; Shaded area: 95% CI; Points: Observations"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey70", fill = NA)
  )

ggplot() +
  geom_point(data = long_df, aes(x = Time, y = Luminescence), alpha = 0.5, size = 1.5, color = "grey20") +
  geom_ribbon(data = harmonic_results, aes(x = Time, ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.3) +
  geom_line(data = harmonic_results, aes(x = Time, y = Fit), color = "lightblue", linewidth = 1) +
  facet_wrap(~Insulin) +  # Removes free scaling to have same y-axis for all panels
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +
  theme_pubclean(base_size = 14) +
  labs(
    x = "Time (Hours)",
    y = "% Change in Luminescence",
    title = "Harmonic Regression (Period constrained 20–26 hrs)",
    subtitle = "Solid lines: Harmonic fit; Shaded area: 95% CI; Points: Observations"
  ) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major = element_line(color = "grey95"),
    panel.border = element_rect(color = "grey70", fill = NA)
  )


