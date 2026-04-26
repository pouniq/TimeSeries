#########################################################################
# Copper Price Analysis with GARCH Model
# Using adjusted column for smooth, dividend-adjusted prices
# Created by: Amirmohammad Abdollahpour 2/dey/1404
#########################################################################

# Load required libraries.
library(rugarch)
library(tsibble)
library(quantmod)
library(fpp2)
library(writexl)
library(FinTS)
library(xts)
library(tseries)
library(ggplot2)

#########################################################################
# 1. Data Import and Preparation.
#########################################################################

# Download copper PPI data from FRED
getSymbols("PCOPPUSDM", src = "FRED")
copper_ppi <- PCOPPUSDM

class(copper_ppi)
# Visualize full dataset
chartSeries(copper_ppi)

# Subset data to 2015-2025
cp <- copper_ppi['2015/2025']
chartSeries(cp, type = "line",
            theme = chartTheme("white", up.col = '#b87133'  ))

# Log transformation
lncp <- log(cp)
chartSeries(lncp)

# First difference of log prices (returns)
copper_returns <- diff(lncp)
chartSeries(copper_returns, type = "line",
            theme = chartTheme("white", up.col = '#b87133'))

# Remove NA values
copper_returns <- na.omit(copper_returns)

# Test for stationarity
adf_result <- adf.test(copper_returns)
print(adf_result)

#########################################################################
# 2. Exploratory Analysis ----
#########################################################################

# Assign to shorter variable name
copper <- copper_returns

# Check ACF and PACF
acf(copper)
pacf(copper)


#########################################################################
# 3. ARMA Model Estimation ----
#########################################################################

# Fit ARMA(4,2) model
copper_arma <- arima(copper, order = c(3, 0, 2))
print(summary(copper_arma))

est <- copper_arma$coef
se  <- sqrt(diag(copper_arma$var.coef))
z_vals <- est / se
p_vals <- 2 * (1 - pnorm(abs(z_vals)))
results <- data.frame(
  Estimate = est,
  StdError = se,
  z_value  = z_vals,
  p_value  = p_vals
)

results$p_value <- formatC(results$p_value, format = "f", digits = 6)
results


# Extract fitted values and residuals
copper_fitted <- xts(fitted(copper_arma), order.by = index(copper))
copper_residuals <- copper - copper_fitted

# Plot actual vs fitted
plot(copper, col = "black", lwd = 2,
     main = "ARMA(4,2): Actual vs Fitted Values",
     ylab = "Log Returns", xlab = "Date")
lines(copper_fitted, col = "#red", lwd = 2)
legend("topleft",
       legend = c("Actual", "Fitted"),
       col = c("black", "#red"),
       lwd = 2.5)

# Ljung-Box test for residual autocorrelation
lb_pvalues <- sapply(1:50, function(l) {
  Box.test(copper_residuals, lag = l, type = "Ljung-Box")$p.value
})
plot(lb_pvalues, type = "h", ylim = c(0, 1),
     main = "Ljung-Box Test P-values",
     xlab = "Lag", ylab = "P-value")
abline(h = 0.05, col = "red", lty = 2)

# ACF of residuals
acf(copper_residuals, main = "ACF of ARMA Residuals")
acf(copper_residuals^2, main = "ACF of ARMA Residuals(squared)")

chartSeries(copper_residuals, type = "line",
            theme = chartTheme("white", up.col = '#b87133'))

# ARCH-LM test
# my arch test accept the h0 meaning there is not need for Garch and ARCH
arch_test <- ArchTest(copper_residuals)
print(arch_test)


#########################################################################
# 4. GARCH Model Estimation ----
#########################################################################

# Specify GARCH(0,1) with ARMA(4,2) mean equation
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(0, 1)),
  mean.model = list(armaOrder = c(4, 2), include.mean = TRUE),
  distribution.model = "std" # Student's t-distribution
)

# Fit the model
garch_fit <- ugarchfit(garch_spec, data = copper)
print(garch_fit)

# Extract fitted values
garch_fitted <- fitted(garch_fit)

# Plot actual vs fitted
plot(copper, col = "black", lwd = 1,
     main = "GARCH: Actual vs Fitted Values",
     ylab = "Log Returns", xlab = "Date")
lines(garch_fitted, col = "#red", lwd = 2)
legend("topleft",
       legend = c("Actual", "Fitted"),
       col = c("black", "#red"),
       lwd = 3)

# Calculate R-squared
r_squared <- cor(as.numeric(copper), as.numeric(garch_fitted))^2
cat("R-squared:", round(r_squared, 4), "\n")

#########################################################################
# 5. Model Diagnostics ----
#########################################################################

# Extract standardized residuals
std_residuals <- residuals(garch_fit, standardize = TRUE)
garch_resid <- residuals(garch_fit)
acf(garch_resid)
acf(garch_resid^2)

sapply(1:50, function(l) {
  Box.test(garch_resid^2, lag = l, type = "Box-Pierce")$p.value
})
# Ljung-Box test on standardized residuals
lb_std_pvalues <- sapply(1:50, function(l) {
  Box.test(std_residuals, lag = l, type = "Ljung-Box")$p.value
})
plot(lb_std_pvalues, type = "h", ylim = c(0, 1),
     main = "Ljung-Box Test P-values (Standardized Residuals)",
     xlab = "Lag", ylab = "P-value")
abline(h = 0.05, col = "red", lty = 2)

# Ljung-Box test on squared standardized residuals
lb_sq_pvalues <- sapply(1:50, function(l) {
  Box.test(std_residuals^2, lag = l, type = "Ljung-Box")$p.value
})
plot(lb_sq_pvalues, type = "h", ylim = c(0, 1),
     main = "Ljung-Box Test P-values (Squared Std. Residuals)",
     xlab = "Lag", ylab = "P-value")
abline(h = 0.05, col = "red", lty = 2)

# ACF of standardized residuals
acf(std_residuals, main = "ACF of Standardized Residuals")

# Plot conditional volatility
plot(sigma(garch_fit), col = "red", lwd = 2,
     main = "Conditional Volatility",
     ylab = "Sigma", xlab = "Date")

ArchTest(garch_resid)

#########################################################################
# 6. in sample fit
#########################################################################



# CORRECT way: Convert each fitted return back to price
# Start with the second price (since first return uses first price)
fitted_prices <- numeric(length(garch_fitted))

for(i in 1:length(garch_fitted)) {
  # Each fitted price = previous ACTUAL price * exp(fitted return)
  fitted_prices[i] <- as.numeric(cp[i]) * exp(as.numeric(garch_fitted[i]))
}

# Create xts object
fitted_original_scale <- xts(fitted_prices, 
                             order.by = index(cp)[(1:length(fitted_prices)) + 1])

# Check the fit
head(cbind(cp, fitted_original_scale), 10)
tail(cbind(cp, fitted_original_scale), 10)

plot(cp, main = "Actual vs Fitted Prices", ylab = "Price", col = "black")
lines(fitted_original_scale, col = "red", lwd = 1.5)
legend("topleft", legend = c("Actual", "Fitted"), 
       col = c("black", "red"), lty = 1, lwd = c(1, 1.5))

# Calculate fit quality
mse <- mean((as.numeric(cp[2:length(cp)]) - fitted_prices)^2)
rmse <- sqrt(mse)
cat("RMSE:", rmse, "\n")
#########################################################################
# 7. Forecast ----
#########################################################################

forecast_horizon <- 10
fc <- ugarchforecast(garch_fit, n.ahead = forecast_horizon)
forecasted_returns <- fitted(fc) 

forecasted_returns <- as.numeric(fc@forecast$seriesFor)
last_actual_price <- as.numeric(tail(cp, 1))

forecasted_prices <- last_actual_price * exp(cumsum(forecasted_returns))

last_date <- index(tail(cp, 1))
forecast_dates <- seq(last_date, by = "month", length.out = forecast_horizon + 1)[-1]
forecasted_prices_xts <- xts(forecasted_prices, order.by = forecast_dates)

plot(cp, col = "black", lwd = 1,
     main = "GARCH: Actual vs Fitted Values",
     ylab = "Log Returns", xlab = "Date")
lines(forecasted_prices_xts, col = "red", lwd = 2)
legend("topleft",
       legend = c("Actual", "Fitted"),
       col = c("black", "red"),
       lwd = 3)






forecast_summary <- data.frame(
  Date = format(forecast_dates, "%Y-%m"),
  Forecasted_Price = round(forecasted_prices, 2),
  Lower_95 = round(lower_prices, 2),
  Upper_95 = round(upper_prices, 2)
)

print(forecast_summary)



