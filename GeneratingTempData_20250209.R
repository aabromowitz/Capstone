# Demoing
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
plotLim = 3000
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = (-2.901738-0.630408) + weather$VPmax..mbar.*0.947056 + weather$Tdew..degC.*0.116543
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "2 Simulations")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew"), 
       col = c("red", "blue", "green"), lwd = 2)
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = (-2.901738-0.630408) + weather$VPmax..mbar.*0.947056 + weather$Tdew..degC.*0.116543
weather$sim3 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "3 Simulations")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
points(weather$sim3[1:plotLim], col = "cyan", pch = 19, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew", "VPMax non-lin"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

################################################################################

# Pulling data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
temp <- weather$T..degC.

# Get correlation values
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.")]
cor_values <- sapply(weather_sub, function(col) cor(col, temp, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted) # VPmax..mbar was highest at 0.97

# Get the linear model for T..degC and VPmax..mbar.
mod <- lm(T..degC. ~ VPmax..mbar., data = weather)
summary(mod) # VPmax..mbar. = 0.947056, Intercept = -2.901738

# Get the best linear model for the residuals
resid = residuals(mod)
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.", "VPmax..mbar.")]
cor_values <- sapply(weather_sub, function(col) cor(col, resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted) # Tdew..degC. was highest at 0.37, not amazing

# Get the linear model for Tdew..degC.
weather$resid = resid
mod <- lm(resid ~ Tdew..degC., data = weather)
summary(mod) # Tdew..degC. = 0.116543, Intercept = -0.630408

# Plot the 3 relationships
plotLim = 3000
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = (-2.901738-0.630408) + weather$VPmax..mbar.*0.947056 + weather$Tdew..degC.*0.116543
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "2 Simulations")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew"), 
       col = c("red", "blue", "green"), lwd = 2)

# Try with the formula for celcius using vpmax
sum((resid(mod))^2) # 160267.2
weather$sim3 = (243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
mod <- lm(T..degC. ~ sim3, data = weather)
summary(mod) # 1.354e+00, intercept = 4.509e+01

# With the non-linear one
plotLim = 3000
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = (-2.901738-0.630408) + weather$VPmax..mbar.*0.947056 + weather$Tdew..degC.*0.116543
weather$sim3 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "3 Simulations")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
points(weather$sim3[1:plotLim], col = "cyan", pch = 19, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+Tdew", "VPMax non-lin"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)
