# Demoing
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather$sim1 = 5.509857 + 0.981397*weather$Tdew..degC.
weather$sim2 = (5.509857-1.689e+00) + 0.981397*weather$Tdew..degC. + 6.500e-03*weather$PAR...mol.m..s.
weather$sim3 = -25.46842 + 15.78015*log(weather$VPact..mbar. + 1)
weather$sim4 = (-25.46842-1.685e+00) + 15.78015*log(weather$VPact..mbar. + 1) + 6.544e-06*weather$p..mbar.*weather$PAR...mol.m..s.

plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # Sim1 = 4.648065
sqrt(sum((weather$T..degC. - weather$sim2)^2)/length(weather$T..degC.)) # Sim2 = 3.74863
sqrt(sum((weather$T..degC. - weather$sim3)^2)/length(weather$T..degC.)) # Sim3 = 4.64174
sqrt(sum((weather$T..degC. - weather$sim4)^2)/length(weather$T..degC.)) # Sim4 = 3.743787

################################################################################

# Pull data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)

# Use Tdew as first relationship
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.",  "VPmax..mbar.", "Tpot..K.", "Tlog..degC.")]
weather$sim1 = 5.509857 + 0.981397*weather$Tdew..degC.
weather$resid = weather$T..degC. - weather$sim1
cor_values <- sapply(weather_sub, function(col) cor(col, weather$resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted[1:25])

# Get second relationship
mod <- lm(resid ~ PAR...mol.m..s., data = weather)
summary(mod)

# Look at errors for both
weather$sim2 = (5.509857-1.689e+00) + 0.981397*weather$Tdew..degC. + 6.500e-03*weather$PAR...mol.m..s.
sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # Sim1 = 4.648065
sqrt(sum((weather$T..degC. - weather$sim2)^2)/length(weather$T..degC.)) # Sim2 = 3.74863

# Plots
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

# What does PAR look like?
plot(weather$PAR...mol.m..s.)
plot(log(weather$PAR...mol.m..s. + 1))
plot(log(weather$PAR...mol.m..s. + 1)[1:3000]) # It probably becomes 0 at night

# Try for non-linear
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.",  "VPmax..mbar.", "Tpot..K.", "Tlog..degC.")]
cols <- names(weather_sub)
for (pair in combn(cols, 2, simplify = FALSE)) { # interactions
  new_col_name <- paste(pair, collapse = "_x_")
  weather_sub[[new_col_name]] <- weather_sub[[pair[1]]] * weather_sub[[pair[2]]]
}
for (col in cols) { # squares
  new_col_name <- paste(col, "_2", sep="")
  weather_sub[[new_col_name]] <- weather_sub[[col]] * weather_sub[[col]]
}
for (col in cols) { # logs
  new_col_name <- paste(col, "_log", sep="")
  weather_sub[[new_col_name]] <- log(weather_sub[[col]]+1) 
}
cor_values <- sapply(weather_sub, function(col) cor(col, weather$T..degC., use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted[1:25])

# USe VPact..mbar._log
weather$VPact..mbar._log = log(weather$VPact..mbar. + 1)
mod <- lm(T..degC. ~ VPact..mbar._log, data = weather)
summary(mod) # -25.46842, 15.78015
weather$sim3 = -25.46842 + 15.78015*log(weather$VPact..mbar. + 1)
sqrt(sum((weather$T..degC. - weather$sim3)^2)/length(weather$T..degC.)) # Sim1 = 4.64174

# Next variable
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.",  "VPmax..mbar.", "Tpot..K.", "Tlog..degC.", "rh....", "VPdef..mbar.", "sh..g.kg.")]
cols <- names(weather_sub)
for (pair in combn(cols, 2, simplify = FALSE)) { # interactions
  new_col_name <- paste(pair, collapse = "_x_")
  weather_sub[[new_col_name]] <- weather_sub[[pair[1]]] * weather_sub[[pair[2]]]
}
for (col in cols) { # squares
  new_col_name <- paste(col, "_2", sep="")
  weather_sub[[new_col_name]] <- weather_sub[[col]] * weather_sub[[col]]
}
for (col in cols) { # logs
  new_col_name <- paste(col, "_log", sep="")
  weather_sub[[new_col_name]] <- log(weather_sub[[col]]+1) 
}
weather_sub$sim3 = -25.46842 + 15.78015*log(weather_sub$VPact..mbar. + 1)
weather_sub$resid = weather$T..degC. - weather_sub$sim3
cor_values <- sapply(weather_sub, function(col) cor(col, weather_sub$resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted[1:25]) # p..mbar._x_PAR...mol.m..s. 0.5911682

# Figure out relationship
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather$sim3 = -25.46842 + 15.78015*log(weather$VPact..mbar. + 1)
weather$resid = weather$T..degC. - weather$sim3
weather$p..mbar._x_PAR...mol.m..s. = weather$p..mbar. * weather$PAR...mol.m..s.
mod <- lm(resid ~ p..mbar._x_PAR...mol.m..s., data = weather)
summary(mod) # -1.685e+00, 6.544e-06

# Look at errors for both
weather$sim1 = 5.509857 + 0.981397*weather$Tdew..degC.
weather$sim2 = (5.509857-1.689e+00) + 0.981397*weather$Tdew..degC. + 6.500e-03*weather$PAR...mol.m..s.
weather$sim3 = -25.46842 + 15.78015*log(weather$VPact..mbar. + 1)
weather$sim4 = (-25.46842-1.685e+00) + 15.78015*log(weather$VPact..mbar. + 1) + 6.544e-06*weather$p..mbar.*weather$PAR...mol.m..s.
sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # Sim1 = 4.648065
sqrt(sum((weather$T..degC. - weather$sim2)^2)/length(weather$T..degC.)) # Sim2 = 3.74863
sqrt(sum((weather$T..degC. - weather$sim3)^2)/length(weather$T..degC.)) # Sim3 = 4.64174
sqrt(sum((weather$T..degC. - weather$sim4)^2)/length(weather$T..degC.)) # Sim4 = 3.743787

# Look at plots
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 18, cex = 1.5)
legend("bottom", legend = c("Temp", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green"), lwd = 2)

# All 4 plots
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin + NonLin Simulations")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green", "cyan", "magenta"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin + NonLin Simulations")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)
legend("bottom", legend = c("Temp", "TDew", "TDew + PAR", "LogVPAct", "LogVPAct + p*PAR"), 
       col = c("red", "blue", "green", "cyan", "magenta"), lwd = 2)

# What does the data look like?
plot(weather$VPact..mbar.)
plot(weather$VPact..mbar.[1:2000]) # This is good, know how to simulate
plot(weather$p..mbar.)
plot(weather$p..mbar.[1:2000]) # This might only have a daily 

# Correlations
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather_sub <- weather[, names(weather) %in% c("T..degC.", "Tdew..degC.", "VPact..mbar.", "p..mbar.", "PAR...mol.m..s.")]
cor(weather_sub) # Tdew and VPact are highly correlated, but not in the same prediction
