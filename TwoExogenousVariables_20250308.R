# Pull in data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)

# Get formula for VPMax + Tdew
mod <- lm(T..degC. ~ VPmax..mbar. + Tdew..degC., data = weather)
summary(mod) # int = -2.274677 VPmax..mbar. = 0.815350 Tdew..degC. =  0.236818

# Get formula for VPMax + VPDef
mod <- lm(T..degC. ~ VPmax..mbar. + VPdef..mbar., data = weather)
summary(mod) # int = -4.070479 VPmax..mbar. = 1.119350 VPdef..mbar. -0.275935

weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim0 = (-2.901738-0.630408) + weather$VPmax..mbar.*0.947056 + weather$Tdew..degC.*0.116543
weather$sim2 = -2.274677  + weather$VPmax..mbar.*0.815350 + weather$Tdew..degC.*0.236818
weather$sim3 = -4.070479  + weather$VPmax..mbar.*1.119350 + weather$VPdef..mbar.*-0.275935

sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # 1.877053
sqrt(sum((weather$T..degC. - weather$sim0 )^2)/length(weather$T..degC.)) # 1.743963
sqrt(sum((weather$T..degC. - weather$sim2 )^2)/length(weather$T..degC.)) # 1.595012
sqrt(sum((weather$T..degC. - weather$sim3 )^2)/length(weather$T..degC.)) # 1.707784

# Now try non-linear
weather$sim4 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
sqrt(sum((weather$T..degC. - weather$sim4 )^2)/length(weather$T..degC.)) # 0.04509353

# Only limited to Tdew and VPDef
weather <- read.csv(file_path, header = TRUE)
weather <- weather[,names(weather) %in% c('date','T..degC.','VPmax..mbar.','Tdew..degC.','VPdef..mbar.')]
weather$sim4 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
weather$resid = weather$T..degC. - weather$sim4

# Looking at corrs 
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.", "resid", "sim4")]
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
cor_values <- sapply(weather_sub, function(col) cor(col, weather$resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted) # VPdef..mbar._2 0.586855865 Tdew..degC._x_VPmax..mbar. 0.242519918

# New cols
weather$VPdef..mbar._2 = weather$VPdef..mbar. * weather$VPdef..mbar.
weather$Tdew..degC._x_VPmax..mbar. = weather$Tdew..degC. * weather$VPmax..mbar.
weather$VPmax..mbar._NonLin = 243.5*log(weather$VPmax..mbar./100)/(17.62-log(weather$VPmax..mbar./100))

# Formula for non-lin VPdef
mod <- lm(T..degC. ~ VPmax..mbar._NonLin + VPdef..mbar._2, data = weather)
summary(mod) # int = 4.495e+01 VPmax..mbar._NonLin = 1.349e+00 VPdef..mbar._2 =  3.345e-04

# Formula for non-lin Tdew
mod <- lm(T..degC. ~ VPmax..mbar._NonLin + Tdew..degC._x_VPmax..mbar., data = weather)
summary(mod) # int = 4.491e+01 VPmax..mbar._NonLin = 1.348e+00 Tdew..degC._x_VPmax..mbar. =  2.793e-04

# Errors
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = -2.274677  + weather$VPmax..mbar.*0.815350 + weather$Tdew..degC.*0.236818
weather$sim3 = -4.070479  + weather$VPmax..mbar.*1.119350 + weather$VPdef..mbar.*-0.275935
weather$sim4 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
weather$sim5 = 44.95  + 1.349*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100)) + 3.345e-04*weather$VPdef..mbar.*weather$VPdef..mbar.
weather$sim6 = 44.91  + 1.348*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100)) + 2.793e-04*weather$VPmax..mbar.*weather$Tdew..degC.

sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # 1.877053
sqrt(sum((weather$T..degC. - weather$sim2 )^2)/length(weather$T..degC.)) # 1.595012
sqrt(sum((weather$T..degC. - weather$sim3 )^2)/length(weather$T..degC.)) # 1.707784
sqrt(sum((weather$T..degC. - weather$sim4 )^2)/length(weather$T..degC.)) # 0.04509353
sqrt(sum((weather$T..degC. - weather$sim5 )^2)/length(weather$T..degC.)) # 0.03074466
sqrt(sum((weather$T..degC. - weather$sim6 )^2)/length(weather$T..degC.)) # 0.04068494

(1.877053-1.595012) / (1.877053-1.707784) # 1.666229
(0.04509353-0.03074466) / (0.04509353-0.04068494) # 3.254753, let's choose VPdef

# Plots
weather$sim1 = -2.901738 + weather$VPmax..mbar.*0.947056
weather$sim2 = -4.070479  + weather$VPmax..mbar.*1.119350 + weather$VPdef..mbar.*-0.275935
weather$sim3 = 45.09 + 1.354*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100))
weather$sim4 = 44.95  + 1.349*(243.5*log(weather$VPmax..mbar./100))/(17.62-log(weather$VPmax..mbar./100)) + 3.345e-04*weather$VPdef..mbar.*weather$VPdef..mbar.

plotLim = 2000
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "VPDef Lin")
points(weather$sim1[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+VPDef"), 
       col = c("red", "blue", "green"), lwd = 2)

plotLim = 20
plot(weather$T..degC.[1:plotLim], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "VPDef NonLin")
points(weather$sim3[1:plotLim], col = "blue", pch = 17, cex = 1.5)
points(weather$sim4[1:plotLim], col = "green", pch = 18, cex = 1.5)
legend("topleft", legend = c("Temp", "VPMax", "VPMax+VPDef"), 
       col = c("red", "blue", "green"), lwd = 2)