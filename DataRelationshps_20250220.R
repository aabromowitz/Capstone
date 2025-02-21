# Dr Sadler thought that the relationship I found between Temp and VPMax was almost too exact.
# So try to look at non-linear relationships that aren't quite that good.

# Demo
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
weather$sim1 = 5.509857 + 0.981397*weather$Tdew..degC.
weather$sim2 = (5.509857-3.675153) + 0.981397*weather$Tdew..degC. + 0.764030*weather$VPdef..mbar.
weather$sim3 = 3.238e+00 + 2.920e-02*weather$rh.... * weather$VPdef..mbar.
weather$sim4 = (3.238e+00-3.218634) + 2.920e-02*weather$rh.... * weather$VPdef..mbar. + 0.008003*weather$Tdew..degC.*weather$rh....
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin + NonLin Simulations")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)
legend("bottom", legend = c("Original Temp", "Tdew Regression", "Tdew+VPdef Regression", "rhxVPdef Regression","rhxVPdef+Tdewxrh Regression"), 
       col = c("red", "blue", "green", "cyan", "magenta"), lwd = 2)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin + NonLin Simulations")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)
legend("bottom", legend = c("Original Temp", "Tdew Regression", "Tdew+VPdef Regression", "rhxVPdef Regression","rhxVPdef+Tdewxrh Regression"), 
       col = c("red", "blue", "green", "cyan", "magenta"), lwd = 2)

sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # Sim1 (1 lin) = 4.648065
sqrt(sum((weather$T..degC. - weather$sim2)^2)/length(weather$T..degC.)) # Sim2 (2 lin) = 1.921759
sqrt(sum((weather$T..degC. - weather$sim3)^2)/length(weather$T..degC.)) # Sim3 (1 nonlin) = 4.001809
sqrt(sum((weather$T..degC. - weather$sim4)^2)/length(weather$T..degC.)) # Sim4 (2 nonlin) = 1.491257

file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
n = 6*24*366
freq=c(f1,f2)

b0 = mean(weather$VPmax..mbar.)
A1 = -3.552730
phi1 = 5.356089 
A2 = -7.932237
phi2 = -0.379505
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
resid = weather$VPmax..mbar. - mean(weather$VPmax..mbar.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
coef=c(A1,A2) # -3.552730 -7.932237
est=est.ar.wge(resid,p=16,method='burg')
phi = est$phi
vara = est$avar
ts1 = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)

b0 = mean(weather$Tdew..degC.)
A1 = 0.124165
A2 = 6.557257
phi1 = -1.181011 + pi # 1.960582
phi2 = 15.013296 - 4*pi # 2.446925
resid = weather$Tdew..degC. - mean(weather$Tdew..degC.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
plotLim = 6*24*366 
plotStart = 1
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$Tdew..degC.)
est=est.ar.wge(resid,p=12,method='burg')
coef=c(A1,A2) # -3.552730 -7.932237
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
phi = est$phi
vara = est$avar
ts2 = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)

cor(weather$VPmax..mbar.,weather$Tdew..degC. ) # 0.7126456
cor(ts1,ts2) # 0.5716898: a little bit lower, but honestly more similar than I thought it would be

alpha <- abs(cor(ts1, ts2)-cor(weather$VPmax..mbar.,weather$Tdew..degC. ))*0.75
ts2_adj <- (1 - alpha) * ts2 + alpha * ts1
ts1_adj <- (1 - alpha) * ts1 + alpha * ts2_adj
cor(ts1, ts2)
cor(ts1_adj, ts2_adj)

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-20,60))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-10,15))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 6*24*366 
plotStart = 1
plot(ts1[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Simulated Data Comparison", ylim=c(-20,60))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax Sim", "VPMax Sim Corr", "Tdew Sim", "Tdew Sim Corr"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(ts1[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Simulated Data Comparison", ylim=c(-10,15))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax Sim", "VPMax Sim Corr", "Tdew Sim", "Tdew Sim Corr"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

################################################################################

# Pull data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
library(tswge)

# Get the first non-linear variable
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

# Find relationship
weather$rh...._x_VPdef..mbar. = weather$rh.... * weather$VPdef..mbar.
mod <- lm(T..degC. ~ rh...._x_VPdef..mbar., data = weather)
summary(mod) # coef = 3.238e+00, weight = 2.920e-02

# Look at plot comparison
weather$sim1 = 3.238e+00 + 2.920e-02*weather$rh.... * weather$VPdef..mbar.
plotLim = 3000 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
# Hm, the simulated data doesn't seem to go negative

# Look at residulas
weather$resid = weather$T..degC. - weather$sim1
plot(weather$resid)
cor_values <- sapply(weather_sub, function(col) cor(col, weather$resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted[1:25])

# Add simulated column
weather$Tdew..degC._x_rh.... = weather$Tdew..degC. * weather$rh....
mod <- lm(resid ~ Tdew..degC._x_rh...., data = weather)
summary(mod) # coef = -3.218634, weight = 0.008003

# Plot second simulated column
weather$sim2 = (3.238e+00-3.218634) + 2.920e-02*weather$rh.... * weather$VPdef..mbar. + 0.008003*weather$Tdew..degC.*weather$rh....
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "NonLin Simulation")
points(weather$sim1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$sim2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)

# Try a linear one, but without certain columns
weather <- read.csv(file_path, header = TRUE)
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.",  "VPmax..mbar.", "Tpot..K.", "Tlog..degC.")]
cor_values <- sapply(weather_sub, function(col) cor(col, weather$T..degC., use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted)

# Correlations in general
weather_sub <- weather[, !names(weather) %in% c("date", "T..degC.")]
cor_values <- sapply(weather_sub, function(col) cor(col, weather$T..degC., use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted) # it looks like Tpot..K., Tlog..decC., and VPmax..mbar. are the only ones above 0.9

# Linear relationship
mod <- lm(T..degC. ~ Tdew..degC., data = weather)
summary(mod) # intercept = 5.509857, weight = 0.981397

weather$sim3 = 5.509857 + 0.981397*weather$Tdew..degC.
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)

plotLim = 3000 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)

# Residuals
weather$resid = weather$T..degC. - weather$sim3
cor_values <- sapply(weather_sub, function(col) cor(col, weather$resid, use = "complete.obs"))
cor_values_sorted <- sort(cor_values, decreasing = TRUE)
print(cor_values_sorted[1:25])

# Linear relationship
mod <- lm(resid ~ VPdef..mbar., data = weather)
summary(mod) # intercept = -3.675153, weight = 0.764030

weather$sim4 = (5.509857-3.675153) + 0.981397*weather$Tdew..degC. + 0.764030*weather$VPdef..mbar.
plotLim = 6*24*366 
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)

plotLim = 3000
plotStart = 1
plot(weather$T..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Lin Simulation")
points(weather$sim3[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
points(weather$sim4[(plotStart):(plotStart+plotLim-1)], col = "magenta", pch = 17, cex = 1.5)

# RMSEs
sqrt(sum((weather$T..degC. - weather$sim1)^2)/length(weather$T..degC.)) # Sim1 = 4.001809
sqrt(sum((weather$T..degC. - weather$sim2)^2)/length(weather$T..degC.)) # Sim2 = 1.491257
sqrt(sum((weather$T..degC. - weather$sim3)^2)/length(weather$T..degC.)) # Sim3 = 4.648065
sqrt(sum((weather$T..degC. - weather$sim4)^2)/length(weather$T..degC.)) # Sim4 = 1.921759

# What does rh look like?
plot(weather$rh....)
plot(weather$rh....[1:3000])

# Making sure it has the daily and annual frequency
y = weather$rh.... - mean(weather$rh....)
t = 1:(24*6*366)
fft_result <- fft(y)
power_spectrum <- Mod(fft_result)^2
freqs <- (0:(length(t)-1)) / length(t)  # Normalize frequencies
dominant_indices <- order(power_spectrum, decreasing = TRUE)[1:5]  # Get top 5 frequencies
freqs[dominant_indices[1]] # 0.006944444 / 144 hours / 1 day
freqs[dominant_indices[2]] # 0.9928848
freqs[dominant_indices[3]] # 1.897389e-05 / 52704.01 hours / 366 days
# I think I can do the same approach as before

# Get the amplitude and shift
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
y = weather$rh.... - mean(weather$rh....) # 14.48794
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
summary(fit) # A1 = 13.242671 A2 = 12.539368 phi1 = -0.943365 phi2 = 0.624285

# Getting residuals so far
A1 = 13.242671 
A2 = 12.539368 
phi1 = -0.942353 + 2*pi # 5.340832
phi2 = 0.170724
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$rh....) 
resid = weather$rh.... - y_pred
plot(resid)
plotLim = 6*24*366 
plotStart = 1
plot(weather$rh....[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_pred[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)

# Figure out p
aic5.ar.wge(resid,p=1:30,type='bic',method='burg') # p = 12  0.2826507
est=est.ar.wge(resid,p=12,method='burg')

# Try a generation
n = 6*24*366
b0 = mean(weather$rh....)
coef=c(A1,A2) 
freq=c(f1,f2)
psi=c(phi1,phi2) 
phi = est$phi
vara = est$avar
y_gen = pmin(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),100)
plotLim = 6*24*366 
plotStart = 1
plot(weather$rh....[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$rh....[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)])

################################################################################

# Is there a way to simulate two datasets
# Try using code from chatGPT

file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
n = 6*24*366
freq=c(f1,f2)

b0 = mean(weather$VPmax..mbar.)
A1 = -3.552730
phi1 = 5.356089 
A2 = -7.932237
phi2 = -0.379505
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
resid = weather$VPmax..mbar. - mean(weather$VPmax..mbar.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
coef=c(A1,A2) # -3.552730 -7.932237
est=est.ar.wge(resid,p=16,method='burg')
phi = est$phi
vara = est$avar
ts1 = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)

b0 = mean(weather$Tdew..degC.)
A1 = 0.124165
A2 = 6.557257
phi1 = -1.181011 + pi # 1.960582
phi2 = 15.013296 - 4*pi # 2.446925
resid = weather$Tdew..degC. - mean(weather$Tdew..degC.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
plotLim = 6*24*366 
plotStart = 1
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$Tdew..degC.)
est=est.ar.wge(resid,p=12,method='burg')
coef=c(A1,A2) # -3.552730 -7.932237
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
phi = est$phi
vara = est$avar
ts2 = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)

# How similar are the correlation structures?
cor(weather$VPmax..mbar.,weather$Tdew..degC. ) # 0.7126456
cor(ts1,ts2) # 0.5716898: a little bit lower, but honestly more similar than I thought it would be 

# How do they all compare?
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Non Correlated", ylim=c(-15,55))
points(ts1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Non Correlated", ylim=c(-10,20))
points(ts1[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

# Now try correlating them
# target_corr <- cor(weather$VPmax..mbar.,weather$Tdew..degC. ) 
target_corr <- cor(weather$VPmax..mbar.,weather$Tdew..degC. ) - 0.0
Sigma <- matrix(c(1, target_corr, target_corr, 1), nrow=2)  # Covariance matrix
chol_decomp <- chol(Sigma)  # Cholesky factorization
ts_matrix <- scale(cbind(ts1, ts2))
correlated_ts <- ts_matrix %*% chol_decomp
ts1_corr <- correlated_ts[,1]
ts2_corr <- correlated_ts[,2]
cor(ts1_corr,ts2_corr) # hm, that seems really high now

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-20,60))
points(ts1_corr[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_corr[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-10,15))
points(ts1_corr[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_corr[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

# They gave me another approach
empirical_corr <- cor(ts1, ts2)
lambda <- target_corr / empirical_corr  # Scaling factor
ts2_adj <- lambda * ts2 + (1 - lambda) * rnorm(length(ts2), mean=0, sd=sd(ts2))
final_corr <- cor(ts1, ts2_adj)
print(final_corr)

# Another approach to slightly nudge the data
alpha <- abs(cor(ts1, ts2)-cor(weather$VPmax..mbar.,weather$Tdew..degC. ))*0.75
ts2_adj <- (1 - alpha) * ts2 + alpha * ts1
ts1_adj <- (1 - alpha) * ts1 + alpha * ts2_adj
cor(ts1, ts2)
cor(ts1_adj, ts2_adj)

plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-20,60))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Correlated", ylim=c(-10,15))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax", "VPMax Sim", "Tdew", "Tdew Sim"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 6*24*366 
plotStart = 1
plot(ts1[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Simulated Data Comparison", ylim=c(-20,60))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax Sim", "VPMax Sim Corr", "Tdew Sim", "Tdew Sim Corr"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)

plotLim = 2000
plotStart = 1
plot(ts1[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, 
     ylab = "Temp", main = "Simulated Data Comparison", ylim=c(-10,15))
points(ts1_adj[(plotStart):(plotStart+plotLim-1)], col = "blue", pch = 17, cex = 1.5)
points(ts2[(plotStart):(plotStart+plotLim-1)], col = "green", pch = 17, cex = 1.5)
points(ts2_adj[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 17, cex = 1.5)
legend("bottom", legend = c("VPMax Sim", "VPMax Sim Corr", "Tdew Sim", "Tdew Sim Corr"), 
       col = c("red", "blue", "green", "cyan"), lwd = 2)