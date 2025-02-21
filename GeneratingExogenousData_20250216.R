# Demo
library(tswge)
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
y_gen = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Saturation Vapor Presure')
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Saturation Vapor Presure')
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)], main = 'Saturation Vapor Presure')

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
y_gen = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Dew Temp')
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Dew Temp')
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)], main = 'Dew Temp')

b0 = mean(log(weather$VPdef..mbar.+1))
A1 = 0.55900
A2 = 0.642318
phi1 = -0.942353 + pi # 2.199243
phi2 = 0.170724 + pi # 3.312313
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(log(weather$VPdef..mbar.+1)) 
resid = log(weather$VPdef..mbar.+1) - y_pred
est=est.ar.wge(resid,p=12,method='burg')
coef=c(A1,A2) 
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar
y_gen = pmax(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Vapor Pressure Deficit')
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5, main = 'Vapor Pressure Deficit')
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], main = 'Vapor Pressure Deficit')

################################################################################

# Pulling in data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)

# Looking at variables of interest: VPmax..mbar., Tdew..degC., VPdef..mbar.
plotts.wge(weather$VPmax..mbar.)
plotts.wge(weather$T..degC.)
plotts.wge(weather$VPmax..mbar.[1:3000])
plotts.wge(weather$T..degC.[1:3000])
# weather and VPmax look pretty similar

plotts.wge(weather$Tdew..degC.)
plotts.wge(weather$T..degC.)
plotts.wge(weather$Tdew..degC.[1:3000])
plotts.wge(weather$T..degC.[1:3000])

plotts.wge(weather$VPdef..mbar.)
plotts.wge(weather$T..degC.)
plotts.wge(weather$VPdef..mbar.[1:3000])
plotts.wge(weather$T..degC.[1:3000])
# All 3/4 of them have the same general pattern of 2 frequencies, not much variability on the small scale, and a clear autocorrelation
# The differences would be in the amplitudes, AR(p) terms, and the variances

################################################################################

# Get the interpolated data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)

# What does it look like when you use a low pass / average filter
VPmax = ma.smooth.wge(weather$VPmax..mbar.,order=1001)
plotts.wge(VPmax$smooth)
plotts.wge(VPmax$smooth[1:3000])

# What does the nls function do?
# x(t)=coef[1]*cos(2*pi*freq[1]*t+psi[1])+coef[2]*cos(2*pi*freq[2]*t+psi[2])+a(t)
f1 = 1/(6*24)
f2 = 1/(6*24*366)
y = weather$VPmax..mbar. - mean(weather$VPmax..mbar.) # 14.48794
t = 1:(6*24*366)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
summary(fit) # A1 = -3.552730 A2 = -7.932237 phi1 = 5.356089 phi2 = -0.379505

# Shift values
# A*cos(2*pi*f*t + phi) = -A * (2*pi*f*t + phi + pi)
A1 = -3.552730
phi1 = 5.356089 
A2 = -7.932237
phi2 = -0.379505
resid = weather$VPmax..mbar. - mean(weather$VPmax..mbar.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
plot(resid)
plot(resid[1:3000])
plotLim = 6*24*366 # 500
plotStart = 1
y=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$VPmax..mbar.)
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
# This looks pretty bad, like the 6*24 frequency is off

# Does FFT give better frequency values?
y = weather$VPmax..mbar. - mean(weather$VPmax..mbar.)
fft_result <- fft(y)
power_spectrum <- Mod(fft_result)^2
freqs <- (0:(length(t)-1)) / length(t)  # Normalize frequencies
dominant_indices <- order(power_spectrum, decreasing = TRUE)[1:5]  # Get top 5 frequencies
freqs[dominant_indices[1]] # 1.897389e-05 / 52704.01 hours / 366 days
freqs[dominant_indices[2]] # 0.999981
freqs[dominant_indices[3]] # 0.006944444 / 144 hours / 1 day
# sigh, it's not the frequencies

# What's the AIC for this?
aic5.ar.wge(resid,p=16:30,type='bic',method='burg') # p=16  -2.802359
est=est.ar.wge(resid,p=16,method='burg') # burg is way quicker

# Let's try some generalizations
n = 6*24*366
b0 = mean(weather$VPmax..mbar.)
coef=c(A1,A2) # -3.552730 -7.932237
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
phi = est$phi
vara = est$avar
y_gen = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 500
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)])
# This doesn't look like an exact match, but it doesn't look terrible either

# It's all over the place on close inspection, how about p = 1?
est=est.ar.wge(resid,p=1,method='burg')
n = 6*24*366
b0 = mean(weather$VPmax..mbar.)
coef=c(A1,A2) # -3.552730 -7.932237
freq=c(1/(6*24),1/(6*24*366))
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
phi = est$phi
vara = est$avar
y_gen = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 5000 
plotStart = 1
plot(weather$VPmax..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)])

A1 = 0.55900
A2 = 0.642318
phi1 = -0.942353 + pi # 2.199243
phi2 = 0.170724 + pi # 3.312313
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(log(weather$VPdef..mbar.+1)) 
resid = log(weather$VPdef..mbar.+1) - y_pred
plot(resid)
plot(resid[1:3000]) # blarg
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_pred)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)

# Figure out p
aic5.ar.wge(resid,p=1:30,type='bic',method='burg') # p = 12  -3.259619
est=est.ar.wge(resid,p=12,method='burg')

# Try a generation
n = 6*24*366
b0 = mean(log(weather$VPdef..mbar.+1))
coef=c(A1,A2) 
freq=c(1/(6*24),1/(6*24*366))
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar
y_gen = pmax(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)])

################################################################################
# Do the same thing, but for dew temp

# Plot variable
plot(weather$Tdew..degC.) # looks similar
plot(weather$Tdew..degC.[1:3000]) # there definitely appears to be a period of some sort

# First, confirm that the periods match
y = weather$Tdew..degC. - mean(weather$Tdew..degC.)
fft_result <- fft(y)
power_spectrum <- Mod(fft_result)^2
freqs <- (0:(length(t)-1)) / length(t)  # Normalize frequencies
dominant_indices <- order(power_spectrum, decreasing = TRUE)[1:5]  # Get top 5 frequencies
freqs[dominant_indices[1]] # 1.897389e-05 / 52704.01 hours / 366 days
freqs[dominant_indices[2]] # 0.999981
freqs[dominant_indices[3]] # 0.0001138434 / 8783.996 hours / 61 days, hmm
freqs[dominant_indices[4]] # 0.9998862
freqs[dominant_indices[5]] # 3.794778e-05, wait is there really no 144 period?

# I guess proceed with assuming that there is the normal behavior
f1 = 1/(6*24)
f2 = 1/(6*24*366)
y = weather$Tdew..degC. - mean(weather$Tdew..degC.) # 5.409322
t = 1:(6*24*366)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
summary(fit) # A1 = -0.124165 A2 = 6.557257 phi1 = -1.181011 phi2 = 15.013296

# Look at the residuals
A1 = 0.124165
A2 = 6.557257
phi1 = -1.181011 + pi # 1.960582
phi2 = 15.013296 - 4*pi # 2.446925
resid = weather$Tdew..degC. - mean(weather$Tdew..degC.) - A1 * cos(2 * pi * f1 * t + phi1) - A2 * cos(2 * pi * f2 * t + phi2)
plot(resid)
plot(resid[1:3000]) # blarg
plotLim = 6*24*366 
plotStart = 1
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$Tdew..degC.)
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_pred[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)

# Figure out p
aic5.ar.wge(resid,p=1:30,type='bic',method='burg') # p = 12  -3.259619
est=est.ar.wge(resid,p=12,method='burg')

# Try a generation
n = 6*24*366
b0 = mean(weather$Tdew..degC.)
coef=c(A1,A2) # -3.552730 -7.932237
freq=c(1/(6*24),1/(6*24*366))
psi=c(phi1,phi2+2*pi) # 5.356089 -0.379505
phi = est$phi
vara = est$avar
y_gen = gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$Tdew..degC.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)])

################################################################################
# Do the same thing, but for vapor pressure deficit
# This is VPmax - VPactual

# Plot variable
plot(weather$VPdef..mbar.) # it seems to go up, but it almost looks like it has more variance as the numbers are larger
plot(weather$VPdef..mbar.[1:3000]) # there definitely appears to be a period of some sort

# First, confirm that the periods match
y = weather$VPdef..mbar. - mean(weather$VPdef..mbar.) # 4.810997
fft_result <- fft(y)
power_spectrum <- Mod(fft_result)^2
freqs <- (0:(length(t)-1)) / length(t)  # Normalize frequencies
dominant_indices <- order(power_spectrum, decreasing = TRUE)[1:5]  # Get top 5 frequencies
freqs[dominant_indices[1]] # 1.897389e-05 / 52704.01 hours / 366 days
freqs[dominant_indices[2]] # 0.999981
freqs[dominant_indices[3]] # 0.006944444 / 144 hours / 1 day

# Fit the sine curves
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
summary(fit) # A1 = 3.530766 A2 = -3.938670 phi1 = -4.063334 phi2 = -0.071566

# Look at the residuals
A1 = 3.530766
A2 = 3.938670
phi1 = -4.063334 + 2*pi # 2.219851
phi2 = -0.071566 + pi # 3.070027
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(weather$VPdef..mbar.) # it looks like it can't be negative
resid = y - y_pred
plot(resid)
plot(resid[1:3000]) # blarg
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_pred[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)

# Figure out p
aic5.ar.wge(resid,p=1:30,type='bic',method='burg') # p = 12  -3.259619
est=est.ar.wge(resid,p=26,method='burg')

# Try a generation
n = 6*24*366
b0 = mean(weather$VPdef..mbar.)
coef=c(A1,A2) 
freq=c(1/(6*24),1/(6*24*366))
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar
y_gen = pmax(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points(y_gen[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge(y_gen[(plotStart):(plotStart+plotLim-1)])
# Sigh, these plots aren't really as good
# I think I would have to figure out a way to make the variance go up as the value goes up, and to remove negatives

################################################################################

# What about logs?
# Plot variable
plot(log(weather$VPdef..mbar.+1)) # it seems to go up, but it almost looks like it has more variance as the numbers are larger
plot(log(weather$VPdef..mbar.+1)[1:3000])

# Fit the sine curves
f1 = 1/(6*24)
f2 = 1/(6*24*366)
t = 1:(6*24*366)
y = log(weather$VPdef..mbar.+1) - mean(log(weather$VPdef..mbar.+1)) 
fit <- nls(y ~ A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2), 
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0),
           data = data.frame(t,y))
summary(fit) # A1 = -0.55900 A2 = -0.642318 phi1 = -0.942353 phi2 = 0.170724

# Look at the residuals
A1 = 0.55900
A2 = 0.642318
phi1 = -0.942353 + pi # 2.199243
phi2 = 0.170724 + pi # 3.312313
y_pred=A1 * cos(2 * pi * f1 * t + phi1) + A2 * cos(2 * pi * f2 * t + phi2) + mean(log(weather$VPdef..mbar.+1)) 
resid = log(weather$VPdef..mbar.+1) - y_pred
plot(resid)
plot(resid[1:3000]) # blarg
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_pred)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)

# Figure out p
aic5.ar.wge(resid,p=1:30,type='bic',method='burg') # p = 12  -3.259619
est=est.ar.wge(resid,p=12,method='burg')

# Try a generation
n = 6*24*366
b0 = mean(log(weather$VPdef..mbar.+1))
coef=c(A1,A2) 
freq=c(1/(6*24),1/(6*24*366))
psi=c(phi1,phi2+2*pi) 
phi = est$phi
vara = est$avar
y_gen = pmax(gen.sigplusnoise.wge(n,b0=b0,b1=0,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=0),0)
plotLim = 6*24*366 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotLim = 2000 
plotStart = 1
plot(weather$VPdef..mbar.[(plotStart):(plotStart+plotLim-1)], type = "p", col = "red", pch = 16, cex = 1.5)
points((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)], col = "cyan", pch = 19, cex = 1.5)
plotts.wge((exp(y_gen)-1)[(plotStart):(plotStart+plotLim-1)])
# That's better