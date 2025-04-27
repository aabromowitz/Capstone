# Pull in library
library(tswge)

# Pull in files
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_TDew.csv"
vecs = read.csv(file_path, header = TRUE)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)

# Create vectors and variables used later
meanVec = vecs$mean
varVec = vecs$var
psi1vec = vecs$psi1
psi1vec = psi1vec %% (2*pi)
psi1vec2 = numeric(366)
for (day in 1:366){
  if (psi1vec[day] < pi){
    psi1vec2[day] = psi1vec[day]
  } else { # psi1vec[day] > pi
    psi1vec2[day] = psi1vec[day] - pi
  }
}
psi1mean = mean(psi1vec2) 
A1vec = vecs$A1
AVec = abs(vecs$A1)
lAvec = log(AVec)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
# phi = c(phi1mean,phi2mean)
phi = 0.9
f=1/366
b0A = -6.55756
b0Psi = 5.57984
b0B = 5.40932
# b0_phi = 0.7
b0_phi = 0.9
b0_avar = 5.48
AA = -0.29247
APsi = 1.15203
AB = 1.54916
# Amp_phi = 0.1603404
Amp_phi = 0.7
Amp_avar = 0.9459506
varA = -0.89346    
varPsi = 0.18136
varB = -3.75244
# Var_phi = 0.2812496
Var_phi = 0.7
Var_avar = 0.4903821
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
numInDay = 6*24
overlap = 20
t=1:366
numDays = 366

################################################################################
# First, make one curve and then plot to show it looks reasonable

# Set seed
set.seed(1)

# Model daily b0 as a sine curve, parameter values (like b0A were obtained through analysis in other workbooks)
b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
b0sin = b0sin_true + gen.arma.wge(n=366, phi=b0_phi, theta=0, sn=0, vara =b0_avar, plot = FALSE)

# Model daily amplitude as a mean value + noise
Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
Asin = Asin_true + gen.arma.wge(n=366, phi=Amp_phi, theta=0, sn=0, vara =Amp_avar, plot = FALSE)

# Model daily variance as the log of a sine curve
logVar_sin_true = varA * cos(2 * pi * f * t + varPsi) + varB
logVar_sin = logVar_sin_true + gen.arma.wge(n=366, phi=Var_phi, theta=0, sn=0, vara =Var_avar, plot = FALSE)
var_sin = exp(logVar_sin)
var_sin_true = exp(logVar_sin_true)

# Generate daily sine curves for TDew, without any noise in the daily means and amplitudes (for comparison)
y_gen_tdew_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew_true[t_day] = st
}

# Daily TDew values can be separated, so blend them together
lastOfDay = seq(numInDay, length(y_gen_tdew_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew_true)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew_true)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew_true[beforeDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew_true[afterDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[afterDay] * ii) / (overlap + 1)
}

# Regenerate daily TDew values, but with more noisy daily means and amplitudes
y_gen_tdew = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin[day]
  b0 = b0sin[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew[t_day] = st
}

# Blend the daily sine curves together
lastOfDay = seq(numInDay, length(y_gen_tdew)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew[beforeDay] = ((y_gen_tdew[lastOfDay] + y_gen_tdew[firstOfDay]) * factor/2 + y_gen_tdew[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew[afterDay] = ((y_gen_tdew[lastOfDay] + y_gen_tdew[firstOfDay]) * factor/2 + y_gen_tdew[afterDay] * ii) / (overlap + 1)
}

# Generate daily variance
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # plus one because of connecting days
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays) {
  zt_int[, day] = gen.arma.wge(n = (6*24+1), phi = phi, theta = 0, vara = var_sin[day], sn = 0, plot=FALSE)
}

# Make sure that the daily variances are connected so the AR aspect of it is maintained
zt_d[1:(6*24)] = zt_int[1:(6*24), 1]
for (day in 2:numDays) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24+1), (day-1)] - zt_int[(6*24), (day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t_day] = sim + diff
}

# Add the signal and the noise together
y_gen_tdew_noise = y_gen_tdew[1:(6*24*366)] + zt_d
# y_gen_tdew_noise[y_gen_tdew_noise < (min(weather$Tdew..degC., na.rm=TRUE)*0.5)] = min(weather$Tdew..degC., na.rm=TRUE)*0.5

# Plot the full realization
numDays = 366
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - 366 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_true[1:(numDays*6*24)],col = "blue")
points(y_gen_tdew_noise[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

# Plot the first 10 days zoomed in, to show low variance and amplitude example
numDays = 10
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_noise[1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - First 10 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_true[1:(numDays*6*24)],col = "blue")
points(y_gen_tdew_noise[1:(numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

# Plot the middle 10 days zoomed in, to show high variance and amplitude example
startDay = 366/2
numDays = 10
ymin = min(c(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_tdew_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_tdew_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - Middle 10 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_tdew_true[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_tdew_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

################################################################################
# Now, generate multiple realizations to confirm that ACF and Spectral Density plots compare well to actual VPMax

# Generate multiple realizations using same approach as above
numSims = 20
y_gen_tdew_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  t = 1:366
  b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
  b0sin = b0sin_true + gen.arma.wge(n=366, phi=0.7, theta=0, sn=0, vara =5.48, plot = FALSE)
  Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
  Asin = Asin_true + gen.arma.wge(n=366, phi=Amp_phi, theta=0, sn=0, vara =Amp_avar, plot = FALSE)
  logVar_sin_true = varA * cos(2 * pi * f * t + varPsi) + varB
  logVar_sin = logVar_sin_true + gen.arma.wge(n=366, phi=Var_phi, theta=0, sn=0, vara =Var_avar, plot = FALSE)
  var_sin = exp(logVar_sin)
  var_sin_true = exp(logVar_sin_true)
  y_gen_tdew = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = Asin[day] 
    b0 = b0sin[day] 
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_tdew[t] = st
  }
  lastOfDay = seq(numInDay, length(y_gen_tdew)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_tdew)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_tdew)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_tdew)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_tdew[beforeDay] = ((y_gen_tdew[lastOfDay]+y_gen_tdew[firstOfDay])*factor/2 + y_gen_tdew[beforeDay]*ii)/(overlap+1)
    y_gen_tdew[afterDay] = ((y_gen_tdew[lastOfDay]+y_gen_tdew[firstOfDay])*factor/2 + y_gen_tdew[afterDay]*ii)/(overlap+1)
  }
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # plus one because of connecting days
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays) {
    zt_int[, day] = gen.arma.wge(n = (6*24+1), phi = phi, theta = 0, vara = var_sin[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1]
  for (day in 2:numDays) {
    t_day = ((day-1)*6*24+1):(day*6*24)
    diff = zt_int[(6*24+1), (day-1)] - zt_int[(6*24), (day-1)]
    sim = zt_int[1:(6*24), day]
    zt_d[t_day] = sim + diff
  }
  y_gen_tdew_noise = y_gen_tdew[1:(6*24*366)] + zt_d
  # y_gen_tdew_noise[y_gen_tdew_noise < (min(weather$Tdew..degC., na.rm=TRUE)*0.5)] = min(weather$Tdew..degC., na.rm=TRUE)*0.5
  y_gen_tdew_mult[, seed] = y_gen_tdew_noise
}

# Compare the ACF from the realizations to the ACF of the actual data and the ACF from the signal (which doesn't change)
ACF = acf(weather$Tdew..degC., plot = "FALSE")
ACF2 = acf(y_gen_tdew_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6, 
     main='ACF Comparison', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_tdew_mult[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

# Compare the spectral density of the first days (low amplitude and variance) from the realizations to the spectral density of the actual data and the spectral density from the signal
numDays = 20
SpecDen = parzen.wge(weather$Tdew..degC.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_tdew_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

# Compare the spectral density of the middle days (low amplitude and variance) from the realizations to the spectral density of the actual data and the spectral density from the signal
startDay = 180
SpecDen = parzen.wge(weather$Tdew..degC.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_tdew_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

################################################################################

# Ug, it looks like I have to calculate Amp_avar and Amp_phi
Amp = abs(vecs$A1)
y = Amp
Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
Amp_res = y-Asin_true
est = est.ar.wge(Amp_res,p=1,factor=TRUE,method="mle")
est$avar # 0.9459506
est$phi # 0.1603404
