# Script to generate new VPMax data and test how similar it is to the original VPMax

# Pull in library
library(tswge)

# Pull in files
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_VPMax.csv"
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
phi = c(phi1mean,phi2mean)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
varAmp = abs(-1.63190)
varPsi = -0.12543 + pi
varB = -3.83464
varPhi = 0.4
varVar = 0.7
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
b0sin_true = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift
b0sin = b0sin_true + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE) # use gen.arma since errors had relationship to prior values
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9 # make sure values don't get too low (i.e. negative)

# Model daily amplitude as the log of a sine curve
Asin_true = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec)) # use gen.arma since errors had relationship to prior values
Asin_error = exp(gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
ASin = Asin_true * Asin_error

# Model daily variance as the log of a sine curve
varSin_true = exp(varAmp * cos(2 * pi * f * t + varPsi) + varB) # use gen.arma since errors had relationship to prior values
varSin_error = exp(gen.arma.wge(n=366,phi=varPhi,theta = 0, vara = varVar, sn = 0, plot=FALSE))
varSin = varSin_true * varSin_error

# Generate daily sine curves for VPMax, without any noise in the daily means and amplitudes (for comparison)
y_gen_vpmax_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax_true[t] = st
}

# Daily VPMax values can be separated, so blend them together
lastOfDay = seq(numInDay, length(y_gen_vpmax_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax_true)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax_true)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax_true[beforeDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax_true[afterDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[afterDay]*ii)/(overlap+1)
}

# Regenerate daily VPMax values, but with more noisy daily means and amplitudes
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = ASin[day]
  b0 = b0sin[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}

# Blend the daily sine curves together
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}

# Generate daily variance
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
}

# Make sure that the daily variances are connected so the AR aspect of it is maintained
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}

# Add the signal and the noise together
y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5

# Plot the full realization
numDays = 366
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - 366 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

# Plot the first 10 days zoomed in, to show low variance and amplitude example
numDays = 10
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_noise[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - First 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

# Plot the middle 10 days zoomed in, to show high variance and amplitude example
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - Middle 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_vpmax_true[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)

################################################################################
# Now, generate multiple realizations to confirm that ACF and Spectral Density plots compare well to actual VPMax

# Generate multiple realizations using same approach as above
numSims = 20
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  t = 1:366
  b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
  b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
  ASin = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec) + gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
  varSin = exp(varAmp * cos(2 * pi * f * t + varPsi) + varB + gen.arma.wge(n=366,phi=varPhi,theta = 0, vara = varVar, sn = 0, plot=FALSE))
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = ASin[day] 
    b0 = b0sin[day] 
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays){
    zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varSin[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1] 
  for (day in 2:numDays){
    t = ((day-1)*6*24+1):(day*6*24)
    diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
    sim = zt_int[1:(6*24), day]
    zt_d[t] = sim + diff
  }
  y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
  y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}

# Compare the ACF from the realizations to the ACF of the actual data and the ACF from the signal (which doesn't change)
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

# Compare the spectral density of the first days (low amplitude and variance) from the realizations to the spectral density of the actual data and the spectral density from the signal
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

# Compare the spectral density of the middle days (low amplitude and variance) from the realizations to the spectral density of the actual data and the spectral density from the signal
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
