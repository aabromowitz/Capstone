# There's that issue where you can get negative temp values in the summer
library(tswge)
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
f=1/366
freq1 = 1/144
dayStart = 1
dayStop = 366
numInDay = 6*24
overlap = 20
numDays = 366
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_VPMax.csv"
vecs_VPMAx = read.csv(file_path, header = TRUE)

# Split points up into 12 chunks, roughly 1 per month.
numSplits = 16
numPerSplit = 366*144/numSplits

# Output the min VPMax for each month
minsPerSplit = numeric(numSplits)
for (split in 1:numSplits){
  t = ((split-1)*numPerSplit+1):(split*numPerSplit)
  min = min(weather$VPmax..mbar.[t])
  minsPerSplit[split] = min
  print(min)
}

# Does that match with the plot?
# plot(weather$VPmax..mbar.)

# Try applying the min
meanVec = vecs_VPMAx$mean
varVec = vecs_VPMAx$var
psi1vec = vecs_VPMAx$psi1
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
A1vec = vecs_VPMAx$A1
AVec = abs(vecs_VPMAx$A1)
lAvec = log(AVec)
pvec = vecs_VPMAx$p
phicoeff1vec = vecs_VPMAx$phi1
phicoeff2vec = vecs_VPMAx$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
# phiShift = 0
phiShift = 0.1
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
if (phiShift > 0){
  b0phi = b0phi + (1-b0phi)*phiShift
} else {
  b0phi = b0phi - b0phi*phiShift
}
# b0var = 6
b0var = 4
f=1/366
meanShift = 3
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
if (phiShift > 0){
  ampPhi = ampPhi + (1-ampPhi)*phiShift
} else {
  ampPhi = ampPhi - ampPhi*phiShift
}
# ampVar = 0.4
ampVar = 0.3
# varAmp = abs(-1.63190)
varAmp = 1.5
varPsi = -0.12543 + pi
# varB = -3.83464
varB = -4
varPhi = 0.4
if (phiShift > 0){
  varPhi = varPhi + (1-varPhi)*phiShift
} else {
  varPhi = varPhi - varPhi*phiShift
}
# varVar = 0.7
varVar = 0.6
# minMultVPMax = 0.6
minMultVPMax = 0.9
t=1:366
numDays = 366
b0sin_true = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift
Asin_true = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec))
y_gen_vpmax_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  psi1 = psi1mean
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax_true[t] = st
}
lastOfDay = seq(numInDay, length(y_gen_vpmax_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax_true)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax_true)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax_true[beforeDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax_true[afterDay] = ((y_gen_vpmax_true[lastOfDay]+y_gen_vpmax_true[firstOfDay])*factor/2 + y_gen_vpmax_true[afterDay]*ii)/(overlap+1)
}
numSims = 20
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  if (seed %% 10 == 0){
    print(seed)
  }
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
    psi1 = psi1mean
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
  # y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*minVPMax)] = min(weather$VPmax..mbar.)*minVPMax
  
  # Apply min per month
  for (split in 1:numSplits){
    t = ((split-1)*numPerSplit+1):(split*numPerSplit)
    min = minsPerSplit[split]
    threshold = min*minMultVPMax
    indices_to_change = t[y_gen_vpmax_noise[t] < threshold]
    y_gen_vpmax_noise[indices_to_change] = threshold
  }
  y_gen_vpmax_mult[, seed] = y_gen_vpmax_noise
}
y_gen_vpmax_adj = y_gen_vpmax_mult
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - VPMax', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
numDays = 366
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - 366 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - First 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][1:(numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - Middle 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_vpmax_true[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)

# Those look alright, so what happens when we add in TDew and Temp
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/vectors_TDew.csv"
vecs_TDew = read.csv(file_path, header = TRUE)
meanVec = vecs_TDew$mean
varVec = vecs_TDew$var
psi1vec = vecs_TDew$psi1
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
A1vec = vecs_TDew$A1
AVec = abs(vecs_TDew$A1)
lAvec = log(AVec)
pvec = vecs_TDew$p
phicoeff1vec = vecs_TDew$phi1
phicoeff2vec = vecs_TDew$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
# phi = c(phi1mean,phi2mean)
# phiShift = 0
phiShift = 0.4
phi = 0.7
# phi = 0.9
if (phiShift > 0){
  phi = phi + (1-phi)*phiShift
} else {
  phi = phi - phi*phiShift
}
f=1/366
b0A = -6.55756
b0Psi = 5.57984
b0B = 5.40932
b0_phi = 0.7
# b0_phi = 0.9
if (phiShift > 0){
  b0_phi = b0_phi + (1-b0_phi)*phiShift
} else {
  b0_phi = b0_phi - b0_phi*phiShift
}
b0_avar = 5.48
AA = -0.29247
APsi = 1.15203
AB = 1.54916
Amp_phi = 0.1603404
# Amp_phi = 0.7
if (phiShift > 0){
  Amp_phi = Amp_phi + (1-Amp_phi)*phiShift
} else {
  Amp_phi = Amp_phi - Amp_phi*phiShift
}
Amp_avar = 0.9459506
varA = -0.89346    
varPsi = 0.18136
varB = -3.75244
Var_phi = 0.2812496
# Var_phi = 0.7
if (phiShift > 0){
  Var_phi = Var_phi + (1-Var_phi)*phiShift
} else {
  Var_phi = Var_phi - Var_phi*phiShift
}
Var_avar = 0.4903821
t=1:366
b0sin_true = b0A * cos(2 * pi * f * t + b0Psi) + b0B
Asin_true = AA * cos(2 * pi * f * t + APsi) + AB
y_gen_tdew_true = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop) {
  t_day = ((day-1)*6*24+1):(day*6*24)
  A1 = Asin_true[day]
  b0 = b0sin_true[day]
  st = b0 + A1 * cos(2 * pi * freq1 * (1:length(t_day)) + psi1)
  y_gen_tdew_true[t_day] = st
}
lastOfDay = seq(numInDay, length(y_gen_tdew_true)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_tdew_true)-numInDay, by = numInDay)
for (ii in 1:overlap) {
  beforeDay = seq(numInDay+1-ii, length(y_gen_tdew_true)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_tdew_true)-numInDay+overlap, by = numInDay)
  factor = overlap + 1 - ii
  y_gen_tdew_true[beforeDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[beforeDay] * ii) / (overlap + 1)
  y_gen_tdew_true[afterDay] = ((y_gen_tdew_true[lastOfDay] + y_gen_tdew_true[firstOfDay]) * factor/2 + y_gen_tdew_true[afterDay] * ii) / (overlap + 1)
}
numSims = 20
y_gen_tdew_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_tdew_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  if (seed %% 10 == 0){
    print(seed)
  }
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
y_gen_vpmax_adj = y_gen_vpmax_mult  
y_gen_tdew_adj = y_gen_tdew_mult 
numAhead = 1
# gapAhead = 0
gapAhead = 40 / numAhead
aheadFlag = 0
iterFlag = 1 # no nudges
lag0BeforeFlag = 0
lag0AfterFlag = 1
numIter = 1
for (sim in 1:numSims){
  if (sim %% 10 == 0){
    print(sim)
  }
  
  # Initialize
  y_gen_vpmax_adj[sim] = y_gen_vpmax_mult[sim]
  y_gen_tdew_adj[sim] = y_gen_tdew_mult[sim]
  
  # Try doing the nudging multiple times
  if (iterFlag == 1){
    for (iter in 1:numIter){
      
      # No lag
      if (lag0BeforeFlag == 1){
        corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
        diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
        alpha <- diff * (1-abs(diff))^2
        y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
        y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
      }
      
      # lag
      if (aheadFlag == 1){
        for (aheadCounter in numAhead:1){
          ahead = aheadCounter * gapAhead
          
          # Positive lag
          corr = cor(weather$VPmax..mbar.[(ahead+1):(6*24*366)],weather$Tdew..degC.[1:(6*24*366-ahead)])
          diff = corr - cor(y_gen_vpmax_adj[(ahead+1):(6*24*366),sim],y_gen_tdew_adj[1:(6*24*366-ahead),sim])
          alpha <- diff * (1-abs(diff))^2
          y_gen_tdew_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_tdew_mult[1:(6*24*366-ahead),sim] + alpha * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim]
          y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_vpmax_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_tdew_adj[1:(6*24*366-ahead),sim]
          
          # Negative lag
          corr = cor(weather$VPmax..mbar.[1:(6*24*366-ahead)],weather$Tdew..degC.[(ahead+1):(6*24*366)])
          diff = corr - cor(y_gen_vpmax_adj[1:(6*24*366-ahead),sim],y_gen_tdew_adj[(ahead+1):(6*24*366),sim])
          alpha <- diff * (1-abs(diff))^2
          y_gen_tdew_adj[(ahead+1):(6*24*366),sim] <- (1 - alpha) * y_gen_tdew_adj[(ahead+1):(6*24*366),sim] + alpha * y_gen_vpmax_adj[1:(6*24*366-ahead),sim]
          y_gen_vpmax_adj[1:(6*24*366-ahead),sim] <- (1 - alpha) * y_gen_vpmax_adj[1:(6*24*366-ahead),sim] + alpha * y_gen_tdew_adj[(ahead+1):(6*24*366),sim]
        }
      }
      
      # No lag
      if (lag0AfterFlag == 1){
        corr = cor(weather$VPmax..mbar.,weather$Tdew..degC.)
        diff = corr - cor(y_gen_vpmax_adj[sim],y_gen_tdew_adj[sim])
        alpha <- diff * (1-abs(diff))^2
        y_gen_tdew_adj[,sim] <- (1 - alpha) * y_gen_tdew_adj[sim] + alpha * y_gen_vpmax_adj[sim]
        y_gen_vpmax_adj[,sim] <- (1 - alpha) * y_gen_vpmax_adj[sim] + alpha * y_gen_tdew_adj[,sim]
      }
      
      # Make sure it isn't too low
      for (split in 1:numSplits){
        t = ((split-1)*numPerSplit+1):(split*numPerSplit)
        min = minsPerSplit[split]
        threshold = min*minMultVPMax
        indices_to_change = t[y_gen_vpmax_noise[t] < threshold]
        y_gen_vpmax_noise[indices_to_change] = threshold
      }
      
    }
  }
  
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
ACF2 = acf(y_gen_vpmax_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - VPMax', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_vpmax_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
ACF = acf(weather$Tdew..degC., plot = "FALSE")
ACF2 = acf(y_gen_tdew_true, plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6,
     main='ACF Comparison - TDew', xlab = 'Lag',ylab = 'ACF')
lines(ACF2$lag, ACF2$acf, lwd = 6, col = "blue")
for (sim in 1:numSims){
  ACF3 = acf(y_gen_tdew_adj[,sim], plot = "FALSE")
  lines(ACF3$lag, ACF3$acf, lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
CCF = ccf(weather$VPmax..mbar.,weather$Tdew..degC.,plot=FALSE)
plot(CCF$lag, CCF$acf , type = "l", lwd = 6, ylim=c(0.71,0.72), xlim=c(-1,1),
     main='CCF Comparison - VPMax and TDew', xlab = 'Lag',ylab = 'ACF')
for (sim in 1:numSims){
  CCF2 = ccf(y_gen_vpmax_adj[,sim], y_gen_tdew_adj[,sim],plot=FALSE)
  lines(CCF2$lag, CCF2$acf, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data"), 
       col = c("black", "red"),pch=1) 
numDays = 366
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - 366 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax_adj[,1][1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - First 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_vpmax_true[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][1:(numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of VPMax - Middle 10 Days', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_vpmax_true[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 366
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_adj[,1][1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_adj[,1][1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - 366 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_true[1:(numDays*6*24)],col = "blue")
points(y_gen_tdew_adj[,1][1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_adj[,1][1:(numDays*6*24)])
ymax = max(weather$Tdew..degC.[1:(numDays*6*24)],y_gen_tdew_adj[,1][1:(numDays*6*24)])
plot(weather$Tdew..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - First 10 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_tdew_true[1:(numDays*6*24)],col = "blue")
points(y_gen_tdew_adj[,1][1:(numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_tdew_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_tdew_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$Tdew..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Realization of TDew - Middle 10 Days', xlab = 'Days',ylab = 'Dew Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_tdew_true[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_tdew_adj[,1][(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('top',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
trunc = 1000
lenPlot = 35
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - First 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_vpmax_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='VPMax Spectral Density Comparison - Middle 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_vpmax_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
numDays = 20
SpecDen = parzen.wge(weather$Tdew..degC.[1:(6*24*numDays)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_tdew_true[1:(6*24*numDays)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='TDew Spectral Density Comparison - First 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_adj[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
startDay = 180
SpecDen = parzen.wge(weather$Tdew..degC.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
SpecDen2 = parzen.wge(y_gen_tdew_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
ymax = max(SpecDen$pzgram,SpecDen2$pzgram)
ymin = min(SpecDen$pzgram,SpecDen2$pzgram)
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='TDew Spectral Density Comparison - Middle 20 Days', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen3$freq,SpecDen3$pzgram, lwd = 2, col = "red")
}
legend('topright',legend = c("Real Data", "Simulated Data - No Error", "Simulated Data - With Error"), 
       col = c("black", "blue", "red"),pch=1)
trunc = 1000
lenPlot = 50
SpecDen = parzen.wge(weather$Tdew..degC.[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_tdew_true[1:(6*24*numDays)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='TDew Spectral Density Comparison - First 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_adj[1:(6*24*numDays),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
SpecDen = parzen.wge(weather$Tdew..degC.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
SpecDen2 = parzen.wge(y_gen_tdew_true[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE", trunc=trunc)
ymax = max(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
ymin = min(SpecDen$pzgram,SpecDen2$pzgram[1:lenPlot])
plot(SpecDen$freq[1:lenPlot],SpecDen$pzgram[1:lenPlot], type = "l", lwd = 6, ylim=c(ymin,ymax),
     main='TDew Spectral Density Comparison - Middle 20 Days - First Peak', xlab = 'Frequency',ylab = 'Pseudo-z Periodogram')
lines(SpecDen2$freq[1:lenPlot],SpecDen2$pzgram[1:lenPlot], lwd = 6, col = "blue")
for (sim in 1:numSims){
  SpecDen3 = parzen.wge(y_gen_tdew_adj[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE", trunc=trunc)
  lines(SpecDen3$freq[1:lenPlot],SpecDen3$pzgram[1:lenPlot], lwd = 2, col = "red")
}
legend('bottomleft',legend = c("Real Data", "Simulated Data - No Noise", "Simulated Data - With Noise"), 
       col = c("black", "blue", "red"),pch=1)
# Those seem alright after 1 nudge

# How about the Temp plots
y_gen_TLin_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_TLin_adj) <- paste0("sim", 1:numSims)
y_gen_TNonlin_adj = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_TNonlin_adj) <- paste0("sim", 1:numSims)
varTemp = 0.1
minTemp = 0.7
numSims = 20
for (sim in 1:numSims){
  # Oh shoot, I might need to add some error variance as well
  set.seed(sim)
  varLin = rnorm(n=length(y_gen_TLin_adj[sim]),mean=0,sd=sqrt(varTemp))
  varNonlin = rnorm(n=length(y_gen_TNonlin_adj[sim]),mean=0,sd=sqrt(varTemp))
  
  y_gen_TLin_adj[sim] = -2.274795  + y_gen_vpmax_mult[sim]*0.815397 + y_gen_tdew_mult[sim]*0.236779 + varLin
  y_gen_TNonlin_adj[sim] = 44.06+243.5*log(y_gen_vpmax_mult[sim]/100)/(17.62-log(y_gen_vpmax_mult[sim]/100))*1.326+y_gen_vpmax_mult[sim]*2.102e-02-y_gen_tdew_mult[sim]*1.000e-03+y_gen_vpmax_mult[sim]*y_gen_tdew_mult[sim]*6.250e-05 + varNonlin
  
}
numDays = 366
ymin = min(weather$T..degC.[1:(numDays*6*24)],y_gen_TLin_adj[1:(numDays*6*24),1],y_gen_TNonlin_adj[1:(numDays*6*24),1])
ymax = max(weather$T..degC.[1:(numDays*6*24)],y_gen_TLin_adj[1:(numDays*6*24),1],y_gen_TNonlin_adj[1:(numDays*6*24),1])
plot(weather$T..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Temp Simulations - 366 Days', xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_TLin_adj[1:(numDays*6*24),1],col = "blue")
points(y_gen_TNonlin_adj[1:(numDays*6*24),1],col = "green")
legend('topleft',legend = c("Real Data", "VPMax + TDew - Linear", "VPMax + TDew - Nonlinear"), 
       col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(weather$T..degC.[1:(numDays*6*24)],y_gen_TLin_adj[1:(numDays*6*24),1],y_gen_TNonlin_adj[1:(numDays*6*24),1])
ymax = max(weather$T..degC.[1:(numDays*6*24)],y_gen_TLin_adj[1:(numDays*6*24),1],y_gen_TNonlin_adj[1:(numDays*6*24),1])
plot(weather$T..degC.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Temp Simulations - First 10 Days', xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0))
points(y_gen_TLin_adj[1:(numDays*6*24),1],col = "blue")
points(y_gen_TNonlin_adj[1:(numDays*6*24),1],col = "green")
legend('top',legend = c("Real Data", "VPMax + TDew - Linear", "VPMax + TDew - Nonlinear"), 
       col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$T..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_TLin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],y_gen_TNonlin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1]))
ymax = max(c(weather$T..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_TLin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],y_gen_TNonlin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1]))
plot(weather$T..degC.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax),
     main='Temp Simulations - Middle 10 Days', xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
ticks <- axTicks(1)
axis(1, at = ticks, labels = round(ticks / 144,0) + startDay)
points(y_gen_TLin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "blue")
points(y_gen_TNonlin_adj[(1+6*24*startDay):(6*24*startDay+numDays*6*24),1],col = "green")
legend('top',legend = c("Real Data", "VPMax + TDew - Linear", "VPMax + TDew - Nonlinear"), 
       col = c("red", "blue", "green"),pch=1)