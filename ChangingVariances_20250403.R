# Plot vars
library(tswge)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/vectors4.csv"
vecs = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
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
vara = mean(varVec)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
plot(varVec)

# One day of one var
vara = mean(varVec)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
set.seed(1)
zt1_1 = gen.arma.wge(n = 144*1, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
zt1_2 = gen.arma.wge(n = 144*1, phi = phi, theta = 0, vara = varVec[1], sn = 0, plot=FALSE)
ymin = min(c(zt1_1,zt1_2))
ymax = max(c(zt1_1,zt1_2))
plot(zt1_1,col = "red",ylim=c(ymin,ymax))
points(zt1_2,col = "blue")
legend('topleft',legend = c("Sim - Mean Var", "Sim - Daily Vars"), col = c("red", "blue"),pch=1) # It definitely is less variable

# Try combining 2 days
numDays = 2
zt_2 = numeric(numDays*6*24)
set.seed(1)
for (day in 1:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  zt_2[t] = gen.arma.wge(n = 6*24*1, phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_1 = gen.arma.wge(n = 6*24*numDays, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
plot(zt_1,col = "red",ylim=c(ymin,ymax))
points(zt_2,col = "blue")
legend('topleft',legend = c("Sim - Mean Var", "Sim - Daily Vars"), col = c("red", "blue"),pch=1) # It definitely is less variable
# The issue is going to be that there could be a jump in the values from one day to the next

# Try with 10 days
numDays = 10
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("sim", 1:numDays)
zt_d = numeric(numDays*6*24)
set.seed(1)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
zt_m = gen.arma.wge(n = 6*24*numDays, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
ymin = min(c(zt_m,zt_d))
ymax = max(c(zt_m,zt_d))
plot(zt_m,col = "red",ylim=c(ymin,ymax))
points(zt_d,col = "blue")
legend('topleft',legend = c("Sim - Mean Var", "Sim - Daily Vars"), col = c("red", "blue"),pch=1) # This looks pretty good

# Try for all 366 days
numDays = 366
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("sim", 1:numDays)
zt_d = numeric(numDays*6*24)
set.seed(1)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
zt_m = gen.arma.wge(n = 6*24*numDays, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
ymin = min(c(zt_m,zt_d))
ymax = max(c(zt_m,zt_d))
plot(zt_m,col = "red",ylim=c(ymin,ymax))
points(zt_d,col = "blue")
legend('topleft',legend = c("Sim - Mean Var", "Sim - Daily Vars"), col = c("red", "blue"),pch=1) 
# This looks like it's working, it is definitely more variable in the middle

# What happens when we try using the real b0 values
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = abs(A1vec[day])
  b0 = meanVec[day] # use non-simulated values
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
vara = mean(vecs$var)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
# zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, plot=FALSE)
numDays = 366
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
set.seed(1)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt_d[1:(numDays*6*24)],col = "green")
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt_d[1:(numDays*6*24)],col = "green")
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt_d[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
# It definitely has similar behavior

# Try ACFs and Spectral Densities
meanVec = vecs$mean
dayStart = 1
dayStop = 366
numSims = 20
psi1 = psi1mean
freq1 = 1/144
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = abs(A1vec[day])
    b0 = meanVec[day] # use data
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  # zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays){
    zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
  for (day in 2:numDays){
    t = ((day-1)*6*24+1):(day*6*24)
    diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
    sim = zt_int[1:(6*24), day]
    zt_d[t] = sim + diff
  }
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt_d
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
# This doesn't look perfect, but it doesn't look bad either.
# I'll show this to Dr Sadler to see if he's fine with it, so we can move on.

################################################################################
# Look at b0 with a sin approximation

# plot original b0
plot(meanVec)

# Try to get a sin curve for it
f = 1/366
t = 1:366
y = meanVec
fit <- nls(y ~ A * cos(2 * pi * f * t + psi), 
           start = list(A = 1, psi = 0),
           data = data.frame(t,y))
summary(fit) # A = -7.9308 psi = -0.3872

# Plot sin curve
A = abs(-7.9308)
psi = -0.3872 + pi
b0sin = A * cos(2 * pi * f * t + psi) + mean(meanVec)
ymin = min(c(meanVec,b0sin))
ymax = max(c(meanVec,b0sin))
plot(meanVec,col = "red",ylim=c(ymin,ymax))
points(b0sin,col = "blue")
legend('topleft',legend = c("Daily Means", "Daily Means - Sin"), col = c("red", "blue"),pch=1) 
# Honestly, this looks to be too low

# Try moving up
meanShift = 1
ampShift = 1
A = abs(-7.9308) + ampShift
psi = -0.3872 + pi
b0sin = A * cos(2 * pi * f * t + psi) + mean(meanVec) + meanShift
ymin = min(c(meanVec,b0sin))
ymax = max(c(meanVec,b0sin))
plot(meanVec,col = "red",ylim=c(ymin,ymax))
points(b0sin,col = "blue")
legend('topleft',legend = c("Daily Means", "Daily Means - Sin"), col = c("red", "blue"),pch=1) 
# I can see why it does it, because even though it does a better job near the max, it does a worse job near the mins

# What do the residuals look like?
A = abs(-7.9308)
psi = -0.3872 + pi
b0sin = A * cos(2 * pi * f * t + psi) + mean(meanVec)
res = meanVec - b0sin
plot(res)

# Figure out AR phi for residuals
est = est.ar.wge(res,p=1,factor=TRUE,method="mle") # phi = 0.7142715, avar = 6.071388

# Plot the sin curve plus the noise
b0A = abs(-7.9308)
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
meanShift = 1 # So it has a harder time going negative
set.seed(1)
b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
ymin = min(c(meanVec,b0sin))
ymax = max(c(meanVec,b0sin))
plot(meanVec,col = "red",ylim=c(ymin,ymax))
points(b0sin,col = "blue")
legend('topleft',legend = c("Daily Means", "Daily Means - Sin"), col = c("red", "blue"),pch=1) 

# Now try with a full realization
set.seed(1)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
t=1:366
meanShift = 3 # So it has a harder time going negative
b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  A1 = abs(A1vec[day])
  # b0 = meanVec[day] # use non-simulated values
  b0 = b0sin[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
vara = mean(vecs$var)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
numDays = 366
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
# This doesn't look unreasonable

# What do the ACFs and Spectral densities look like
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.3 # to drop ACF plot
b0var = 6
f=1/366
meanShift = 3 
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
numSims = 20
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:10)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  t=1:366
  b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
  b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = abs(A1vec[day])
    b0 = b0sin[day] # use simulated b0
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
  lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
  firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
  for (ii in 1:overlap){
    beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
    afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
    factor = overlap+1-ii
    y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
    y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
  }
  # zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  numDays = 366
  zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
  colnames(zt_int) <- paste0("day", 1:numDays)
  zt_d = numeric(numDays*6*24)
  for (day in 1:numDays){
    zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
  }
  zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
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
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
# This looks about the same as before, so pretty good.  I had to drop the phi value from 0.7 to 0.3 so the ACFs would match though.

################################################################################

# Now try with generating the A1 values

# Plot the amplitudes
plot(abs(vecs$A1))
AVec = abs(vecs$A1)

# Find the sin curve
f = 1/366
t = 1:366
y = AVec
fit <- nls(y ~ A * cos(2 * pi * f * t + psi), 
           start = list(A = 1, psi = 0),
           data = data.frame(t,y))
summary(fit) # A = -3.1011 phi = -0.2016

# Plot sin curve
A = abs(-3.1011)
psi = -0.2016 + pi
ASin = A * cos(2 * pi * f * t + psi) + mean(AVec)
ymin = min(c(AVec,ASin))
ymax = max(c(AVec,ASin))
plot(AVec,col = "red",ylim=c(ymin,ymax))
points(ASin,col = "blue")
legend('topleft',legend = c("Daily Means", "Daily Means - Sin"), col = c("red", "blue"),pch=1) 
# This doesn't look terrible

# What do the residuals look like?
Aamp = abs(-3.1011)
psiAmp = -0.2016 + pi
b0sin = Aamp * cos(2 * pi * f * t + psiAmp)
res = AVec - ASin
plot(res) # Right, these don't really look like they're being generated by a white noise or an AR process

# What happens when you take the log of AR?
plot(log(AVec)) # That looks a little bit more normal

# Get the A and psi values
f = 1/366
t = 1:366
y = log(AVec)
fit <- nls(y ~ A * cos(2 * pi * f * t + psi), 
           start = list(A = 1, psi = 0),
           data = data.frame(t,y))
summary(fit) # A = -1.00422 psi = -0.07425

# Now how about the plot
lAvec = log(AVec)
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ASin = Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec)
ymin = min(c(lAvec,ASin))
ymax = max(c(lAvec,ASin))
plot(lAvec,col = "red",ylim=c(ymin,ymax))
points(ASin,col = "blue")
legend('topleft',legend = c("Daily Means", "Daily Means - Sin"), col = c("red", "blue"),pch=1) 

# Now how to the residuals look like
ASin = Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec)
res = lAvec - ASin
plot(res) 

# Figure out AR phi for residuals
est = est.ar.wge(res,p=1,factor=TRUE,method="mle") # phi = 0.3537245, avar = 0.4929509

# Plot the sin curve plus the noise
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.3
set.seed(1)
ASin = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec) + gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
ymin = min(c(AVec,ASin))
ymax = max(c(AVec,ASin))
plot(AVec,col = "red",ylim=c(ymin,ymax))
points(ASin,col = "blue")
legend('topleft',legend = c("Daily Amplitudes", "Daily Amplitudes - Sin"), col = c("red", "blue"),pch=1) 
# This actually looks pretty similar

# Now try with a full realization
set.seed(1)
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.7
b0var = 6
f=1/366
t=1:366
meanShift = 3 # So it has a harder time going negative
b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
Aamp = abs(-1.00422)
psiAmp = -0.07425 + pi
ampPhi = 0.4
ampVar = 0.4
ASin = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec) + gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24+1)
  # A1 = abs(A1vec[day])
  A1 = ASin[day]
  b0 = b0sin[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
numInDay = 144
overlap = 20
lastOfDay = seq(numInDay, length(y_gen_vpmax)-numInDay, by = numInDay)
firstOfDay = seq(numInDay+1, length(y_gen_vpmax)-numInDay, by = numInDay)
for (ii in 1:overlap){
  beforeDay = seq(numInDay+1-ii, length(y_gen_vpmax)-numInDay, by = numInDay)
  afterDay = seq(numInDay+ii, length(y_gen_vpmax)-numInDay+overlap, by = numInDay)
  factor = overlap+1-ii
  y_gen_vpmax[beforeDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[beforeDay]*ii)/(overlap+1)
  y_gen_vpmax[afterDay] = ((y_gen_vpmax[lastOfDay]+y_gen_vpmax[firstOfDay])*factor/2 + y_gen_vpmax[afterDay]*ii)/(overlap+1)
}
vara = mean(vecs$var)
pvec = vecs$p
phicoeff1vec = vecs$phi1
phicoeff2vec = vecs$phi2
phi1mean = mean(phicoeff1vec[pvec==2]) 
phi2mean = mean(phicoeff2vec[pvec==2])
phi = c(phi1mean,phi2mean)
numDays = 366
zt_int = data.frame(matrix(0, nrow = (6*24+1), ncol = numDays))  # add an extra value, because we want to see what the difference is
colnames(zt_int) <- paste0("day", 1:numDays)
zt_d = numeric(numDays*6*24)
for (day in 1:numDays){
  zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
}
zt_d[1:(6*24)] = zt_int[1:(6*24), 1] # Set the first day
for (day in 2:numDays){
  t = ((day-1)*6*24+1):(day*6*24)
  diff = zt_int[(6*24*1+1),(day-1)] - zt_int[(6*24*1),(day-1)]
  sim = zt_int[1:(6*24), day]
  zt_d[t] = sim + diff
}
y_gen_vpmax_noise = y_gen_vpmax[1:(6*24*366)] + zt_d
y_gen_vpmax_noise[y_gen_vpmax_noise < (min(weather$VPmax..mbar.)*0.5)] = min(weather$VPmax..mbar.)*0.5
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[1:(numDays*6*24)],col = "green")
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax_noise[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
# It looks pretty similar, but does kind of have lower values

# What do the ACFs and Spectral densities look like
ampShift = 1.5
b0A = abs(-7.9308) + ampShift
b0psi = -0.3872 + pi
b0phi = 0.4 
b0var = 6
f=1/366
t=1:366
meanShift = 3 
freq1 = 1/144
psi1 = psi1mean
dayStart = 1
dayStop = 366
numSims = 20
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  b0sin = b0A * cos(2 * pi * f * t + b0psi) + mean(meanVec) + meanShift + gen.arma.wge(n=366,phi=b0phi,theta = 0, vara = b0var, sn = 0, plot=FALSE)
  b0sin[b0sin < (min(meanVec)*0.9)] = min(meanVec)*0.9
  Aamp = abs(-1.00422)
  psiAmp = -0.07425 + pi
  ampPhi = 0.5
  ampVar = 0.3
  ASin = exp(Aamp * cos(2 * pi * f * t + psiAmp) + mean(lAvec) + gen.arma.wge(n=366,phi=ampPhi,theta = 0, vara = ampVar, sn = 0, plot=FALSE))
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    A1 = ASin[day] # use simulated A1
    b0 = b0sin[day] # use simulated b0
    st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
    y_gen_vpmax[t] = st
  }
  numInDay = 144
  overlap = 20
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
    zt_int[, day] = gen.arma.wge(n = (6*24*1+1), phi = phi, theta = 0, vara = varVec[day], sn = 0, plot=FALSE)
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
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
}
numDays = 20
SpecDen = parzen.wge(weather$VPmax..mbar.[1:(6*24*numDays)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[1:(6*24*numDays),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}
startDay = 180
SpecDen = parzen.wge(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*numDays+6*24*startDay)], plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram , type = "l", lwd = 6)
for (sim in 1:numSims){
  SpecDen2 = parzen.wge(y_gen_vpmax_mult[(1+6*24*startDay):(6*24*numDays+6*24*startDay),sim], plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}