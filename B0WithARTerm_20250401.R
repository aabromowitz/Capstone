# I'm going to see what happens if instead of me using a white noise error term for B0, 
#   I make it have an AR(1) term

# Plot mean values
library(tswge)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/vectors4.csv"
vecs = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
meanVec = vecs$mean
plot(meanVec)

# Plot triangle with residuals
t = (1:366)
y = 24.15488 - abs(((366/2)-t))*0.1050146
res = y - meanVec
plot(res)

# What kind of AR(1) phi does it give?
est = est.ar.wge(res, p = 1)
# phi is about 0.8, so I'll just try that

# What do the mean values look like
set.seed(1)
maxRes = 8
meanVar = var(res[abs(res)<maxRes]) 
y = 25.2 - abs(((366/2)-t))*0.096
# b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
b0sim = y + gen.arma.wge(n=366, phi=0.8,vara=meanVar, plot=FALSE)
plot(meanVec,col = "red")
points(b0sim,col = "blue") # This looks alright

# Try a full simulation without the overlap
dayStart = 1
dayStop = 366
y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
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
absA1vec = abs(A1vec)
freq1 = 1/144
for (day in dayStart:dayStop){
  t = ((day-1)*6*24+1):(day*6*24)
  psi1 = psi1mean
  A1 = abs(A1vec[day])
  # b0 = meanVec[day] 
  b0 = b0sim[day]
  st = b0 + A1 * cos(2 * pi * freq1 * t + psi1)
  y_gen_vpmax[t] = st
}
plot(y_gen_vpmax[1:2000])

# Now try adding the overlap part
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
zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
y_gen_vpmax = y_gen_vpmax + zt
plot(y_gen_vpmax[1:1000])

# Check the plots
numDays = 366
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
numDays = 10
ymin = min(c(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],y_gen_vpmax[1:(numDays*6*24)])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[1:(numDays*6*24)],col = "blue")
points(y_gen_vpmax[1:(numDays*6*24)]+zt[1:(numDays*6*24)],col = "green")
legend('topleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
startDay = 366/2
numDays = 10
ymin = min(c(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]))
ymax = max(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)])
plot(weather$VPmax..mbar.[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "red",ylim=c(ymin,ymax))
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "blue")
points(y_gen_vpmax[(1+6*24*startDay):(6*24*startDay+numDays*6*24)]+zt[(1+6*24*startDay):(6*24*startDay+numDays*6*24)],col = "green")
legend('bottomleft',legend = c("VPMax", "Sim - no noise", "Sim - with noise"), col = c("red", "blue", "green"),pch=1)
# This looks somewhat reasonable, but a little low.
# Maybe try making it higher, and figuring out a way to prevent it from going negative next
# Also, try to plot the ACFs, just to make sure this helped at all.

# See if ACFs are better
numSims = 10
y_gen_vpmax_mult = data.frame(matrix(0, nrow = 6*24*366, ncol = numSims))  
colnames(y_gen_vpmax_mult) <- paste0("sim", 1:numSims)
phiForb0 = 0.8
for (seed in 1:numSims){
  print(seed)
  set.seed(seed)
  # b0sim = y + rnorm(n=366,mean=0,sd=sqrt(meanVar))
  b0sim = y + gen.arma.wge(n=366, phi=phiForb0,vara=meanVar, plot=FALSE)
  dayStart = 1
  dayStop = 366
  y_gen_vpmax = numeric(6*24*(dayStop-dayStart+1))
  for (day in dayStart:dayStop){
    t = ((day-1)*6*24+1):(day*6*24)
    psi1 = psi1mean
    A1 = abs(A1vec[day])
    b0 = b0sim[day] # use simulated b0
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
  zt = gen.arma.wge(n = 144*366, phi = phi, theta = 0, vara = vara, sn = 0, plot=FALSE)
  y_gen_vpmax_mult[, seed] = y_gen_vpmax + zt
}
ACF = acf(weather$VPmax..mbar., plot = "FALSE")
plot(ACF$lag, ACF$acf , type = "l", lwd = 6)
for (sim in 1:numSims){
  ACF2 = acf(y_gen_vpmax_mult[,sim], plot = "FALSE")
  lines(ACF2$lag, ACF2$acf, lwd = 2, col = "red")
} # This is a little better.  It has more variety of ACFs, but that variety does contain the "real" ACF

