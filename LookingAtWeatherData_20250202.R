# Demoing
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
temp <- weather$T..degC.
library(tswge)

n = 366*144
vara = .1
phi = .998
freq = c(1/144,1/(144*366))
psi = c(-3.926991, -3.787508) 
sn = 10 # A number I've never tried before, so hopefully it looks alright
# coef = c(1.836553, 6.411447) 
coef = c(1.8, 6.5) 
gendat = gen.sigplusnoise.wge(n=n,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn) + mean(temp) + 4 # 10.81848
plotts.wge(temp)
plotts.wge(gendat)
max(gendat)
max(temp)
min(gendat)
min(temp)
plotts.wge(temp[1:2000])
plotts.wge(gendat[1:2000])

################################################################################

# Pull in weather.csv data
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Datasets/weather.csv"
weather <- read.csv(file_path, header = TRUE)
temp <- weather$T..degC.
library(tswge)
plotts.wge(temp)
# It is only a year's worth of data, so there isn't annual periodicity

# Look at the ACF and parzen window of the data
plotts.sample.wge(temp)
# Ha, this was taking forever, so I'll try a subset of the data

# Try looking at data subsets
temp10000 <- temp[1:10000]
plotts.sample.wge(temp10000)
temp100 <- temp[1:100]
plotts.wge(temp100) # oh yeah, it's not every hour, it's every 10 minutes
# so an hour would be 6 datapoints and a day would be 6*24 = 144 datapoints
# The full year would be 144*365 = 52,560, which is the size of the whole dataset
temp1000 <- temp[1:1000]
plotts.wge(temp1000) # now the periodicity is a bit more obvious
# 1/144 = 0.006944444

# Look at the parzen window for the 0.007
parzen.wge(temp1000) # 0.007 isn't visible
parzen.wge(temp1000,trunc = 500) # There are a lot of bumps here

# I remember there was a way to remove really high periods from the data, more than 12
p = 144
phi1 = 2*cos(pi/(p/2))
phi1 # 1.998096
phi2 = -1
f=1/(2*pi)*acos(phi1/(2*sqrt(-phi2)))
print(f) # 0.006944444, which matches the 1/144

# Try differencing out the 144 period
d10000 = artrans.wge(temp10000,phi.tr=c(phi1,-1))
dev.off()
acf(d10000,lag.max = 1000)
dfull = artrans.wge(temp,phi.tr=c(2*cos(pi/(p/2)),-1))
plotts.wge(dfull)

# I wonder if this passes one of those white noise tests?
ljung.wge(dfull,K=24) # default is K=24, want p>0.05 to show that it is white noise
ljung.wge(dfull,K=48) # p = 0 for both, so I guess that means it isn't white noise
plotts.wge(d10000)
ljung.wge(d10000,K=24) 
ljung.wge(d10000,K=48) # p = 0 for smaller dataset
d1000 = artrans.wge(temp1000,phi.tr=c(phi1,-1))
plotts.wge(d1000)
ljung.wge(d1000,K=24) 
ljung.wge(d1000,K=48) # p = 0 for 1000
plotts.sample.wge(d1000)

# Looking at aic
aic5.wge(dfull,p=0:10,q=0:4,type='aic') # taking too long
est.ar.wge(dfull,p=1) # -0.2362896
est.ar.wge(temp,p=1) # 0.9994192

# I wonder how useful this even is.  We know that the data would have two periods, 144 and ~52,560
max(temp) # 34.8
which.max(temp) # 31766
weather$date[31766] # "2020-08-08 15:40:00"
temp[31766-144*366/2+144/2] # 0.56, 7 away from min
min(temp) # -6.44
which.min(temp) # 3065
weather$date[3065] # "2020-01-22 06:50:00", not quite 6 months apart, but pretty close
temp[3065+144*366/2-144/2] # 20.04, over 14 away from max
# I feel like that means 8/8 is closer to a "typical" max than 1/22 would be from a "typical" min

################################################################################
# What are the inputs to the function to generate sine waves again?
?gen.sigplusnoise.wge # x(t)=coef[1]*cos(2*pi*freq[1]*t+psi[1])+coef[2]*cos(2*pi*freq[2]*t+psi[2])+a(t)
# n = 52696, to match the weather dataset
# b0 = mean(temp) = 10.81824
# b1 = 0, shouldn't really have an uphill or downhill slope, just a couple of sine curves
# coef is the numerb multiplied by 1 to get the range
# freq = 1/144 and 1/(144*365)
# phi = 0?  I feel like there's probably a better value to choose for this
# vara = variance of the data when I difference out the 144 period = var(dfull) = 0.05812857

# What time is the hottest temperature of the day?
num_days = length(temp)/144 # huh, I guess this year is a leap year, also maybe a day is missing somewhere in here
max_inds = numeric(floor(num_days))
for (ii in 1:length(max_inds))
{
  sub_temp = temp[((ii-1)*144+1):(ii*144)]
  max_inds[ii] = which.max(sub_temp)
}
median(max_inds) # 85, which is around 2 pm, that's believable

# See if there is a missing time
dates = weather$date
datediffs = numeric(length(dates)-1)
for (ii in 1:length(datediffs))
{
  time1 <- as.POSIXct(dates[ii])
  time2 <- as.POSIXct(dates[ii+1])
  datediffs[ii] <- as.numeric(difftime(time2, time1, units = "mins"))
} # that didn't really work.  there were some weird conversion problems that messed up this approach.  sigh

# It looks like on 5/29, it jumps from 9:30 to 11:10
# There was also an issue on 5/12 where 6:00 was repeated, but I just deleted the extra time
# This is day 150
tempNo529 = c(temp[1:(149*144)],temp[(150*144-8):length(temp)])
length(tempNo529)/144 # 365, yay
num_days = 365 # remove 5/29
max_inds = numeric(num_days) 
for (ii in 1:num_days) {
  sub_temp = tempNo529[((ii-1)*144+1):(ii*144)]
  max_inds[ii] = which.max(sub_temp)
}
median(max_inds) # 90, that looks like 2:50 PM, which makes some sense

# Now get the average max temp
max_temps = numeric(num_days) 
for (ii in 1:num_days) {
  sub_temp = tempNo529[((ii-1)*144+1):(ii*144)]
  max_temps[ii] = max(sub_temp)
}
mean(max_temps) # 15.33167, that's lower than I'd expect, hm

# Try the same thing out for min dates
min_inds = numeric(num_days) 
for (ii in 1:num_days) {
  sub_temp = tempNo529[((ii-1)*144+1):(ii*144)]
  min_inds[ii] = which.min(sub_temp)
}
median(min_inds) # 35, that looks like 5:50 AM, which makes some sense

# Now get the average min temp
min_temps = numeric(num_days) 
for (ii in 1:num_days) {
  sub_temp = tempNo529[((ii-1)*144+1):(ii*144)]
  min_temps[ii] = min(sub_temp)
}
mean(min_temps) # 6.148904, really, the temperature only changes by 9.182767 C or 16.52898 F on average per day?
mean(max_temps) - mean(min_temps) # 9.182767
coef1 = (mean(max_temps) - mean(min_temps))/2

# I guess for the annual range, we could do max - min for the whole year to get the range
# Then subtract off the daily variation of 9.18 to get the annual variation
max(temp) # 34.8
which.max(temp) # 31765
min(temp) # -6.44
which.min(temp) # 3065
coef2 = (max(temp) - min(temp) - coef1*2)/2 # 16.02862 C, or 28.85151 F, so close to 60 F

# For psi, we know that cos(0) = 1, which is the max
# so 2*pi*1/144*t+psi1 = 0 at the max
# I believe that means we want for -#days to 8/8 to be our phi for the year
# This would be 31 + 29 + 31 + 30 + 31 + 30 + 31 + 7 = 220 days
# So 220 * 144 + 90 = 31770, which is pretty close to the actual max
psi1 = -2*pi/(144/90) # 3.926991
psi2 = -2*pi/(144*366/31770) # 3.787508

# Figure out coef
diffs72 = numeric(length(temp)-144/2)
for (ii in 1:length(diffs72))
{  
  diffs72[ii] = abs(temp[ii]-temp[ii+72])
}
diffs72[length(diffs72)] # it's a number
mean(diffs72) # 4.785391, I'm not sure if I trust this

n = 366*144
coef = c(coef1,coef2) # 4.591384 16.028616
freq = c(1/144,1/(144*366))
psi = c(psi1,psi2) # -3.926991 -3.787508
# dfullNo529 = artrans.wge(tempNo529,phi.tr=c(2*cos(pi/(144/2)),-1))
# vara = var(dfullNo529) # 0.05787472
vara = var(temp) # 55.77919
sn = 1
gendat = gen.sigplusnoise.wge(n=n,coef=coef,freq=freq,psi=psi,vara=vara,sn=1) + mean(temp) # 10.81848
plotts.wge(temp)
plotts.wge(gendat) # I feel like one is too much variance, and the other is too little variance

# Trying with less variation, and some phi values
vara = sqrt(var(temp)) # 7.468546
phi = .9
sn = 2
gendat = gen.sigplusnoise.wge(n=n,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn) + mean(temp) # 10.81848
plotts.wge(temp)
plotts.wge(gendat)

# This actually looks pretty similar, but with a little too much range
vara = 1
phi = .99
sn = 3
gendat = gen.sigplusnoise.wge(n=n,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn) + mean(temp) # 10.81848
plotts.wge(temp)
plotts.wge(gendat)

# This actually looks pretty close
vara = .3
phi = .993
sn = 4
coef = c(coef1*.4,coef2*.4) # 1.836553 6.411447
gendat = gen.sigplusnoise.wge(n=n,coef=coef,freq=freq,psi=psi,phi=phi,vara=vara,sn=sn) + mean(temp) + 4 # 10.81848
plotts.wge(temp)
plotts.wge(gendat)
max(gendat)
max(temp)
min(gendat)
min(temp)

################################################################################
# What's the correlation like between variables?
cor_matrix <- cor(weather) # ERROR: Apparently it has to be numeric
print(cor_matrix)

# Try to get only numeric variables
library(tidyverse)
weather_numeric = select_if(weather, is.numeric) # apparently the only non-numeric thing were the datetimes
cor_matrix <- cor(weather_numeric) 
print(cor_matrix)

# As expected, a lot of these are really high correlation.  Maybe a few can be removed?
threshold = 0.9
high_cor_pairs <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)
high_cor_pairs <- as.data.frame(high_cor_pairs)
high_cor_pairs <- high_cor_pairs[high_cor_pairs$row < high_cor_pairs$col, ]  # Remove duplicates
high_cor_pairs$Var1 <- rownames(cor_matrix)[high_cor_pairs$row]
high_cor_pairs$Var2 <- colnames(cor_matrix)[high_cor_pairs$col]
high_cor_pairs$Correlation <- cor_matrix[cbind(high_cor_pairs$row, high_cor_pairs$col)]
high_cor_pairs <- high_cor_pairs[, c("Var1", "Var2", "Correlation")]
print(high_cor_pairs)
# T..degC., Tpot..K.: makes sense, I think these are just temps in C vs K
# T..degC., VPmax..mbar.: this is kinda interesting.  I think mbar is a unit of pressure.  
#   Looking up VPMax, it says it's something called saturation vapor pressure
#   Apparently there is a formula relating VPmax and T, but it's something like V = a*exp(b*T/(T+c))
# sh is specific humidity
# VPact is similar to VPmax, but instead of the max vapor pressure the air can hold, 
#   Vpact is the actual vapor pressure at the moment
# Both VPact and sh measure water vapor in the air, but in different ways
