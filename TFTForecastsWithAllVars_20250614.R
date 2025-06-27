# Common variables
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_1_to_50_varTemp001.csv"
sims_temp <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/vpmax_1_to_50_varTemp001.csv"
sims_vpmax <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/tdew_1_to_50_varTemp001.csv"
sims_tdew <- read.csv(file_path, header = TRUE)
len_total = length(sims[,1])
h = 720
vlines = c(-720,0, 96, 192, 336)
b = h*2

# Plot all 3 for sim 1 for TFT
sim_num = 1
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_1_1_1_forecast.csv"
tft_1_fores = read.csv(file_path, header = TRUE)

# Temp
sim = sims_temp[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_1_fores$forecast))
y_max = max(c(sim,tft_1_fores$forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('Temp TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_1_fores$forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))

# VPMax
sim = sims_vpmax[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_1_fores$vpmax_forecast))
y_max = max(c(sim,tft_1_fores$vpmax_forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('VPMax TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Maximum Vapor Pressure (mbar)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_1_fores$vpmax_forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))

# TDew
sim = sims_tdew[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_1_fores$tdew_forecast))
y_max = max(c(sim,tft_1_fores$tdew_forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('TDew TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Dew Temperature (c)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_1_fores$tdew_forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))


# Plot all 3 for sim 19 for TFT
# b=h*2
b=h*20
sim_num = 19
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_19_1_1_forecast.csv"
tft_19_fores = read.csv(file_path, header = TRUE)

# Temp
sim = sims_temp[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_19_fores$forecast))
y_max = max(c(sim,tft_19_fores$forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('Temp TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_19_fores$forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))

# VPMax
sim = sims_vpmax[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_19_fores$vpmax_forecast))
y_max = max(c(sim,tft_19_fores$vpmax_forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('VPMax TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Maximum Vapor Pressure (mbar)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_19_fores$vpmax_forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))

# TDew
sim = sims_tdew[(len_total-h-b+1):len_total,sim_num]
y_min = min(c(sim,tft_19_fores$tdew_forecast))
y_max = max(c(sim,tft_19_fores$tdew_forecast))
plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
     main=paste('TDew TFT Forecast - Simulation', sim_num), 
     xlab = 'Days',ylab = 'Dew Temperature (c)',xaxt='n')
for (iLine in 1:length(vlines)){
  vline = vlines[iLine] + b
  if (vline > 0){
    abline(v=vline,col="black",lty=2,lwd=0.5)
  }
}
if(-720+b > 0){
  abline(h=sim[-720+b],col="black",lty=2,lwd=0.5)
}
abline(h=sim[b],col="black",lty=2,lwd=0.5)
points(seq(1+b,h+b,1),tft_19_fores$tdew_forecast,col = "magenta")
legend('top',legend = c("Simulation", "Forecast"), 
       col = c("blue", "magenta"),pch=1)
ticks <- axTicks(1)
axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))

################################################################################

# Looking at the correlation for a simulation
library(tswge)
sim_num = 1
train_len = len_total - h
season = 6*24
vpmax_sim = sims_vpmax[,sim_num]
tdew_sim = sims_tdew[,sim_num]
temp_sim = sims_temp[,sim_num]
temp_trim <- temp_sim[1:train_len]
vpmax_trim <- vpmax_sim[1:train_len]
tdew_trim <- tdew_sim[1:train_len]
fit = arima(temp_trim, order = c(4, 0, 1), xreg = cbind(vpmax_trim, tdew_trim)) 
sum_fit <- summary(fit)
coefs <- sum_fit$coef[c('vpmax_trim','tdew_trim')] 
ses <- c(sqrt(sum_fit$var.coef['vpmax_trim','vpmax_trim']), 
         sqrt(sum_fit$var.coef['tdew_trim','tdew_trim']))
z_vals <- coefs / ses
p_vals <- 2 * (1 - pnorm(abs(z_vals)))

sim_num = 2
train_len = len_total - h
season = 6*24
vpmax_sim = sims_vpmax[,sim_num]
tdew_sim = sims_tdew[,sim_num]
temp_sim = sims_temp[,sim_num]
temp_trim <- temp_sim[1:train_len]
vpmax_trim <- vpmax_sim[1:train_len]
tdew_trim <- tdew_sim[1:train_len]
fit = arima(temp_trim, order = c(4, 0, 1), xreg = cbind(vpmax_trim, tdew_trim)) 
summary(fit)
sum_fit <- summary(fit)
coefs <- sum_fit$coef[c('vpmax_trim','tdew_trim')] 
ses <- c(sqrt(sum_fit$var.coef['vpmax_trim','vpmax_trim']), 
         sqrt(sum_fit$var.coef['tdew_trim','tdew_trim']))
z_vals <- coefs / ses
z_vals
p_vals <- 2 * (1 - pnorm(abs(z_vals)))