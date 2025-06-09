# First, let's put all the ARIMA forecasts into one file.

baseDir = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_"
numSims = 50
forecasts = data.frame(matrix(0, nrow = 720, ncol = numSims))  
colnames(forecasts) <- paste0("sim", 1:numSims)
for (iSim in 1:numSims){
  file_path = paste(baseDir,iSim,"_1_1.csv",sep="")
  f <- read.csv(file_path, header = TRUE)
  forecasts[,iSim] = f$forecast
}
write.csv(forecasts,"C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_ARIMA_combined.csv")

################################################################################

# Now let's plot some of the interesting forecasts

# Common variables
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_1_to_50_attempt2.csv"
sims <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/Signal.csv"
sig_tot = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/ARIMA_Forecasts_1_to_50_varTemp001.csv"
arima_fores = read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_ARIMA_combined.csv"
tft_fores = read.csv(file_path, header = TRUE)
len_total = length(sims[,1])
h = 720
vlines = c(-720,0, 96, 192, 336)
b = h*2
sig = sig_tot[(len_total-h-b+1):len_total,2]
interesting_plots = c(1,2,3,6,9,10,11,12,16,19,20,23,24,32,38,43,44,45,47,49)
num_plots = length(interesting_plots)

# Loop through plots
for (iPlot in 1:num_plots){
  sim_num = interesting_plots[iPlot]
  arima_fore = arima_fores[,sim_num]
  tft_fore = tft_fores[,sim_num]
  sim = sims[(len_total-h-b+1):len_total,sim_num]
  y_min = min(c(sim,arima_fore,tft_fore,sig))
  y_max = max(c(sim,arima_fore,tft_fore,sig))
  plot(seq(1,h+b,1),sim, col='blue',ylim=c(y_min,y_max),
       main=paste('Forecast Comparison - Simulation', sim_num), xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
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
  points(seq(1,h+b,1),sig, col='green')
  points(seq(1+b,h+b,1),arima_fore,col = "red")
  points(seq(1+b,h+b,1),tft_fore,col = "magenta")
  legend('top',legend = c("Simulation", "Signal", "ARIMA Forecast", "ARIMA Forecast"), 
         col = c("blue", "green", "red", "magenta"),pch=1)
  ticks <- axTicks(1)
  axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))
}

################################################################################

# Look at the distributions of errors

# Pull errors
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_1_to_50_varTemp001.csv"
arima_errors <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_1_to_50_varTemp001.csv"
tft_errors <- read.csv(file_path, header = TRUE)

# Histograms comparing MSE errors
err_cut = 10
arimaMse96 = arima_errors$mse_96
arimaMse192 = arima_errors$mse_192
arimaMse336 = arima_errors$mse_336
arimaMse720 = arima_errors$mse_720
tftMse96 = tft_errors$MSE.96
tftMse192 = tft_errors$MSE.192
tftMse336 = tft_errors$MSE.336
tftMse720 = tft_errors$MSE.720
diffMse96 = tftMse96 - arimaMse96
diffMse192 = tftMse192 - arimaMse192
diffMse336 = tftMse336 - arimaMse336
diffMse720 = tftMse720 - arimaMse720
max_err = ceiling(max(arimaMse96,tftMse96)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse96,breaks=cuts,col=rgb(1,0,0,0.5),main='MSE Errors for 96 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse96,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMse192,tftMse192)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse192,breaks=cuts,col=rgb(1,0,0,0.5),main='MSE Errors for 192 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse192,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMse336,tftMse336)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse336,breaks=cuts,col=rgb(1,0,0,0.5),main='MSE Errors for 336 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse336,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMse720,tftMse720)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse720,breaks=cuts,col=rgb(1,0,0,0.5),main='MSE Errors for 720 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse720,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
err_cut = 5
cuts = seq(from = floor(min(diffMse96)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMse96)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMse96,breaks=cuts,col='green',main='MSE Error Differences between ARIMA and ARIMA for 96 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMse192)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMse192)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMse192,breaks=cuts,col='green',main='MSE Error Differences between ARIMA and ARIMA for 192 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMse336)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMse336)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMse336,breaks=cuts,col='green',main='MSE Error Differences between ARIMA and ARIMA for 336 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMse720)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMse720)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMse720,breaks=cuts,col='green',main='MSE Error Differences between ARIMA and ARIMA for 720 Horizon',xlab="Difference in MSE")

# Tests for comparing MSEs
wilcox.test(arimaMse96, tftMse96) # p-value = 1.63e-06
wilcox.test(arimaMse192, tftMse192) # p-value = 1.63e-06
wilcox.test(arimaMse336, tftMse336) # p-value = 1.109e-05
wilcox.test(arimaMse720, tftMse720) # p-value = 9.681e-05
wilcox.test(arimaMse96, tftMse96, paired=TRUE) # p-value = 4.971e-06
wilcox.test(arimaMse192, tftMse192, paired=TRUE) # p-value = 1.53e-06
wilcox.test(arimaMse336, tftMse336, paired=TRUE) # p-value = 2.14e-06
wilcox.test(arimaMse720, tftMse720, paired=TRUE) # p-value = 1.281e-05

# Histograms comparing MSEs to sim vs signal for ARIMA
err_cut = 10
tftMseSig96 = tft_errors$MSE.SIG.96
tftMseSig192 = tft_errors$MSE.SIG.192
tftMseSig336 = tft_errors$MSE.SIG.336
tftMseSig720 = tft_errors$MSE.SIG.720
diffMseSigTft96 = tftMse96 - tftMseSig96
diffMseSigTft192 = tftMse192 - tftMseSig192
diffMseSigTft336 = tftMse336 - tftMseSig336
diffMseSigTft720 = tftMse720 - tftMseSig720
max_err = ceiling(max(tftMseSig96,tftMse96)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig96,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 96 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse96,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Signal", "Sim"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(tftMseSig192,tftMse192)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig192,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 192 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse192,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Signal", "Sim"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(tftMseSig336,tftMse336)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig336,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 336 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse336,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Signal", "Sim"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(tftMseSig720,tftMse720)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig720,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 720 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(tftMse720,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Signal", "Sim"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
err_cut = 5
cuts = seq(from = floor(min(diffMseSigTft96)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigTft96)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigTft96,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 96 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigTft192)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigTft192)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigTft192,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 192 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigTft336)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigTft336)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigTft336,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 336 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigTft720)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigTft720)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigTft720,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 720 Horizon',xlab="Difference in MSE")

# Tests for comparing ARIMA MSEs
wilcox.test(tftMseSig96, tftMse96, paired=TRUE) # p-value = 0.04167
wilcox.test(tftMseSig192, tftMse192, paired=TRUE) # p-value = 0.01499
wilcox.test(tftMseSig336, tftMse336, paired=TRUE) # p-value = 0.00236
wilcox.test(tftMseSig720, tftMse720, paired=TRUE) # p-value = 0.000132

# Histograms comparing MSEs to sim vs signal for ARIMA
err_cut = 10
arimaMseSig96 = arima_errors$mse_sig_96
arimaMseSig192 = arima_errors$mse_sig_192
arimaMseSig336 = arima_errors$mse_sig_336
arimaMseSig720 = arima_errors$mse_sig_720
diffMseSig96 = tftMseSig96 - arimaMseSig96
diffMseSig192 = tftMseSig192 - arimaMseSig192
diffMseSig336 = tftMseSig336 - arimaMseSig336
diffMseSig720 = tftMseSig720 - arimaMseSig720
max_err = ceiling(max(arimaMseSig96,tftMseSig96)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig96,breaks=cuts,col=rgb(1,0,0,0.5),main='Signal MSE Errors for 96 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig96,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig192,tftMseSig192)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig192,breaks=cuts,col=rgb(1,0,0,0.5),main='Signal MSE Errors for 192 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig192,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig336,tftMseSig336)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig336,breaks=cuts,col=rgb(1,0,0,0.5),main='Signal MSE Errors for 336 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig336,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig720,tftMseSig720)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(tftMseSig720,breaks=cuts,col=rgb(1,0,0,0.5),main='Signal MSE Errors for 720 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig720,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("ARIMA", "ARIMA"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
err_cut = 5
cuts = seq(from = floor(min(diffMseSig96)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSig96)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSig96,breaks=cuts,col='green',main='Signal MSE Error Differences between ARIMA and ARIMA for 96 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSig192)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSig192)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSig192,breaks=cuts,col='green',main='Signal MSE Error Differences between ARIMA and ARIMA for 192 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSig336)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSig336)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSig336,breaks=cuts,col='green',main='Signal MSE Error Differences between ARIMA and ARIMA for 336 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSig720)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSig720)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSig720,breaks=cuts,col='green',main='Signal MSE Error Differences between ARIMA and ARIMA for 720 Horizon',xlab="Difference in MSE")

# Tests for comparing Signal MSEs
wilcox.test(tftMseSig96, arimaMseSig96, paired=TRUE) # p-value = 0.569
wilcox.test(tftMseSig192, arimaMseSig192, paired=TRUE) # p-value = 0.569
wilcox.test(tftMseSig336, arimaMseSig336, paired=TRUE) # p-value = 0.6852
wilcox.test(tftMseSig720, arimaMseSig720, paired=TRUE) # p-value = 0.5054

# Histograms comparing MSEs to sim vs signal for ARIMA
err_cut = 10
diffMseSigArima96 = arimaMse96 - arimaMseSig96
diffMseSigArima192 = arimaMse192 - arimaMseSig192
diffMseSigArima336 = arimaMse336 - arimaMseSig336
diffMseSigArima720 = arimaMse720 - arimaMseSig720
max_err = ceiling(max(arimaMseSig96,arimaMse96)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse96,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 96 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig96,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Sim", "Signal"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig192,arimaMse192)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse192,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 192 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig192,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Sim", "Signal"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig336,arimaMse336)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse336,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 336 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig336,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Sim", "Signal"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
max_err = ceiling(max(arimaMseSig720,arimaMse720)+err_cut)
cuts = seq(from = 0, to = max_err, by = err_cut)
hist(arimaMse720,breaks=cuts,col=rgb(1,0,0,0.5),main='ARIMA MSE Errors for 720 Horizon',xlab="MSE",xlim=c(0,max_err))
hist(arimaMseSig720,breaks=cuts,col=rgb(0,1,0,0.5),add=TRUE)
legend("topright", legend = c("Sim", "Signal"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), pch=16)
err_cut = 5
cuts = seq(from = floor(min(diffMseSigArima96)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigArima96)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigArima96,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 96 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigArima192)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigArima192)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigArima192,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 192 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigArima336)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigArima336)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigArima336,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 336 Horizon',xlab="Difference in MSE")
cuts = seq(from = floor(min(diffMseSigArima720)/err_cut)*err_cut-err_cut/2, to = ceiling(max(diffMseSigArima720)/err_cut)*err_cut+err_cut/2, by = err_cut)
hist(diffMseSigArima720,breaks=cuts,col='green',main='ARIMA MSE Error Differences between Signal and Simulation for 720 Horizon',xlab="Difference in MSE")

# Tests for comparing ARIMA MSEs
wilcox.test(arimaMse96, arimaMseSig96, paired=TRUE) # p-value = 0.0007033
wilcox.test(arimaMse192, arimaMseSig192, paired=TRUE) # p-value = 0.001263
wilcox.test(arimaMse336, arimaMseSig336, paired=TRUE) # p-value = 0.01143
wilcox.test(arimaMse720, arimaMseSig720, paired=TRUE) # p-value = 0.556

################################################################################

# Redo some plots

# Pull in weather data and simulations
file_path = "https://raw.githubusercontent.com/aabromowitz/Capstone/refs/heads/main/weather_interpolation.csv"
weather <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/vpmax_1_to_50_varTemp001.csv"
vpmax_sims <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_1_to_50_varTemp001.csv"
temp_sims <- read.csv(file_path, header = TRUE)

# VPMax plot
dayStop = 366
dayStart = 1
numDays = dayStop - dayStart + 1
ymin = min(weather$VPmax..mbar.[1:(numDays*6*24)],vpmax_sims[1:(numDays*6*24),1])
ymax = max(weather$VPmax..mbar.[1:(numDays*6*24)],vpmax_sims[1:(numDays*6*24),1])
plot(weather$VPmax..mbar.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Realization of Maximum Vapor Pressure - Full Year', xlab = 'Days',ylab = 'Max Vapor Pressure (mbar)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(vpmax_sims[1:(numDays*6*24),1],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)

# Temp plot
dayStop = 366
dayStart = 1
numDays = dayStop - dayStart + 1
ymin = min(weather$T..degC.[1:(numDays*6*24)],temp_sims[1:(numDays*6*24),1])
ymax = max(weather$T..degC.[1:(numDays*6*24)],temp_sims[1:(numDays*6*24),1])
plot(weather$T..degC.[1:(numDays*6*24)],col = "blue",ylim=c(ymin,ymax),
     main='Air Temperature Simulation', xlab = 'Days',ylab = 'Air Temperature (C)',xaxt='n')
ticks <- axTicks(1)  # get default x-axis tick positions
axis(1, at = ticks, labels = round(ticks / 144,0))
points(temp_sims[1:(numDays*6*24),1],col = "red")
legend('topleft',legend = c("Real Data", "Simulated Data"), 
       col = c("blue", "red"),pch=1)