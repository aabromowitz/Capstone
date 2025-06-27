# We wanted to see what happens if move the TFT forecasts to line up with the last forecast

# Pull in sims
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_1_to_50_varTemp001.csv"
sims <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/Signal.csv"
sig <- read.csv(file_path, header = TRUE)

# Create metric dataframe
numSims = 50
metric_names <- c("MSE", "MAE", "MSE.MA", "MAE.MA", "MSE.SIG", "MAE.SIG", "MSE.SIG.MA", "MAE.SIG.MA")
horizons <- c(96, 192, 336, 720)
columns <- c("Run.Number", 
             as.vector(sapply(metric_names, function(metric) paste0(metric, ".", horizons))))
metrics_df <- data.frame(matrix(nrow = numSims, ncol = length(columns)))
colnames(metrics_df) <- columns

# Loop through all of the sims
totLen = length(sims[,1])
lastTrain = totLen - 720
firstTest = lastTrain + 1
for (iSim in 1:numSims){
  file_path = paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_", iSim, "_1_1.csv",sep="")
  fore = read.csv(file_path, header = TRUE)
  sim = sims[,iSim]
  for (iH in 1:length(horizons)){
    h = horizons[iH]
    metrics_df$Run.Number[iSim] = iSim
    # MSE
    shift = sim[lastTrain] - fore$forecast[1]
    newFore = fore$forecast[1:h] + shift
    colname = paste0("MSE.",h)
    metrics_df[[colname]][iSim] = mean((sim[firstTest:(firstTest + h - 1)] - newFore)^2)
    # MAE
    colname = paste0("MAE.",h)
    metrics_df[[colname]][iSim] = mean(abs(sim[firstTest:(firstTest + h - 1)] - newFore))
    # MSE.MA
    shift = sim[lastTrain] - fore$smoothed_forcast[1]
    newFore = fore$smoothed_forcast[1:h] + shift
    colname = paste0("MSE.MA.",h)
    metrics_df[[colname]][iSim] = mean((sim[firstTest:(firstTest + h - 1)] - newFore)^2)
    # MAE.MA
    colname = paste0("MAE.MA.",h)
    metrics_df[[colname]][iSim] = mean(abs(sim[firstTest:(firstTest + h - 1)] - newFore))
    # MSE.SIG
    shift = sim[lastTrain] - fore$forecast[1]
    newFore = fore$forecast[1:h] + shift
    colname = paste0("MSE.SIG.",h)
    metrics_df[[colname]][iSim] = mean((sig$T..degC.[firstTest:(firstTest + h - 1)] - newFore)^2)
    # MAE.SIG
    colname = paste0("MAE.SIG.",h)
    metrics_df[[colname]][iSim] = mean(abs(sig$T..degC.[firstTest:(firstTest + h - 1)] - newFore))
    # MSE.SIG.MA
    shift = sim[lastTrain] - fore$smoothed_forcast[1]
    newFore = fore$smoothed_forcast[1:h] + shift
    colname = paste0("MSE.SIG.MA.",h)
    metrics_df[[colname]][iSim] = mean((sig$T..degC.[firstTest:(firstTest + h - 1)] - newFore)^2)
    # MAE.SIG.MA
    colname = paste0("MAE.SIG.MA.",h)
    metrics_df[[colname]][iSim] = mean(abs(sig$T..degC.[firstTest:(firstTest + h - 1)] - newFore))
  }
}

# Save off dataframe
write.csv(metrics_df, "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/TFT_Shifted_Errors_1_to_50_varTemp001.csv", row.names = FALSE)

################################################################################
# Compare errors

# Pull in ARIMA errors
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_With_Trimmed_1_to_50_varTemp001.csv"
arima_errors <- read.csv(file_path, header = TRUE)

# Tests for comparing MSEs
mean(arima_errors$mse_96) # 3.565827
mean(metrics_df$MSE.96) # 5.989477
mean(metrics_df$MSE.MA.96) # 2.765151
wilcox.test(arima_errors$mse_96, metrics_df$MSE.96, paired=TRUE) # p-value = 0.009954
wilcox.test(arima_errors$mse_96, metrics_df$MSE.MA.96, paired=TRUE) # p-value = 0.2883

mean(arima_errors$mse_192) # 4.278645
mean(metrics_df$MSE.192) # 7.660755
mean(metrics_df$MSE.MA.192) # 4.421161
wilcox.test(arima_errors$mse_192, metrics_df$MSE.192, paired=TRUE) # p-value = 0.0006554
wilcox.test(arima_errors$mse_192, metrics_df$MSE.MA.192, paired=TRUE) # p-value = 0.3745

mean(arima_errors$mse_336) # 5.572522
mean(metrics_df$MSE.336) # 10.1456
mean(metrics_df$MSE.MA.336) # 6.555292
wilcox.test(arima_errors$mse_336, metrics_df$MSE.336, paired=TRUE) # p-value = 0.0001001
wilcox.test(arima_errors$mse_336, metrics_df$MSE.MA.336, paired=TRUE) # p-value = 0.1225

mean(arima_errors$mse_720) # 8.979426
mean(metrics_df$MSE.720) # 14.98075
mean(metrics_df$MSE.MA.720) # 11.2277
wilcox.test(arima_errors$mse_336, metrics_df$MSE.720, paired=TRUE) # p-value = 6.985e-07
wilcox.test(arima_errors$mse_336, metrics_df$MSE.MA.720, paired=TRUE) # p-value = 1.983e-05

################################################################################
# Plots
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/ARIMA_Forecasts_With_Trimmed_1_to_50_varTemp001.csv"
arima_fores = read.csv(file_path, header = TRUE)
num_plots = numSims
# num_plots = 1
len_total = length(sims[,1])
h = 720
lastTrain = len_total - h
vlines = c(-720,0, 96, 192, 336)
b = h
dest.dir <- "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Plots/Jun7WithShiftedTFT"
for (iPlot in 1:num_plots){
  sim_num = iPlot
  arima_fore = arima_fores[,sim_num]
  file_path = paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/1_to_50_varTemp001_", sim_num, "_1_1.csv",sep="")
  fore = read.csv(file_path, header = TRUE)
  shift = sims[lastTrain,sim_num] - fore$forecast[1]
  tft_fore = fore$forecast + shift
  shift = sims[lastTrain,sim_num] - fore$smoothed_forcast[1]
  tft_fore_smooth = fore$smoothed_forcast + shift
  sim = sims[(len_total-h-b+1):len_total,sim_num]
  y_min = min(c(sim,arima_fore,tft_fore,tft_fore_smooth))
  y_max = max(c(sim,arima_fore,tft_fore,tft_fore_smooth))
  out.path <- file.path(dest.dir, paste0("forecast_sim_", sim_num, ".png"))
  png(filename = out.path, width = 800, height = 600)
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
  points(seq(1+b,h+b,1),arima_fore,col = "red")
  points(seq(1+b,h+b,1),tft_fore,col = "magenta")
  points(seq(1+b,h+b,1),tft_fore_smooth, col='green')
  legend('top',legend = c("Simulation", "ARIMA Forecast", "TFT Forecast", "TFT Forecast - Smoothed"), 
         col = c("blue", "red", "magenta", "green"),pch=1)
  ticks <- axTicks(1)
  axis(1, at = ticks, labels = round((ticks + len_total-h-b) / 144,0))
  dev.off()
}
