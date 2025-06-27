# Packages
library(tswge)
library(readr)

# Parameters
numSims = 1 # set to 50 for real sims
startSim = 1
toPlot = TRUE # set this to TRUE for testing, and FALSE when doing all 50 sims
fileStr = '1_to_50_varTemp001'

# Constants
season = 6*24 
horizons = c(96, 192, 336, 720)
h_max = max(horizons)
train_len = 6*24*366 - h_max 

# Pull in all sims
file_path = paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/vpmax_",fileStr,".csv",sep="")
vpmax_sims <- read.csv(file_path, header = TRUE)
file_path = paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/tdew_",fileStr,".csv",sep="")
tdew_sims <- read.csv(file_path, header = TRUE)
file_path = paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/temp_",fileStr,".csv",sep="")
temp_sims <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Sims/Signal.csv"
sig <- read.csv(file_path, header = TRUE)
sig_tail_max = tail(sig$T..degC.,h_max)

# Initialize dataframes
error_df = data.frame(sim = 1:numSims)
for (iHor in 1:length(horizons)){
  h = horizons[iHor]
  colname = paste('mse_',h,sep="")
  error_df[[colname]] = rep(0,numSims)
}
for (iHor in 1:length(horizons)){
  h = horizons[iHor]
  colname = paste('mae_',h,sep="")
  error_df[[colname]] = rep(0,numSims)
}
for (iHor in 1:length(horizons)){
  h = horizons[iHor]
  colname = paste('mse_sig_',h,sep="")
  error_df[[colname]] = rep(0,numSims)
}
for (iHor in 1:length(horizons)){
  h = horizons[iHor]
  colname = paste('mae_sig_',h,sep="")
  error_df[[colname]] = rep(0,numSims)
}
fore_df = data.frame(matrix(0, nrow = h_max, ncol = numSims)) 
colnames(fore_df) <- paste0("sim", 1:numSims)

# Loop over sims
start_time <- Sys.time()
for (iSim in startSim:numSims){
  
  # Output time
  elapsed <- Sys.time() - start_time
  current_time <- format(Sys.time(), "%I:%M:%S %p")
  print(paste("Sim", iSim, "- Time elapsed:", round(elapsed, 2), "- Current Time:", current_time))
  
  # Get loop sims
  vpmax_sim = vpmax_sims[,iSim]
  tdew_sim = tdew_sims[,iSim]
  temp_sim = temp_sims[,iSim]
  
  # Get training data
  temp_trim <- temp_sim[1:train_len]
  vpmax_trim <- vpmax_sim[1:train_len]
  tdew_trim <- tdew_sim[1:train_len]
  
  # VPMax
  # vpmax_d = diff(vpmax_sim,season)
  vpmax_d = diff(vpmax_trim,season)
  invisible(capture.output({
    vpmax_est <- est.arma.wge(vpmax_d, p=9, q=5) # NOTE: In future, used trimmed data here
  }))  
  vpmax_phi <- vpmax_est$phi
  vpmax_theta <- vpmax_est$theta
  invisible(capture.output({
    # vpmax_fore_max <- fore.arima.wge(vpmax_sim, phi= vpmax_phi, theta = vpmax_theta, s=season, 
    #                                  lastn=TRUE, n.ahead = h_max, plot = FALSE)$f # NOTE: In future, used trimmed data here
    vpmax_fore_max <- fore.arima.wge(vpmax_trim, phi= vpmax_phi, theta = vpmax_theta, s=season, 
                                     lastn=FALSE, n.ahead = h_max, plot = FALSE)$f
  }))
  vpmax_tail_max <- tail(vpmax_sim, h_max)
  if (toPlot == TRUE){
    plot(1:h_max, vpmax_tail_max, type = "l", col = "black", lwd = 2,
         ylim = range(c(vpmax_tail_max, vpmax_fore_max)),
         xlab = "Time Step", ylab = "VPMax",
         main = paste("ARIMA(9,5 with s = 144): Last 720 Observed + Forecasted - VPMax - Sim",iSim))
    lines(1:h_max, vpmax_fore_max, col = "blue", lwd = 2)
    legend("topleft", legend = c("Observed", "Forecasted"),
           col = c("black", "blue"), lty = 1, lwd = 2)
  }
  
  # TDew
  # tdew_d = diff(tdew_sim,season)
  tdew_d = diff(tdew_trim,season)
  invisible(capture.output({
    tdew_est <- est.arma.wge(tdew_d, p=7, q=5) # NOTE: In future, used trimmed data here
  }))
  tdew_phi <- tdew_est$phi
  tdew_theta <- tdew_est$theta
  invisible(capture.output({
    # tdew_fore_max <- fore.arima.wge(tdew_sim, phi = tdew_phi, theta = tdew_theta, s=season, 
    #                                 lastn = TRUE, n.ahead = h_max, plot = FALSE)$f # NOTE: In future, used trimmed data here
    tdew_fore_max <- fore.arima.wge(tdew_trim, phi = tdew_phi, theta = tdew_theta, s=season, 
                                    lastn = FALSE, n.ahead = h_max, plot = FALSE)$f
  }))
  tdew_tail_max <- tail(tdew_sim, h_max)
  if (toPlot == TRUE){
    plot(1:h_max, tdew_tail_max, type = "l", col = "black", lwd = 2,
         ylim = range(c(tdew_tail_max, tdew_fore_max)),
         xlab = "Time Step", ylab = "TDew",
         main = paste("ARIMA(4,4 with s = 144): Last 720 Observed + Forecasted - TDew - Sim",iSim))
    lines(1:h_max, tdew_fore_max, col = "blue", lwd = 2)
    legend("topleft", legend = c("Observed", "Forecasted"),
           col = c("black", "blue"), lty = 1, lwd = 2)
  }
  
  # Temp
  fit = arima(temp_trim, order = c(4, 0, 1), xreg = cbind(vpmax_trim, tdew_trim)) 
  xreg_fore_max <- cbind(vpmax_fore_max, tdew_fore_max)
  temp_fore_max <- predict(fit, newxreg = xreg_fore_max)$pred
  temp_tail_max <- tail(temp_sim, h_max)
  if (toPlot == TRUE){
    plot(1:h_max, temp_tail_max, type = "l", col = "black", lwd = 2,
         ylim = range(c(temp_tail_max, temp_fore_max)),
         xlab = "Time Step", ylab = "Temperature (°C)",
         main = paste("MLR + ARMA(4,1) Forecast with Forecasted Exogenous Variables - Sim",iSim))
    lines(1:h_max, temp_fore_max, col = "blue", lwd = 2)
    legend("topleft", legend = c("Observed", "Forecasted"), col = c("black", "blue"), lty = 1, lwd = 2)
  }
  fore_df[,iSim] = temp_fore_max
  
  # Error metrics
  for (iHor in 1:length(horizons)){
    h = horizons[iHor]
    colname = paste('mse_',h,sep="")
    err = mean((temp_tail_max[1:h]-temp_fore_max[1:h])^2)
    error_df[iSim,colname] = err
    colname = paste('mae_',h,sep="")
    err = mean(abs(temp_tail_max[1:h]-temp_fore_max[1:h]))
    error_df[iSim,colname] = err
    colname = paste('mse_sig_',h,sep="")
    err = mean((sig_tail_max[1:h]-temp_fore_max[1:h])^2)
    error_df[iSim,colname] = err
    colname = paste('mae_sig_',h,sep="")
    err = mean(abs(sig_tail_max[1:h]-temp_fore_max[1:h]))
    error_df[iSim,colname] = err
  }
  
}

# Save off error metrics and forecasts
write_csv(error_df, paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_With_Trimmed_",fileStr,".csv",sep=""))
write_csv(fore_df, paste("C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Forecasts/ARIMA_Forecasts_With_Trimmed_",fileStr,".csv",sep=""))
