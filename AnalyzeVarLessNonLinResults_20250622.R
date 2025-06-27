################################################################################
# Look at less variance results

# Pull in ARIMA results
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_1_to_50_varTemp001varFac05.csv"
arima_results_lessVar <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_Errors_With_Trimmed_1_to_50_varTemp001.csv"
arima_results_orig <- read.csv(file_path, header = TRUE)

# Look at means
mean(arima_results_orig$mse_96) # 3.565827
mean(arima_results_orig$mse_192) # 4.278645
mean(arima_results_orig$mse_336) # 5.572522
mean(arima_results_orig$mse_720) # 8.979426
mean(arima_results_lessVar$mse_96) # 1.855132
mean(arima_results_lessVar$mse_192) # 2.238425
mean(arima_results_lessVar$mse_336) # 2.919909
mean(arima_results_lessVar$mse_720) # 4.72879
mean(arima_results_orig$mse_96) / mean(arima_results_lessVar$mse_96) # 1.922142
mean(arima_results_orig$mse_192) / mean(arima_results_lessVar$mse_192) # 1.911453
mean(arima_results_orig$mse_336) / mean(arima_results_lessVar$mse_336) # 1.908458
mean(arima_results_orig$mse_720) / mean(arima_results_lessVar$mse_720) # 1.898885

# Wilson tests
wilcox.test(arima_results_orig$mse_96, arima_results_lessVar$mse_96, paired=TRUE) # p-value = 4.919e-09
wilcox.test(arima_results_orig$mse_192, arima_results_lessVar$mse_192, paired=TRUE) # p-value = 1.513e-09
wilcox.test(arima_results_orig$mse_336, arima_results_lessVar$mse_336, paired=TRUE) # p-value = 7.79e-10
wilcox.test(arima_results_orig$mse_720, arima_results_lessVar$mse_720, paired=TRUE) # p-value = 7.79e-10

# Pull in TFT results
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/TFT_Errors_1_to_50_varTemp001varFac05.csv"
tft_results_lessVar <- read.csv(file_path, header = TRUE)
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/TFT_Errors_1_to_50_varTemp001.csv"
tft_results_orig <- read.csv(file_path, header = TRUE)

# Look at means
mean(tft_results_orig$MSE.96) # 15.95684
mean(tft_results_orig$MSE.192) # 17.61195
mean(tft_results_orig$MSE.336) # 19.81339
mean(tft_results_orig$MSE.720) # 22.94062
mean(tft_results_lessVar$MSE.96) # 8.699824
mean(tft_results_lessVar$MSE.192) # 8.822522
mean(tft_results_lessVar$MSE.336) # 9.186374
mean(tft_results_lessVar$MSE.720) # 11.01318
mean(tft_results_orig$MSE.96) / mean(tft_results_lessVar$MSE.96) # 1.834157
mean(tft_results_orig$MSE.192) / mean(tft_results_lessVar$MSE.192) # 1.996249
mean(tft_results_orig$MSE.336) / mean(tft_results_lessVar$MSE.336) # 2.156824
mean(tft_results_orig$MSE.720) / mean(tft_results_lessVar$MSE.720) # 2.083015

# Wilson tests
wilcox.test(tft_results_orig$MSE.96, tft_results_lessVar$MSE.96, paired=TRUE) # p-value = 0.001705
wilcox.test(tft_results_orig$MSE.192, tft_results_lessVar$MSE.192, paired=TRUE) # p-value = 7.263e-05
wilcox.test(tft_results_orig$MSE.336, tft_results_lessVar$MSE.336, paired=TRUE) # p-value = 1.339e-05
wilcox.test(tft_results_orig$MSE.720, tft_results_lessVar$MSE.720, paired=TRUE) # p-value = 1.767e-06

# Hm, it's interesting how different the p-values are, even with the ratios in means is similar
# Look at number of times one value is less than another
sum(arima_results_orig$mse_96 < arima_results_lessVar$mse_96) # 2
sum(arima_results_orig$mse_192 < arima_results_lessVar$mse_192) # 1
sum(arima_results_orig$mse_336 < arima_results_lessVar$mse_336) # 0
sum(arima_results_orig$mse_720 < arima_results_lessVar$mse_720) # 0
sum(tft_results_orig$MSE.96 < tft_results_lessVar$MSE.96) # 14
sum(tft_results_orig$MSE.192 < tft_results_lessVar$MSE.192) # 10
sum(tft_results_orig$MSE.336 < tft_results_lessVar$MSE.336) # 10
sum(tft_results_orig$MSE.720 < tft_results_lessVar$MSE.720) # 7

# Wilson tests for ARIMA vs TFT
wilcox.test(arima_results_lessVar$mse_96, tft_results_lessVar$MSE.96, paired=TRUE) # p-value = 1.983e-05
wilcox.test(arima_results_lessVar$mse_192, tft_results_lessVar$MSE.192, paired=TRUE) # p-value = 5.687e-05
wilcox.test(arima_results_lessVar$mse_336, tft_results_lessVar$MSE.336, paired=TRUE) # p-value = 2.567e-05
wilcox.test(arima_results_lessVar$mse_720, tft_results_lessVar$MSE.720, paired=TRUE) # p-value = 0.0002442

################################################################################
# Look at non-linear simulation results

# Pull in ARIMA results
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/ARIMA_NonLin_Errors_1_to_50_varTemp001.csv"
arima_results_nonlin <- read.csv(file_path, header = TRUE)

# Look at means
mean(arima_results_orig$mse_96) # 3.565827
mean(arima_results_orig$mse_192) # 4.278645
mean(arima_results_orig$mse_336) # 5.572522
mean(arima_results_orig$mse_720) # 8.979426
mean(arima_results_nonlin$mse_96) # 7.463704
mean(arima_results_nonlin$mse_192) # 10.16487
mean(arima_results_nonlin$mse_336) # 14.07899
mean(arima_results_nonlin$mse_720) # 21.82759
mean(arima_results_nonlin$mse_96) / mean(arima_results_orig$mse_96) # 2.09312
mean(arima_results_nonlin$mse_192) / mean(arima_results_orig$mse_192) # 2.375721
mean(arima_results_nonlin$mse_336) / mean(arima_results_orig$mse_336) # 2.526502
mean(arima_results_nonlin$mse_720) / mean(arima_results_orig$mse_720) # 2.430845

# Number of results that are less
sum(arima_results_nonlin$mse_96 < arima_results_orig$mse_96) # 11
sum(arima_results_nonlin$mse_192 < arima_results_orig$mse_192) # 7
sum(arima_results_nonlin$mse_336 < arima_results_orig$mse_336) # 5
sum(arima_results_nonlin$mse_720 < arima_results_orig$mse_720) # 5

# Pull in TFT results
file_path = "C:/Users/aabro/OneDrive/Desktop/SMU Program/Capstone/Results/TFT_NonLin_Errors_1_to_50_varTemp001.csv"
tft_results_nonlin <- read.csv(file_path, header = TRUE)

# Look at means
mean(tft_results_orig$MSE.96) # 15.95684
mean(tft_results_orig$MSE.192) # 17.61195
mean(tft_results_orig$MSE.336) # 19.81339
mean(tft_results_orig$MSE.720) # 22.94062
mean(tft_results_nonlin$MSE.96) # 29.40285
mean(tft_results_nonlin$MSE.192) # 31.74817
mean(tft_results_nonlin$MSE.336) # 36.26499
mean(tft_results_nonlin$MSE.720) # 43.94629
mean(tft_results_nonlin$MSE.96) / mean(tft_results_orig$MSE.96) # 1.842648
mean(tft_results_nonlin$MSE.192) / mean(tft_results_orig$MSE.192) # 1.802649
mean(tft_results_nonlin$MSE.336) / mean(tft_results_orig$MSE.336) # 1.830327
mean(tft_results_nonlin$MSE.720) / mean(tft_results_orig$MSE.720) # 1.915654

# Wilson tests
wilcox.test(tft_results_orig$MSE.96, tft_results_nonlin$MSE.96, paired=TRUE) # p-value = 0.01949
wilcox.test(tft_results_orig$MSE.192, tft_results_nonlin$MSE.192, paired=TRUE) # p-value = 0.006115
wilcox.test(tft_results_orig$MSE.336, tft_results_nonlin$MSE.336, paired=TRUE) # p-value = 0.001103
wilcox.test(tft_results_orig$MSE.720, tft_results_nonlin$MSE.720, paired=TRUE) # p-value = 0.0001733

# How many are less
sum(tft_results_nonlin$MSE.96 < tft_results_orig$MSE.96) # 20
sum(tft_results_nonlin$MSE.192 < tft_results_orig$MSE.192) # 15
sum(tft_results_nonlin$MSE.336 < tft_results_orig$MSE.336) # 13
sum(tft_results_nonlin$MSE.720 < tft_results_orig$MSE.720) # 9

# Wilson tests for ARIMA vs TFT
wilcox.test(arima_results_nonlin$mse_96, tft_results_nonlin$MSE.96, paired=TRUE) # p-value = 4.628e-05
wilcox.test(arima_results_nonlin$mse_192, tft_results_nonlin$MSE.192, paired=TRUE) # p-value = 2.843e-06
wilcox.test(arima_results_nonlin$mse_336, tft_results_nonlin$MSE.336, paired=TRUE) # p-value = 3.454e-05
wilcox.test(arima_results_nonlin$mse_720, tft_results_nonlin$MSE.720, paired=TRUE) # p-value = 0.0001269
