# (Optional) Clear memory space
# ls() # List items
# rm() # Remove items

# Load all libs -----------------------------------------------------------------------------------------

library(spectrolab)
library(dplyr)
library(tidyverse)
library(readxl)
library(prospect)
library(pracma) #function(fmincon) used for inversion
library(writexl)

spectra2text <- function(inputdir, metadata, sheet, outputdir) {
  setwd(inputdir) 
  files = Sys.glob('*.asd')
  specs = read_spectra(path=files, type = "target_reflectance", format="asd", extract_metadata=TRUE)
  
  splice_bands = c(1001, 1801) #defined by manufacturer
  matched = match_sensors(x = specs, splice_at = splice_bands,
                          interpolate_wvl = c(5, 1))
  
  metadata_table = read_excel(metadata, sheet = sheet)
  metadata_table$Leaf <- as.character(metadata_table$Leaf)
  
  summary = tibble(filename=files,
                   type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(metadata_table)),
                   Position = rep(metadata_table$Position, each=20),
                   tree.ID= rep(metadata_table$tree.ID, each=20),
                   Site = rep(metadata_table$Site, each=20),
                   Leaf = rep(metadata_table$Leaf, each=20),
                   Comment = rep(metadata_table$Comment, each=20))
  df <- cbind(summary, matched)
  CR = matrix(, nrow = nrow(df)/20, ncol = ncol(df)-8)
  #bis hier ohne jump correction 
  
  n = 1
  for (i in seq(1,nrow(df),20)){
    for (j in 9:ncol(df)){
      WR = mean(df[(i+1):(i+4),j])
      WRL = mean(df[(i+6):(i+9),j])
      BR = mean(df[(i+11):(i+14),j])
      BRL = mean(df[(i+16):(i+19),j])
      CR[n,j-8] = (BRL*WR - WRL*BR)/(WR-BR)
    }
    n = n+1
  }
  R.A = cbind(metadata_table,CR) #bind metadata to df 
  start_col <- 7 
  end_col <- 2157
  current_names <- colnames(R.A)
  new_names <- c(350:2500)
  colnames(R.A)[start_col:end_col] <- new_names #rename col names
  write.csv2(R.A, file = outputdir) #write csv
  print("The .csv has been exported in the defined output directory")
}
#Call Function

#arg1 = input directory (asd) (str)
#arg2 = input directory (metadata.xlsx) (str)
#arg3 = sheet of metadata.xlsx (str)
#arg4 = output directory (str)

spectra2text("/Users/Dave/Desktop/All_Data/SpectroscopyData/Curated/AfterDroughtData/20230906_B_C_D_I_J_K/Block K",
             "/Users/Dave/Desktop/All_Data/SpectroscopyData/Curated/AfterDroughtData/All_Metadata/Metadata.xlsx",
             "Block K",
             "/Users/Dave/Desktop/All_Data/SpectroscopyData/CuratedSpectraCSV/AfterDrought/Block_B_C_D_I_J_K/Block_K.csv")

# Calculation of measurement uncertainties (Absolute and Relative) (Cheng et al. 2023) (x) -----------------------------

setwd("/Users/Dave/Desktop/All_Data/SpectroscopyData/Curated/AfterDroughtData/20230906_B_C_D_I_J_K/Block J")
files = Sys.glob('*.asd')
raw_spectra = read_spectra(path=files, type = "target_reflectance", format="asd", extract_metadata=TRUE)
splice_bands = c(1001, 1801)
raw_spectra_matched = match_sensors(x = raw_spectra, splice_at = splice_bands,
                                    interpolate_wvl = c(5, 1))
metadata_table = read_excel("/Users/Dave/Desktop/All_Data/SpectroscopyData/Curated/PreDroughtData/All_Metadata/Metadata.xlsx")
metadata_table$Leaf <- as.character(metadata_table$Leaf)

summary = tibble(filename=files,
                 type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(metadata_table)),
                 Position = rep(metadata_table$Position, each=20),
                 tree.ID= rep(metadata_table$tree.ID, each=20),
                 Site = rep(metadata_table$Site, each=20),
                 Leaf = rep(metadata_table$Leaf, each=20),
                 Comment = rep(metadata_table$Comment, each=20))

df <- cbind(summary, raw_spectra_matched)
CR = matrix(, nrow = nrow(df)/20, ncol = ncol(df)-8)
AU = matrix(, nrow = nrow(df)/20, ncol = ncol(df)-8)
RU = matrix(, nrow = nrow(df)/20, ncol = ncol(df)-8)
n = 1
for (i in seq(1,nrow(df),20)){
  for (j in 9:ncol(df)){
    WR = mean(df[(i+1):(i+4),j])
    WRL = mean(df[(i+6):(i+9),j])
    BR = mean(df[(i+11):(i+14),j])
    BRL = mean(df[(i+16):(i+19),j])
    STD_WR = sd(df[(i+1):(i+4),j])
    STD_WRL = sd(df[(i+6):(i+9),j])
    STD_BR = sd(df[(i+11):(i+14),j])
    STD_BRL = sd(df[(i+16):(i+19),j])
    CR[n,j-8] = (BRL*WR - WRL*BR)/(WR-BR)
    AU[n,j-8] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/sqrt(5))^2 + (BR/(WR-BR))^2 * (STD_WRL/sqrt(5))^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/sqrt(5))^2 + (WR/(WR-BR))^2 * (STD_BRL/sqrt(5))^2)
    RU[n,j-8] = ((sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/sqrt(5))^2 + (BR/(WR-BR))^2 * (STD_WRL/sqrt(5))^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/sqrt(5))^2 + (WR/(WR-BR))^2 * (STD_BRL/sqrt(5))^2)) / ((BRL*WR - WRL*BR)/(WR-BR)))*100
  }
  n = n+1
}

result_AU = cbind(metadata_table, AU)
result_RU = cbind(metadata_table, RU)
start_col <- 7 
end_col <- 2157
current_names <- colnames(result_AU)
current_names <- colnames(result_RU)
new_names <- c(350:2500)
colnames(result_AU)[start_col:end_col] <- new_names
colnames(result_RU)[start_col:end_col] <- new_names
write.csv2(result_AU, file = "/Users/Dave/Desktop/All_Data/Statistical_Analysis_All/Uncertainties/Temp_AU.csv")
write.csv2(result_RU, file = "/Users/Dave/Desktop/All_Data/Statistical_Analysis_All/Uncertainties/Temp_RU.csv" )
print("The .csv has been exported in the defined output directory")

# Prospect Forwards (Tryout, can be ignored, only for visualisation) (x) --------------------------------------------------
# Title: PROSPECT modelling on leaf spectroscopy data during common garden experiment 2023
# Description: This script serves to model and analyse the spectral data from ASD FieldSpec 4 during the common garden experiment 2023
# Source Modelling: Feret J, de Boissieu F (2023). prospect: PROSPECT leaf radiative transfer model and inversion routines. 
# R package version 1.3.0, https://gitlab.com/jbferet/prospect.
# Author: Dave Kurath
# Data of Creation: 13.09.2023
# Data of last modification: 27.09.2023, Dave Kurath

# Function PROSPECT runs model for individual samples

# Variables:
# SpecPROSPECT: DF including refractive index and specific absorption coefficients for a given spectral range (default: 400nm - 2500nm)

# Important: Calling prospect-d requires LMA to be set with >0.00 as it is composed of CBC + PROT, running prospect-pro
# requires LMA to be set to 0.00

# Biochemical and biophysical input variables (model parameters):
# N (default = 1.5) [1-4]
# CHL Cab (default = 40.0 𝜇𝑔.𝑐𝑚2) [0-100] 
# CAR (default = 8.0 𝜇𝑔.𝑐𝑚2) [0-40]
# ANT (default = 0.0 𝜇𝑔.𝑐𝑚2) [0-40] 
# BROWN (default = 0.0 arbitrary units) 
# EWT (default = 0.01 𝑔.𝑐𝑚2) [0.05]
# LMA, Dry Matter (default = 0.008 𝑔.𝑐𝑚2) [0-0.05] Only in Prospect-D
# PROT (default = 0.0 𝑔.𝑐𝑚2)  [0-0.03] Only in Prospect-Pro
# CBC (default = 0.0 𝑔.𝑐𝑚2) [0-0.02] Only in Prospect-Pro
# alpha (default = 40.0 degrees)



# Running PROSPECT in default over default range:
# LRT_default <- PROSPECT(SpecPROSPECT)

# Running PROSPECT with defined set of of parameter over default range:

LRT_VSWIR <- PROSPECT(SpecPROSPECT,N = 1.4,CHL = 30,CAR = 6,EWT = 0.02,LMA = 0.01)
x_values <- LRT_VSWIR$wvl
y_values <- LRT_VSWIR$Reflectance
plot(x_values, y_values, type = "b", pch = 10, col = "black", cex = 0.1,
     xlab = "Wavelength (nm)", ylab = "Reflectance",
     main = "Tryout Forward Modelling Whole Spectral Domain")
lines(x_values, y_values, col = "black")

# Running PROSPECT in defined spectral domain (e.g. 400nm-1000nm)

wvlRange = list()
wvlRange$lb <- which(abs(SpecPROSPECT$lambda-400)==min(abs(SpecPROSPECT$lambda-400)))
wvlRange$ub <- which(abs(SpecPROSPECT$lambda-1000)==min(abs(SpecPROSPECT$lambda-1000)))
SpecPROSPECT_VNIR <- SpecPROSPECT[wvlRange$lb:wvlRange$ub,]
LRT_VNIR <- PROSPECT(SpecPROSPECT_VNIR,N = 1.4,CHL = 30,CAR = 6,EWT = 0.02,LMA = 0.01)
x_values <- LRT_VNIR$wvl
y_values <- LRT_VNIR$Reflectance
plot(x_values, y_values, type = "b", pch = 10, col = "black", cex = 0.1,
     xlab = "Wavelength (nm)", ylab = "Reflectance",
     main = "Tryout Forward Modelling VIS-NIR Domain")
lines(x_values, y_values, col = "black")



# Prospect Inversed (Validation, All Data) (x, done for both optical and whole) --------------------------------

# Variables:
# SpecPROSPECT: DF including refractive index and specific absorption coefficients for a given spectral range (default: 400nm - 2500nm)
# Refl: numeric: individual leaf reflectance corresponding to the spectral domain defined in SpecPROSPECT. Set to NULL if inversion on transmittance only
# Tran: numeric: individual leaf transmittance corresponding to the spectral domain defined in SpecPROSPECT. Set to NULL if inversion on reflectance only
# --> Since the data is reflectance. Set Tran: NULL
# Parms2Estimate list. Parameters to estimate. Set to ‘ALL’ by default. (Keep in mind which Prospect Model to use)
# InitValues data frame including initial values of PROSPECT input parameters. During optimization, they are used either as initialization values for parameters to estimate,
# or as fix values for other parameters. Parameters not compatible with PROSPECT_version are not taken into account.
# PROSPECT_version character. corresponds to the PROSPECT version. Following versions: ‘5’, ‘5B’, ‘D’, ‘DB’, ‘PRO’, ‘PROB’
# --> either D or PRO
# MeritFunction character. name of the function to be used as merit function with given criterion to minimize (default = RMSE)
# --> One could opt for different optimization strategies!
# xlub data.frame. Boundaries of the parameters to estimate. The data.frame must have columns corresponding to first line being the lower boundaries and second
# line the upper boundaries.
# alphaEst boolean.
# --> most results with alpha as default

# Output: List containing estimated values of Input variables

# Run Prospect-D inversion over full domain (Testsample)

# Simulate leaf optical properties
CHL <- 45;      CAR <- 10;      ANT <- 0.2
EWT <- 0.012;   LMA <- 0.010;   N   <- 1.3
LRT_D <- PROSPECT(SpecPROSPECT,CHL=CHL,CAR=CAR,ANT=ANT,EWT=EWT,LMA=LMA,N=N)

# Define set of parameters to be estimated
Parms2Estimate  <- 'ALL'
# Define initial values for the inversion (should not impact final results)
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
# Invert PROSPECT with simulated leaf optical properties
OutPROSPECT_Test <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                    Refl = LRT_D$Reflectance, Tran = LRT_D$Transmittance,
                                    Parms2Estimate = 'ALL', PROSPECT_version = 'D')

# Run Prospect-D inversion over full domain (Tryout for one validation sample only)
# Import Data
spectrasample <- read_excel("/Users/Dave/Desktop/All_Data/Prospect Modelling/Pruned For Inverse Model/OneSample.xlsx")
spectra_list <- list(wvl = (spectrasample$wvl), reflectance = matrix(as.numeric(spectrasample$reflectance), ncol =1))

Parms2Estimate  <- 'ALL'
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
# invert PROSPECT with simulated leaf optical properties
OutPROSPECT_OneSample <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                         Refl = spectra_list$reflectance, Tran = NULL,
                                         Parms2Estimate = 'ALL', PROSPECT_version = 'D')
print(OutPROSPECT_OneSample)

# Inverse Mode: Iterative optimization over full domain no prior estimation of N

# Load Data
spectrasample <- read_excel("/Users/Dave/Desktop/All_Data/Prospect Modelling/Pruned For Inverse Model/Validation/Validation_Transformed.xlsx")

#Create For-Loop later to iterate through all the data (TBC)
spectra_list_wvl <- list(wvl = (spectrasample$wvl))
spectra_list <- list()
spectra_list <- append(spectra_list, list(reflectanceW1 = matrix(as.numeric(spectrasample$W1), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW2 = matrix(as.numeric(spectrasample$W2), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW3 = matrix(as.numeric(spectrasample$W3), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW4 = matrix(as.numeric(spectrasample$W4), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW5 = matrix(as.numeric(spectrasample$W5), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW6 = matrix(as.numeric(spectrasample$W6), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW7 = matrix(as.numeric(spectrasample$W7), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW8 = matrix(as.numeric(spectrasample$W8), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW9 = matrix(as.numeric(spectrasample$W9), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW10 = matrix(as.numeric(spectrasample$W10), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW11 = matrix(as.numeric(spectrasample$W11), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW12 = matrix(as.numeric(spectrasample$W12), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW13 = matrix(as.numeric(spectrasample$W13), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW14 = matrix(as.numeric(spectrasample$W14), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW15 = matrix(as.numeric(spectrasample$W15), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW16 = matrix(as.numeric(spectrasample$W16), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW17 = matrix(as.numeric(spectrasample$W17), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW18 = matrix(as.numeric(spectrasample$W18), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW19 = matrix(as.numeric(spectrasample$W19), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW20 = matrix(as.numeric(spectrasample$W20), ncol =1)))

# PROSPECT-D full spectra domain
Parms2Estimate  <- 'ALL'
InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
df_D <- data_frame()
for (i in spectra_list){
  OutPROSPECT <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                 Refl = i, Tran = NULL,
                                 Parms2Estimate = 'ALL', PROSPECT_version = 'D')
  df_D <- rbind(df_D, OutPROSPECT)
  print(OutPROSPECT)
}

# Inversion of PROSPECT-D using the optimal configuration for the estimation of leaf constituents (Féret et al., 2019 & 2020)
# Needs specification (!)

Parms2Estimate  = c('CHL','CAR','ANT','EWT','LMA')
InitValues <- data.frame(CHL=40, CAR=8, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=1.5)
# call Invert_PROSPECT_OPT
df_D_opt <- data_frame()

for (i in spectra_list){
  OutPROSPECT_Opt <- Invert_PROSPECT_OPT(SpecPROSPECT=SpecPROSPECT, lambda=spectra_list_wvl$wvl, 
                                         Refl = i, Tran = NULL,
                                         PROSPECT_version = 'D',
                                         Parms2Estimate = Parms2Estimate,
                                         InitValues = InitValues)
  df_D_opt <- rbind(df_D_opt, OutPROSPECT_Opt)
  print(OutPROSPECT_Opt)
}



# Inverse Mode: Testing Estimation with R only and prior estimation of N (Validation HPLC Data) 

# Prior N-Estimation based on Spafford et al. 2021
# Best results based on Spafford et al. 2021 decision tree:
# Full domain reflectance with prior N estimation: CAR, EWT
# Optimal subdomain reflectance with prior N estimation: All, Cab, LMA
# ANT is not considered


#Load Data
spectrasample <- read_excel("/Users/Dave/Desktop/All_Data/Prospect Modelling/Pruned For Inverse Model/Validation/Validation_Transformed.xlsx")

spectra_list_wvl <- list(wvl = (spectrasample$wvl))
spectra_list <- list()
spectra_list <- append(spectra_list, list(reflectanceW1 = matrix(as.numeric(spectrasample$W1), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW2 = matrix(as.numeric(spectrasample$W2), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW3 = matrix(as.numeric(spectrasample$W3), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW4 = matrix(as.numeric(spectrasample$W4), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW5 = matrix(as.numeric(spectrasample$W5), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW6 = matrix(as.numeric(spectrasample$W6), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW7 = matrix(as.numeric(spectrasample$W7), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW8 = matrix(as.numeric(spectrasample$W8), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW9 = matrix(as.numeric(spectrasample$W9), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW10 = matrix(as.numeric(spectrasample$W10), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW11 = matrix(as.numeric(spectrasample$W11), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW12 = matrix(as.numeric(spectrasample$W12), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW13 = matrix(as.numeric(spectrasample$W13), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW14 = matrix(as.numeric(spectrasample$W14), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW15 = matrix(as.numeric(spectrasample$W15), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW16 = matrix(as.numeric(spectrasample$W16), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW17 = matrix(as.numeric(spectrasample$W17), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW18 = matrix(as.numeric(spectrasample$W18), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW19 = matrix(as.numeric(spectrasample$W19), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceW20 = matrix(as.numeric(spectrasample$W20), ncol =1)))
#Get N prior estimation based on Qiu et al. (2018)

lambda <- spectra_list_wvl$wvl
Nprior_lys = list()
for (i in spectra_list){
  Nprior_R = Get_Nprior(SpecPROSPECT,lambda,Refl=i)
  
  Nprior_lys <- append(Nprior_lys, Nprior_R)
  print(Nprior_R)
}

# PROSPECT-D full spectra domain but with prior N estimation
Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_D_nprior <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                 Refl = i, Tran = NULL,
                                 Parms2Estimate = Parms2Estimate, PROSPECT_version = 'D', InitValues = InitValues)
  df_D_nprior <- rbind (df_D_nprior, OutPROSPECT)
  print(OutPROSPECT)
  
}

# PROSPECT-D, prior estimation of N and optimal spectral domain for each constituent

Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_D_nprior_opt <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=8, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, alpha=40, PROT=0, CBC=0, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT_OPT <- Invert_PROSPECT_OPT(SpecPROSPECT=SpecPROSPECT, lambda =spectra_list_wvl$wvl,
                                         Refl = i, Tran = NULL,
                                         Parms2Estimate = Parms2Estimate , PROSPECT_version = 'D', InitValues = InitValues)
  df_D_nprior_opt <- rbind (df_D_nprior_opt, OutPROSPECT_OPT)
  print(OutPROSPECT_OPT)
  
}
# Inverse Mode: Testing Estimation with R only and prior estimation of N (Validation Whitebag)
spectrasample <- read_excel("/Users/Dave/Desktop/All_Data/Prospect Modelling/Pruned For Inverse Model/Validation_LMA_EWT/Validation_Transformed_Whitebag.xlsx")

spectra_list_wvl <- list(wvl = (spectrasample$wvl))
spectra_list <- list()
spectra_list <- append(spectra_list, list(reflectanceA1 = matrix(as.numeric(spectrasample$A1), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA2 = matrix(as.numeric(spectrasample$A2), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA3 = matrix(as.numeric(spectrasample$A3), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA4 = matrix(as.numeric(spectrasample$A4), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA5 = matrix(as.numeric(spectrasample$A5), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA6 = matrix(as.numeric(spectrasample$A6), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA7 = matrix(as.numeric(spectrasample$A7), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA8 = matrix(as.numeric(spectrasample$A8), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA9 = matrix(as.numeric(spectrasample$A9), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA10 = matrix(as.numeric(spectrasample$A10), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA11 = matrix(as.numeric(spectrasample$A11), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA12 = matrix(as.numeric(spectrasample$A12), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA13 = matrix(as.numeric(spectrasample$A13), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA14 = matrix(as.numeric(spectrasample$A14), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA15 = matrix(as.numeric(spectrasample$A15), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA16 = matrix(as.numeric(spectrasample$A16), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA17 = matrix(as.numeric(spectrasample$A17), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceA18 = matrix(as.numeric(spectrasample$A18), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB1 = matrix(as.numeric(spectrasample$B1), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB2 = matrix(as.numeric(spectrasample$B2), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB3 = matrix(as.numeric(spectrasample$B3), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB4 = matrix(as.numeric(spectrasample$B4), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB5 = matrix(as.numeric(spectrasample$B5), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB6 = matrix(as.numeric(spectrasample$B6), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB7 = matrix(as.numeric(spectrasample$B7), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB8 = matrix(as.numeric(spectrasample$B8), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB9 = matrix(as.numeric(spectrasample$B9), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB10 = matrix(as.numeric(spectrasample$B10), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB11 = matrix(as.numeric(spectrasample$B11), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB12 = matrix(as.numeric(spectrasample$B12), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB13 = matrix(as.numeric(spectrasample$B13), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB14 = matrix(as.numeric(spectrasample$B14), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB15 = matrix(as.numeric(spectrasample$B15), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB16 = matrix(as.numeric(spectrasample$B16), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB17 = matrix(as.numeric(spectrasample$B17), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceB18 = matrix(as.numeric(spectrasample$B18), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC1 = matrix(as.numeric(spectrasample$C1), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC2 = matrix(as.numeric(spectrasample$C2), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC3 = matrix(as.numeric(spectrasample$C3), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC4 = matrix(as.numeric(spectrasample$C4), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC5 = matrix(as.numeric(spectrasample$C5), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC6 = matrix(as.numeric(spectrasample$C6), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC7 = matrix(as.numeric(spectrasample$C7), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC8 = matrix(as.numeric(spectrasample$C8), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC9 = matrix(as.numeric(spectrasample$C9), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC10 = matrix(as.numeric(spectrasample$C10), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC11 = matrix(as.numeric(spectrasample$C11), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC12 = matrix(as.numeric(spectrasample$C12), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC13 = matrix(as.numeric(spectrasample$C13), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC14 = matrix(as.numeric(spectrasample$C14), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC16 = matrix(as.numeric(spectrasample$C16), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC17 = matrix(as.numeric(spectrasample$C17), ncol =1)))
spectra_list <- append(spectra_list, list(reflectanceC18 = matrix(as.numeric(spectrasample$C18), ncol =1)))

#Get N prior estimation based on Qiu et al. (2018)
lambda <- spectra_list_wvl$wvl
Nprior_lys = list()
for (i in spectra_list){
  Nprior_R = Get_Nprior(SpecPROSPECT,lambda,Refl=i)
  
  Nprior_lys <- append(Nprior_lys, Nprior_R)
  print(Nprior_R)
}

# PROSPECT-D full spectra domain but with prior N estimation
Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_D_nprior <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                 Refl = i, Tran = NULL,
                                 Parms2Estimate = Parms2Estimate, PROSPECT_version = 'D', InitValues = InitValues)
  df_D_nprior <- rbind (df_D_nprior, OutPROSPECT)
  print(OutPROSPECT)
  
}

# PROSPECT-D, prior estimation of N and optimal spectral domain for each constituent

Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_D_nprior_opt <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=8, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, alpha=40, PROT=0, CBC=0, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT_OPT <- Invert_PROSPECT_OPT(SpecPROSPECT=SpecPROSPECT, lambda =spectra_list_wvl$wvl,
                                         Refl = i, Tran = NULL,
                                         Parms2Estimate = Parms2Estimate , PROSPECT_version = 'D', InitValues = InitValues)
  df_D_nprior_opt <- rbind (df_D_nprior_opt, OutPROSPECT_OPT)
  print(OutPROSPECT_OPT)
  
}

# PROSPECT-D, all samples, prior N estimation, whole spectra and optical subdomains (FOR ALL DATA - IMPORTANT)
spectrasample <- read_excel("/Users/Dave/Desktop/All_Data/Prospect Modelling/Pruned For Inverse Model/All_Transformed.xlsx")

spectra_list_wvl <- list(wvl = (spectrasample$wvl))
spectra_list <- list()

for (i in 1:2634) {
  col_name <- paste0("reflectanceW", i)
  spectra_list[[col_name]] <- matrix(as.numeric(spectrasample[[paste0("W", i)]]), ncol = 1)
  print(i)
}

lambda <- spectra_list_wvl$wvl
Nprior_lys = list()
c <- 0
for (i in spectra_list){
  c <- c + 1
  Nprior_R = Get_Nprior(SpecPROSPECT,lambda,Refl=i)
  
  Nprior_lys <- append(Nprior_lys, Nprior_R)
  print(Nprior_R)
  print(c)
}

df_N <- data_frame(Nprior_lys)
hist(as.numeric(df_N$Nprior_lys), breaks = 30)

# PROSPECT-D inverse full spectra domain
Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_all_full <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=10, ANT=0.1, BROWN=0, EWT=0.01, LMA=0.01, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT <- Invert_PROSPECT(SpecPROSPECT=SpecPROSPECT,
                                 Refl = i, Tran = NULL,
                                 Parms2Estimate = Parms2Estimate, PROSPECT_version = 'D', InitValues = InitValues)
  df_all_full <- rbind (df_all_full, OutPROSPECT)
  print(OutPROSPECT)
  
}

# PROSPECT-D inverse optical domain

Parms2Estimate  <- c("CHL", "CAR", "ANT", "EWT", "LMA")
df_all_opt <- data_frame()
c <- 0
for (i in spectra_list){
  c <- c + 1
  print (c)
  InitValues <- data.frame(CHL=40, CAR=8, ANT=0, BROWN=0, EWT=0.01, LMA=0.01, alpha=40, PROT=0, CBC=0, N=Nprior_lys[c])
  print(Nprior_lys[i])
  OutPROSPECT_OPT <- Invert_PROSPECT_OPT(SpecPROSPECT=SpecPROSPECT, lambda =spectra_list_wvl$wvl,
                                         Refl = i, Tran = NULL,
                                         Parms2Estimate = Parms2Estimate , PROSPECT_version = 'D', InitValues = InitValues)
  df_all_opt <- rbind (df_all_opt, OutPROSPECT_OPT)
  print(OutPROSPECT_OPT)
  
}

write_xlsx(df_all_full,"/Users/Dave/Desktop/All_Data/Statistical_Analysis_All/all_constituents_full.xlsx")
write_xlsx(df_all_opt,"/Users/Dave/Desktop/All_Data/Statistical_Analysis_All/all_constituents_opt.xlsx")
write_xlsx(data.frame(Nprior_lys),"/Users/Dave/Desktop/All_Data/Statistical_Analysis_All/all_n_prior.xlsx")


