################################################################################################################################################################
########################   Analysis Code  ######################################################################################################################
################################################################################################################################################################

# Title: Data script for "Leaf spectroscopy reveals drought response variation in Fagus sylvatica saplings from across the species' range"
# Description: This script has the codes written for the manuscript
# Author: xx, xx
# Data of Creation: 16.10.2023
# Data of last modification: 10.10.2024, xx

# (Optional) Clear memory space ----
# ls() # List items
# rm() # Remove items

# load all libs ----

library(spectrolab)
library(dplyr)
library(tidyverse)
library(readxl)
library(prospect)
library(pracma) #function(fmincon) used for inversion
library(writexl)
library(openxlsx)
library(ggplot2)
library(showtext)
library(lme4)
library(stats)
library(ggpubr)
library(broom)
library(plotrix)
library(lemon)
library(gridExtra)
library(lme4)
library(ggsci)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(car)
library(afex)
library(vioplot)
library(agricolae)

# conversion .asd to .txt (Spectrolab) ----

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

spectra2text("arg1,arg2,arg3,arg4")



# PROSPECT inversion ----
# Title: PROSPECT modelling on leaf spectroscopy data during common garden experiment 2023
# Description: This script serves to model and analyse the spectral data from ASD FieldSpec 4 during the common garden experiment 2023
# Source Modelling: Feret J, de Boissieu F (2023). prospect: PROSPECT leaf radiative transfer model and inversion routines. 
# R package version 1.3.0, https://gitlab.com/jbferet/prospect.
# Author: xx
# Data of Creation: 13.09.2023
# Data of last modification: 27.09.2023, xx

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

# SpecPROSPECT: DF including refractive index and specific absorption coefficients for a given spectral range (default: 400nm - 2500nm)
# Refl: numeric: individual leaf reflectance corresponding to the spectral domain defined in SpecPROSPECT. Set to NULL if inversion on transmittance only
# Tran: numeric: individual leaf transmittance corresponding to the spectral domain defined in SpecPROSPECT. Set to NULL if inversion on reflectance only
# data: reflectance. Set Tran: NULL
# MeritFunction character. name of the function to be used as merit function with given criterion to minimize (default = RMSE)
# xlub data.frame. Boundaries of the parameters to estimate. The data.frame must have columns corresponding to first line being the lower boundaries and second
# line the upper boundaries.
# alphaEst boolean.


# PROSPECT-D, all samples, prior N estimation, whole spectral domain (FOR ALL DATA)
spectrasample <- read_excel("")

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
write_xlsx(df_all_full,"")



# statistical analysis ----

df <- read.xlsx("") 

df$Length <- as.numeric(df$Length)
df$Wilting <- as.numeric(df$Wilting)
df$RGR <- as.numeric(df$RGR)
df$ARI <- as.numeric(df$ARI)
df$CARI <- as.numeric(df$CARI)
df$NDSR <- as.numeric(df$NDSR)
df$CAI <- as.numeric(df$CAI)
df$CCI <- as.numeric(df$CCI)
df$CIRE <- as.numeric(df$CIRE)
df$NDLMA <- as.numeric(df$NDLMA)
df$NDNI <- as.numeric(df$NDNI)
df$ChlNDI <- as.numeric(df$ChlNDI)
df$NDLI <- as.numeric(df$NDLI)
df$NDWI <- as.numeric(df$NDWI)
df$PRI <- as.numeric(df$PRI)
df$EWI <- as.numeric(df$EWI)
df$MSI <- as.numeric(df$MSI)
df$WABI <- as.numeric(df$WABI)
df$WI <- as.numeric(df$WI)
df$N_whole <- as.numeric(df$N_whole)
df$Ccab_whole <- as.numeric(df$Ccab_whole)
df$Ccar_whole <- as.numeric(df$Ccar_whole)
df$Cant_whole <- as.numeric(df$Cant_whole)
df$EWT_whole <- as.numeric(df$EWT_whole)
df$LMA_whole <- as.numeric(df$LMA_whole)
df$N_prior <- as.numeric(df$N_prior)
df$ratio_chl_car <- as.numeric(df$ratio_chl_car)

Treatment_Subset <- df[df$Group == "1Treatment" & 
                         df$Block %in% c("A", "C", "E") & 
                         df$DoT > 12, ]

Control_Subset <- df[df$Group == "1Control" & 
                       df$Block %in% c("B", "D", "F") & 
                       df$DoT > 12, ]

RestControl_Subset <- df[df$Group == "1Control" & 
                           df$Block %in% c("J", "H"),]

RestTreatment_Subset <- df[df$Group == "1Treatment" & 
                             df$Block %in% c("G", "I"),]

Combined <- rbind(Treatment_Subset, Control_Subset, RestControl_Subset, RestTreatment_Subset)

filteredwilting_combined <- Combined[Combined$Wilting <= 1, ]
filtered_combined <- filteredwilting_combined[filteredwilting_combined$NDWI >= 0, ]


tree.means <- filtered_combined %>% #do means out of leaf 1 and leaf 2
  group_by(tree.ID, Date, Position, Block, Group, DoT, Comment, Site) %>%
  summarise_all(mean)
tree.means_df <- as.data.frame(tree.means)
tree.means_df$Cluster <- as.factor(tree.means_df$Cluster.23)

lys <- c("RGR", "ARI", "CARI", "NDSR", "CAI", "CCI", "CIRE", "NDLMA", "NDNI",
         "ChlNDI", "NDLI", "NDWI", "PRI", "EWI", "MSI", "WABI", "WI", "EWT_whole", "Ccab_whole", "Ccar_whole", "LMA_whole", "Cant_whole", "ratio_chl_car")

# Test qq
quartz()
par(mfrow = c(5,5), mar = c(2, 2, 1, 1))
for (i in lys) {
  qqnorm(tree.means_df[[i]], main = paste("qq of", i), ylab = "", xlab = "")
  qqline(tree.means_df[[i]])
}

# Test histograms
quartz()
par(mfrow = c(5,5), mar = c(2, 2, 1, 1))
for (i in lys) {
  hist(tree.means_df[[i]], main = paste("hist of", i), ylab = "", xlab = "")
}

#ANOVA model

lm_model_avg <- function(variable_name, data, output_folder){
  
  result_shapiro <- shapiro.test(data[[variable_name]])
  
  lm2 <- lm(get(variable_name) ~ Length + Group * Cluster, data = data)
  print(anova(lm2))
  
  output_file <- file.path(output_folder, paste(variable_name, "_analysis.txt", sep = ""))
  sink(output_file, append = FALSE)
  
  
  cat("Variable Name: ", variable_name, "\n")
  cat(paste("Shapiro-Wilk test for normality p-value:", result_shapiro$p.value), "\n")
  cat("\n Summary of ANOVA Model with Length:\n")
  print(variable_name)
  print(anova(lm2))
  
  sink()
  
  pdf_file <- file.path(output_folder, paste(variable_name, "_plots.pdf", sep = ""))
  pdf(pdf_file)
  par(mfrow = c(2, 2))
  plot(lm2)
  dev.off()
}

output_folder <- ""
for (i in lys) {
  print(i)
  lm_model_avg(i, tree.means_df, output_folder)
  
}

lm_model_avg(ratio_chl_car, tree.means_df, output_folder)
lm2 <- lm(ratio_chl_car ~ Length + Group * Cluster, data = tree.means_df)
print(anova(lm2))

# Tukey HSD
lm2 <- lm(NDLI ~ Length + Group * Cluster, data = tree.means_df)
print(anova(lm2))
posthoc_result <- HSD.test(lm2, "Cluster")
print(summary(posthoc_result))
print(posthoc_result$groups)
print(posthoc_result$pvalues)
print(posthoc_result$means)
print(posthoc_result$parameters)

lm2 <- lm(NDNI ~ Length + Group * Cluster, data = tree.means_df)
print(anova(lm2))
posthoc_result <- HSD.test(lm2, "Cluster")
print(summary(posthoc_result))
print(posthoc_result$groups)
print(posthoc_result$pvalues)
print(posthoc_result$means)
print(posthoc_result$parameters)

# stomatal conductance
df_sc <- read.xlsx("") 

lm_model_sc <- function(variable_name, data, output_folder){
  
  result_shapiro <- shapiro.test(data[[variable_name]])
  
  lm2 <- lm(get(variable_name) ~ length + group * cluster.23, data = data)
  print(anova(lm2))
  
  output_file <- file.path(output_folder, paste(variable_name, "_analysis.txt", sep = ""))
  sink(output_file, append = FALSE)
  
  
  cat("Variable Name: ", variable_name, "\n")
  cat(paste("Shapiro-Wilk test for normality p-value:", result_shapiro$p.value), "\n")
  cat("\n Summary of ANOVA Model with Length:\n")
  print(variable_name)
  print(anova(lm2))
  
  sink()
  
  pdf_file <- file.path(output_folder, paste(variable_name, "_plots.pdf", sep = ""))
  pdf(pdf_file)
  par(mfrow = c(2, 2))
  plot(lm2)
  dev.off()
}

lm_model_sc("stomcondtop", df_sc, "")



# calculate LRR ----

#create desired order
desired_order <- c("Cluster 3", "Cluster 2", "Cluster 1")
tree.means_df$Cluster <- factor(tree.means_df$Cluster, levels = desired_order)

#calculate LRR
C = 1
log_response_ratios <- tree.means_df %>%
  group_by(Cluster) %>%
  summarize(
    LRR_RGR = log(mean(RGR[Group == "1Treatment"]+C) / mean(RGR[Group == "1Control"]+C)),
    LRR_ARI = log(mean(ARI[Group == "1Treatment"]+3.33) / mean(ARI[Group == "1Control"]+3.33)), #one issue: max neg value is 2.4 -- constant does not make it pos.
    LRR_CAI = log(mean(CAI[Group == "1Treatment"]+C) / mean(CAI[Group == "1Control"]+C)),
    LRR_CCI = log(mean(CCI[Group == "1Treatment"]+C) / mean(CCI[Group == "1Control"]+C)),
    LRR_CARI = log(mean(CARI[Group == "1Treatment"]+C) / mean(CARI[Group == "1Control"]+C)),
    LRR_NDSR = log(mean(NDSR[Group == "1Treatment"]+C) / mean(NDSR[Group == "1Control"]+C)),
    LRR_CIRE = log(mean(CIRE[Group == "1Treatment"]+C) / mean(CIRE[Group == "1Control"]+C)),
    LRR_NDLMA = log(mean(NDLMA[Group == "1Treatment"]+C) / mean(NDLMA[Group == "1Control"]+C)),
    LRR_NDNI = log(mean(NDNI[Group == "1Treatment"]+C) / mean(NDNI[Group == "1Control"]+C)),
    LRR_NDLI = log(mean(NDLI[Group == "1Treatment"]+C) / mean(NDLI[Group == "1Control"]+C)),
    LRR_ChlNDI = log(mean(ChlNDI[Group == "1Treatment"]+C) / mean(ChlNDI[Group == "1Control"]+C)),
    LRR_NDWI = log(mean(NDWI[Group == "1Treatment"]+C) / mean(NDWI[Group == "1Control"]+C)),
    LRR_PRI = log(mean(PRI[Group == "1Treatment"]+C) / mean(PRI[Group == "1Control"]+C)),
    LRR_EWI = log(mean(EWI[Group == "1Treatment"]+C) / mean(EWI[Group == "1Control"]+C)),
    LRR_MSI = log(mean(MSI[Group == "1Treatment"]+C) / mean(MSI[Group == "1Control"]+C)),
    LRR_WABI = log(mean(WABI[Group == "1Treatment"]+C) / mean(WABI[Group == "1Control"]+C)),
    LRR_WI = log(mean(WI[Group == "1Treatment"]+C) / mean(WI[Group == "1Control"]+C)),
    SE_RGR = sqrt( ((sd(RGR[Group == "1Treatment"]+C))**2 / (length(RGR[Group == "1Treatment"]+C) * (mean(RGR[Group == "1Treatment"]+C))**2)) + ((sd(RGR[Group == "1Control"]+C))**2 / (length(RGR[Group == "1Control"]+C) * (mean(RGR[Group == "1Control"]+C))**2)) ),
    SE_ARI = sqrt( ((sd(ARI[Group == "1Treatment"]+2.5))**2 / (length(ARI[Group == "1Treatment"]+2.5) * (mean(ARI[Group == "1Treatment"]+2.5))**2)) + ((sd(ARI[Group == "1Control"]+2.5))**2 / (length(ARI[Group == "1Control"]+2.5) * (mean(ARI[Group == "1Control"]+2.5))**2)) ),  
    SE_CAI = sqrt( ((sd(CAI[Group == "1Treatment"]+C))**2 / (length(CAI[Group == "1Treatment"]+C) * (mean(CAI[Group == "1Treatment"]+C))**2)) + ((sd(CAI[Group == "1Control"]+C))**2 / (length(CAI[Group == "1Control"]+C) * (mean(CAI[Group == "1Control"]+C))**2)) ),
    SE_CCI = sqrt( ((sd(CCI[Group == "1Treatment"]+C))**2 / (length(CCI[Group == "1Treatment"]+C) * (mean(CCI[Group == "1Treatment"]+C))**2)) + ((sd(CCI[Group == "1Control"]+C))**2 / (length(CCI[Group == "1Control"]+C) * (mean(CCI[Group == "1Control"]+C))**2)) ),
    SE_CARI = sqrt( ((sd(CARI[Group == "1Treatment"]+C))**2 / (length(CARI[Group == "1Treatment"]+C) * (mean(CARI[Group == "1Treatment"]+C))**2)) + ((sd(CARI[Group == "1Control"]+C))**2 / (length(CARI[Group == "1Control"]+C) * (mean(CARI[Group == "1Control"]+C))**2)) ),
    SE_NDSR = sqrt( ((sd(NDSR[Group == "1Treatment"]+C))**2 / (length(NDSR[Group == "1Treatment"]+C) * (mean(NDSR[Group == "1Treatment"]+C))**2)) + ((sd(NDSR[Group == "1Control"]+C))**2 / (length(NDSR[Group == "1Control"]+C) * (mean(NDSR[Group == "1Control"]+C))**2)) ),
    SE_CIRE = sqrt( ((sd(CIRE[Group == "1Treatment"]+C))**2 / (length(CIRE[Group == "1Treatment"]+C) * (mean(CIRE[Group == "1Treatment"]+C))**2)) + ((sd(CIRE[Group == "1Control"]+C))**2 / (length(CIRE[Group == "1Control"]+C) * (mean(CIRE[Group == "1Control"]+C))**2)) ),
    SE_NDLMA = sqrt( ((sd(NDLMA[Group == "1Treatment"]+C))**2 / (length(NDLMA[Group == "1Treatment"]+C) * (mean(NDLMA[Group == "1Treatment"]+C))**2)) + ((sd(NDLMA[Group == "1Control"]+C))**2 / (length(NDLMA[Group == "1Control"]+C) * (mean(NDLMA[Group == "1Control"]+C))**2)) ),
    SE_NDNI = sqrt( ((sd(NDNI[Group == "1Treatment"]+C))**2 / (length(NDNI[Group == "1Treatment"]+C) * (mean(NDNI[Group == "1Treatment"]+C))**2)) + ((sd(NDNI[Group == "1Control"]+C))**2 / (length(NDNI[Group == "1Control"]+C) * (mean(NDNI[Group == "1Control"]+C))**2)) ),
    SE_NDLI = sqrt( ((sd(NDLI[Group == "1Treatment"]+C))**2 / (length(NDLI[Group == "1Treatment"]+C) * (mean(NDLI[Group == "1Treatment"]+C))**2)) + ((sd(NDLI[Group == "1Control"]+C))**2 / (length(NDLI[Group == "1Control"]+C) * (mean(NDLI[Group == "1Control"]+C))**2)) ),
    SE_ChlNDI = sqrt( ((sd(ChlNDI[Group == "1Treatment"]+C))**2 / (length(ChlNDI[Group == "1Treatment"]+C) * (mean(ChlNDI[Group == "1Treatment"]+C))**2)) + ((sd(ChlNDI[Group == "1Control"]+C))**2 / (length(ChlNDI[Group == "1Control"]+C) * (mean(ChlNDI[Group == "1Control"]+C))**2)) ),
    SE_NDWI = sqrt( ((sd(NDWI[Group == "1Treatment"]+C))**2 / (length(NDWI[Group == "1Treatment"]+C) * (mean(NDWI[Group == "1Treatment"]+C))**2)) + ((sd(NDWI[Group == "1Control"]+C))**2 / (length(NDWI[Group == "1Control"]+C) * (mean(NDWI[Group == "1Control"]+C))**2)) ),
    SE_PRI = sqrt( ((sd(PRI[Group == "1Treatment"]+C))**2 / (length(PRI[Group == "1Treatment"]+C) * (mean(PRI[Group == "1Treatment"]+C))**2)) + ((sd(PRI[Group == "1Control"]+C))**2 / (length(PRI[Group == "1Control"]+C) * (mean(PRI[Group == "1Control"]+C))**2)) ),
    SE_EWI = sqrt( ((sd(EWI[Group == "1Treatment"]+C))**2 / (length(EWI[Group == "1Treatment"]+C) * (mean(EWI[Group == "1Treatment"]+C))**2)) + ((sd(EWI[Group == "1Control"]+C))**2 / (length(EWI[Group == "1Control"]+C) * (mean(EWI[Group == "1Control"]+C))**2)) ),
    SE_MSI = sqrt( ((sd(MSI[Group == "1Treatment"]+C))**2 / (length(MSI[Group == "1Treatment"]+C) * (mean(MSI[Group == "1Treatment"]+C))**2)) + ((sd(MSI[Group == "1Control"]+C))**2 / (length(MSI[Group == "1Control"]+C) * (mean(MSI[Group == "1Control"]+C))**2)) ),
    SE_WABI = sqrt( ((sd(WABI[Group == "1Treatment"]+C))**2 / (length(WABI[Group == "1Treatment"]+C) * (mean(WABI[Group == "1Treatment"]+C))**2)) + ((sd(WABI[Group == "1Control"]+C))**2 / (length(WABI[Group == "1Control"]+C) * (mean(WABI[Group == "1Control"]+C))**2)) ),
    SE_WI = sqrt( ((sd(WI[Group == "1Treatment"]+C))**2 / (length(WI[Group == "1Treatment"]+C) * (mean(WI[Group == "1Treatment"]+C))**2)) + ((sd(WI[Group == "1Control"]+C))**2 / (length(WI[Group == "1Control"]+C) * (mean(WI[Group == "1Control"]+C))**2)) ),
    LRR_Ccab_opt = log(mean(Ccab_opt[Group == "1Treatment"]) / mean(Ccab_opt[Group == "1Control"])),
    LRR_LMA_opt = log(mean(LMA_opt[Group == "1Treatment"]) / mean(LMA_opt[Group == "1Control"])),
    LRR_Cant_opt = log(mean(Cant_opt[Group == "1Treatment"]) / mean(Cant_opt[Group == "1Control"])),
    LRR_EWT_opt = log(mean(EWT_opt[Group == "1Treatment"]) / mean(EWT_opt[Group == "1Control"])),
    LRR_Ccar_opt = log(mean(Ccar_opt[Group == "1Treatment"]) / mean(Ccar_opt[Group == "1Control"])),
    LRR_CCR_opt = log(mean(CCR_opt[Group == "1Treatment"]) / mean(CCR_opt[Group == "1Control"])),
    LRR_Ccab_whole = log(mean(Ccab_whole[Group == "1Treatment"]) / mean(Ccab_whole[Group == "1Control"])),
    LRR_LMA_whole = log(mean(LMA_whole[Group == "1Treatment"]) / mean(LMA_whole[Group == "1Control"])),
    LRR_Cant_whole = log(mean(Cant_whole[Group == "1Treatment"]) / mean(Cant_whole[Group == "1Control"])),
    LRR_EWT_whole = log(mean(EWT_whole[Group == "1Treatment"]) / mean(EWT_whole[Group == "1Control"])),
    LRR_Ccar_whole = log(mean(Ccar_whole[Group == "1Treatment"]) / mean(Ccar_whole[Group == "1Control"])),
    LRR_CCR_whole = log(mean(CCR_whole[Group == "1Treatment"]) / mean(CCR_whole[Group == "1Control"])),
    SE_Ccab_opt = sqrt( ((sd(Ccab_opt[Group == "1Treatment"]))**2 / (length(Ccab_opt[Group == "1Treatment"]) * (mean(Ccab_opt[Group == "1Treatment"]))**2)) + ((sd(Ccab_opt[Group == "1Control"]))**2 / (length(Ccab_opt[Group == "1Control"]) * (mean(Ccab_opt[Group == "1Control"]))**2)) ),
    SE_LMA_opt = sqrt( ((sd(LMA_opt[Group == "1Treatment"]))**2 / (length(LMA_opt[Group == "1Treatment"]) * (mean(LMA_opt[Group == "1Treatment"]))**2)) + ((sd(LMA_opt[Group == "1Control"]))**2 / (length(LMA_opt[Group == "1Control"]) * (mean(LMA_opt[Group == "1Control"]))**2)) ),  
    SE_Cant_opt = sqrt( ((sd(Cant_opt[Group == "1Treatment"]))**2 / (length(Cant_opt[Group == "1Treatment"]) * (mean(Cant_opt[Group == "1Treatment"]))**2)) + ((sd(Cant_opt[Group == "1Control"]))**2 / (length(Cant_opt[Group == "1Control"]) * (mean(Cant_opt[Group == "1Control"]))**2)) ), 
    SE_EWW_opt = sqrt( ((sd(EWT_opt[Group == "1Treatment"]))**2 / (length(EWT_opt[Group == "1Treatment"]) * (mean(EWT_opt[Group == "1Treatment"]))**2)) + ((sd(EWT_opt[Group == "1Control"]))**2 / (length(EWT_opt[Group == "1Control"]) * (mean(EWT_opt[Group == "1Control"]))**2)) ), 
    SE_Ccar_opt = sqrt( ((sd(Ccar_opt[Group == "1Treatment"]))**2 / (length(Ccar_opt[Group == "1Treatment"]) * (mean(Ccar_opt[Group == "1Treatment"]))**2)) + ((sd(Ccar_opt[Group == "1Control"]))**2 / (length(Ccar_opt[Group == "1Control"]) * (mean(Ccar_opt[Group == "1Control"]))**2)) ),
    SE_CCR_opt = sqrt( ((sd(CCR_opt[Group == "1Treatment"]))**2 / (length(CCR_opt[Group == "1Treatment"]) * (mean(CCR_opt[Group == "1Treatment"]))**2)) + ((sd(CCR_opt[Group == "1Control"]))**2 / (length(CCR_opt[Group == "1Control"]) * (mean(CCR_opt[Group == "1Control"]))**2)) ),
    SE_Ccab_whole = sqrt( ((sd(Ccab_whole[Group == "1Treatment"]))**2 / (length(Ccab_whole[Group == "1Treatment"]) * (mean(Ccab_whole[Group == "1Treatment"]))**2)) + ((sd(Ccab_whole[Group == "1Control"]))**2 / (length(Ccab_whole[Group == "1Control"]) * (mean(Ccab_whole[Group == "1Control"]))**2)) ),
    SE_LMA_whole = sqrt( ((sd(LMA_whole[Group == "1Treatment"]))**2 / (length(LMA_whole[Group == "1Treatment"]) * (mean(LMA_whole[Group == "1Treatment"]))**2)) + ((sd(LMA_whole[Group == "1Control"]))**2 / (length(LMA_whole[Group == "1Control"]) * (mean(LMA_whole[Group == "1Control"]))**2)) ),  
    SE_Cant_whole = sqrt( ((sd(Cant_whole[Group == "1Treatment"]))**2 / (length(Cant_whole[Group == "1Treatment"]) * (mean(Cant_whole[Group == "1Treatment"]))**2)) + ((sd(Cant_whole[Group == "1Control"]))**2 / (length(Cant_whole[Group == "1Control"]) * (mean(Cant_whole[Group == "1Control"]))**2)) ), 
    SE_EWT_whole = sqrt( ((sd(EWT_whole[Group == "1Treatment"]))**2 / (length(EWT_whole[Group == "1Treatment"]) * (mean(EWT_whole[Group == "1Treatment"]))**2)) + ((sd(EWT_whole[Group == "1Control"]))**2 / (length(EWT_whole[Group == "1Control"]) * (mean(EWT_whole[Group == "1Control"]))**2)) ), 
    SE_Ccar_whole = sqrt( ((sd(Ccar_whole[Group == "1Treatment"]))**2 / (length(Ccar_whole[Group == "1Treatment"]) * (mean(Ccar_whole[Group == "1Treatment"]))**2)) + ((sd(Ccar_whole[Group == "1Control"]))**2 / (length(Ccar_whole[Group == "1Control"]) * (mean(Ccar_opt[Group == "1Control"]))**2)) ),
    SE_CCR_whole = sqrt( ((sd(CCR_whole[Group == "1Treatment"]))**2 / (length(CCR_whole[Group == "1Treatment"]) * (mean(CCR_whole[Group == "1Treatment"]))**2)) + ((sd(CCR_whole[Group == "1Control"]))**2 / (length(CCR_whole[Group == "1Control"]) * (mean(CCR_whole[Group == "1Control"]))**2)) )
  )

#stomatal conductance LRR

#create desired order
desired_order <- c("Cluster 3", "Cluster 2", "Cluster 1")
df_sc$cluster.23 <- factor(df_sc$cluster.23, levels = desired_order)

lrr_sc <- df_sc %>%
  group_by(cluster.23) %>%
  summarize(
    lrr = log(mean(stomcond[group == "treatment"]) / mean(stomcond[group == "control"])),
    se = sqrt( ((sd(stomcond[group == "treatment"]))**2 / (length(stomcond[group == "treatment"]) * (mean(stomcond[group == "treatment"]))**2)) + ((sd(stomcond[group == "control"]))**2 / (length(stomcond[group == "control"]) * (mean(stomcond[group == "control"]))**2)) ),
    lrr_top = log(mean(stomcondtop[group == "treatment"]) / mean(stomcondtop[group == "control"])),
    se_top = sqrt( ((sd(stomcondtop[group == "treatment"]))**2 / (length(stomcondtop[group == "treatment"]) * (mean(stomcondtop[group == "treatment"]))**2)) + ((sd(stomcondtop[group == "control"]))**2 / (length(stomcondtop[group == "control"]) * (mean(stomcondtop[group == "control"]))**2)) )
  )

#visualisation (fig.3)
#change MSI to inverse values
lrr_new <- log_response_ratios %>%
  mutate(LRR_MSI = LRR_MSI * (-1), SE_MSI = SE_MSI * (-1))

#chlorophyll
plot_chl <- ggplot(lrr_new, aes(x = LRR_NDSR, y = Cluster)) +
  geom_point(position = position_nudge(y = 0.2), color = "#66c2a4") +
  geom_errorbarh(aes(xmin = LRR_NDSR - SE_NDSR, xmax = LRR_NDSR + SE_NDSR), height = 0.1, position = position_nudge(y = 0.2), color = "#66c2a4") + 
  geom_point(data = lrr_new, aes(x = LRR_CIRE, y = Cluster), color = "#238b45") +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_CIRE - SE_CIRE, xmax = LRR_CIRE + SE_CIRE), height = 0.1, color = "#238b45") +
  geom_point(data = lrr_new, aes(x = LRR_Ccab_whole, y = Cluster), color = "#00441b", position = position_nudge(y = -0.2), shape = 17) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_Ccab_whole - SE_Ccab_whole, xmax = LRR_Ccab_whole + SE_Ccab_whole), position = position_nudge(y = -0.2), height = 0.1, color = "#00441b") +
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") ) +
  
  quartz()
print(plot_chl)

#carotenoids

plot_car <- ggplot(lrr_new, aes(x = LRR_CCI, y = Cluster)) +
  geom_point(position = position_nudge(y = 0.2), color = "#fe9929") +
  geom_errorbarh(aes(xmin = LRR_CCI - SE_CCI, xmax = LRR_CCI + SE_CCI), height = 0.1, position = position_nudge(y = 0.2), color = "#fe9929") + 
  geom_point(data = lrr_new, aes(x = LRR_PRI, y = Cluster), color = "#cc4c02") +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_PRI - SE_PRI, xmax = LRR_PRI + SE_PRI), height = 0.1, color = "#cc4c02") +
  geom_point(data = lrr_new, aes(x = LRR_Ccar_whole, y = Cluster), color = "#662506", position = position_nudge(y = -0.2), shape = 17) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_Ccar_whole - SE_Ccar_whole, xmax = LRR_Ccar_whole + SE_Ccar_whole), position = position_nudge(y = -0.2), height = 0.1, color = "#662506") +
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") ) +
  
  quartz()
print(plot_car)

#antocyanins

plot_ant <- ggplot(lrr_new, aes(x = LRR_RGR, y = Cluster)) +
  geom_point(position = position_nudge(y = 0.2), color = "#df65b0") +
  geom_errorbarh(aes(xmin = LRR_RGR - SE_RGR, xmax = LRR_RGR + SE_RGR), height = 0.1, position = position_nudge(y = 0.2), color = "#df65b0") + 
  geom_point(data = lrr_new, aes(x = LRR_ARI, y = Cluster), color = "#ce1256") +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_ARI - SE_ARI, xmax = LRR_ARI + SE_ARI), height = 0.1, color = "#ce1256") +
  geom_point(data = lrr_new, aes(x = LRR_Cant_whole, y = Cluster), color = "#67001f", position = position_nudge(y = -0.2), shape = 17) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_Cant_whole - SE_Cant_whole, xmax = LRR_Cant_whole + SE_Cant_whole), position = position_nudge(y = -0.2), height = 0.1, color = "#67001f") +
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") ) +
  
  quartz()
print(plot_ant)

#water 

plot_water <- ggplot(lrr_new, aes(x = LRR_NDWI, y = Cluster)) +
  geom_point(position = position_nudge(y = 0.3), color = "#9ecae1") +
  geom_errorbarh(aes(xmin = LRR_NDWI - SE_NDWI, xmax = LRR_NDWI + SE_NDWI), height = 0.1, position = position_nudge(y = 0.3), color = "#9ecae1") + 
  geom_point(data = lrr_new, aes(x = LRR_MSI, y = Cluster), color = "#6baed6", position = position_nudge(y = 0.1)) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_MSI - SE_MSI, xmax = LRR_MSI + SE_MSI), height = 0.1, position = position_nudge(y=0.1), color = "#6baed6") +
  geom_point(data = lrr_new, aes(x = LRR_WABI, y = Cluster), color = "#2171b5", position = position_nudge(y = -0.1)) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_WABI - SE_WABI, xmax = LRR_WABI + SE_WABI), height = 0.1, position = position_nudge(y=-0.1), color = "#2171b5") +
  geom_point(data = lrr_new, aes(x = LRR_EWT_whole, y = Cluster), color = "#08306b", position = position_nudge(y = -0.3), shape = 17) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_EWT_whole - SE_EWT_whole, xmax = LRR_EWT_whole + SE_EWT_whole), position = position_nudge(y = -0.3), height = 0.1, color = "#08306b") +
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") ) +
  
  quartz()
print(plot_water)

#rest

plot_rest <- ggplot(lrr_new, aes(x = LRR_NDNI, y = Cluster)) +
  geom_point(position = position_nudge(y = 0.3), color = "#dadaeb") +
  geom_errorbarh(aes(xmin = LRR_NDNI - SE_NDNI, xmax = LRR_NDNI + SE_NDNI), height = 0.1, position = position_nudge(y = 0.3), color = "#dadaeb") + 
  geom_point(data = lrr_new, aes(x = LRR_NDLI, y = Cluster), color = "#9e9ac8", position = position_nudge(y = 0.1)) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_NDLI - SE_NDLI, xmax = LRR_NDLI + SE_NDLI), height = 0.1, position = position_nudge(y=0.1), color = "#9e9ac8") +
  geom_point(data = lrr_new, aes(x = LRR_EWI, y = Cluster), color = "#6a51a3", position = position_nudge(y = -0.1)) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_EWI - SE_EWI, xmax = LRR_EWI + SE_EWI), height = 0.1, position = position_nudge(y=-0.1), color = "#6a51a3") +
  geom_point(data = lrr_new, aes(x = LRR_LMA_whole, y = Cluster), color = "#3f007d", position = position_nudge(y = -0.3), shape = 17) +
  geom_errorbarh(data = lrr_new, aes(xmin = LRR_LMA_whole - SE_LMA_whole, xmax = LRR_LMA_whole + SE_LMA_whole), position = position_nudge(y = -0.3), height = 0.1, color = "#3f007d") +
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") ) +
  quartz()
print(plot_rest)

#sc top
plot_sc_top <- ggplot(lrr_sc, aes(x = lrr_top, y = cluster.23)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lrr_top - se_top, xmax = lrr_top + se_top), height = 0.1, position = position_dodge(width = 0.4), color = "black") +
  #geom_jitter(data = filtered_all, aes(x = NDSR, y = Site, color = Group), width = 0.1, height = 0, alpha = 0.9, size = 0.8) +
  #scale_color_manual(values = c("1Control" = "#5BBCD6", "1Treatment" = "#CB7A5C")) + 
  theme_minimal() +
  xlab(NULL) +  
  ylab(NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.y = element_blank()) +
  ggtitle(NULL) +
  theme(plot.title = element_text(size = 10, family = "Arial"), axis.text.x = element_text(size = 6, family = "Arial"), axis.text.y = element_text(size = 7, family = "Arial") )

quartz()
grid.arrange(plot_sc_top, ncol = 1, nrow = 1)

# raw spectra plot (fig2) ---- 

df_spectra_first <- read.xlsx("") 
#average leaf spectra (did not happen in analysis, the traits were averaged but not the spectra)
df_spectra_first <- df_spectra_first %>%
  group_by(tree.ID, Date, Position, Group, Comment, Site, Cluster.23) %>%
  summarise_all(mean)
df_spectra_first$Cluster.23 <- as.factor(df_spectra_first$Cluster.23)
df_spectra_first$Group <- as.factor(df_spectra_first$Group)

#subsetting
df_spectra_first <- subset(df_spectra_first, select = -c(tree.ID,Date,Comment,Site,Leaf,SPAD))
# and deleting columns
df_spectra_first <- df_spectra_first[, -c(4:54)]

#separate in the three clusters

cluster_1 <- df_spectra_first %>% filter(Cluster.23 == 1)
cluster_2 <- df_spectra_first %>% filter(Cluster.23 == 2)
cluster_3 <- df_spectra_first %>% filter(Cluster.23 == 3)

color_mapping = c("1Treatment" = "#56B4E9" , "1Control" = "#E69F00")
ylim_range <- c(0.0, 0.7)

spectra1 <- as_spectra(cluster_1, name_idx = 1, meta_idx = 2:3)
quartz()
plot(spectra1, lwd = 0.75, lty = 1, main = "1", cex.main=1.4, cex.lab=1.2,
     xlab = "Wavelength [nm]", ylab = "Reflectance [%]", col = color_mapping[cluster_1$Group], ylim = ylim_range)
line_types <- c("Treatment" = 1, "Control" = 1)  # num for line type
#legend("topright", legend = levels(cluster_1$Group), col = color_mapping, lty = line_types)
#color_mapping <- scale_color_manual(values = c("Treatment" = "blue", "Control" = "red"))

spectra2 <- as_spectra(cluster_2, name_idx = 1, meta_idx = 2:3)
quartz()
plot(spectra2, lwd = 0.75, lty = 1, main = "2", cex.main=1.4, cex.lab=1.2,
     xlab = "Wavelength [nm]", ylab = "Reflectance [%]", col = color_mapping[cluster_2$Group], ylim = ylim_range)
line_types <- c("Treatment" = 1, "Control" = 1)  # num for line type
#legend("topright", legend = levels(cluster_2$Group), col = color_mapping, lty = line_types)
#color_mapping <- scale_color_manual(values = c("Treatment" = "blue", "Control" = "red"))

spectra3 <- as_spectra(cluster_3, name_idx = 1, meta_idx = 2:3)
quartz()
plot(spectra3, lwd = 0.75, lty = 1, main = "3", cex.main=1.4, cex.lab=1.2,
     xlab = "Wavelength [nm]", ylab = "Reflectance [%]", col = color_mapping[cluster_3$Group], ylim = ylim_range)
line_types <- c("Treatment" = 1, "Control" = 1)  # num for line type
#legend("topright", legend = levels(cluster_3$Group), col = color_mapping, lty = line_types, inset = c(0.05, 0.05))
#color_mapping <- scale_color_manual(values = c("Treatment" = "blue", "Control" = "red"))
# 
# prospect inversion validation (fig4) ---- 

df_integration_area <- read.xlsx("") #note:"2" means lambda=665 chromatgram
keep <- c("name" ,"chla2","chladev2", "pheoa2", "pheoadev2", "chlb2", "chlbdev2", "pheob2", "pheobdev2", "caradev", "carb", "neo", "neodev", "lut", "lutdev", "millingweight", "dryweight", "leafarea", "group")
df_integration_area = df_integration_area[keep]


# calculate concentration based on calibration lines
df_integration_area <- df_integration_area %>%
  mutate(
    
    chla2 = ifelse(!is.na(chla2), chla2 /104.92 , chla2), #chl-a standards
    chladev2 = ifelse(!is.na(chladev2), chladev2 /104.92, chladev2),
    pheoa2 = ifelse(!is.na(pheoa2), pheoa2 /104.92, pheoa2),
    pheoadev2 = ifelse(!is.na(pheoadev2), pheoadev2 /104.92, pheoadev2),
    
    chlb2 = ifelse(!is.na(chlb2), chlb2/31.502 , chlb2), #chl-b standards
    chlbdev2 = ifelse(!is.na(chlbdev2), chlbdev2/31.502 , chlbdev2),
    pheob2 = ifelse(!is.na(pheob2), pheob2/31.502 ,  pheob2),
    pheobdev2 = ifelse(!is.na(pheobdev2), pheobdev2/31.502 ,  pheobdev2),
    
    caradev = ifelse(!is.na(caradev), caradev/405.37, caradev), #alpha-car standards
    
    carb = ifelse(!is.na(carb), (carb-285.43)/204.98, carb), #beta-car standards
    
    neo = ifelse(!is.na(neo), neo/204.98, neo),  #neoxanthin standards
    neodev = ifelse(!is.na(neodev), neodev/204.98, neodev),
    
    lut = ifelse(!is.na(lut), lut/366.08, lut), #lut (xant) standard
    lutdev = ifelse(!is.na(lutdev), lut/366.08, lutdev),
    
    #calculate concentration per area
    chla2 = chla2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,  
    chladev2 = chladev2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    pheoa2 = pheoa2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    pheoadev2 = pheoadev2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    
    chlb2 = chlb2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    chlbdev2 = chlbdev2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    pheoa2 = pheoa2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    pheoadev2 = pheoadev2*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    
    carb = carb*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    caradev = caradev*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    
    neo = neo*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    neodev = neodev*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    
    lut = lut*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    lutdev = lutdev*0.3*(1/millingweight)*dryweight*(1/leafarea)*1000,
    
  )

df_integration_area <- df_integration_area %>%
  mutate(
    name = name,
    group = group,
    chl = chla2+chlb2+chladev2+chlbdev2+pheoa2+pheob2+pheoadev2+pheobdev2,
    car = caradev+carb+neo+neodev+lut+lutdev
    
  )

df_validation_prospect <- read.xlsx("")

df_combined <- data.frame(chl = df_integration_area$chl, chl_prospect = df_validation_prospect$chl_prospect, group = df_integration_area$group)
correlation_coefficient <- cor(df_combined$chl, df_combined$chl_prospect)
quartz()
plot(df_combined$chl, df_combined$chl_prospect, 
     xlab = "measured chl", ylab = "prospect chl",
     main = paste("Scatterplot chl (Correlation =", round(correlation_coefficient, 2), ")"),
     col = ifelse(df_combined$group == "t", "brown", "blue"), pch = 16)

abline(lm(df_combined$chl_prospect ~ df_combined$chl), col = "red")
abline(a = 0, b = 1, col = "black")
legend("bottomright", legend = unique(df_combined$group), col = c("brown","blue"), pch = 16)


df_combined <- data.frame(car = df_integration_area$car, car_prospect = df_validation_prospect$car_prospect, group = df_integration_area$group)
correlation_coefficient <- cor(df_combined$car, df_combined$car_prospect)
quartz()
plot(df_combined$car, df_combined$car_prospect, 
     xlab = "measured car", ylab = "prospect car",
     main = paste("Scatterplot car (Correlation =", round(correlation_coefficient, 2), ")"),
     col = ifelse(df_combined$group == "t", "brown", "blue"), pch = 16)

abline(lm(df_combined$car_prospect ~ df_combined$car), col = "red")
abline(a = 0, b = 1, col = "black")
legend("bottomright", legend = unique(df_combined$group), col = c("brown","blue"), pch = 16)


#ewt_lma
df_validation_ewt <- read.xlsx("")

correlation_coefficient <- cor(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect)
quartz()
plot(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect, 
     xlab = "measured ewt", ylab = "prospect ewt",
     main = paste("Scatterplot ewt (Correlation =", round(correlation_coefficient, 2), ")"),
     col = ifelse(df_validation_ewt$group == "t", "brown", "blue"), pch = 16)

abline(lm(df_validation_ewt$ewt_prospect ~ df_validation_ewt$ewt_measured), col = "red")
abline(a = 0, b = 1, col = "black")
legend("bottomright", legend = unique(df_validation_ewt$group), col = c("brown","blue"), pch = 16)

correlation_coefficient <- cor(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect)
quartz()
plot(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect, 
     xlab = "measured lma", ylab = "prospect lma",
     main = paste("Scatterplot lma (Correlation =", round(correlation_coefficient, 2), ")"),
     col = ifelse(df_validation_ewt$group == "t", "brown", "blue"), pch = 16)

abline(lm(df_validation_ewt$lma_prospect ~ df_validation_ewt$lma_measured), col = "red")
abline(a = 0, b = 1, col = "black")
legend("bottomright", legend = unique(df_validation_ewt$group), col = c("brown","blue"), pch = 16)


df_combined_chl <- data.frame(chl = df_integration_area$chl, chl_prospect = df_validation_prospect$chl_prospect, group = df_integration_area$group)
df_combined_car <- data.frame(car = df_integration_area$car, car_prospect = df_validation_prospect$car_prospect, group = df_integration_area$group)

# change values to mg.cm^-1
df_validation_ewt <- df_validation_ewt %>%
  mutate(across(-group, ~ . * 1000))

# Open a graphics device
quartz()

# Set up a 2x2 panel layout
par(mfrow = c(2, 2))

# Plot 1: chl
quartz()
correlation_coefficient_chl <- cor(df_combined_chl$chl, df_combined_chl$chl_prospect)
plot(df_combined_chl$chl, df_combined_chl$chl_prospect, 
     xlab = "measured chl", ylab = "prospect chl",
     main = NULL,
     col = ifelse(df_combined_chl$group == "t", "#E69F00", "#56B4E9"), pch = 16)
abline(lm(df_combined_chl$chl_prospect ~ df_combined_chl$chl), col = "black")
abline(a = 0, b = 1, col = "black", lty = "dashed")

# Plot 2: car
quartz()
correlation_coefficient_car <- cor(df_combined_car$car, df_combined_car$car_prospect)
plot(df_combined_car$car, df_combined_car$car_prospect, 
     xlab = "measured car", ylab = "prospect car",
     main = NULL,
     col = ifelse(df_combined_car$group == "t", "#E69F00", "#56B4E9"), pch = 16)
abline(lm(df_combined_car$car_prospect ~ df_combined_car$car), col = "black")
abline(a = 0, b = 1, col = "black", lty = "dashed")

# Plot 3: ewt
quartz()
correlation_coefficient_ewt <- cor(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect)
plot(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect, 
     xlab = "measured ewt", ylab = "prospect ewt",
     main = NULL,
     col = ifelse(df_validation_ewt$group == "t", "#E69F00", "#56B4E9"), pch = 16)
abline(lm(df_validation_ewt$ewt_prospect ~ df_validation_ewt$ewt_measured), col = "black")
abline(a = 0, b = 1, col = "black", lty = "dashed")

# Plot 4: lma
quartz()
correlation_coefficient_lma <- cor(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect)
plot(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect, 
     xlab = "measured lma", ylab = "prospect lma",
     main = NULL,
     col = ifelse(df_validation_ewt$group == "t", "#E69F00", "#56B4E9"), pch = 16)
abline(lm(df_validation_ewt$lma_prospect ~ df_validation_ewt$lma_measured), col = "black")
abline(a = 0, b = 1, col = "black", lty = "dashed")
legend("bottomright", legend = unique(df_validation_ewt$group), col = c("#E69F00","#56B4E9"), pch = 16, cex = 1, inset = c(0.05, 0.05))




# Plot 1: chl
correlation_coefficient_chl <- cor(df_combined_chl$chl, df_combined_chl$chl_prospect)
chl_lm <- lm(df_combined_chl$chl_prospect ~ df_combined_chl$chl)
r_squared_chl <- summary(chl_lm)$r.squared
cat("R-squared value for chl:", r_squared_chl, "\n")
correlation_test <- cor.test(df_combined_chl$chl, df_combined_chl$chl_prospect)
correlation_coefficient_chl <- correlation_test$estimate
p_value_chl <- correlation_test$p.value
chl_lm <- lm(df_combined_chl$chl_prospect ~ df_combined_chl$chl)
r_squared_chl <- summary(chl_lm)$r.squared
cat("R-squared value for chl:", r_squared_chl, "\n")
cat("P-value for correlation:", p_value_chl, "\n")


# Plot 2: car
correlation_coefficient_car <- cor(df_combined_car$car, df_combined_car$car_prospect)
car_lm <- lm(df_combined_car$car_prospect ~ df_combined_car$car)
r_squared_car <- summary(car_lm)$r.squared
cat("R-squared value for car:", r_squared_car, "\n")
correlation_test_car <- cor.test(df_combined_car$car, df_combined_car$car_prospect)
correlation_coefficient_car <- correlation_test_car$estimate
p_value_car <- correlation_test_car$p.value
car_lm <- lm(df_combined_car$car_prospect ~ df_combined_car$car)
r_squared_car <- summary(car_lm)$r.squared
cat("R-squared value for car:", r_squared_car, "\n")
cat("P-value for correlation:", p_value_car, "\n")



# Plot 3: ewt
correlation_coefficient_ewt <- cor(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect)
ewt_lm <- lm(df_validation_ewt$ewt_prospect ~ df_validation_ewt$ewt_measured)
r_squared_ewt <- summary(ewt_lm)$r.squared
cat("R-squared value for ewt:", r_squared_ewt, "\n")
correlation_test_ewt <- cor.test(df_validation_ewt$ewt_measured, df_validation_ewt$ewt_prospect)
correlation_coefficient_ewt <- correlation_test_ewt$estimate
p_value_ewt <- correlation_test_ewt$p.value
ewt_lm <- lm(df_validation_ewt$ewt_prospect ~ df_validation_ewt$ewt_measured)
r_squared_ewt <- summary(ewt_lm)$r.squared
cat("R-squared value for ewt:", r_squared_ewt, "\n")
cat("P-value for correlation:", p_value_ewt, "\n")

# Plot 4: lma
correlation_coefficient_lma <- cor(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect)
lma_lm <- lm(df_validation_ewt$lma_prospect ~ df_validation_ewt$lma_measured)
r_squared_lma <- summary(lma_lm)$r.squared
cat("R-squared value for lma:", r_squared_lma, "\n")
correlation_test_lma <- cor.test(df_validation_ewt$lma_measured, df_validation_ewt$lma_prospect)
correlation_coefficient_lma <- correlation_test_lma$estimate
p_value_lma <- correlation_test_lma$p.value
lma_lm <- lm(df_validation_ewt$lma_prospect ~ df_validation_ewt$lma_measured)
r_squared_lma <- summary(lma_lm)$r.squared
cat("R-squared value for lma:", r_squared_lma, "\n")
cat("P-value for correlation:", p_value_lma, "\n")

#

# tms and weather plot ----
df_weather <- read.xlsx("") 

x <- df_weather$number

quartz()
layout(matrix(1:2, nrow = 2), heights = c(2, 2))

par(mar=c(2,2,2,5))
ylim <- c(500, 2500)
plot(x, df_weather$sm_t_avg, type="l", lty =1, ylim = ylim, xlab= "", ylab= "", xaxt="n", family="Courier", cex.axis=0.7, col = "#ff7f00")
lines(x, df_weather$sm_c_avg, type = "l", lty = 1, col = "#5aae61")
abline(v = c(1777, 3313), col="black", lty=2)
abline(v = c(5329, 6961), col="black", lty=2)
par(family="Courier")

par(new=TRUE)
ylim <- c(0, 40)
plot(x, df_weather$precip, yaxt="n", type="h", lwd=3, col="#386cb0", xlab= "", ylab= "", axes=FALSE, ylim = ylim)
axis(side=4,col.axis = "#386cb0", cex.axis=0.7, family="Courier",col="#386cb0")
par(family="Courier")
#legend("topright", legend=c("Treatment", "Control", "Precipitation"),
#col=c("#ff7f00", "#5aae61", "transparent"), lty=c(1, 1, 0),
#pch=c(NA, NA, 15), fill=c(NA, NA, "#386cb0"), border=c(NA, NA, NA),
#cex=0.5, box.lty=0)


par(mar=c(2,2,0,5))
ylim <- c(5, 45)
plot(x, df_weather$temp_t_avg, type="l", lty =1, ylim = ylim, xlab= "", ylab= "", family="Courier", cex.axis=0.7, col = "#9970ab", yaxt="n")
axis(side=2, col.axis = "#9970ab", cex.axis=0.7, family="Courier",col="#9970ab")
abline(v = c(1777, 3313), col="black", lty=2)
abline(v = c(5329, 6961), col="black", lty=2)
#lines(x, df_weather$temp_c_avg, type = "l", lty = 1, col = "grey") no need for it because values are essentially the same

par(new=TRUE)
ylim <- c(5, 45)
plot(x, df_weather$airtemp, yaxt="n", type="p", lty=1, cex=0.7, lwd=1, col="#f0027f", xlab= "", ylab= "", axes=FALSE, ylim = ylim)
axis(side=4,col.axis = "#f0027f", cex.axis=0.7, family="Courier",col="#f0027f")

par(new=TRUE)
ylim <- c(-200, 100)
plot(x, df_weather$rh, yaxt="n", type="p", cex=0.7, pch=16, lwd=0.5, col="#01665e", xlab= "", ylab= "", axes=FALSE, ylim = ylim)
axis(side=4,col.axis = "#01665e", cex.axis=0.7, family="Courier", col="#01665e", line=3)
par(family="Courier")
#legend("bottom", legend=c("Temperature TMS-4", "Mean temperature ATMOS41", "Relative humidity"),
#col=c("#9970ab", "#f0027f", "#01665e", "black"), lty=c(1, 0, 0, 2), pch=c(NA, 1, 16, NA),
#cex=0.5, box.lty=0)


# length, soil moisture, stomatal conductance ----
df_sm_length <- read.xlsx("") 

p <- ggplot(df_sm_length, aes(x = Length, y = Mean, color = Group)) +
  geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = Group, fill = Group)) +
  theme_minimal()

p <- p + labs(
  x = "Length [cm]",
  y = "Mean raw uncalibrated TDT values [-]"
)

custom_colors <- c("Control" = "#35978F", "Treatment" = "#BF812D")

p <- p + theme(
  text = element_text(family = "Arial")
) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) 

quartz()
print(p)

df_sm_length$Group <- factor(df_sm_length$Group)
group_names <- levels(df_sm_length$Group)

for (group in group_names) {
  subset_data <- subset(df_sm_length, Group == group)
  lm_model <- lm(Mean ~ Length, data = subset_data)
  summary_lm <- summary(lm_model)
  r_squared <- summary_lm$r.squared
  
  cat(sprintf("R-squared for Group '%s': %.3f\n", group, r_squared))
}

for (group in group_names) {
  subset_data <- subset(df_sm_length, Group == group)
  lm_model <- lm(Mean ~ Length, data = subset_data)
  summary_lm <- summary(lm_model)
  r_squared <- summary_lm$r.squared
  p_value <- summary_lm$coefficients[2, 4]  # Extracting the p-value for the Length coefficient
  
  cat(sprintf("R-squared for Group '%s': %.3f\n", group, r_squared))
  cat(sprintf("P-value for Group '%s': %.3f\n", group, p_value))
}

# Set up the custom colors
custom_colors <- c("Control" = "#56B4E9", "Treatment" = "#E69F00")

# Create an empty plot with custom labels and no points
quartz()
plot(df_sm_length$Length, df_sm_length$Mean, type = "n",
     xlab = "Length [cm]",
     ylab = "Mean raw uncalibrated TDT values [-]",
     main = "",
     family = "Arial")

# Add points with colors according to the group
points(df_sm_length$Length[df_sm_length$Group == "Control"], 
       df_sm_length$Mean[df_sm_length$Group == "Control"], 
       col = custom_colors["Control"], pch = 19)

points(df_sm_length$Length[df_sm_length$Group == "Treatment"], 
       df_sm_length$Mean[df_sm_length$Group == "Treatment"], 
       col = custom_colors["Treatment"], pch = 19)

# Add linear regression lines and confidence intervals for each group
for (group in unique(df_sm_length$Group)) {
  subset_data <- subset(df_sm_length, Group == group)
  
  # Fit linear model
  fit <- lm(Mean ~ Length, data = subset_data)
  
  # Get prediction data
  pred <- predict(fit, interval = "confidence")
  
  # Add regression line
  lines(subset_data$Length, pred[, "fit"], col = custom_colors[group], lwd = 2)
  
  # Add confidence intervals
  polygon(c(subset_data$Length, rev(subset_data$Length)),
          c(pred[, "lwr"], rev(pred[, "upr"])),
          col = adjustcolor(custom_colors[group], alpha.f = 0.2),
          border = NA)
}

# Add a legend
legend("topright", legend = names(custom_colors), col = custom_colors, pch = 19, bty = "n")

print(d)

df_sc_length <- read.xlsx("") 

slength <- ggplot(df_sc_length, aes(x = length, y = sc, color = Group)) +
  geom_point(aes(color = Group)) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, aes(color = Group, fill = Group)) +
  theme_minimal()

p <- slength + labs(
  x = "Length [cm]",
  y = expression(Stomatal~conductance~"["~mmol/m^2~s~"]"))

custom_colors <- c("Control" = "#35978F", "Treatment" = "#BF812D")

p <- p + theme(
  text = element_text(family = "Arial")
) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) 

quartz()
print(p)

df_sc_length$Group <- factor(df_sc_length$Group)
group_names <- levels(df_sc_length$Group)

for (group in group_names) {
  subset_data <- subset(df_sc_length, Group == group)
  lm_model <- lm(sc ~ length, data = subset_data)
  summary_lm <- summary(lm_model)
  r_squared <- summary_lm$r.squared
  
  cat(sprintf("R-squared for Group '%s': %.3f\n", group, r_squared))
}

for (group in group_names) {
  subset_data <- subset(df_sc_length, Group == group)
  lm_model <- lm(sc ~ length, data = subset_data)
  summary_lm <- summary(lm_model)
  r_squared <- summary_lm$r.squared
  p_value <- summary_lm$coefficients[2, 4]  # Extracting the p-value for the Length coefficient
  
  cat(sprintf("R-squared for Group '%s': %.3f\n", group, r_squared))
  cat(sprintf("P-value for Group '%s': %.3f\n", group, p_value))
}


#lines(x, df_weather$sm_t_upper, type = "l", lty = 1, col = "black")
#lines(x, df_weather$sm_t_lower, type = "l", lty = 1, col = "black")
#lines(x, df_weather$sm_c_upper, type = "l", lty = 1, col = "red")
#lines(x, df_weather$sm_c_lower, type = "l", lty = 1, col = "red")
#polygon(c(x, rev(x)), c(df_weather$sm_c_upper, rev(df_weather$sm_c_lower)), col = adjustcolor("green", 0.5), , border = NA)
#polygon(c(x, rev(x)), c(df_weather$sm_t_upper, rev(df_weather$sm_t_lower)), col = adjustcolor("blue", 0.5), , border = NA)

# HPLC visualisation ----

# reference chroms
f_chroms <- read.xlsx("") 
df_chroms$time <- as.numeric(df_chroms$time)
df_chroms$W3_450 <- as.numeric(df_chroms$W3_450)
df_chroms$W3_550 <- as.numeric(df_chroms$W3_550)
df_chroms$W3_665 <- as.numeric(df_chroms$W3_665)
df_chroms$W7_450 <- as.numeric(df_chroms$W7_450)
df_chroms$W7_665 <- as.numeric(df_chroms$W7_665)
df_chroms$W13_450 <- as.numeric(df_chroms$W13_450)
df_chroms$W13_665 <- as.numeric(df_chroms$W13_665)
df_chroms_grouped <- df_chroms %>%
  group_by(chunk = rep(1:(n() %/% 10), each = 10, length.out = n())) 
df_chroms_averaged <- df_chroms_grouped %>%
  summarise(across(everything(), mean, na.rm = TRUE))
df_chroms_averaged <- select(df_chroms_averaged, -chunk)
write_xlsx(df_chroms_averaged, "")

#W3
plot3_1 <-ggplot(df_chroms_averaged, aes(x = time, y = W3_450)) +
  geom_line(color = "#feb24c", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))

plot3_2 <- ggplot(df_chroms_averaged, aes(x = time, y = W3_550)) +
  geom_line(color = "#7fcdbb", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))

plot3_3 <- ggplot(df_chroms_averaged, aes(x = time, y = W3_665)) +
  geom_line(color = "#2c7fb8", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))

#W7
plot7_1 <-ggplot(df_chroms_averaged, aes(x = time, y = W7_450)) +
  geom_line(color = "#feb24c", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))

plot7_3 <- ggplot(df_chroms_averaged, aes(x = time, y = W7_665)) +
  geom_line(color = "#2c7fb8", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))


#W13
plot13_1 <-ggplot(df_chroms_averaged, aes(x = time, y = W13_450)) +
  geom_line(color = "#feb24c", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))

plot13_3 <- ggplot(df_chroms_averaged, aes(x = time, y = W13_665)) +
  geom_line(color = "#2c7fb8", size=0.5) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(family = "Courier"))


quartz()
grid.arrange(plot3_1, plot3_2, plot3_3, ncol = 1)
quartz()
grid.arrange(plot7_1, plot7_3, ncol = 1)
quartz()
grid.arrange(plot13_1, plot13_3, ncol = 1)

# calibration curves
quartz()
x <- c(194.0, 116.4, 77.6, 38.8, 9.7, 3.9, 0.4)
y <- c(19935.6, 12493.4, 8627.8, 4332.0, 1078.9, 449.7, 46.3)

lm_model <- lm(y ~ x)
plot(x, y, main="Linear Regression", xlab="Concentration [µg/mL]", ylab="Absorption Units [mAU]",
     family="Courier", font.lab=2)
abline(lm_model, col="red")
coef_a <- coef(lm_model)[1]
coef_b <- coef(lm_model)[2]
r_squared <- summary(lm_model)$r.squared
r_squared_rounded <- round(r_squared, 5)
equation <- paste("Y =", round(coef_a, 2), "+", round(coef_b, 2), "X")
r_squared_text <- paste("R-squared =", r_squared_rounded)
text(50, 15000, equation, pos=4, family="Courier")
text(50, 14000, paste("R-squared =", r_squared_rounded), pos=4, family="Courier")

# Chlorophyll A standard calibration curve
x_chla <- c(194.0, 116.4, 77.6, 38.8, 9.7, 3.9, 0.4) 
y_chla <- c(19935.6, 12493.4, 8627.8, 4332.0, 1078.9, 449.7, 46.3)

# Chlorophyll B standard calibration curve
x_chlb <- c(225.0, 90.0, 45.0, 28.8, 11.2, 4.5, 0.45) 
y_chlb <- c(7094.7, 2804.7, 1432.9, 956.0, 349.1, 148.2, 14.6)

# Cyanidine chloride standard calibration curve
x_cyan <- c(233.0, 93.1, 46.6, 11.7, 4.7) 
y_cyan <- c(39568.4, 15485.1, 8019.0, 1520.74, 578.0)

# Pelargonidine chloride standard calibration curve
x_pel <- c(207.8, 135.1, 83.2, 41.5, 10.4, 4.2, 0.4) 
y_pel <- c(32257.0, 21477.0, 12763.0, 6282.0, 1198.0, 442.8, 54.4)

# Delphinidine chloride standard calibration curve
x_del <- c(183.5, 110.1, 73.4, 36.7, 9.2, 3.7, 0.4) 
y_del <- c(34689.0, 20875.0, 13982.0, 6641.0, 1627.0, 551.4, 17.34)

# Alpha carotene standard calibration curve
x_acar <- c(88.0, 52.80, 35.2, 17.6, 4.4, 1.76, 0.44) 
y_acar <- c(28041.3, 17146.3, 11469.5, 5837.7, 1448.9, 578.1, 146.2)

# Beta carotene standard calibration curve
x_bcar <- c(187.0, 112.4, 75.0, 37.0, 9.0, 3.75, 0.9) 
y_bcar <- c(36047.0, 20856.0, 14424.0, 7275.0, 2052.0, 806.5, 188.7)

# Xantophyll (Marigold) standard calibration curve
x_xant <- c(166.67, 100.0, 66.67, 33.33, 11.11, 3.7, 1.01) 
y_xant <- c(58848.1, 37539.1, 27240.8, 14083.4, 5361.8, 1784.9, 529.9)

# Lutein standard calibration curve
x_lut <- c(224.00, 134.4, 89.6, 44.80, 11.20, 4.48, 1.12) 
y_lut <- c(64957.0, 37712.0, 24629.4, 12660.0, 3236.7, 1358.4, 311.4)

# Neoxanthin standard calibration curve
x_neo <- c(203.33, 122.0, 81.3, 40.67, 10.17, 4.07, 1.02) 
y_neo <- c(52399.0, 33459.0, 22446.6, 11503.6, 2729.0, 1072.0, 245.4)

# Zeaxanthin standard calibration curve
x_zea <- c(139.33, 69.7, 46.44, 27.87, 6.97, 2.79, 0.7) 
y_zea <- c(54576.3, 29061.2, 19764.4, 12321.0, 2976.2, 1255.2, 272.7)

standard_name <- c('Chlorophyll A in Ace (665nm)', "Chlorophyll B in Ace (665nm)", "Cyanidin in MeOH (540nm)", "Pelargonidin in MeOH (526nm)",
                   'Delphinidine in MeOH (546nm)', "Alpha-Carotene in EtOAc (480nm)", "Beta-Carotene in Ace (478nm)", "Xantophyll in Ace (450nm)", 
                   "Lutein in Ace (472nm)", "Neoxanthin in Ace (450nm)", "Zeaxanthin in Ace (450nm)")
eluent_name <- c('Acetone', 'Acetone', 'MeOH', 'MeOH', 'MeOH', "Ethylacetate", "Acetone", "Acetone", "Acetone", "Acetone", "Acetone")


x_data_list <- list(x_chla, x_chlb, x_cyan, x_pel, x_del, x_acar, x_bcar, x_xant, x_lut, x_neo, x_zea)
y_data_list <- list(y_chla, y_chlb, y_cyan, y_pel, y_del, y_acar, y_bcar, y_xant, y_lut, y_neo, y_zea)

quartz()
par(mfrow=c(4, 3), mar=c(3, 2, 1, 1))
for (i in 1:11) {
  x <- x_data_list[[i]]
  y <- y_data_list[[i]]
  name <- standard_name[[i]]
  eluentname <- eluent_name[[i]]
  lm_model <- lm(y ~ x)
  plot(x, y, main=name, xlab="", ylab="", cex.axis=0.6, cey.axis = 0.6,
       family="Courier", font.lab=1, cex.main=0.8)
  abline(lm_model, col="red")
  coef_a <- coef(lm_model)[1]
  coef_b <- coef(lm_model)[2]
  r_squared <- summary(lm_model)$r.squared
  r_squared_rounded <- round(r_squared, 5)
  equation <- paste("Y =", round(coef_a, 2), "+", round(coef_b, 2), "X")
  axis(1, at=pretty(x), labels=NA, tick=FALSE)
  mtext(side=1, text=paste("R-squared =", r_squared_rounded), line=2.1, at=mean(range(x)), cex=0.4,family="Courier")
  mtext(side=1, text=equation, line=1.6, at=mean(range(x)), cex=0.4, family="Courier")
}

#spectra of standards and phaeophytin 

#standards
df_spectra_std <- read.xlsx("") 
quartz()
par(mfrow = c(4, 3), mar = c(2, 2, 1, 1))  # 4 rows and 3 columns for 12 plots, adjust margins
for (i in 2:ncol(df_spectra_std)) {
  plot(df_spectra_std$wvl, df_spectra_std[, i], type = "l", col = i - 1, xlab = "", ylab = "")
  title(main = colnames(df_spectra_std)[i], font.main = 4, cex.main = 1.2, family = "Courier")
}

#phaeophytins
df_spectra_phae <- read.xlsx("") 
quartz()
par(mfrow = c(1, 2), mar = c(2, 2, 1, 1))
for (i in 2:ncol(df_spectra_phae)) {
  plot(df_spectra_phae$wvl, df_spectra_phae[, i], type = "l", col = i - 1, xlab = "", ylab = "")
  title(main = colnames(df_spectra_phae)[i], font.main = 4, cex.main = 1.2, family = "Courier")
}
