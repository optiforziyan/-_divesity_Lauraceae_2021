################################################################################ .
# Structural equation models                                                   #
# Ziyan Liao - 2021-05-31                                                      #
# Email: ziyan.liao@wsl.ch                                                     #
################################################################################ .
# ------------------------------------------------------------------------------ ----  
# (0)   Help function for drawing correlation plot                               ----
chart.Correlation <- function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                                               "spearman"), ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = "pearson", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel, ...)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
}

# ------------------------------------------------------------------------------ ---- 
# (1)   Variable Definition, Data preparation, standardized                      ----

# Abiotic drivers  
# MAT:	Mean Annual Temperature [¡æ]
# MAP:	Mean Annual Precipitation [mm]
# T A:	Temperature Anomaly [¡æ]
# P A:	Precipitation Anomaly [mm]
# TSN:	Temperature Seasonality [¡æ]
# PSN:  Precipitation Seasonality [mm]

# Biotic drivers
# NVF:	Number of sub-vegetation Types
# PDS:	Proportion of Deciduous Species
# CWE:	Corrected Weighted Endemisim
# LSD:	Local Schoener¡¯s D

# Migration filter
# GMD:	Geographically Migration Distance [km]

# SR_LSA Species Richness of Local Species Assemblage

# Remove objects 
# rm(list = ls())

# loading 'lavaan' package
library(lavaan)

# data for runing SEMs, including: 
# (a) Complete
# (b) Tropical monsoon forest, rain forest 
# (c) Subtropical evergreen broad-leaved forest
# (d) Warm temperate deciduous broad-leaved forest data sets of Lauraceae
# (e) Alpine vegetation on the Qinghai-Tibet Plateau 

data <-read.csv('SEM.csv')

# The names of all variables (see Main text for more details)
names(data)

# Overview of the dataset
summary(data)

# Figure 3a:
data1 <- data
# Figure 3b: 
data2 <- data[data$Type=="Tropical rain forest",]
# Figure 3c: 
data3 <- data[data$Type=="subtropical evergreen broad-leaved forest",]
# Figure 3d: 
data4 <- data[data$Type=="Warm temperate deciduous broad-leaved forest",]
# Figure 3e:  
data5 <- data[data$Type=="Alpine vegetation on the Qinghai-Tibet Plateau",]

## variables with much larger variances can dominate the estimation algorithm's 
## search for the best-fitting parameter estimates, paying less attention to 
## important discrepancies among variables with less variance
## Standardized by dividing variances

data1 <-as.data.frame(apply(data1[,3:17], 2, function(x) x/sd(x[which(!is.na(x))])))
data2 <-as.data.frame(apply(data2[,3:17], 2, function(x) x/sd(x[which(!is.na(x))])))
data3 <-as.data.frame(apply(data3[,3:17], 2, function(x) x/sd(x[which(!is.na(x))])))
data4 <-as.data.frame(apply(data4[,3:17], 2, function(x) x/sd(x[which(!is.na(x))])))
data5 <-as.data.frame(apply(data5[,3:17], 2, function(x) x/sd(x[which(!is.na(x))])))

# ------------------------------------------------------------------------------ ----  
# (2a)  Developed SEM-a                                                          ----

# Overview of the standardized data
summary(data1)
names(data1)
# Before developing a SEM, let us see the correlation of the predictors 
library("PerformanceAnalytics")
chart.Correlation(data1[,c(2:15)],method='pearson',pch='*',col='grey')
# However, for structural equation models, 
# covariance is generally not a problem for building a SEM

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  SR_LSA ~ MAT + MAP + TSN + TA + PA + NVT +  PDS + CWE + LSD  + GMD 
  ¦Âsne ~ MAT + MAP + TSN + PA + NVT + PDS  + CWE + GMD + SR_LSA  
  ¦Âsim ~ MAT + MAP + TA + PA  + CWE + LSD  + GMD + SR_LSA
  # covariances
  ¦Âsne ~~ ¦Âsim
'
# Run the model
fit <- sem(model, data=data1, auto.var = TRUE)

# Summary information about the variables
varTable(fit)

# Show the results
summary(fit, standardized=TRUE, rsquare=TRUE)
 
# Fit Measures for a Latent Variable Model
fitmeasures(fit,c("chisq","df","pvalue","rmsea","gfi","aic","cfi"))

# Plot the results using 'lavaanPlot' package
library(lavaanPlot)
lavaanPlot(model=fit,node_options = list(shape="box", fontname="Helvetica"),
           edge_options=list(color="grey"),coefs=TRUE,covs=TRUE,stars=c("regress"))


# Assessing the need to consider some missing relationships 
# between variables through MI values

mf <- modificationindices(fit)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)

# ------------------------------------------------------------------------------ ----  
# (2b)  Variance partition for general zones                                     ----

library(vegan)
showvarparts(3, bg=2:4,Xnames=c("Abiotic",
                                "Biotic",
                                "Dispersal"))
# for SR_LSA
attach(data1)
SR_LSA_res_part1 <- varpart(SR_LSA,~ MAT + MAP + TA + PA + PSN + TSN,
                    ~ NVT + PDS + CWE + LSD, GMD, data=data1)
# for ¦Âsor
¦Âsor_res_part1 <- varpart(¦Âsor,~ MAT + MAP + TA + PA + PSN + TSN,
                            ~ NVT + PDS + CWE + LSD, GMD, data=data1)
# for ¦Âsne
¦Âsne_res_part1 <- varpart(¦Âsne,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data1)
# for ¦Âsim
¦Âsim_res_part1 <- varpart(¦Âsim,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data1)

plot(SR_LSA_res_part1 , Xnames=c("Abiotic",
                         "Biotic",
                         "Dispersal"), bg=2:4, digits=2, 
     main="SR_LSA_res_part1")
plot(¦Âsor_res_part1 , Xnames=c("Abiotic",
                                 "Biotic",
                                 "Dispersal"), bg=2:4, digits=2,
      main="¦Âsor_res_part1")
plot(¦Âsne_res_part1 , Xnames=c("Abiotic",
                         "Biotic",
                         "Dispersal"), bg=2:4, digits=2,
      main="¦Âsne_res_part1")
plot(¦Âsim_res_part1 , Xnames=c("Abiotic",
                         "Biotic",
                         "Dispersal"), bg=2:4, digits=2,
      main="¦Âsim_res_part1")

part_result_zone1 <- rbind(SR_LSA_res_part1[["part"]][["indfract"]],
                     ¦Âsor_res_part1[["part"]][["indfract"]],
                     ¦Âsne_res_part1[["part"]][["indfract"]],
                     ¦Âsim_res_part1[["part"]][["indfract"]]) 
write.csv(part_result_zone1,"part1_result.csv")
# ------------------------------------------------------------------------------ ----
# (3a)  Developed SEM-b                                                          ----

# Before developing a SEM, let us see the correlation of the predictors 
library("PerformanceAnalytics")
chart.Correlation(data2[,c(2,3,7)],method='pearson',pch='.',col='black')
# However, for structural equation models, 
# covariance is generally not a problem

#Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  SR_LSA ~ MAT  +  TA  + NVT
  ¦Âsne ~ MAT  + TA + NVT + SR_LSA
  ¦Âsim ~  MAT + TSN  + TA 
  # covariances
  ¦Âsne ~~ ¦Âsim
'
# Run the model
fit <- sem(model, data=data2, auto.var = TRUE)

# Summary information about the variables
varTable(fit)

# Show the results
summary(fit, standardized=TRUE, rsquare=TRUE)

# Fit Measures for a Latent Variable Model
fitmeasures(fit,c("chisq","df","pvalue","rmsea","gfi","aic","cfi"))

# Plot the results using 'lavaanPlot' package
library(lavaanPlot)
lavaanPlot(model=fit,node_options = list(shape="box", fontname="Helvetica"),
           edge_options=list(color="grey"),coefs=TRUE,covs=TRUE,stars=c("regress"))

# Assessing the need to consider some missing relationships 
# between variables through MI values

mf <- modificationindices(fit)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)
# ------------------------------------------------------------------------------ ----  
# (3b)  Variance partition for Tropical rain forest                              ----

library(vegan)
showvarparts(3, bg=2:4,Xnames=c("Abiotic",
                                "Biotic",
                                "Dispersal"))
# for SR_LSA
attach(data2)
SR_LSA_res_part2 <- varpart(SR_LSA,~ MAT + MAP + TA + PA + PSN + TSN,
                            ~ NVT + PDS + CWE + LSD, GMD, data=data2)
# for ¦Âsor
¦Âsor_res_part2 <- varpart(¦Âsor,~ MAT + MAP + TA + PA + PSN + TSN,
                          ~ NVT + PDS + CWE + LSD, GMD, data=data2)
# for ¦Âsne
¦Âsne_res_part2 <- varpart(¦Âsne,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data2)
# for ¦Âsim
¦Âsim_res_part2 <- varpart(¦Âsim,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data2)

plot(SR_LSA_res_part2 , Xnames=c("Abiotic",
                                 "Biotic",
                                 "Dispersal"), bg=2:4, digits=2, 
     main="SR_LSA_res_part2")
plot(¦Âsor_res_part2 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsor_res_part2")
plot(¦Âsne_res_part2 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsne_res_part2")
plot(¦Âsim_res_part2 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsim_res_part2")

part_result_zone2 <- rbind(SR_LSA_res_part2[["part"]][["indfract"]],
                           ¦Âsor_res_part2[["part"]][["indfract"]],
                           ¦Âsne_res_part2[["part"]][["indfract"]],
                           ¦Âsim_res_part2[["part"]][["indfract"]]) 
write.csv(part_result_zone2,"part2_result.csv")
# ------------------------------------------------------------------------------ ----
# (4a)  Developed SEM-c                                                          ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  SR_LSA ~ MAT + MAP + PSN + TA + PA + PDS +  LSD + GMD
  ¦Âsne ~ MAT + MAP  + PSN + NVT + PDS  + LSD + GMD 
  ¦Âsim ~   MAP + PSN + TA  + NVT  + PDS + CWE + GMD + SR_LSA 
  # covariances
  ¦Âsne ~~ ¦Âsim
'
# Run the model
fit <- sem(model, data=data3, auto.var = TRUE)

# Summary information about the variables
varTable(fit)

# Show the results
summary(fit, standardized=TRUE, rsquare=TRUE)

# Fit Measures for a Latent Variable Model
fitmeasures(fit,c("chisq","df","pvalue","rmsea","gfi","aic","cfi"))

# Plot the results using 'lavaanPlot' package
library(lavaanPlot)
lavaanPlot(model=fit,node_options = list(shape="box", fontname="Helvetica"),
           edge_options=list(color="grey"),coefs=TRUE,covs=TRUE,stars=c("regress"))


# Assessing the need to consider some missing relationships 
# between variables through MI values

mf <- modificationindices(fit)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)
# ------------------------------------------------------------------------------ ----  
# (4b)  Variance partition for Subtropical evergreen broad-leaved forest         ----

library(vegan)
showvarparts(3, bg=2:4,Xnames=c("Abiotic",
                                "Biotic",
                                "Dispersal"))
# for SR_LSA
attach(data3)
SR_LSA_res_part3 <- varpart(SR_LSA,~ MAT + MAP + TA + PA + PSN + TSN,
                            ~ NVT + PDS + CWE + LSD, GMD, data=data3)
# for ¦Âsor
¦Âsor_res_part3 <- varpart(¦Âsor,~ MAT + MAP + TA + PA + PSN + TSN,
                          ~ NVT + PDS + CWE + LSD, GMD, data=data3)
# for ¦Âsne
¦Âsne_res_part3 <- varpart(¦Âsne,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data3)
# for ¦Âsim
¦Âsim_res_part3 <- varpart(¦Âsim,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data3)

plot(SR_LSA_res_part3 , Xnames=c("Abiotic",
                                 "Biotic",
                                 "Dispersal"), bg=2:4, digits=2, 
     main="SR_LSA_res_part3")
plot(¦Âsor_res_part3 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsor_res_part3")
plot(¦Âsne_res_part3 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsne_res_part3")
plot(¦Âsim_res_part3 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsim_res_part3")

part_result_zone3 <- rbind(SR_LSA_res_part3[["part"]][["indfract"]],
                           ¦Âsor_res_part3[["part"]][["indfract"]],
                           ¦Âsne_res_part3[["part"]][["indfract"]],
                           ¦Âsim_res_part3[["part"]][["indfract"]]) 
write.csv(part_result_zone3,"part3_result.csv")
# ------------------------------------------------------------------------------ ----
# (5a)  Developed SEM-d                                                          ----

# Structural (Here, the final structural equations were showed)
model <- '
 # regressions
  SR_LSA ~  PSN + TA + PDS + LSD 
  ¦Âsne ~ MAT + MAP + PSN + PA + LSD + GMD + SR_LSA 
  ¦Âsim ~ MAP + TSN +PSN + TA + PA + PDS  + LSD + GMD 
  # covariances
  ¦Âsne ~~ ¦Âsim
'
# Run the model
fit <- sem(model, data=data4, auto.var = TRUE)

# Summary information about the variables
varTable(fit)

# Show the results
summary(fit, standardized=TRUE, rsquare=TRUE)

# Fit Measures for a Latent Variable Model
fitmeasures(fit,c("chisq","df","pvalue","rmsea","gfi","aic","cfi"))

# Plot the results using 'lavaanPlot' package
library(lavaanPlot)
lavaanPlot(model=fit,node_options = list(shape="box", fontname="Helvetica"),
           edge_options=list(color="grey"),coefs=TRUE,covs=TRUE,stars=c("regress"))


# Assessing the need to consider some missing relationships 
# between variables through MI values

mf <- modificationindices(fit)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)


# Before developing a SEM, let us see the correlation of the predictors 
library("PerformanceAnalytics")
chart.Correlation(data4[,c(2:3,5:10,11:15)],method='pearson',pch='.',col='grey')
chart.Correlation(data4[,c(13,2:3)],method='pearson',pch='.',col='grey')

# ------------------------------------------------------------------------------ ----  
# (5b)  Variance partition for Warm temperate deciduous broad-leaved forest      ----

library(vegan)
showvarparts(3, bg=2:4,Xnames=c("Abiotic",
                                "Biotic",
                                "Dispersal"))
# for SR_LSA
attach(data4)
SR_LSA_res_part4 <- varpart(SR_LSA,~ MAT + MAP + TA + PA + PSN + TSN,
                            ~ NVT + PDS + CWE + LSD, GMD, data=data4)
# for ¦Âsor
¦Âsor_res_part4 <- varpart(¦Âsor,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data4)
# for ¦Âsne
¦Âsne_res_part4 <- varpart(¦Âsne,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data4)
# for ¦Âsim
¦Âsim_res_part4 <- varpart(¦Âsim,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data4)

plot(SR_LSA_res_part4 , Xnames=c("Abiotic",
                                 "Biotic",
                                 "Dispersal"), bg=2:4, digits=2, 
     main="SR_LSA_res_part4")
plot(¦Âsor_res_part4, Xnames=c("Abiotic",
                              "Biotic",
                              "Dispersal"), bg=2:4, digits=2,
      main="¦Âsor_res_part4")
plot(¦Âsne_res_part4 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsne_res_part4")
plot(¦Âsim_res_part4 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsim_res_part4")
part_result_zone4 <- rbind(SR_LSA_res_part4[["part"]][["indfract"]],
                           ¦Âsor_res_part4[["part"]][["indfract"]],
                           ¦Âsne_res_part4[["part"]][["indfract"]],
                           ¦Âsim_res_part4[["part"]][["indfract"]]) 
write.csv(part_result_zone4,"part4_result.csv")
# ------------------------------------------------------------------------------ ----
# (6a)  Developed SEM-e                                                          ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  SR_LSA ~  MAT + MAP + TA + PA + PDS + CWE + GMD
  ¦Âsne ~   MAP + PSN   + SR_LSA 
  ¦Âsim ~  TA  + LSD + GMD + SR_LSA
  # covariances
  ¦Âsne ~~ ¦Âsim
'
# Run the model
fit <- sem(model, data=data5, auto.var = TRUE)

# Summary information about the variables
varTable(fit)

# Show the results
summary(fit, standardized=TRUE, rsquare=TRUE)

# Fit Measures for a Latent Variable Model
fitmeasures(fit,c("chisq","df","pvalue","rmsea","gfi","aic","cfi"))

# Plot the results using 'lavaanPlot' package
library(lavaanPlot)
lavaanPlot(model=fit,node_options = list(shape="box", fontname="Helvetica"),
           edge_options=list(color="grey"),coefs=TRUE,covs=TRUE,stars=c("regress"))


# Assessing the need to consider some missing relationships 
# between variables through MI values

mf <- modificationindices(fit)
mf <- mf[order(mf$mi, decreasing = TRUE), ]
head(mf)
# ------------------------------------------------------------------------------ ----  
# (6b)  Variance partition for Alpine vegetation on the Qinghai-Tibet Plateau    ----

library(vegan)
showvarparts(3, bg=2:4,Xnames=c("Abiotic",
                                "Biotic",
                                "Dispersal"))
# for SR_LSA
attach(data5)
SR_LSA_res_part5 <- varpart(SR_LSA,~ MAT + MAP + TA + PA + PSN + TSN,
                            ~ NVT + PDS + CWE + LSD , GMD, data=data5)
# for ¦Âsor
¦Âsor_res_part5 <- varpart(¦Âsor,~ MAT + MAP + TA + PA + PSN + TSN,
                          ~ NVT + PDS + CWE + LSD , GMD, data=data5)
# for ¦Âsne
¦Âsne_res_part5 <- varpart(¦Âsne,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD , GMD, data=data5)
# for ¦Âsim
¦Âsim_res_part5 <- varpart(¦Âsim,~ MAT + MAP + TA + PA + PSN + TSN,
                           ~ NVT + PDS + CWE + LSD, GMD, data=data5)

plot(SR_LSA_res_part5 , Xnames=c("Abiotic",
                                 "Biotic",
                                 "Dispersal"), bg=2:4, digits=2, 
     main="SR_LSA_res_part4")
plot(¦Âsor_res_part5 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsor_res_part5")
plot(¦Âsne_res_part5 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsne_res_part5")
plot(¦Âsim_res_part5 , Xnames=c("Abiotic",
                               "Biotic",
                               "Dispersal"), bg=2:4, digits=2,
      main="¦Âsim_res_part5")

part_result_zone5 <- rbind(SR_LSA_res_part5[["part"]][["indfract"]],
                           ¦Âsor_res_part5[["part"]][["indfract"]],
                           ¦Âsne_res_part5[["part"]][["indfract"]],
                           ¦Âsim_res_part5[["part"]][["indfract"]]) 
write.csv(part_result_zone5,"part5_result.csv")

# ------------------------------------------------------------------------------ ----

# EOF                                                                            ----
################################################################################ .



