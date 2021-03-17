################################################################################ .
# Structural equation models                                                   #
# ZY - 2021-03-13                                                              #
################################################################################ .
# ------------------------------------------------------------------------------ ----  
# (1)  Variable Definition,Data preparation, standardized                        ----
# MAT	Mean annual temperature [¡æ]
# MAP	Mean annual precipitation [mm]
# TA	Temperature anomaly [¡æ]
# PA	Precipitation anomaly [mm]
# TSN	Temperature seasonality [¡æ]
# PSN	Precipitation seasonality [mm]
# NVF	Number of sub-vegetation types
# PCP	Potential colonization probability 
# PDC	Proportion of deciduous composition
# CWE	Corrected-weighted endemisim
# LSD	Local Schoener¡¯s D
# PSA	Proportion of suitable area 
# GMD	Geographically migration distance 
# LSP Local species pool

# Remove objects 
rm(list = ls())

# loading 'lavaan' package
library(lavaan)

# data for runing SEMs, including: 
# (a) Complete
# (b) Subtropical evergreen broad-leaved forest
# (c) Tropical monsoon forest, rain forest 
# (d) Alpine vegetation on the Qinghai-Tibet Plateau 
# (e) Warm temperate deciduous broad-leaved forest data sets of Lauraceae
data <-read.csv('SEM.csv')

# The names of all variables (see Main text for more details)
names(data)

# Overview of the dataset
summary(data)

# Figure 3a:
data1 <- data
# Figure 3b: 
data2 <-data[data$Type=="subtropical evergreen broad-leaved forest",]
# Figure 3c: 
data3 <-data[data$Type=="Tropical rain forest",]
# Figure 3d: 
data4 <-data[data$Type=="Alpine vegetation on the Qinghai-Tibet Plateau",]
# Figure 3e:  
data5 <-data[data$Type=="Warm temperate deciduous broad-leaved forest",]

## variables with much larger variances can dominate the estimation algorithm's 
## search for the best-fitting parameter estimates, paying less attention to 
## important discrepancies among variables with less variance
## For instance, when selecting the data1
## Standardized by dividing variances
data1 <-as.data.frame(apply(data1[,3:19], 2, function(x) x/sd(x[which(!is.na(x))])))
data2 <-as.data.frame(apply(data2[,3:19], 2, function(x) x/sd(x[which(!is.na(x))])))
data3 <-as.data.frame(apply(data3[,3:19], 2, function(x) x/sd(x[which(!is.na(x))])))
data4 <-as.data.frame(apply(data4[,3:19], 2, function(x) x/sd(x[which(!is.na(x))])))
data5 <-as.data.frame(apply(data5[,3:19], 2, function(x) x/sd(x[which(!is.na(x))])))

# Overview of the standardized data
summary(data1)

# Before developing a SEM, the covariance of the predictor variables 
# also needs to be considered using 'PerformanceAnalytics' package

library("PerformanceAnalytics")
# For instance, check for the abiotic variables
chart.Correlation(data1[,5:10],method='pearson')

# ------------------------------------------------------------------------------ ----  
# (2)  Developed a SEM                                                           ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  LSP ~ MAT + MAP + TA + PA + PDC + CWE + LSD + PSA + GMD
  beta.sne ~ MAT + MAP + TA + TSN + PSN + NVT + PCP + PDC + LSD + PSA + GMD + LSP
  beta.sim ~ TA + PA + PSN + PCP + PSA + GMD + LSP
  # covariances
  beta.sne ~~ beta.sim
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
# (3)  Developed b SEM                                                           ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  LSP ~ MAT + MAP + TA + PA  + PSN + NVT + PDC + CWE + LSD + PSA + GMD
  beta.sne ~ MAT + MAP + PA  + PSN + NVT + PCP + PDC  + LSD + GMD + LSP
  beta.sim ~ MAT + TA + TSN + PSN + NVT + PCP + PDC + CWE + GMD + LSP
  # covariances
  beta.sne ~~ beta.sim
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
# (4)  Developed c SEM                                                           ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  LSP ~  MAT + MAP + TSN + PSN + NVT + PCP
  beta.sne ~ MAT + MAP + TA + TSN + NVT + PCP  + LSP
  beta.sim ~  TA + TSN   + PSA 
  # covariances
  beta.sne ~~ beta.sim
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
# (5)  Developed d SEM                                                           ----

# Structural (Here, the final structural equations were showed)
model <- '
  # regressions
  LSP ~  MAT + PA + PCP  + CWE  + GMD
  beta.sne ~  PSN  + CWE  + PSA + LSP
  beta.sim ~ MAP + TSN + PCP  + LSD + PSA
  # covariances
  beta.sne ~~ beta.sim
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
# ------------------------------------------------------------------------------ ----
# (6)  Developed e SEM                                                           ----

# Structural (Here, the final structural equations were showed)
model <- '
 # regressions
  LSP ~  MAT + TA + PA + PSN + PCP + PDC + LSD + GMD
  beta.sne ~ MAT + MAP + PA  +  PSN + LSD +  LSP + GMD
  beta.sim ~ MAP + TA + PA + TSN +PSN + PDC + LSD +GMD
  # covariances
  beta.sne ~~ beta.sim
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
# EOF                                                                            ----
################################################################################ .
# ------------------------------------------------------------------------------ ----
