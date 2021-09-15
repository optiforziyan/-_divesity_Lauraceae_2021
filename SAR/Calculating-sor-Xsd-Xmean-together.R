################################################################################ .
# simultaneous autoregressive (SARerr) models                                  #
################################################################################ .
# (1)  Tuning SARerr model                                                       ----
rm(list=ls())

# Read the selected LSAs
discrete_LSAs <- read.csv("./input/grids_xy.csv", header = T)

# because the grid cells in Hainan and Taiwan are isolated, for spatial analysis we removed.
# the FID is 38, 222, 374
discrete_LSAs <- discrete_LSAs[-c(1,11,25),]

# Read the current climate and past climate change predictors
library(readxl)
SEM_new_data <- read_excel("./input/SAR_data.xlsx")

# Final data for analysis
final_data <- merge(SEM_new_data, discrete_LSAs, by.x = "FID", by.y = "FID")
names(final_data)
test_data <- final_data[,c("Lon","Lat","beta.sor","beta.sim","beta.sne", 
                           "TA_sd", "PA_sd", "TSN_sd", "PSN_sd", 
                           "MAT_sd","MAP_sd","TCV_sd","PCV_sd",
                           "TA_mean", "PA_mean", "TSN_mean", "PSN_mean", 
                           "MAT_mean","MAP_mean","TCV_mean","PCV_mean")]
names(test_data)
if(!require(DataExplorer)){
  install.packages("DataExplorer")
  require(DataExplorer)
}
library(tidyverse)
library(funModeling)
plot_histogram(test_data)

library("PerformanceAnalytics")
names(test_data)

if(!dir.exists("./output")){dir.create("./output")}
# chart.Correlation(final_data[,c(5,8:30)],method='pearson',pch='*',col='grey')
x <- cor(test_data[,c("TA_mean", "TA_sd", "PA_mean", "PA_sd", 
                      "TCV_mean","TCV_sd","PCV_mean","PCV_sd",
                      "MAT_mean","MAT_sd","MAP_mean","MAP_sd",
                      "TSN_mean","TSN_sd","PSN_mean","PSN_sd" )],
                  method = "pearson")
write.csv(x,"./output/correlation_pearson.csv")

# Z-score 
data <- as.data.frame(cbind(test_data[,1:2],
                            apply(test_data[,3:21], 2, function(x)scale(x))))
plot_histogram(data)

#coordinates
coords <- as.matrix(data[,1:2])

# lm model
respv <- "beta.sor"
# because of the multicollinearity,
# we excluded the MAP_mean 
expvp <- c("TA_mean", "TA_sd", "PA_mean", "PA_sd", 
           "TCV_mean","TCV_sd","PCV_mean","PCV_sd",
           "MAT_mean","MAT_sd",           "MAP_sd",
           "TSN_mean","TSN_sd","PSN_mean","PSN_sd" ) 
# lm function
lm_build <- function(data, respv, expv) {
  # calibrate liner regression model
  formula <- as.formula(paste0(respv, " ~ ", paste0(expv, collapse = " + ")))
  m <- lm(formula, data = data)
  #show(summary(m))
  return(m)
}
m <- lm_build(data, respv= respv, expv=expvp)
summary(m)

library(car)
vif(m)

# distance
d <- seq(150,6000,150)
sw <- c("W","B","S")

# 
df = merge(d,sw)
recorded.res <- matrix(NA,dim(df)[1],6)
colnames(recorded.res) <- c("Distance(km)","Weights",'minRSA','maxI', 'adjR2', 'AIC')

# tuning function
minRSA <- function(Model, ND, Wstyle=c("W","B","S")){
  #Load packages necessary for analyses below
  library(spdep) # SAR
  library(ncf)   # SAR
  library(tripack)
  library(maptools)
  library(foreign)
  library(sp)
  library(SparseM)
  library(boot)
  library(spatialreg)
  library(spData)
  library(sf)
  library(psych)
  library(performance) #calcalating AICc, however, for SARerr seems not work...
  # Define neighbourhood 
  nb1.5 <- dnearneigh(coords,0,ND,longlat=TRUE)
  # Spatial weights
  nb1.5.w <- nb2listw(nb1.5, glist=NULL, style=Wstyle, zero.policy=FALSE)
  # SARerr model
  sem.nb1.5.w <- errorsarlm(Model, listw=nb1.5.w,zero.policy=TRUE, tol.solve=1e-15)
  sem.nb1.5.w.res <- summary(sem.nb1.5.w, Nagelkerke=TRUE, Hausman=TRUE, adj.se=TRUE)
  # Correlograms
  cor.sem.nb1.5.w<-correlog(coords[,1],
                            coords[,2],
                            z=residuals(sem.nb1.5.w), 
                            na.rm=T, increment=1, resamp=1)
  # minRSA is a measure of the autocorrelation of model residuals in space 
  # and is the sum of the absolute value of Moran¡¯s I at the first 20 distance classes
  minRSA <- sum(abs(cor.sem.nb1.5.w$correlation[1:20]))
  maxI <- max(abs(cor.sem.nb1.5.w$correlation[1:20]))
  aic <- AIC(sem.nb1.5.w)
  # build and return results
  res <- list(sem.nb1.5.w.res, minRSA, sem.nb1.5.w.res[["NK"]],aic, maxI )
  names(res) <- c('SARerr','minRSA', 'adjR2', 'AIC','maxI')
  return(res)
}

# test for the function
minRSA(Model=m, ND =150, Wstyle = "W")

# loop
for (i in 1:dim(df)[1]){
  cat("tuning:",  i ," times","\n")
  aaa <- minRSA(Model = m, ND = df[i,1], Wstyle = df[i,2])
  recorded.res[i,] <- c(df[i,1], df[i,2],round(aaa$minRSA,3),round(aaa$maxI,3),round(aaa$adjR2,3),round(aaa$AIC,3))
}

recorded.res <-as.data.frame(recorded.res)
saveRDS(recorded.res,"./output/beta.sor_Xall_SARerr.Rds")

# ------------------------------------------------------------------------------ ----
# (2)  Reducing variables                                                        ----
# choosing minRSA and corresponding distance and weights
recorded.res <- as.data.frame(readRDS("./output/beta.sor_Xall_SARerr.Rds"))
recorded.res <- recorded.res[order(recorded.res$minRSA),]
head(recorded.res)


expvp <- c("TA_mean", "TA_sd", "PA_mean", "PA_sd", 
           "TCV_mean","TCV_sd","PCV_mean","PCV_sd",
           "MAT_mean","MAT_sd",           "MAP_sd",
           "TSN_mean","TSN_sd","PSN_mean","PSN_sd" ) 

aaa <- minRSA(Model=m, ND =recorded.res$`Distance(km)`[1], Wstyle = recorded.res$Weights[1])
aaa 
# step by step delete by P value
# 1st PA_mean 0.822071
# 2nd PSN_sd 0.7389903
# 3rd TCV_sd 0.7306621
# 4th PCV_mean 0.6870814
# 5th TSN_sd   0.6676129
# 6th MAP_sd  0.5460780

expvp <- c("TA_mean", "TA_sd", "PA_sd", 
           "TCV_mean","PCV_sd",
           "MAT_mean","MAT_sd",
           "TSN_mean","PSN_mean")
m0 <- lm_build (data, respv="beta.sor", expv= expvp)
minRSA(Model=m0, ND =recorded.res$`Distance(km)`[1], Wstyle = recorded.res$Weights[1])


nn <- choose(9,1)+choose(9,2)+choose(9,3)+choose(9,4)+choose(9,5)+
      choose(9,6)+choose(9,7)+choose(9,8)+choose(9,9) 

variable_test <-data.frame(variable_combination = NA,adjR2 =NA,AIC=NA)

for (j in 1:9) {
  cat(j, "variables begin","\n")
  comb <- combn(expvp,j)
  for (i in 1:dim(comb)[2]) {
  cat("tuning:",  i ," times","\n")
  m <- lm_build (data, respv ="beta.sor", expv = comb[,i])
  aaa <- minRSA(Model=m, ND =recorded.res$`Distance(km)`[1], Wstyle = recorded.res$Weights[1])
  var <- paste(comb[,i][1],comb[,i][2],comb[,i][3],comb[,i][4],comb[,i][5],comb[,i][6],comb[,i][7],comb[,i][8],comb[,i][9],sep='+')
  variable_test <- rbind(variable_test,c(var,round(aaa$adjR2,3),round(aaa$AIC,3)))
  }
}
variable_test <- variable_test[-1,]
write.csv(variable_test,"./output/beta.sor_Xsd_SARerr_AIC.csv",row.names = F)
varb_test <- read.csv("./output/beta.sor_Xsd_SARerr_AIC.csv")
varb_test <- varb_test[order(varb_test$AIC, decreasing = F),]
head(varb_test)

# the minimize AIC
minAIC <- varb_test$AIC[1]
# We then choose a combination of these variables within 2 units of the difference between the AIC and the minimum AIC, judging by adj R2
varb_test <- subset(varb_test, varb_test$AIC < (minAIC + 2),)
dim(varb_test) # remaining 26 rows

# ranking by adj R2
varb_test <- varb_test[order(varb_test$adjR2, decreasing = T),]
head(varb_test)
varb_test$variable_combination[1]

# final SARerr model
m <- lm_build (data, respv="beta.sor", 
               expv= c("MAT_mean","TSN_mean","PSN_mean",
                       "TA_mean", "TA_sd", "PA_sd", 
                       "TCV_mean","PCV_sd"))
summary(m)
r1 <- minRSA(m, 150, "S")
r1
#stepwise 

# ------------------------------------------------------------------------------ ----
# EOF
