#---------------------------------------------------------------------------------------#
# Date: July 10 2019
# Description code: Estimation of indirect social genetic effects for skin lesion count 
#                   in group-housed pigs by quantifying behavioral interactions
#---------------------------------------------------------------------------------------#
rm(list = ls())
setwd("~/Documents/Belcy_MSU_HCPP/2_Files_Matrix_Estimation_Varcomponents_Reml/Files_github_SGE/")

# 1. load library and files and matrix 
library(regress)
library(gwaR)
library(Matrix)
library(MASS)
library(dplyr)

# 1.1. G matrix (VanRaden,2008), R object: G, class: matrix, Dimension 1079
load("G_matrix.Rdata")
# 1.2. Phenotypes: Lesion score log(y+1) file, finisher post-mixing, R object: lcount.data
# class: data.frame
load("lesion_Count_Data.Rdata")

# 1.3 Matrices with time of interaction (seconds) in aggressive behaviors Reciprocal fight, 
#     Attacks, between pairs of individuals and uniform social interaction

# 1.3.1. Matrix in Reciprocal fight behavior, R object:Zc.RF.Mix, class: matrix
load("Zc_Undirectional_Mix.Rdata")
# 1.3.2. Matrix in Attacks beharior, R object:Zc.ATB.Mix, class: matrix
load("Zc_Directional_Mix.Rdata")
# 1.3.3. Matrix uniform social interaction, R object:Zc.Unif.Big, class: matrix
load("Zc_Uniform_undir_Mix.Rdata")

#---------------------------------------------------------------------------------------------
# 2. Fitting Direct Genetic Effect model (DGE) including interaction times 
#     as covariate for each trait
#---------------------------------------------------------------------------------------------
#-----
# 2.1. Fitting DGE model including interaction times in Reciprocal Fight as covariate
#    for each trait
#-------
# 2.1.1.  Trait: Anterior lesion count

gbfl.covRF<-gblup(rsp = "Finisher_Post_Front_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Front_ls + RFight.std +
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))
# 2.1.2. Trait: Central lesion count
gbml.covRF<-gblup(rsp = "Finisher_Post_Middle_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Middle_ls + RFight.std+
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))


# 2.1.3. Trait: Caudal lesion count 

gbrl.covRF<-gblup(rsp = "Finisher_Post_Rear_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Rear_ls + RFight.std+
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))
save(gbfl.covRF,gbml.covRF,gbrl.covRF, file = "DGE_covartime_RF.Rdata")
#------
# 2.2. Fitting DGE model including interaction times in Receives attacks as covariate
#    for each trait
#-----
# 2.2.1. Trait: Anterior lesion count 
gbfl.covAT<-gblup(rsp = "Finisher_Post_Front_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Front_ls + Receive_Time.std+
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))

# 2.2.2. Trait: Central lesion count

gbml.covAT<-gblup(rsp = "Finisher_Post_Middle_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Middle_ls + Receive_Time.std+
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))

# 2.2.3. Trait: Caudal lesion count 

gbrl.covAT<-gblup(rsp = "Finisher_Post_Rear_ls", data = lcount.data,
                  design = c(y ~ Sex + Rep + Finisher_Pre_Rear_ls + Receive_Time.std +
                               Finisher_Pre_Obs:Finisher_Post_Obs+Finisher_Wt,
                             ~ Finisher_Pen), G=G, pos= c(T,T,T))

save(gbfl.covAT, gbml.covAT,gbrl.covAT, file = "DGE_covartime_AT.Rdata")


#---------------------------------------------------------------------------------------------
# 3. Fitting Traditional Social Genetic model (TSGE) including interaction times 
#     as covariate for each trait
#---------------------------------------------------------------------------------------------

# load set of Functions: Implementation of the REML estimates of (co)variance components 
#                 through the EM algorithm and the asymptotic variances of the estimates

#source("Functions_Reml_EM.R")
source("~/Documents/Belcy_MSU_HCPP/2_Files_Matrix_Estimation_Varcomponents_Reml/Files_github/Functions_Reml_EM.R")

# 3.1. Standardaized interaction social matrix
Zunif.st<-st.Zcmat(Zc.Unif.Big) 

#-----
# 3.2. Fitting TSGE-RF model including interaction times in Reciprocal Fight as covariate
#    for each trait
#-------

# 3.2.1. Trait: Anterior lesion count 
flcov.tsgeRF<-igest(gbfl.covRF,Zunif.st,tol = 10^-5,k_iter = 300)
# 3.2.1.1. Variance-covariance matrix of REML estimates and standard error  
fltsgeRF.sd<-invImat(flcov.tsgeRF)
# 3.2.1.2. Heritability and standard error
hfltsgeRF<-varcompreml(flcov.tsgeRF)
# 3.2.1.3. Correlation between direct and social genetic effects and its standard error
fltsgeRF.rds<-varcovrml(flcov.tsgeRF)

save(fltsgeRF.sd,hfltsgeRF,fltsgeRF.rds,file = "TSGE_RF_Frontls_Geneparameter.Rdata")

# 3.2.2. Trait: Central lesion count 
mlcov.tsgeRF<-igest(gbml.covRF,Zunif.st,tol = 10^-5,k_iter = 300)
# 3.2.2.1. Variance-covariance matrix of REML estimates and standard error 
mltsgeRF.sd<-invImat(mlcov.tsgeRF)
# 3.2.2.2. Heritability and standard error
hmltsgeRF<-varcompreml(mlcov.tsgeRF)
# 3.2.2.3. Correlation between direct and social genetic effects and its standard error
mltsgeRF.rds<-varcovrml(mlcov.tsgeRF)

save(mltsgeRF.sd,hmltsgeRF,mltsgeRF.rds,file = "TSGE_RF_Middlels_Geneparameter.Rdata")

# 3.2.3. Trait: Caudal lesion count 
rlcov.tsgeRF<-igest(gbrl.covRF,Zunif.st,tol = 10^-5,k_iter = 300)
# 3.2.3.1. Variance-covariance matrix of REML estimates and standard error
rltsgeRF.sd<-invImat(rlcov.tsgeRF)
# 3.2.3.2 Heritability and standard error
hrltsgeRF<-varcompreml(rlcov.tsgeRF)
# 3.2.3.3. Correlation between direct and social genetic effects and its standard error
rltsge.RF.rds<-varcovrml(rlcov.tsgeRF)
save(rltsgeRF.sd,hrltsgeRF,rltsge.RF.rds,file = "TSGE_RF_Rearls_Geneparameter.Rdata")
save(flcov.tsgeRF, mlcov.tsgeRF,rlcov.tsgeRF, file = "TSGE_RF_covartime.Rdata")

#-----
# 3.3. Fitting TSGE-AT model including interaction times in Receives attacks as covariate
#    for each trait
#-------

# 3.3.1. Trait: Anterior lesion count 
flcov.tsgeAT<-igest(gbfl.covAT,Zunif.st,tol = 10^-5,k_iter = 300)
# 3.3.1.1. Variance-covariance matrix of REML estimates and standard error  
fltsgeAT.sd<-invImat(flcov.tsgeAT)
# 3.3.1.2. Heritability and standard error
hfltsgeAT<-varcompreml(flcov.tsgeAT)
# 3.3.1.3. Correlation between direct and social genetic effects and its standard error
fltsgeAT.rds<-varcovrml(flcov.tsgeAT)
save(fltsgeAT.sd,hfltsgeAT,fltsgeAT.rds,file = "TSGE_AT_Frontls_Geneparameter.Rdata")

# 3.3.2. Trait: Central lesion count 
mlcov.tsgeAT<-igest(gbml.covAT,Zunif.st,tol = 10^-5,k_iter = 560)
# 3.3.2.1. Variance-covariance matrix of REML estimates and standard error 
mltsgeAT.sd<-invImat(mlcov.tsgeAT)
# 3.3.2.2. Heritability and standard error
hmltsgeAT<-varcompreml(mlcov.tsgeAT)
# 3.3.2.3. Correlation between direct and social genetic effects and its standard error
mltsgeAT.rds<-varcovrml(mlcov.tsgeAT)
save(mltsgeAT.sd,hmltsgeAT,mltsgeAT.rds,file = "TSGE_AT_Middlels_Geneparameter.Rdata")

# 3.3.3. Trait: Caudal lesion count
rlcov.tsgeAT<-igest(gbrl.covAT,Zunif.st,tol = 10^-5,k_iter = 500)
# 3.3.3.1. Variance-covariance matrix of REML estimates and standard error
rltsgeAT.sd<-invImat(rlcov.tsgeAT)
# 3.3.3.2 Heritability and standard error
hrltsgeAT<-varcompreml(rlcov.tsgeAT)
# 3.3.3.3. Correlation between direct and social genetic effects and its standard error
rltsge.AT.rds<-varcovrml(rlcov.tsgeAT)
save(rltsgeAT.sd,hrltsgeAT,rltsge.AT.rds,file = "TSGE_AT_Rearls_Geneparameter.Rdata")

save(flcov.tsgeAT, mlcov.tsgeAT,rlcov.tsgeAT, file = "TSGE_AT_covartime.Rdata")


#---------------------------------------------------------------------------------------------
# 4. Fitting Interaction-base Social Genetic Effect model with Reciprocal Figth behavior
#     (ISGE-RF) including interaction time in Reciprocal Fight behavior as covariate  
#---------------------------------------------------------------------------------------------

# 4.1. Standardaized interaction social matrix
Zcst.RF<-st.Zcmat(Zc.RF.Mix) 

# 4.2.Trait: Anterior lesion count 
flcov.isgeRF<-igest(gbfl.covRF,Zcst.RF,tol = 10^-5,k_iter = 400)
# 4.2.1. Variance-covariance matrix of REML estimates and standard error
flisgeRF.sd<-invImat(flcov.isgeRF)
# 4.2.2. Heritability and standard error
hflisgeRF<-varcompreml(flcov.isgeRF)
# 4.2.3. Correlation between direct and social genetic effects and its standard error
flisgeRF.rds<-varcovrml(flcov.isgeRF)
save(flisgeRF.sd,hflisgeRF,flisgeRF.rds,file = "ISGE_RF_Frontls_Geneparameter.Rdata")

# 4.3. Trait: Central lesion count 
mlcov.isgeRF<-igest(gbml.covRF,Zcst.RF,tol = 10^-5,k_iter = 400)
# 4.3.1. Variance-covariance matrix of REML estimates and standard error 
mlisgeRF.sd<-invImat(mlcov.isgeRF)
# 4.3.2. Heritability and standard error
hmlisgeRF<-varcompreml(mlcov.isgeRF)
# 4.3.3. Correlation between direct and social genetic effects and its standard error
mlisgeRF.rds<-varcovrml(mlcov.isgeRF)
save(mlisgeRF.sd,hmlisgeRF,mlisgeRF.rds,file = "ISGE_RF_Middlels_Geneparameter.Rdata")

# 4.4. Trait: Caudal lesion count 
rlcov.isgeRF<-igest(gbrl.covRF,Zcst.RF,tol = 10^-5,k_iter = 400)
# 4.4.1. Variance-covariance matrix of REML estimates and standard error 
rlisgeRF.sd<-invImat(rlcov.isgeRF)
# 4.4.2. Heritability and standard error
hrlisgeRF<-varcompreml(rlcov.isgeRF)
# 4.4.3. Correlation between direct and social genetic effects and its standard error
rlisge.RF.rds<-varcovrml(rlcov.isgeRF)
save(rlisgeRF.sd,hrlisgeRF,rlisge.RF.rds,file = "ISGE_RF_Rearls_Geneparameter.Rdata")
save(flcov.isgeRF, mlcov.isgeRF,rlcov.isgeRF, file = "ISGE_RF_covartime.Rdata")

#---------------------------------------------------------------------------------------------
# 5. Fitting Interaction-based Social Genetic Effect model with Attack behavior
#     (ISGE-AT) 
#---------------------------------------------------------------------------------------------

# 5.1. Standardaized interaction social matrix
Zcst.AT<-st.Zcmat(Zc.ATB.Mix) 

# 5.2.Trait: Anterior lesion count 
flcov.isgeAT<-igest(gbfl.covAT,Zcst.AT,tol = 10^-5,k_iter = 300)
# 5.2.1. Variance-covariance matrix of REML estimates and standard error
flisgeAT.sd<-invImat(flcov.isgeAT)
# 5.2.2. Heritability and standard error
hflisgeAT<-varcompreml(flcov.isgeAT)
# 5.2.3. Correlation between direct and social genetic effects and its standard error
flisgeAT.rds<-varcovrml(flcov.isgeAT)
save(flisgeAT.sd,hflisgeAT,flisgeAT.rds,file = "ISGE_AT_Frontls_Geneparameter.Rdata")

# 5.3.Trait: Central lesion count 
mlcov.isgeAT<-igest(gbml.covAT,Zcst.AT,tol = 10^-5,k_iter = 300)
# 5.3.1. Variance-covariance matrix of REML estimates and standard error 
mlisgeAT.sd<-invImat(mlcov.isgeAT)
# 5.3.2. Heritability and standard error
hmlisgeAT<-varcompreml(mlcov.isgeAT)
# 5.3.3. Correlation between direct and social genetic effects and its standard error
mlisgeAT.rds<-varcovrml(mlcov.isgeAT)
save(mlisgeAT.sd,hmlisgeAT,mlisgeAT.rds,file = "ISGE_AT_Middlels_Geneparameter.Rdata")

# 5.4. Trait: Caudal lesion count 
rlcov.isgeAT<-igest(gbrl.covAT,Zcst.AT,tol = 10^-5,k_iter = 300)
# 5.4.1. Variance-covariance matrix of REML estimates and standard error 
rlisgeAT.sd<-invImat(rlcov.isgeAT)
# 5.4.2. Heritability and standard error
hrlisgeAT<-varcompreml(rlcov.isgeAT)
# 5.4.3. Correlation between direct and social genetic effects and its standard error
rlisge.AT.rds<-varcovrml(rlcov.isgeAT)
save(flcov.isgeAT, mlcov.isgeAT,rlcov.isgeAT, file = "ISGE_AT_covartime.Rdata")

#---------------------------------------------------------------------------------------------
# 6. Likelihood Ratio Test (LRT) for Social genetic variance and 
#    covariance of direct and social effect in Social Genetic Effects models
#---------------------------------------------------------------------------------------------

# 6.1. Likelihood Ratio Test DGE-RF model vs TSGE-RF model

# 6.1.1. Trait: Anterior lesion count 
Lrt(gbfl.covRF,flcov.tsgeRF,2)

# 6.1.2. Trait: Central lesion count 
Lrt(gbml.covRF,mlcov.tsgeRF,2)

# 6.1.3. Trait: Caudal lesion count
Lrt(gbrl.covRF,rlcov.tsgeRF,2)

#---------------------------------------------------------------------------------------------
# 6.2. Likelihood Ratio Test DGE-AT model vs TSGE-AT model

# 6.2.1. Trait: Anterior lesion count 
Lrt(gbfl.covAT,flcov.tsgeAT,2)

# 6.2.2. Trait: Central lesion count 
Lrt(gbml.covAT,mlcov.tsgeAT,2)

# 6.2.3. Trait: Caudal lesion count
Lrt(gbrl.covAT,rlcov.tsgeAT,2)

#---------------------------------------------------------------------------------------------
# 6.3. Likelihood Ratio Test DGE-RF model vs ISGE-RF model

# 6.3.1. Trait: Anterior lesion count 
Lrt(gbfl.covRF,flcov.isgeRF,2)

# 6.3.2. Trait: Central lesion count 
Lrt(gbml.covRF,mlcov.isgeRF,2)

# 6.3.3. Trait: Caudal lesion count
Lrt(gbrl.covRF,rlcov.isgeRF,2)

#---------------------------------------------------------------------------------------------
# 6.4. Likelihood Ratio Test DGE-AT model vs ISGE-AT model

# 6.4.1. Trait: Anterior lesion count 
Lrt(gbfl.covAT,flcov.isgeAT,2)

# 6.4.2. Trait: Central lesion count 
Lrt(gbml.covAT,mlcov.isgeAT,2)

# 6.4.3. Trait: Caudal lesion count
Lrt(gbrl.covAT,rlcov.isgeAT,2)




