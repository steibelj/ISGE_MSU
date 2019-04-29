#---------------------------------------------------------------------------------------#
# Date: April 29 2019
# Description code: Estimation of indirect social genetic effects for skin lesion count 
#                   in group-housed pigs by quantifying behavioral interactions
#---------------------------------------------------------------------------------------#

rm(list = ls())

#set this folder to the location where all data and scripts are stored:
setwd("~/Documents/Belcy_MSU_HCPP/2_Files_Matrix_Estimation_Varcomponents_Reml/Files_github/")

# 1. load library and files and matrix 
library(regress)
library(gwaR)
library(Matrix)
library(MASS)

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
# 2. Fitting Direct Genetic Effect model (DGE) for each trait
#---------------------------------------------------------------------------------------------

# 2.1.  Trait: Anterior lesion count

gb.Anterior.lc<-gblup(rsp = "Finisher_Post_Front_ls", data = lcount.data,
                       design = c(y ~ Sex + Rep + Finisher_Pre_Front_ls +
                                    Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                  ~ Finisher_Pen), G=G, pos= c(T,T,T))
save(gb.Anterior.lc, file = "gbF.Rdata")
# 2.2. Trait: Central lesion count

gb.Central.lc<-gblup(rsp = "Finisher_Post_Middle_ls", data = lcount.data,
                        design = c(y ~ Sex + Rep + Finisher_Pre_Middle_ls +
                                     Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                   ~ Finisher_Pen), G=G, pos= c(T,T,T))
save(gb.Central.lc, file = "gbM.Rdata")
# 2.3. Trait: Caudal lesion count 

gb.Caudal.lc<-gblup(rsp = "Finisher_Post_Rear_ls", data = lcount.data,
                      design = c(y ~ Sex + Rep + Finisher_Pre_Rear_ls +
                                   Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                 ~ Finisher_Pen), G=G, pos= c(T,T,T))

save(gb.Caudal.lc, file = "gbC.Rdata")

#---------------------------------------------------------------------------------------------
# 3. Fitting Traditional Social Genetic model (TSGE) 
#---------------------------------------------------------------------------------------------

# load set of Functions: Implementation of the REML estimates of (co)variance components 
#                 through the EM algorithm and the asymptotic variances of the estimates

source("Functions_Reml_EM.R")

# 3.1. Standardaized interaction social matrix
ZsU<-st.Zcmat(Zc.Unif.Big) 

# 3.2.Trait: Anterior lesion count (300 iter)
u.Alc.reml<-igest(gb.Anterior.lc,ZsU,tol = 10^-4,k_iter = 300)
# 3.2.1. Variance-covariance matrix of REML estimates and standard error  
u.vc.Alc.reml<-invImat(u.Alc.reml)
# 3.2.2. Heritability and standard error
u.h2.Alc.reml<-varcompreml(u.Alc.reml)
# 3.2.3. Correlation between direct and social genetic effects and its standard error
u.r.Alc.reml<-varcovrml(u.Alc.reml)
save(u.Alc.reml,u.vc.Alc.reml,u.h2.Alc.reml, u.r.Alc.reml, file = "SgbA_Unif.Rdata")

# 3.3. Trait: Central lesion count (300 iter)
u.centlc.reml<-igest(gb.Central.lc,ZsU,tol = 10^-4,k_iter = 300)
# 3.3.1. Variance-covariance matrix of REML estimates and standard error 
u.vc.centlc.reml<-invImat(u.centlc.reml)
# 3.3.2. Heritability and standard error
u.h2.centlc.reml<-varcompreml(u.centlc.reml)
# 3.3.3. Correlation between direct and social genetic effects and its standard error
u.r.centlc.reml<-varcovrml(u.centlc.reml)
save(u.centlc.reml,u.vc.centlc.reml,u.h2.centlc.reml,u.r.centlc.reml, file = "SgbM_Unif.Rdata")

# 3.4. Trait: Caudal lesion count (300 iter)
u.caudlc.reml<-igest(gb.Caudal.lc,ZsU,tol = 10^-4,k_iter = 300)
# 3.4.1. Variance-covariance matrix of REML estimates and standard error
u.vc.caudlc.reml<-invImat(u.caudlc.reml)
# 3.4.2. Heritability and standard error
u.h2.caudlc.reml<-varcompreml(u.caudlc.reml)
# 3.4.3. Correlation between direct and social genetic effects and its standard error
u.r.caudlc.reml<-varcovrml(u.caudlc.reml)
save(u.caudlc.reml,u.vc.caudlc.reml,u.h2.caudlc.reml, u.r.caudlc.reml, file = "SgbC_Unif.Rdata")

#---------------------------------------------------------------------------------------------
# 4. Fitting Interaction-base Social Genetic Effect model with Reciprocal Figth behavior
#     (ISGE-RF) 
#---------------------------------------------------------------------------------------------

# 4.1. Standardaized interaction social matrix
Zs.Rf<-st.Zcmat(Zc.RF.Mix) 

# 4.2.Trait: Anterior lesion count (200 iter)
Alc.reml<-igest(gb.Anterior.lc,Zs.Rf,tol = 10^-3,k_iter = 200)
# 4.2.1. Variance-covariance matrix of REML estimates and standard error
vc.Alc.reml<-invImat(Alc.reml)
# 4.2.2. Heritability and standard error
h2.Alc.reml<-varcompreml(Alc.reml)
# 4.2.3. Correlation between direct and social genetic effects and its standard error
r.Alc.reml<-varcovrml(Alc.reml)
save(Alc.reml,vc.Alc.reml,h2.Alc.reml,r.Alc.reml, file = "SgbA_RF.Rdata")

# 4.3. Trait: Central lesion count (200 iter)
Centlc.reml<-igest(gb.Central.lc,Zs.Rf,tol = 10^-4,k_iter = 200)
# 4.3.1. Variance-covariance matrix of REML estimates and standard error 
vc.Centlc.reml<-invImat(Centlc.reml)
# 4.3.2. Heritability and standard error
h2.Centlc.reml<-varcompreml(Centlc.reml)
# 4.3.3. Correlation between direct and social genetic effects and its standard error
r.Centlc.reml<-varcovrml(Centlc.reml)
save(Centlc.reml,vc.Centlc.reml,h2.Centlc.reml,r.Centlc.reml, file = "SgbM_RF.Rdata")

# 4.4. Trait: Caudal lesion count (400 iter)
Caudlc.reml<-igest(gb.Caudal.lc,Zs.Rf,tol = 10^-4,k_iter = 400)
# 4.4.1. Variance-covariance matrix of REML estimates and standard error 
vc.Caudlc.reml<-invImat(Caudlc.reml)
# 4.4.2. Heritability and standard error
h2.Caudlc.reml<-varcompreml(Caudlc.reml)
# 4.4.3. Correlation between direct and social genetic effects and its standard error
r.Caudlc.reml<-varcovrml(Caudlc.reml)
save(Caudlc.reml,vc.Caudlc.reml,h2.Caudlc.reml,r.Caudlc.reml, file = "SgbC_RF.Rdata")

#---------------------------------------------------------------------------------------------
# 5. Fitting Interaction-based Social Genetic Effect model with Attack behavior
#     (ISGE-AT) 
#---------------------------------------------------------------------------------------------

# 5.1. Standardaized interaction social matrix
Zs.Ab<-st.Zcmat(Zc.ATB.Mix) 

# 5.2.Trait: Anterior lesion count (300 iter)
Ab.alc.reml<-igest(gb.Anterior.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 5.2.1. Variance-covariance matrix of REML estimates and standard error
Ab.vc.alc.reml<-invImat(Ab.alc.reml)
# 5.2.2. Heritability and standard error
Ab.h2.alc.reml<-varcompreml(Ab.alc.reml)
# 5.2.3. Correlation between direct and social genetic effects and its standard error
Ab.r.alc.reml<-varcovrml(Ab.alc.reml)
save(Ab.alc.reml,Ab.vc.alc.reml,Ab.h2.alc.reml, Ab.h2.alc.reml,Ab.r.alc.reml, file = "SgbA_AB.Rdata" )

# 5.3.Trait: Central lesion count (300 iter)
Ab.centlc.reml<-igest(gb.Central.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 5.3.1. Variance-covariance matrix of REML estimates and standard error 
Ab.vc.centlc.reml<-invImat(Ab.centlc.reml)
# 5.3.2. Heritability and standard error
Ab.h2.centlc.reml<-varcompreml(Ab.centlc.reml)
# 5.3.3. Correlation between direct and social genetic effects and its standard error
Ab.r.centlc.reml<-varcovrml(Ab.centlc.reml)
save(Ab.centlc.reml, Ab.vc.centlc.reml, Ab.h2.centlc.reml, Ab.r.centlc.reml,file = "SgbM_AB.Rdata")

# 5.4. Trait: Caudal lesion count (300 iter)
Ab.caudlc.reml<-igest(gb.Caudal.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 5.4.1. Variance-covariance matrix of REML estimates and standard error 
Ab.vc.caudlc.reml<-invImat(Ab.caudlc.reml)
# 5.4.2. Heritability and standard error
Ab.h2.caudlc.reml<-varcompreml(Ab.caudlc.reml)
# 5.4.3. Correlation between direct and social genetic effects and its standard error
Ab.r.caudlc.reml<-varcovrml(Ab.caudlc.reml)
save(Ab.caudlc.reml,Ab.vc.caudlc.reml,Ab.h2.caudlc.reml,Ab.r.caudlc.reml,file = "SgbC_AB.Rdata")

#---------------------------------------------------------------------------------------------
# 6. Likelihood Ratio Test (LRT) for Social genetic variance and 
#    covariance of direct and social effect in Social Genetic Effects models
#---------------------------------------------------------------------------------------------

# 6.1. Likelihood Ratio in TSGE models

# 6.1.1. Trait: Anterior lesion count 
lr1<-Lrt(gb.Anterior.lc,u.Alc.reml,2)

# 6.1.2. Trait: Central lesion count 
lr2<-Lrt(gb.Central.lc,u.centlc.reml,2)

# 6.1.3. Trait: Caudal lesion count
lr3<-Lrt(gb.Caudal.lc,u.caudlc.reml,2)

#---------------------------------------------------------------------------------------------
# 6.2. Likelihood Ratio Test in ISGE-RF model

# 6.2.1. Trait: Anterior lesion count 
lr1.rf<-Lrt(gb.Anterior.lc,Alc.reml,2)

# 6.2.2. Trait: Central lesion count 
lr2.rf<-Lrt(gb.Central.lc,Centlc.reml,2)

# 6.2.3. Trait: Caudal lesion count
lr3.rf<-Lrt(gb.Caudal.lc,Caudlc.reml,2)

#---------------------------------------------------------------------------------------------
# 6.3. Likelihood Ratio Test in  ISGE-AT model

# 6.3.1. Trait: Anterior lesion count 
lr1.at<-Lrt(gb.Anterior.lc,Ab.alc.reml,2)

# 6.3.2. Trait: Central lesion count 
lr2.at<-Lrt(gb.Central.lc,Ab.centlc.reml,2)

# 6.3.3. Trait: Caudal lesion count
lr3.at<-Lrt(gb.Caudal.lc,Ab.caudlc.reml,2)

