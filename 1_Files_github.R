#---------------------------------------------------------------------------------#
# Date: April 22 2019
# Description code: Files,matrices used in the estimation of variance components
#                   of Social Genetics Efects (SGE) on lesions counts in diferent
#                   regions of the body in pigs
#---------------------------------------------------------------------------------#
rm(list = ls())
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

# 1.3 Matrices with time of interaction (seconds) in aggressive behaviors Reciprocal fight, Attacks,
#     between pairs of individuals and uniform social interaction

# 1.3.1. Matrix in Reciprocal fight behavior, R object:Zc.RF.Mix, class: matrix
load("Zc_Undirectional_Mix.Rdata")
# 1.3.2. Matrix in Attacks beharior, R object:Zc.ATB.Mix, class: matrix
load("Zc_Directional_Mix.Rdata")
# 1.3.3. Matrix uniform social interaction, R object:Zc.Unif.Big, class: matrix
load("Zc_Uniform_undir_Mix.Rdata")

#---------------------------------------------------------------------------------------------
# 2. Estimating (co)variance components with Direct Genetic Effect model (DGE) for each trait
#---------------------------------------------------------------------------------------------

# 2.1. (co)variance components Anterior lesion count trait

gb.Anterior.lc<-gblup(rsp = "Finisher_Post_Front_ls", data = lcount.data,
                       design = c(y ~ Sex + Rep + Finisher_Pre_Front_ls +
                                    Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                  ~ Finisher_Pen), G=G, pos= c(T,T,T))
save(gb.Anterior.lc, file = "gbF.Rdata")
# 2.1. (co)variance components Central lesion count trait

gb.Central.lc<-gblup(rsp = "Finisher_Post_Middle_ls", data = lcount.data,
                        design = c(y ~ Sex + Rep + Finisher_Pre_Middle_ls +
                                     Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                   ~ Finisher_Pen), G=G, pos= c(T,T,T))
save(gb.Central.lc, file = "gbM.Rdata")
# 2.1. (co)variance components Caudal lesion count trait

gb.Caudal.lc<-gblup(rsp = "Finisher_Post_Rear_ls", data = lcount.data,
                      design = c(y ~ Sex + Rep + Finisher_Pre_Rear_ls +
                                   Finisher_Pre_Obs:Finisher_Post_Obs + Finisher_Wt,
                                 ~ Finisher_Pen), G=G, pos= c(T,T,T))

save(gb.Caudal.lc, file = "gbC.Rdata")
#---------------------------------------------------------------------------------------------
# 3. Estimating (co)variance components with Interaction-base Social Genetic Effect model (ISGE) for each trait
#---------------------------------------------------------------------------------------------
# load set of Functions
source("Functions_Reml_EM.R")

# 3.1. Estimating (co)variance components with ISGE model Reciprocal Figth behavior

# 1. Standardaized interaction social matrix
Zs.Rf<-st.Zcmat(Zc.RF.Mix) 

# 2.(co)variance components Anterior lesion count trait (200 iter)
Alc.reml<-igest(gb.Anterior.lc,Zs.Rf,tol = 10^-3,k_iter = 200)
# 2.1. Variance-covariance matrix and standard error for estimates 
vc.Alc.reml<-invImat(Alc.reml)
# 2.2. Heritability and standard error
h2.Alc.reml<-varcompreml(Alc.reml)
# 2.3. Correlation genetic direct and Indirect and standard error
r.Alc.reml<-varcovrml(Alc.reml)

save(Alc.reml,vc.Alc.reml,h2.Alc.reml,r.Alc.reml, file = "SgbA_RF.Rdata")

# 3.(co)variance components Central lesion count trait (200 iter)
Centlc.reml<-igest(gb.Central.lc,Zs.Rf,tol = 10^-4,k_iter = 200)
# 3.1. Variance-covariance matrix and standard error for estimates 
vc.Centlc.reml<-invImat(Centlc.reml)
# 3.2. Heritability and standard error
h2.Centlc.reml<-varcompreml(Centlc.reml)
# 3.3. Correlation genetic direct and Indirect and standard error
r.Centlc.reml<-varcovrml(Centlc.reml)
save(Centlc.reml,vc.Centlc.reml,h2.Centlc.reml,r.Centlc.reml, file = "SgbM_RF.Rdata")

# 4.(co)variance components Caudal lesion count trait (400 iter)
Caudlc.reml<-igest(gb.Caudal.lc,Zs.Rf,tol = 10^-4,k_iter = 400)
# 4.1. Variance-covariance matrix and standard error for estimates 
vc.Caudlc.reml<-invImat(Caudlc.reml)
# 4.2. Heritability and standard error
h2.Caudlc.reml<-varcompreml(Caudlc.reml)
# 4.3. Correlation genetic direct and Indirect and standard error
r.Caudlc.reml<-varcovrml(Caudlc.reml)
save(Caudlc.reml,vc.Caudlc.reml,h2.Caudlc.reml,r.Caudlc.reml, file = "SgbC_RF.Rdata")

# 3.2. Estimating (co)variance components with ISGE model Attack behavior

# 1. Standardaized interaction social matrix
Zs.Ab<-st.Zcmat(Zc.ATB.Mix) 

# 2.(co)variance components Anterior lesion count trait  (300 iter)
Ab.alc.reml<-igest(gb.Anterior.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 2.1. Variance-covariance matrix and standard error for estimates 
Ab.vc.alc.reml<-invImat(Ab.alc.reml)
# 2.2. Heritability and standard error
Ab.h2.alc.reml<-varcompreml(Ab.alc.reml)
# 2.3. Correlation genetic direct and Indirect and standard error
Ab.r.alc.reml<-varcovrml(Ab.alc.reml)
save(Ab.alc.reml,Ab.vc.alc.reml,Ab.h2.alc.reml, Ab.h2.alc.reml,Ab.r.alc.reml, file = "SgbA_AB.Rdata" )

# 3.(co)variance components Central lesion count trait (300 iter)
Ab.centlc.reml<-igest(gb.Central.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 3.1. Variance-covariance matrix and standard error for estimates 
Ab.vc.centlc.reml<-invImat(Ab.centlc.reml)
# 3.2. Heritability and standard error
Ab.h2.centlc.reml<-varcompreml(Ab.centlc.reml)
# 3.3. Correlation genetic direct and Indirect and standard error
Ab.r.centlc.reml<-varcovrml(Ab.centlc.reml)
save(Ab.centlc.reml, Ab.vc.centlc.reml, Ab.h2.centlc.reml, Ab.r.centlc.reml,file = "SgbM_AB.Rdata")

# 4.(co)variance components Caudal lesion count trait (300 iter)
Ab.caudlc.reml<-igest(gb.Caudal.lc,Zs.Ab,tol = 10^-4,k_iter = 300)
# 4.1. Variance-covariance matrix and standard error for estimates 
Ab.vc.caudlc.reml<-invImat(Ab.caudlc.reml)
# 4.2. Heritability and standard error
Ab.h2.caudlc.reml<-varcompreml(Ab.caudlc.reml)
# 4.3. Correlation genetic direct and Indirect and standard error
Ab.r.caudlc.reml<-varcovrml(Ab.caudlc.reml)
save(Ab.caudlc.reml,Ab.vc.caudlc.reml,Ab.h2.caudlc.reml,Ab.r.caudlc.reml,file = "SgbC_AB.Rdata")
#---------------------------------------------------------------------------------------------
# 4. Estimating (co)variance components with TSGE model for each trait
#---------------------------------------------------------------------------------------------
# 1. Standardaized interaction social matrix
ZsU<-st.Zcmat(Zc.Unif.Big) 

# 2.(co)variance components Anterior lesion count trait (300 iter)
u.Alc.reml<-igest(gb.Anterior.lc,ZsU,tol = 10^-4,k_iter = 300)
# 2.1. Variance-covariance matrix and standard error for estimates 
u.vc.Alc.reml<-invImat(u.Alc.reml)
# 2.2. Heritability and standard error
u.h2.Alc.reml<-varcompreml(u.Alc.reml)
# 2.3. Correlation genetic direct and Indirect and standard error
u.r.Alc.reml<-varcovrml(u.Alc.reml)
save(u.Alc.reml,u.vc.Alc.reml,u.h2.Alc.reml, u.r.Alc.reml, file = "SgbA_Unif.Rdata")

# 3.(co)variance components Central lesion count trait (300 iter)
u.centlc.reml<-igest(gb.Central.lc,ZsU,tol = 10^-4,k_iter = 300)
# 3.1. Variance-covariance matrix and standard error for estimates 
u.vc.centlc.reml<-invImat(u.centlc.reml)
# 3.2. Heritability and standard error
u.h2.centlc.reml<-varcompreml(u.centlc.reml)
# 3.3. Correlation genetic direct and Indirect and standard error
u.r.centlc.reml<-varcovrml(u.centlc.reml)
save(u.centlc.reml,u.vc.centlc.reml,u.h2.centlc.reml,u.r.centlc.reml, file = "SgbM_Unif.Rdata")

# 4.(co)variance components Caudal lesion count trait (300 iter)
u.caudlc.reml<-igest(gb.Caudal.lc,ZsU,tol = 10^-4,k_iter = 300)
# 4.1. Variance-covariance matrix and standard error for estimates 
u.vc.caudlc.reml<-invImat(u.caudlc.reml)
# 4.2. Heritability and standard error
u.h2.caudlc.reml<-varcompreml(u.caudlc.reml)
# 4.3. Correlation genetic direct and Indirect and standard error
u.r.caudlc.reml<-varcovrml(u.caudlc.reml)
save(u.caudlc.reml,u.vc.caudlc.reml,u.h2.caudlc.reml, u.r.caudlc.reml, file = "SgbC_Unif.Rdata")



