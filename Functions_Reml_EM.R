#----------------------------------------------------------------------------------------#
# Date: April 11 2019
# Description code: Set of functions for estimate REML-EM variance components
# 1. st.Zcmat
# 2. igest
# 3. invImat
# 4. varcompreml
# 5. sevrml
# 6. varcovrml
# 7. getxfromgblup 
# 8. Lrt
#----------------------------------------------------------------------------------------#

library(MASS)
library(regress)
library(gwaR)
library(Matrix)

#-----------------------------------------------------------
# 1. Function for standardized interaction social matrix
# input object: class matrix
# output: class matrix
#------------------------------------------------------------
st.Zcmat<-function(x){
  
  # 1. Check x object be matrix 
  if (class(x)!="matrix") {stop("x must be matrix")}
  
  # 2. Check x object is or is not symetric
  if(isSymmetric(x)!=T){x<-t(x)}else{x<-x}
  
  # 3. Function for standardized matrix
  y<-sweep(x,1,FUN = "/",STATS = sqrt(rowSums(x^2))) 
  
  # 4. output matrix standardized
  return(y)
}

#--------------------------------------------------------------------------------
# 2. Function igest:REML-EM Algorithm for Estimating of the variance components: 
# direct effect, competence effects, covariance direct effec-competence effects,
# pen effect and error variance, using genomic relationship matrix (G)
# in lesions Scores trait in pigs
#
# inputs: 1. dgo:gblup object, 2. Zc: interaction social matrix standardized
#         3. tol: tolerance,  4. k_iter: iterations number
# output: 1. object list, class sgblup, with estimated varcomponents, DBV, SBV,
#           matrices used in mixed model equations, loglikelihood.
#--------------------------------------------------------------------------------

igest <- function(dgo,Zc,tol=1E-6,k_iter=50) {
  
  if(class(dgo)!="gblup") stop("class of dgo object should be gblup")
  if(class(Zc)!="matrix") stop("class of Zc objec should be matrix")
  
  # 1. Construct Y vector and Desing matrix (X,Zd,Zpen,Zc)
  y<-dgo$model$y
  
  X1<-model.matrix(dgo$model$formula,dgo$model[all.vars(dgo$model$formula)])
  rankMatrix(X1)
  ef<-lm(y~X1-1)
  X<-X1[,!is.na(ef$coefficients)]
  if (rankMatrix(X)!=ncol(X)) stop("X should be full rank")
  
  Ig<-diag(ncol(dgo$model$G)) ## Zd
  rownames(Ig)<-colnames(Ig)<-rownames(dgo$model$G)
  
  rf<-update(dgo$model$Vformula,~.-G-In-1) ## Zpen
  Zpen<-model.matrix(rf,dgo$model[all.vars(rf)])
  rownames(Zpen)<-rownames(dgo$model$G)

  iZc<-intersect(rownames(Zc),rownames(dgo$model$G)) ## Zc
  Zc<-Zc[iZc,iZc]
  
  # 2. Construct matrix XtX, XtZpen,
  Xp<-t(X)
  XpX<-Xp%*%X
  XpZc<-Xp%*%Zc
  ZcX<-t(XpZc)
  XpZpen<-Xp%*%Zpen
  ZpX<-t(XpZpen)
  
  ZcZc<-t(Zc)%*%Zc
  ZcZp<-t(Zc)%*%Zpen
  ZpZp<-t(Zpen)%*%Zpen
  
  # 3.Contruct Ginv, initial value of variance covariance components and Sigma matrix
  G<-dgo$model$G
  invG<-solve(G)
  
  s_2_e0<-dgo$sigma[3]
  s_2_u0<-dgo$sigma[1]
  s_uc0<-0.0
  s_2_p0<-dgo$sigma[2] 
  pro_pen<-dgo$sigma[2]/sum(dgo$sigma)
  if (pro_pen<0.1){
    s_2_p0<-0.1*sum(dgo$sigma)  
  }
  s_2_c0<-s_2_p0
  
  
  # 4. Effect fixed Matrix
  FME<-rbind(cbind(XpX,Xp,XpZc,XpZpen),
             cbind(X,Ig,Zc,Zpen),
             cbind(ZcX,t(Zc),ZcZc,ZcZp),
             cbind(ZpX,t(Zpen),t(ZcZp),ZpZp))
  
  # * Check FME be symetric
  if (isSymmetric(FME)!=T) {
    stop("Effec fixed matrix should be symmetric")
  }
  
  # 5. Indicator Variables for Mix Model Equation and right hand side vector
  idx<-c(rep(T, ncol(X)), rep(F, ncol(Ig)), rep(F,ncol(Ig)+ncol(Zpen)))
  idzd<-c(rep(F, ncol(X)),rep(T, ncol(Ig)), rep(F,ncol(Ig)+ncol(Zpen)))
  idzc<-c(rep(F,ncol(X)), rep(F,ncol(Ig)), rep(T, ncol(ZcZc)), rep(F,ncol((Zpen))))
  idzp<-c(rep(F,ncol(X)+2*ncol(Ig)), rep(T,ncol(Zpen)))
  
  # 6. p= fixed effects, q= Number individuals, r= pen effects
  p<-ncol(X)
  q<-ncol(Ig)
  r<-ncol(Zpen)
  n<-length(y)
  
  # 8. Iteration cycle for Estimation variance components and Loglikelihood 
  k_iter<- k_iter
  llod0<-(-200000)
  
  iteration<-vector("numeric")
  s2e_hat<-vector("numeric")
  s2u_hat<-vector("numeric")
  s2c_hat<-vector("numeric")
  suc_hat<-vector("numeric")
  s2p_hat<-vector("numeric")
  llik<-vector("numeric")
  tol_hat<-vector("numeric")
  
  
  for (i in 1:k_iter) {
    
    # 9. G0 matrix = Sigma matrix
    Sm<-matrix(c(s_2_u0,s_uc0,s_uc0,s_2_c0),ncol = 2)
    invSM<-solve(Sm)
    d_invSm<-det(invSM)
    if (det(invSM)<=0.0) {
      stop("determinant of Sigma inverse matrix should be positive, 
           Sigma inverse matrix should be positive definite")
    }
    
    # 10.Random Matrix (G0)^-1
    rmdd<-invSM[1,1]*invG*s_2_e0
    rmdc1<-invSM[1,2]*invG*s_2_e0
    rmc<-invSM[2,2]*invG*s_2_e0
    
    
    # 11. Coefficient matrix of Mix Model Equation (MME) and right hand side vector(RHS)
    MME<-FME
    MME[idzd,idzd]<-MME[idzd,idzd]+rmdd
    MME[idzd,idzc]<-MME[idzd,idzc]+rmdc1
    MME[idzc,idzd]<-MME[idzc,idzd]+rmdc1
    MME[idzc,idzc]<-MME[idzc,idzc]+rmc
    MME[idzp,idzp]<-MME[idzp,idzp]+(diag(ncol(ZpZp))*(s_2_e0/s_2_p0))
    
    RHS<-rbind(Xp%*%y,
               Ig%*%y,
               t(Zc)%*%y,
               t(Zpen)%*%y)
    
    # 12. Inverse Coefficient Matrix of Mix Model Equation (CM)
    CM<-solve(MME)
    rownames(CM)<-colnames(CM)<-rownames(MME)
    
    # 13. MME Solution (Sol_MME) and etimates of parameters (bhat,uhat,Zc_hat,Zp_hat)
    Sol_MME<-CM%*%RHS
  
    bhat<-Sol_MME[idx,]
    uhat<-Sol_MME[idzd,]
    Zc_hat<-Sol_MME[idzc,]
    Zp_hat<-Sol_MME[idzp,]
    
    # 14.Predicted values (y_hat), Marginal Residuals (e_hat0), Conditional Residual (e_hat)
    y_hat<-X%*%bhat+uhat+Zc%*%Zc_hat+Zpen%*%Zp_hat
    e_hat<-y-y_hat
    e_hat0<-y-X%*%bhat   #*******
    V<-G*s_2_u0+(G%*%t(Zc)+Zc%*%G)*s_uc0+Zc%*%G%*%t(Zc)*s_2_c0+diag(1,length(y))*s_2_e0+Zpen%*%t(Zpen)*s_2_p0
    dt_V<-determinant(V) # output value=log|V|
    
    
    # 15.C_bb:fixed effects, C_uu:additive effects,C_cc: competence effects,
    #   C_uc:covariance additive-competence,C_pen: pen effects
    C_bb<-CM[idx,idx]
    C_uu<-CM[idzd,idzd]
    C_cc<-CM[idzc,idzc]
    C_uc<-CM[idzd,idzc]
    C_pen<-CM[idzp,idzp]
    
    # 16.Quadratic form: e'e, u'Au, c'Ac, u'Ac, p'p, loglikelihood
    ete<-sum(e_hat^2)
    ete0<-sum(e_hat0^2) #***
    uAu<-t(uhat)%*%invG%*%uhat
    cAc<-t(Zc_hat)%*%invG%*%Zc_hat
    uAc<-t(uhat)%*%invG%*%Zc_hat
    ptp<-t(Zp_hat)%*%Zp_hat
    
    
    #**** loglikelihood marginal error
    lLog0<- -0.5*((dt_V$modulus[1])+log(det(Xp%*%solve(V)%*%X))+t(e_hat0)%*%solve(V)%*%e_hat0) 
    
    # 17.Traces
    tr_uAu<-sum(diag(invG%*%C_uu))
    tr_cAc<-sum(diag(invG%*%C_cc))
    tr_uAc<-sum(diag(invG%*%C_uc))
    tr_ptp<-sum(diag(C_pen))
    
    # 18.Estimation Variance Compontents : sigma2_e, sigma2_u, sigma2_c, sigma_uc, sigma2_pen
    s<-tr_uAu*invSM[1,1]+ 2*tr_uAc*invSM[1,2]+tr_cAc*invSM[2,2]+tr_ptp/s_2_p0
    s_2_ek<-((ete+(p+2*q+r-s*s_2_e0)*s_2_e0)/n)
    s_2_uk<-((uAu+(tr_uAu*s_2_e0))/q)
    s_2_ck<-((cAc+(tr_cAc*s_2_e0))/q)
    s_uc_k<-((uAc+(tr_uAc*s_2_e0))/q)
    s_2_pk<-((ptp+(tr_ptp*s_2_e0))/r)
    
    
    s_2_e0<-as.numeric(s_2_ek)
    s_2_u0<-s_2_uk[1,1]
    s_2_c0<-s_2_ck[1,1]
    s_uc0<-s_uc_k[1,1]
    s_2_p0<-s_2_pk[1,1]
    
    # 19.Tolerance Criterion
    
    tol_c<-abs((-2*lLog0)-(-2*llod0))/abs(-2*lLog0)
    tol<-tol
    if (tol_c< tol) {
      print ("The Tolerance criterion has been met")
      break
    }
    
    # 20. Updating of variance components estimates values and loglikelihood
    llod0<-lLog0
    
    s2e_hat[i]<-as.numeric(s_2_ek)
    s2u_hat[i]<-s_2_uk
    s2c_hat[i]<-s_2_ck
    suc_hat[i]<-s_uc_k
    s2p_hat[i]<-s_2_pk
    llik[i]<-llod0
    
    tol_hat[i]<-tol_c
    iteration[i]<-i
    
    } # End iteration
  
  # 21.Result for iterations
  data_sigmas<-cbind(iteration,s2e_hat,s2u_hat,s2c_hat,suc_hat,s2p_hat, llik,tol_hat)
  colnames(data_sigmas)<-c("iter","sigma e", "sigma u","sigma c", "sigma uc", "sigma p","llik","tolerance")
  data_sigmas<-as.data.frame(data_sigmas)
  
  sigma.hat<-matrix(c(data_sigmas$`sigma u`[nrow(data_sigmas)],
                      data_sigmas$`sigma c`[nrow(data_sigmas)],
                      data_sigmas$`sigma uc`[nrow(data_sigmas)],
                      data_sigmas$`sigma p`[nrow(data_sigmas)],
                      data_sigmas$`sigma e`[nrow(data_sigmas)]),ncol=1)
  rownames(sigma.hat)<-c("sigma u","sigma c", "sigma uc", "sigma p","sigma e")
  colnames(sigma.hat)<-c("Estimate")
  
  
  llik0<-data_sigmas$llik[nrow(data_sigmas)]
  
  # 22. Mixed model equation and solution with variance component estimated
  # 1* Sigama matrix and inverse 
  Sigma1<-matrix(c(sigma.hat[1,1],sigma.hat[3,1],
                   sigma.hat[3,1],sigma.hat[2,1]),ncol = 2)
  invSigma1<-solve(Sigma1)
  
  # 2*.Random Matrix (G0)^-1
  rmdd1<-invSigma1[1,1]*invG*sigma.hat[5,1]
  rmdc11<-invSigma1[1,2]*invG*sigma.hat[5,1]
  rmc1<-invSigma1[2,2]*invG*sigma.hat[5,1]
  
  
  # 3*.Coefficient matrix of Mix Model Equation (CME) and right hand side vector(RHS)
  CME<-FME
  CME[idzd,idzd]<-CME[idzd,idzd]+rmdd1
  CME[idzd,idzc]<-CME[idzd,idzc]+rmdc11
  CME[idzc,idzd]<-CME[idzc,idzd]+rmdc11
  CME[idzc,idzc]<-CME[idzc,idzc]+rmc1
  CME[idzp,idzp]<-CME[idzp,idzp]+(diag(ncol(ZpZp))*(sigma.hat[5,1]/sigma.hat[4,1]))
  
  RHS1<-rbind(Xp%*%y,
              Ig%*%y,
              t(Zc)%*%y,
              t(Zpen)%*%y)
  
  # 4*.Inverse Coefficient Matrix of Mix Model Equation (Cmatinv)
  Cmatinv<-solve(CME)
  rownames(Cmatinv)<-colnames(Cmatinv)<-rownames(CME)
  
  # 5*. MME Solution (SolMME), etimates of fixed and random effects (bhat,uhat,Zc_hat,Zp_hat)
  # and residuals marginal (ehat) 
  SolMME<-Cmatinv%*%RHS1
  
  bhat1<-SolMME[idx,]
  uhat1<-SolMME[idzd,]
  Zc_hat1<-SolMME[idzc,]
  Zp_hat1<-SolMME[idzp,]
  
  bhat1<-as.matrix(bhat1)
  colnames(bhat1)<-"Estimates"
  
  uhat1<-as.matrix(uhat1)
  colnames(uhat1)<-"uhat"
  
  ac_hat1<-as.matrix(Zc_hat1)
  colnames(ac_hat1)<-"ac.hat"
  
  p_hat1<-as.matrix(Zp_hat1)
  colnames(p_hat1)<-"pen.hat"
  
  e_hat1<-y-X%*%bhat1
  colnames(e_hat1)<-"ehat"
  
  V1<-G*sigma.hat[1,1]+(G%*%t(Zc)+Zc%*%G)*sigma.hat[3,1]+
    Zc%*%G%*%t(Zc)*sigma.hat[2,1]+ diag(1,length(y))*sigma.hat[5,1]+Zpen%*%t(Zpen)*sigma.hat[4,1]
  Vinv<-solve(V1)
  dt_V1<-determinant(V1)
  
  ## loglikelihood
  loglike<--0.5*((dt_V1$modulus[1])+log(det(Xp%*%Vinv%*%X))+t(e_hat1)%*%Vinv%*%e_hat1)
  llik<-as.numeric(loglike)
  
  y<-y
  X<-X
  coefm<-rbind(bhat1,
               sigma.hat)
  iter<-as.integer(rownames(data_sigmas)[nrow(data_sigmas)])
  formula<-dgo$model$formula
  
  #-----------------------#
  # 23. Outputs object    #
  #-----------------------#
  
  result <- list(data_sigmas=data_sigmas, bhat=bhat1,uhat=uhat1,achat=ac_hat1,
                 penhat=p_hat1, ehat = e_hat1, sigma=sigma.hat,llik=llik,coefm=coefm,
                 cycle=iter,V=V1, Vinv=Vinv,X=X,y=y,Zc=Zc,Zpen=Zpen,
                 G=G,invG=invG,invSM=invSigma1,C.inv=Cmatinv,
                 model=list(y=y,formula=formula))
  result$name<-dgo$name
  class(result) <- "sgblup"
  return(result) 
  
}


#--------------------------------------------------------------------------------------
# 3. invImat function for construct Information matrix Reml (Im) second derivates 
#  of Loglikelihood with construction of P matrix
#    input: gb: object class sgblup
#    output:class list with: 1. Im: information matrix (I(θ)), 
#                            2. VC: inverse information matrix [I(θ)]^(-1)
#                            3.varhat: matrix with the variance components and standard errors
#--------------------------------------------------------------------------------------

invImat<-function(gb){
  if(class(gb)!="sgblup") stop("class object should be class gblup")
  
  # 1. Definition of Matrices 
  Vinv<-gb$Vinv
  X<-gb$X
  G<-gb$G
  Zc<-gb$Zc
  Zpen<-gb$Zpen
  
  # 2. Construction P matrix
  P<- Vinv-Vinv%*%X%*%solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv
  
  # 3. Diagonals Elements Information Matrix (Im) for Reml estimation
  Im11<-0.5*(sum(diag(P%*%G%*%P%*%P)))
  Im22<-0.5*(sum(diag(P%*%Zc%*%G%*%t(Zc)%*%P%*%Zc%*%G%*%t(Zc))))
  Im33<-0.5*(sum(diag(P%*%(G%*%t(Zc)+Zc%*%G)%*%P%*%(G%*%t(Zc)+Zc%*%G))))
  Im44<-0.5*(sum(diag(P%*%Zpen%*%t(Zpen)%*%P%*%Zpen%*%t(Zpen))))
  Im55<-0.5*(sum(diag(P%*%P)))
  
  # 4. Off- Diagonals Elements Im
  Im12<-0.5*(sum(diag(P%*%G%*%P%*%Zc%*%G%*%t(Zc))))
  Im13<-0.5*(sum(diag(P%*%G%*%P%*%(G%*%t(Zc)+Zc%*%G))))
  Im14<-0.5*(sum(diag(P%*%G%*%P%*%Zpen%*%t(Zpen))))
  Im15<-0.5*(sum(diag(P%*%G%*%P)))
  Im23<-0.5*(sum(diag(P%*%Zc%*%G%*%t(Zc)%*%P%*%(G%*%t(Zc)+Zc%*%G))))
  Im24<-0.5*(sum(diag(P%*%G%*%P%*%Zpen%*%t(Zpen))))
  Im25<-0.5*(sum(diag(P%*%Zc%*%G%*%t(Zc))))
  Im34<-0.5*(sum(diag(P%*%(G%*%t(Zc)+Zc%*%G)%*%P%*%Zpen%*%t(Zpen))))
  Im35<-0.5*(sum(diag(P%*%(G%*%t(Zc)+Zc%*%G)%*%P)))
  Im45<-0.5*(sum(diag(P%*%Zpen%*%t(Zpen)%*%P)))
  
  # 5. Construct Im
  Im<-as.matrix(rbind(cbind(Im11,Im12,Im13,Im14,Im15),
                      cbind(Im12,Im22,Im23,Im24,Im25),
                      cbind(Im13,Im23,Im33,Im34,Im35),
                      cbind(Im14,Im24,Im34,Im44,Im45),
                      cbind(Im15,Im25,Im35,Im45,Im55)))
  colnames(Im)<-rownames(Im)<-c("sigma u","sigma c", "sigma uc", "sigma p","sigma e")
  
  # 6. Inverse of Im is Asymptotic covariance matrix of REML estimates
  Iminv<-solve(Im)
  # 7. Diagonal Elements of Inverse of Im
  Iminvii<-diag(Iminv)
  
  # 8.Standard Error for variance components estimated with Reml-EM
  se<-as.matrix(sqrt(Iminvii))
  colnames(se)<-"Std Error"
  
  # 9. Matrix with variance components estimated and standard error
  estim<-cbind(gb$sigma,se) 
  
  # 10. Output 
  outp<-list(Im=Im,VC=Iminv,varhat=estim)
  
  class(outp) <- "sgblup"
  return(outp) 
}
#----------------------------------------------------------------------------------------------
# 4.varcompreml function to calculate heritability and standard error for REML-EM estimates
#   input: object: gb class sgblup
#   ouput: class data.frame with heritability and standard error
#----------------------------------------------------------------------------------------------
varcompreml<-function(gb){
  
  re1<-as.matrix(gb$coefm[grep("sigma",rownames(gb$coefm)),])
  colnames(re1)<-"Estimates"
  re<-as.matrix(re1[c(1,4,5),])
  colnames(re)<-"Estimates"
  h2<-re[,1]/sum(re[,1])
  
  se<-sevrml(gb)
  result<-data.frame(re,Std.Error=se$Se_sig,prop.var=h2,se=se$se.prp)
  return(result)
  
}

##----------------------------------------------------------------------------------------------
# 5. sevrml function calculated standard error for heritability estimated
#   input: object: gb class sgblup
#   ouput: class data.frame standard error of heritability
#----------------------------------------------------------------------------------------------
sevrml<-function(gb)
{
  prp1<-gb$sigma[c(1,4,5),]
  propvar<-prp1/sum(prp1)
  Itita<-invImat(gb)
  Se_sig<-Itita$varhat[c(1,4,5),2]
  VC1<-Itita$VC
  idvc<-intersect(names(prp1),rownames(VC1))
  VC<-VC1[idvc,idvc]
  vrs<-diag(VC)
  tv<-sum(vrs)
  cvs<-rowSums(VC)
  T1<-propvar^2
  T2<-vrs/prp1^2
  T3<-tv/sum(prp1)^2
  T4<-2*(cvs/(prp1*sum(prp1)))
  
  se.prp<-sqrt(T1*(T2+T3-T4))
  otpt<-data.frame(Se_sig,se.prp)
  
  return(otpt)
}

#----------------------------------------------------------------------------------------------
# 6. varcovrml, function to calculate genetic correlation direct-social
# for Estimates Reml-EM 
#   input: object: gb class sgblup
#   ouput: class list with genetic correlation direct-social and standard error
#----------------------------------------------------------------------------------------------
varcovrml<-function(gb){
  
  # Correlation direct-Indirect effects
  
  vrs<-gb$sigma
  g1<-vrs[1,1]
  g2<-vrs[2,1]
  gc<-vrs[3,1]
  gr<-gc/sqrt(g1*g2)
  
  # Variance of correlation and Standard error
  
  Itita<-invImat(gb)
  VC<-Itita$VC[1:3,1:3]
  vcg1<-VC[1,1]
  vcg2<-VC[2,2]
  vcgc<-VC[3,3]
  vcg12<-sum(VC[1,2])
  vcg1c<-sum(VC[1,3])
  vcg2c<-sum(VC[2,3])
  
  
  rgse<-(gr^2)*(vcg1/(4*(g1^2)) + vcg2/(4*g2^2) +vcgc/(gc^2) + vcg12/(2*g1*g2)
                - vcg1c/(g1*gc) - vcg2c/(g2*gc))
  
  se<-sqrt(rgse)
  
  return(list(correlation=gr, St.Error=se))
}

#-------------------------------------------------------------------------------------
# 7. getxfromgblup: function for to construct X matrix from object gblup
#   input: object: gb class sgblup
#   ouput: object class matrix
#--------------------------------------------------------------------------------------
getxfromgblup<-function(dgo){
  y<-dgo$model$y
  X1<-model.matrix(dgo$model$formula,dgo$model[all.vars(dgo$model$formula)])
  ef<-lm(y~X1-1)
  X<-X1[,!is.na(ef$coefficients)]
  if (rankMatrix(X)!=ncol(X)) stop("X should be full rank")
  return(X)
}


#--------------------------------------------------------------------------------------
# 8. Lrt: function for to calculate the Loglikelihood Ratio Test 
#         input: object: 1. gb0: class gblup, 2. gb1: class sgblup,
#                        3. df:degree freedom
#         ouput: object class list with 1. Lrtest=statstic LRT, 2. pval=p-value,
#                                       3. llik.full=loglikelihood full model, 
#                                       4. llik.null=loglikelihood null model)
#--------------------------------------------------------------------------------------

Lrt<-function(gb0,gb1,df){
  
  if(class(gb0)!="gblup") stop("class object should be class gblup")
  if(class(gb1)!="sgblup") stop("class of gb1 objec should be class gblup")
  Xg0<-getxfromgblup(gb0)
  
  if(all.equal(gb1$X,Xg0)!= T) {warning("X Matrix are not equals")}
  if(all.equal.character(colnames(Xg0),colnames(gb1$X))!=T) {warning("X Matrix have different fixed effects")}
  
  # 1.loglikelihood null model (gb direct effects)
  
  w0<- determinant(gb0$V) # log|V|
  g0<- determinant(t(Xg0)%*%gb0$Vinv%*%Xg0) # log|X'VinvX|
  e0<- t(gb0$ehat)%*%gb0$Vinv%*%gb0$ehat
  llkm0<- -0.5*((w0$modulus[1])+(g0$modulus[1])+e0)
  
  # 2. loglikelihood complete model (gb direct + social effect)
  
  w1<-determinant(gb1$V)
  a<-determinant(t(gb1$X)%*%gb1$Vinv%*%gb1$X)
  e<-gb1$ehat
  colnames(e)<-NULL
  e1<-t(e)%*%gb1$Vinv%*%(e)
  
  llkm1<- -0.5*((w1$modulus[1])+(a$modulus[1])+e1)
  
  # 3. Likelihood ratio test
  Lrtest<- 2*(llkm1-llkm0)
  colnames(Lrtest)<- c("LRT")
  
  # 4. p-value
  if(df==1){
    pvalue1<- pchisq(Lrtest[1], 1, lower.tail = F)
  }else{
    pvalue1<- (pchisq(Lrtest[1], 1, lower.tail = F)+pchisq(Lrtest[1], 2, lower.tail = F))/2}
  
  if (df==2) {
    pvalue1<- pchisq(Lrtest[1], 2, lower.tail = F)
  }
  
  resul<-list(Lrtest=Lrtest, pval=pvalue1,llik.full=llkm1, llik.null=llkm0)
  return(resul)
}