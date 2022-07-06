## =====================================
## 1, data processing, calculation and plotting ####
## =====================================
library(ggpmisc)
library(ggplot2)
library(tidyverse)
library(Rage)
library(Rcompadre)
library(popbio)
library(rlist)

rm(list=ls())

## ==============================
# 1. Get functions ######
## ==============================
# check if other place in the survival matrix is 0 or not 
# (not the sub-diagonal and nxn)
overall_survival_check <- function(mat){
  kunk <- as.matrix(as.data.frame(mat))
  # diag(kunk[-1,-ncol(kunk)]) <- 0
  n <- ncol(kunk)
  for (i in 1:n-1) {
    kunk[i+1, i] <-0
  }
  kunk[n,n] <- 0
  if(all(round(kunk,5)==0)){
    overall_survival <- 0
  }else{
    overall_survival <- 1
  }
  overall_survival
}

#### For damping paper only #####
# check if MatF only has one non-zero value in the first row
single_fertility <- function(mat){
  kunk <- as.matrix(as.data.frame(mat))[1,]
  if(length(which(round(kunk, 5) != 0))>1){
    single.fertility <- 0
  }else{
    single.fertility <- 1
  }
  single.fertility
}

# check sigularity of (I-U)
sing_check <- function(m) class(try(solve(m),silent=T))=="matrix"

#check if the first row in MatF has more than 1 constant values
constant_fertility <-function(mat){
  kunk <- as.matrix(as.data.frame(mat))[1,]
  non_zero <- c()
  for (col in 1:length(kunk)) {
    if(round(kunk[col],5)!=0){
      non_zero <- c(non_zero,kunk[col])
    }else{}
  }
  if(length(unique(non_zero))>1){
    constant.fertility <- 0
  }else{
    constant.fertility <-1
  }
  constant.fertility
}

matU <- matrix(c(0,1, 1.1, 0.9,
                 0.8, 0.3, 0.1,
                 0.2, 0.3, 0.2), ncol = 3, byrow = T)
# check for stasis
stasis_check <- function(matU){
  prop = diag(matU)/colSums(matU)
  sum = sum(prop[-length(prop)])
  stasis = sum/(ncol(matU)-1)
  return(stasis)
}

## partial life-cycle model
Partial_lc_model <- function(alpha, omega, beta, F_bar, Pj, Pa) {
  n = beta
  M = matrix(0, n, n)
  
  # lower diagonal has survivorship ratios
  for (i in 1:alpha) {
    M[i+1,i] <- Pj
  }
  for (i in (alpha+1):(beta-1)) {
    M[i+1,i] <- Pa
  }
  M[n,n] <- 0 ## what is the sx for open-end interval?
  
  # first row has net maternity contributions
  for(i in alpha:omega) {
    M[1,i] <- F_bar
  }
  M
}

damping_approx_agemodel <-function(age, phi){
  R0 = sum(phi)
  # phi = phi/R0
  # phi.1 <- (lx * mx)/R0
  
  Tc = weighted.mean(age, phi) ## don't add tau, has already changed age
  var = weighted.mean((age-Tc)^2, phi)
  sigma = sqrt(var)
  mu3 = weighted.mean((age-Tc)^3, phi) 
  mu4 = weighted.mean((age-Tc)^4, phi) 
  k3 = mu3
  k4 = mu4 - 3*sigma^4
  a = Tc
  b = -sigma^2/2
  c = k3/6
  d = -k4/24
  
  A = 1/a
  B = -b/a^3
  C = (2*b^2-a*c)/a^5
  D = (5*a*b*c - a^2*d-5*b^2)/a^7
  
  # r0.approx.p1 = log(R0)/Tc
  # r0.approx.p2 = var*(log(R0))^2/(2*Tc^3)
  # r0.approx =  r0.approx.p1 + r0.approx.p2
  # r1.approx = r0.approx - (2*pi^2*var/(Tc)^3)
  r0.approx =  A*log(R0) + B*log(R0)^2
  r0.approx1 = A*log(R0)
  r1.approx = r0.approx - 4*pi^2*B
  
  r0.approx_higher3 = A*log(R0) + B*log(R0)^2 + C*log(R0)^3
  r1.approx_higher3 = r0.approx_higher3 - 4*pi^2*B -12*pi^2*log(R0)*C
  
  r0.approx_higher4 = A*log(R0) + B*log(R0)^2 + C*log(R0)^3 + D*log(R0)^4
  r1.approx_higher4 = r0.approx_higher4 - 4*pi^2*B -12*pi^2*log(R0)*C-
    (24*pi^2*log(R0)^2 - 16*pi^4)*D
  
  # s1.approx = (2*pi/Tc + 2*pi*var*log(R0)/(Tc)^3)
  s1.approx = (2*pi*A + 4*pi*log(R0)*B)
  s1.approx_higher4 = (2*pi*A + 4*pi*log(R0)*B + (6*pi*log(R0)^2-8*pi^3)*C+
                        (8*pi*log(R0)^3 - 32*pi^3*log(R0)))
  damping.approx = r0.approx - r1.approx
  damping.approx_higher4 = r0.approx_higher4 - r1.approx_higher4
  damping.approx_higher3 = r0.approx_higher3 - r1.approx_higher3
  
  mylist <- list(
    "R0" = R0,
    "Tc" =  Tc,
    "var" = var,
    "sigma" = sigma,
    "mu3" = mu3,
    "mu4" = mu4,
    "k3" = k3,
    "k4" = k4,
    "r0.approx" = r0.approx,
    "r0.approx1" = r0.approx1,
    "r1.approx" = r1.approx,
    "r0.approx_higher4" = r0.approx_higher4,
    "r1.approx_higher4" = r1.approx_higher4,
    "r0.approx_higher3" = r0.approx_higher3,
    "r1.approx_higher3" = r1.approx_higher3,
    "s1.approx" = s1.approx,
    "s1.approx_higher4" = s1.approx_higher4,
    "damping.approx" = damping.approx,
    "damping.approx_higher4" = damping.approx_higher4,
    "damping.approx_higher3" = damping.approx_higher3
    )
  return(mylist)
}

damping_approx_stagemodel <-function(F_mat, P, N, tau){
  R0 <- F_mat[1,] %*% N ## from caswell (2001)
  R0 <- R0[1,1]
  
  Tc <- F_mat %*% N %*% N ## Y'(1)/P
  Tc <- Tc[1,1]/R0*tau # add tau
  T2 <-F_mat %*% (2 * N %*% N %*% P %*% N + N %*% N) ## (Y'(1)+ Y''(1))/P
  T2 <- T2[1,1]/R0*tau^2 ## times tau square
  sigma <- sqrt(T2 - Tc^2)
  var <- T2 - Tc^2

  Y1 <- N %*% P %*% N
  Y2 <- 2 * N %*% P %*% Y1
  Y3 <- 3 * N %*% P %*% Y2
  Y4 <- 4 * N %*% P %*% Y3
  
  Tc_check <- (F_mat %*% Y1)[1,1]/R0*tau + 1*tau # same as the simplified version before
  T2_check <- (F_mat %*% (Y1 +Y2))[1,1]/R0*tau^2 + 2*Tc*tau - 1*tau^2
  T3 <- (F_mat %*% (Y1 + 3*Y2 + Y3))[1,1]/R0*tau^3 + 3*T2*tau - 3*Tc*tau^2 + 1*tau^3
  #in the above equation, every term needs to have /R0*tau^2, but for T2, it already includes /R0^tau, so only add *tau.
  T4 <- (F_mat %*% (Y1 + 7*Y2 + 6*Y3 + Y4))[1,1]/R0*tau^4 + 4*T3*tau - 6*T2*tau^2 + 4*Tc*tau^3 - 1*tau^4
  
  mu3 <- T3 - 3*Tc*T2 + 2*Tc^3
  mu4 <- T4 - 4*Tc*T3 + 6*Tc^2*T2 - 3*Tc^4
  
  k3 = mu3
  k4 = mu4 - 3*sigma^4
  a = Tc
  b = -sigma^2/2
  c = k3/6
  d = -k4/24
  
  A = 1/a
  B = -b/a^3
  C = (2*b^2-a*c)/a^5
  D = (5*a*b*c - a^2*d-5*b^2)/a^7
  
  r0.approx =  A*log(R0) + B*log(R0)^2
  r0.approx1 = A*log(R0)
  r1.approx = r0.approx - 4*pi^2*B
  
  r0.approx_higher3 = A*log(R0) + B*log(R0)^2 + C*log(R0)^3
  r1.approx_higher3 = r0.approx_higher3 - 4*pi^2*B -12*pi^2*log(R0)*C
  
  r0.approx_higher4 = A*log(R0) + B*log(R0)^2 + C*log(R0)^3 + D*log(R0)^4
  r1.approx_higher4 = r0.approx_higher4 - 4*pi^2*B -12*pi^2*log(R0)*C-
    (24*pi^2*log(R0)^2 - 16*pi^4)*D
  
  # s1.approx = (2*pi/Tc + 2*pi*var*log(R0)/(Tc)^3)
  s1.approx = (2*pi*A + 4*pi*log(R0)*B)
  s1.approx_higher4 = (2*pi*A + 4*pi*log(R0)*B + (6*pi*log(R0)^2-8*pi^3)*C+
                         (8*pi*log(R0)^3 - 32*pi^3*log(R0)))
  damping.approx = r0.approx - r1.approx
  damping.approx_higher4 = r0.approx_higher4 - r1.approx_higher4
  damping.approx_higher3 = r0.approx_higher3 - r1.approx_higher3
  
  mylist <- list(
    "R0" = R0,
    "Tc" =  Tc,
    "var" = var,
    "sigma" = sigma,
    "mu3" = mu3,
    "mu4" = mu4,
    "k3" = k3,
    "k4" = k4,
    "r0.approx" = r0.approx,
    "r0.approx1" = r0.approx1,
    "r1.approx" = r1.approx,
    "r0.approx_higher4" = r0.approx_higher4,
    "r1.approx_higher4" = r1.approx_higher4,
    "r0.approx_higher3" = r0.approx_higher3,
    "r1.approx_higher3" = r1.approx_higher3,
    "s1.approx" = s1.approx,
    "s1.approx_higher4" = s1.approx_higher4,
    "damping.approx" = damping.approx,
    "damping.approx_higher4" = damping.approx_higher4,
    "damping.approx_higher3" = damping.approx_higher3
  )
  return(mylist)
}

## use simple stage model
Stage_model <- function(sample){
  n.case <- length(sample$matU)
  eigenvalue <- data.frame()
  for (i in 1: (n.case)) {
    
    # if a matrix has 6-month interval, it would be 0.5 since it is 0.5 years.
    # if(sample$ProjectionInterval[[i]] %in% c("NA",0)) next
    Usample <- unlist(sample$matU[i])
    Fsample <- unlist(sample$matF[i])
    tau <- sample$ProjectionInterval[[i]]
    P <- matrix(Usample, ncol = sqrt(length(Usample)), byrow = F)
    F_mat <- matrix(Fsample, ncol = sqrt(length(Fsample)), byrow = F)
    
    ns <- ncol(P) #number stages
    a1 <-  F_mat + P
    e <- eigen(a1)
    
    lam1 <- e$values[1]
    lam1.check <- e$values[Im(e$values) == 0][1]
    if (Mod(lam1) - Mod(lam1.check) > 10^-10){ 
      #choose >10^-10 instead of == 0, 
      #because sometimes they look exactly the same from the result that R returns, only have very small difference 
      stop("warning!!")
    }else{}
    
    lam2 <- e$values[2]
    
    r0.cal <- log(Mod(lam1))/tau # add tau
    r1.cal <- log(Mod(lam2))/tau
    s1.cal <- Arg(lam2)/tau
    damping.cal <- r0.cal - r1.cal
    
    # create a flag for more eigenvalues
    lam3 <- e$values[3]
    if (length(e$values) > 3){
      lam4 <- e$values[4]
    }else{
      lam4 <- NA
    }

    phi <- c()
    I <- diag(1, ns, ns)
    L <- I
    B1.sum <- diag(0, ns, ns)
    
    n.age <- c(1:50)*tau
    
    for (ages in 1:length(n.age)) { 
      junk <- F_mat %*% L
      junk <- junk[1,1]
      phi[ages] <- junk ## equation in the paper
      L <- P %*% L
      B1.junk <- n.age[ages] * F_mat %*% L
      B1.sum <- B1.sum + B1.junk
    }
    
    # ## another way to calculate phi
    # # calculate lx
    # lx <- ageSpecificSurv(P, startLife = 1, N = 150)
    # # calculate mx
    # mx <- ageSpecificRepro(P, F_mat, startLife = 1, N = 150)

    # check singularity here
    if(sing_check(I-P)!= TRUE) next
    N1 <- solve(I - P) %*% (I-P^(length(n.age))) ## geometric finite sum formula
    N <- solve(I - P)## inv(I - P) in matlab
    
    # N <- solve(I - P) %*% (I-P^(length(n.age))) ## geometric finite sum formula
    # N1 <- solve(I - P)## inv(I - P) in matlab
    
    R0 <- F_mat[1,] %*% N ## from caswell (2001)
    R01 <- F_mat[1,] %*% N1 ## from caswell (2001)
    R0 <- R0[1,1]
    R01 <- R01[1,1]
    
    phi <- phi/R0

    ## (1) revised matrix
    Tc <- F_mat %*% N %*% N ## Y'(1)/P
    Tc <- Tc[1,1]/R0*tau # add tau
    T2 <-F_mat %*% (2 * N %*% N %*% P %*% N + N %*% N) ## (Y'(1)+ Y''(1))/P
    T2 <- T2[1,1]/R0*tau^2 ## times tau or tau square?
    sigma <- sqrt(T2 - Tc^2)
    var <- T2 - Tc^2

    ## age at first reproduction by alpha, 
    ## age at last reproduction by omega, 
    ## and maximum longevity by beta
    mx.mat = F_mat[1,]
    age.mat = c(1:ncol(F_mat))*tau
    data.use <- as.data.frame(cbind(age.mat,mx.mat)) 
    if(sum(data.use$mx.mat) != 0){
      alpha = filter(data.use, mx.mat != 0)$age.mat[1]
      omega = tail(filter(data.use, mx.mat != 0)$age.mat, n=1)
      beta = tail(data.use$age.mat, n=1)
      
      alpha.position = match(alpha, data.use$age.mat)
      omega.position = match(omega, data.use$age.mat)
      beta.position = match(beta, data.use$age.mat)
    }else{
      alpha = NA
      omega = NA
      beta = tail(data.use$age, n=1)
      
      alpha.position = NA
      omega.position = NA
      beta.position = match(beta, data.use$age)
      
    }
    
    
    ## ---------------------
    ## approximation
    ## ---------------------
    r0.approx.p1 <- log(R0)/Tc
    r0.approx.p2 <- var*(log(R0)^2)/(2*Tc^3)
    r0.approx <- r0.approx.p1+r0.approx.p2
    r1.approx <- r0.approx - (2*pi^2*var/(Tc)^3)
    s1.approx <- (2*pi/Tc + 2*pi*var*log(R0)/(Tc)^3)
    damping.approx = r0.approx - r1.approx
    
    eigenvalue1 <- data.frame(
      # long.id = sample$id[i], 
      id = paste(sample$SpeciesAccepted[i], sample$Authors[i], sample$Journal[i],
                 sample$YearPublication[i], sample$MatrixStartYear[i], sample$MatrixEndYear[i],
                 sample$MatrixComposite[i], sample$MatrixTreatment[i], sep = "_"),
      MatrixID = sample$MatrixID[i],
      StudiedSex = sample$StudiedSex[[i]],
      Kingdom = sample$Kingdom[i],
      Phylum = sample$Phylum[i],
      Class = sample$Class[i],
      Order = sample$Order[i],
      Family = sample$Family[i],
      SpeciesAccepted = sample$SpeciesAccepted[i],
      Authors = sample$Authors[i],
      YearPublication = sample$YearPublication[i],
      Journal = sample$Journal[i],
      Altitude = sample$Altitude[i],
      CommonName = sample$CommonName[i],
      e = e$values[1],
      R0, Tc, sigma, var,
      r0.cal, r0.approx, 
      r1.cal, r1.approx,
      s1.cal, s1.approx, 
      damping.cal, damping.approx,
      alpha, omega, beta, alpha.position, omega.position, beta.position, tau,
      ProjectionInterval = sample$ProjectionInterval[i], OrganismType = sample$OrganismType[i], 
      lam1, lam2, lam3, lam4)
    
    eigenvalue <-  rbind(eigenvalue, eigenvalue1)
    print(paste(i,"/",n.case))  
  }
  eigenvalue
}

## use simple age model
Age_model <- function(sample2){
  eigenvalue.age <- data.frame()
  n.case <- length(sample2$matU)
  for (j in 1: (n.case)) {
    # if a matrix has 6-month interval, it would be 0.5 since it is 0.5 years.
    # if(sample2$ProjectionInterval[[j]] %in% c("NA",0)) next
    Usample <- unlist(sample2$matU[[j]])
    Fsample <- unlist(sample2$matF[[j]])
    tau <- sample2$ProjectionInterval[[j]]
    
    P <- matrix(Usample, ncol = sqrt(length(Usample)), byrow = F)
    F_mat <- matrix(Fsample, ncol = sqrt(length(Fsample)), byrow = F)
    
    ns <- ncol(P) #number stages
    a1 <-  F_mat + P
    e <- eigen(a1)
    
    # order(Mod(e$values),decreasing=T)[1]
    lam1 <- e$values[1]
    lam1.check <- e$values[Im(e$values) == 0][1]
    if (Mod(lam1) - Mod(lam1.check) > 10^-10){ 
      #choose >10^-10 instead of == 0, 
      #because sometimes they look exactly the same from the result that R returns, only have very small difference 
      stop("warning!!")
    }else{}
    
    lam2 <- e$values[2]
    
    r0.cal <- log(Mod(lam1))/tau # add tau
    r1.cal <- log(Mod(lam2))/tau
    s1.cal <- Arg(lam2)/tau
    damping.cal <- r0.cal - r1.cal
    
    # create a flag for more eigenvalues
    lam3 <- e$values[3]
    if (length(e$values) > 3){
      lam4 <- e$values[4]
    }else{
      lam4 <- NA
    }
    
    phi <- c()
    L.select = 1
    lx <- c()
    mx <- c()
    px <- c()
    for (ages in 1:ns) {
      if (ages < ns){
        F.select <- F_mat[1,ages]
        P.select <- P[ages+1, ages]
      }else{
        F.select <- F_mat[1,ages]
        P.select <- P[ages, ages]
      }
      L.select <- L.select * P.select
      lx[ages] <- L.select
      mx[ages] <- F.select
      px[ages] <- P.select
      phi[ages] <- F.select * L.select
    }
    age <- c(1:ns)*tau #add tau
    data.use <- as.data.frame(cbind(age,lx,mx,px)) 
    
    # ################################
    # ## for partial life-cycle model
    # ###############################
    dominate_e = Mod(e$vectors[1,])
    ## age at first reproduction by alpha, 
    ## age at last reproduction by omega, 
    ## and maximum longevity by beta
    if(sum(data.use$mx) != 0){
      alpha = filter(data.use, mx != 0)$age[1]
      omega = tail(filter(data.use, mx != 0)$age, n=1)
      beta = tail(data.use$age, n=1)
      
      alpha.position = match(alpha, data.use$age)
      omega.position = match(omega, data.use$age)
      beta.position = match(beta, data.use$age)
      
      # change all sx to px
      Fx = data.use$mx * data.use$px ## can change this to Pi
      F_bar_cal = sum(dominate_e[alpha.position:omega.position]*Fx[alpha.position:omega.position])/sum(dominate_e[alpha.position:omega.position])
      Pj_cal = sum(dominate_e[1:alpha.position]*data.use$px[1:alpha.position])/sum(dominate_e[1:alpha.position])## can change this to Pi
      if((alpha+1) < omega){
        Pa_cal = sum(dominate_e[(alpha.position+1):(omega.position-1)]*data.use$px[(alpha.position+1):(omega.position-1)])/sum(dominate_e[(alpha.position+1):(omega.position-1)])
      } else if ((alpha+1) == omega){ ## only one age has mx
        Pa_cal = sum(dominate_e[(alpha.position+1):(omega.position)]*data.use$px[(alpha.position+1):(omega.position)])/sum(dominate_e[(alpha.position+1):(omega.position)])
      }else {
        Pa_cal = sum(dominate_e[(alpha.position):(omega.position)]*data.use$px[(alpha.position):(omega.position)])/sum(dominate_e[(alpha.position):(omega.position)])
      }
    }else{
      alpha = NA
      omega = NA
      beta = tail(data.use$age, n=1)
      
      alpha.position = NA
      omega.position = NA
      beta.position = match(beta, data.use$age)
      
      F_bar_cal = NA
      Pj_cal = NA
      Pa_cal = NA
    }
    
    ## ---------------------
    ## approximation
    ## ---------------------
    R0 = sum(lx * mx)
    phi = phi/R0
    phi.1 <- (lx * mx)/R0
    
    Tc = sum(age*lx*mx)/R0 ## don't add tau, has already changed age
    var = sum((age-Tc)^2*lx*mx)/R0
    sigma = sqrt(var)
    r0.approx.p1 = log(R0)/Tc
    r0.approx.p2 = var*(log(R0))^2/(2*Tc^3)
    r0.approx =  r0.approx.p1 + r0.approx.p2
    r1.approx = r0.approx - (2*pi^2*var/(Tc)^3) 
    s1.approx = (2*pi/Tc + 2*pi*var*log(R0)/(Tc)^3)
    damping.approx = r0.approx - r1.approx
    # damping.approx.p1 = r0.approx.p1 - r1.approx
    
    eigenvalue1 <- data.frame(
      # long.id = sample2$id[[j]], 
      id = paste(sample2$SpeciesAccepted[[j]], sample2$Authors[[j]], sample2$Journal[[j]],
                 sample2$YearPublication[[j]], sample2$MatrixStartYear[[j]], sample2$MatrixEndYear[[j]],
                 sample2$MatrixComposite[[j]], sample2$MatrixTreatment[[j]], sep = "_"),
      MatrixID = sample2$MatrixID[[j]],
      StudiedSex = sample2$StudiedSex[[j]],
      Kingdom = sample2$Kingdom[[j]],
      Phylum = sample2$Phylum[[j]],
      Class = sample2$Class[[j]],
      Order = sample2$Order[[j]],
      Family = sample2$Family[[j]],
      SpeciesAccepted = sample2$SpeciesAccepted[[j]],
      Authors = sample2$Authors[[j]],
      YearPublication = sample2$YearPublication[[j]],
      Journal = sample2$Journal[j],
      Altitude = sample2$Altitude[[j]],
      CommonName = sample2$CommonName[[j]],
      e = e$values[1], R0, Tc, sigma, var,
      r0.cal, r0.approx,
      r1.cal, r1.approx,
      s1.cal, s1.approx,
      damping.cal, damping.approx,
      F_bar_cal, Pj_cal, Pa_cal,
      alpha, omega, beta, alpha.position, omega.position, beta.position, tau,
      ProjectionInterval = sample2$ProjectionInterval[j], OrganismType = sample2$OrganismType[j],
      lam1, lam2, lam3, lam4)
    
    eigenvalue.age <-  rbind(eigenvalue.age, eigenvalue1)
    print(paste(j,"/",n.case))  
    # }
  }
  eigenvalue.age
}

## =========================================
# 2, Read master database and calculate ####
## =========================================
## change to your file path
db <- cdb_fetch(paste0("./tulja_lab_master_com_p_adre.2.1.RData"))%>%
  cdb_unnest
length(unique(db@data$MatrixID))

## (1) Plant by stage ####
db_plant_stage <- db[(db$check_is_plant == TRUE) & (db$check_is_stage == TRUE)]
Compadre_Stage  <- Stage_model(db_plant_stage)%>%
  mutate(db_source = "Compadre_Stage")
length(unique(Compadre_Stage$MatrixID))                                                                                                                                                                         

# check for stasis
db_animal <- db[(db$check_is_animal == TRUE)]
stasis_plant <- sapply(matU(db_plant_stage), stasis_check)
stasis_animal <- sapply(matU(db_animal), stasis_check)

summary(stasis_plant)
summary(stasis_animal)

## (2) Animal by stage ####
db_animal_stage <- db[(db$check_is_animal == TRUE) & (db$check_is_stage == TRUE) &
                        (!db$CommonName %in% "Chinook salmon")]

Comadre_Stage <- Stage_model(db_animal_stage)%>%
  mutate(db_source = "Comadre_Stage")

## (3) Animal from age to stage ####
db$overall_survival <- sapply(matU(db), overall_survival_check)
db$single.fertility <- sapply(matF(db), single_fertility) #check if the fertility matrix only has one value
db_animal_age2stage <- db[(db$check_is_animal == TRUE) & (db$check_is_age2stage == TRUE) &
                        (!db$CommonName %in% "Chinook salmon") & (db$overall_survival == 0)]

Comadre_age2stage <- Stage_model(db_animal_age2stage)%>%
  mutate(db_source = "Comadre_adjust_2Stage")
Comadre_age2stage_phi <- out[[2]]

## (4) Animal by age ####
db_animal_age <- db[(db$check_is_animal == TRUE) & (db$check_is_age == TRUE) &
                      (!db$CommonName %in% "Chinook salmon") & (db$single.fertility == 0)]

Comadre_Age <- Age_model(db_animal_age)%>%
  mutate(db_source = "Comadre_Age")

summary(Comadre_Age$Tc - Comadre_Age2$Tc)
summary(Comadre_Age$sigma - Comadre_Age2$sigma)
summary(Comadre_Age$damping.cal - Comadre_Age2$damping.cal)

## (5) JMG_MO ####
## change to your file path
data <- read.csv("./combined dataset of JMG and MO without duplicates.csv")
id.name <- unique(data$id)

## simple leslie matrix
Leslie <- function(Lf, m, s) {
  n = length(Lf)
  M = matrix(0, n, n)
  
  # lower diagonal has survivorship ratios
  for (i in 1:(n-1)) {
    M[i+1,i] <- s[i]
  }
  M[n,n] <- s[n] ## what is the sx for open-end interval?
  
  # first row has net maternity contributions
  for(i in 1:(n-1)) {
    if(m[i] != 0 | m[i+1] != 0) {
      M[1,i] <- m[i]
    }
  }
  if (m[n] > 0) M[1,n] <-  m[n]
  M
}  

## calculation
eigenvalue <- data.frame()

for (i in 1:length(id.name)) {
  data.use <- filter(data,id == id.name[i])
  
  ## simple Leslie
  M = Leslie(data.use$lx, data.use$mx, data.use$sx)
  e <- eigen(M)
  e$values[1]
  log(e$values[1])
  arv <- abs(Re(e$vectors[,1]))
  stable <- arv/sum(arv)
  
  ## calculte r0 and r1
  lam1 <- e$values[1] # lam1 are all real numbers here
  lam2 <- e$values[2]
  tau <- sum(diff(data.use$age)) / (nrow(data.use)-1)
  ProjectionInterval = 1/tau
  r0.cal <- log(Mod(lam1))/tau
  r1.cal <- log(Mod(lam2))/tau
  s1.cal <- Arg(lam2)/tau
  damping.cal <- r0.cal - r1.cal
  
  ## calculate reproductive period
  dominate_e = Mod(e$vectors[1,])
  ## age at first reproduction by alpha, 
  ## age at last reproduction by omega, 
  ## and maximum longevity by beta
  alpha = filter(data.use, mx != 0)$age[1]
  omega = tail(filter(data.use, mx != 0)$age, n=1)
  beta = tail(data.use$age, n=1)
  
  alpha.position = match(alpha, data.use$age)
  omega.position = match(omega, data.use$age)
  beta.position = match(beta, data.use$age)
  
  tau <- sum(diff(data.use$age)) / (nrow(data.use)-1)

  ## ---------------------
  ## approximation
  ## ---------------------
  mylist <- damping_approx_agemodel(data.use$age, data.use$lx * data.use$mx)
  
eigenvalue1 <- data.frame(
    id = data.use$id[1],
    SpeciesAccepted = data.use$SpeciesAccepted[1],
    Authors = data.use$Authors[1],
    YearPublication = data.use$YearPublication[1],
    Altitude = data.use$Altitude[1],
    group = data.use$group[1],
    db_source = data.use$db_source[1],
    citation = data.use$citation[1],
    specie = data.use$specie[1],
    bodymass = data.use$bodymass[1],
    Order_Madan = data.use$Order_Madan[1],
    Family_Madan = data.use$Family_Madan[1],
    common_name = data.use$common_name[1],
    e = e$values[1], R0 = mylist$R0, Tc = mylist$Tc, sigma = mylist$sigma,
    mu3 = mylist$mu3, mu4 = mylist$mu4,
    k3 = mylist$k3, k4 = mylist$k4, r0.cal, r0.approx = mylist$r0.approx, r0.approx1 = mylist$r0.approx1,
    r0.approx_higher3 = mylist$r0.approx_higher3, r0.approx_higher4 = mylist$r0.approx_higher4,
    r1.cal, r1.approx = mylist$r1.approx, r1.approx_higher3 = mylist$r1.approx_higher3,r1.approx_higher4 = mylist$r1.approx_higher4,
    s1.cal, s1.approx= mylist$s1.approx, s1.approx_higher4= mylist$s1.approx_higher4,
    damping.cal, damping.approx = mylist$damping.approx, damping.approx_higher4 = mylist$damping.approx_higher4,damping.approx_higher3 = mylist$damping.approx_higher3,
    alpha, omega, beta, alpha.position, omega.position, beta.position,
    lam1, lam2,ProjectionInterval)
  
  eigenvalue <-  rbind(eigenvalue, eigenvalue1)
  
  print(paste(i, id.name[i]))
}
JM_MO <-eigenvalue %>%
  mutate(id = as.factor(id)) %>%
  mutate(SpeciesAccepted.original = SpeciesAccepted,
         Class = "Mammalia",
         Phylum = "Mammalia",
         Kingdom = "Animalia")%>%
  mutate(YearPublication = as.factor(YearPublication),
         Altitude = as.factor(Altitude),
         e = as.complex(e))%>%
  dplyr::rename(Order = Order_Madan,
         Family = Family_Madan,
         CommonName = common_name)%>%
  dplyr::select(-SpeciesAccepted.original)

## (6) Combine/filter and output data ####
sigma.filter <- exp(-15)
Tc.filter <- exp(5)
damping.filter <- exp(-15)

animal <- Comadre_Stage %>%
  bind_rows(Comadre_Age) %>%
  bind_rows(JM_MO) %>%
  bind_rows(Comadre_age2stage) %>%
  mutate(db_source1 = case_when(
    db_source %in% c("Choose MO","Choose JMG", "JMG", "MO") ~"GO_Age",
    db_source %in% c("Comadre_Stage", "Comadre_adjust_2Stage") ~ "Comadre_Stage",
    db_source %in% c("Comadre_Age") ~ "Comadre_Age"))%>%
  dplyr::rename(db_source.specifc = db_source,
         db_source = db_source1)%>%
  filter(sigma > sigma.filter,
         damping.cal > damping.filter,
         Tc < Tc.filter)%>%
  filter(!sigma %in% c("NA","Inf"),
         !damping.cal %in% c("NA","Inf"),
         ! Tc %in% c("NA","Inf"))

dublipcate.id.go <- c("Cervus elaphus-Benton_et_al-1995-NA-NA", #8.37, 3.41
                      "Phacochoerus aethiopicus-Rodgers-1984-NA-NA", #4.99, 2.72
                      "Rangifer tarandus-Messier_et_al-1988-NA-NA", #6.86, 3.03
                      "Spermophilus dauricus-Luo&Fox-1990-NA-NA", #2.78, 1.42
                      "Eumetopias jubatus-Calkins&Pitcher-1982-NA-NA", #11.18, 6.90 # the Order is Carnivora
                      "Hemitragus jemlahicus-Caughley-1966-NA-NA" #5.43, 2.42
)
animal <- animal%>%
  filter(!id %in% dublipcate.id.go)
animal$CommonName[animal$CommonName == "Big horn sheep"]<-"Bighorn sheep"

# Returns string without leading or trailing white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
animal$SpeciesAccepted <- trim(animal$SpeciesAccepted)

plant <- Compadre_Stage%>%
  filter(sigma > sigma.filter,
         damping.cal > damping.filter,
         Tc < Tc.filter) %>%
  filter(!sigma %in% c("NA","Inf"),
         !damping.cal %in% c("NA","Inf"),
         ! Tc %in% c("NA","Inf"))%>%
  mutate(db_source.specifc = "Compadre_Stage")
plant$CommonName[plant$SpeciesAccepted == "Fucus vesiculosus"]<-"Bladder wrack"
plant$SpeciesAccepted <- trim(plant$SpeciesAccepted)

# combine animals and plants
all_data <- plant %>%
  mutate(Altitude = as.factor(Altitude))%>%
  bind_rows(animal) %>%
  mutate(db_taxa = case_when(
    db_source %in% c("GO_Age","Comadre_Stage", "Comadre_Age","Comadre_adjust_2Stage") ~"Animal",
    db_source %in% c("Compadre_Stage", "Compadre_Age") ~ "Plant"
  ))%>%
  mutate(db_sep = case_when(
    db_source %in% c("GO_Age","Comadre_Age") ~"Animal by age",
    db_source %in% c("Comadre_Stage") ~"Animal by stage",
    db_source %in% c("Compadre_Stage", "Compadre_Age") ~ "Plant by stage"
  ))%>%
  mutate(damping.time = 1/damping.cal)
write.csv2(all_data, paste("./full animal and plant data v2.csv"),
           quote = T,na = "NA")

## ==============================
### 3, plot #####
## ==============================
## 3.1 read data and count no. of matrices and species #####
all_data <- read.csv2("./full animal and plant data v2.csv")%>%
  mutate(R0_range = case_when(
    R0<=0.9 ~ "R0 <= 0.9",
    R0<=1.1 & R0>0.9~ "0.9 < R0 <= 1.1",
    TRUE ~ "R0 > 1.1"))
all_data$R0_range <- factor(all_data$R0_range,
                            levels = c("R0 <= 0.9", "0.9 < R0 <= 1.1", "R0 > 1.1"))

Animal <- filter(all_data, db_taxa %in% "Animal")
Plant <- filter(all_data, db_taxa %in% "Plant")
Comadre_Stage <- filter(all_data, db_sep %in% "Animal by stage")
Compadre_Stage <- filter(all_data, db_sep %in% "Plant by stage")
Animal_Age <- filter(all_data, db_sep %in% c("Animal by age"))

## count no. of species
all_data %>%
  group_by(db_source.specifc)%>%
  summarise(count = n_distinct(SpeciesAccepted))

all_data %>%
  group_by(db_source)%>%
  summarise(count = n_distinct(SpeciesAccepted))

all_data %>%
  group_by(db_taxa)%>%
  summarise(count = n_distinct(SpeciesAccepted))

## count no. of matrices
all_data %>%
  group_by(db_source.specifc)%>%
  count()

all_data %>%
  group_by(db_source)%>%
  count()

all_data %>%
  group_by(db_taxa)%>%
  count()

#### 3.2 linear regression ####
lm <- lm(log(Animal_Age$sigma)~log(Animal_Age$Tc))
lm <- lm(log(Comadre_Stage$sigma)~log(Comadre_Stage$Tc))
lm <- lm(log(Compadre_Stage$sigma)~log(Compadre_Stage$Tc))
summary(lm)
confint(lm, 'log(Animal_Age$Tc)', level=0.95)
confint(lm, 'log(Comadre_Stage$Tc)', level=0.95)
confint(lm, 'log(Compadre_Stage$Tc)', level=0.95)

lm <- lm(log(1/Animal_Age$damping.cal)~log(Animal_Age$Tc))
lm <- lm(log(1/Comadre_Stage$damping.cal)~log(Comadre_Stage$Tc))
lm <- lm(log(1/Compadre_Stage$damping.cal)~log(Compadre_Stage$Tc))
summary(lm)
confint(lm, 'log(Animal_Age$Tc)', level=0.95)
confint(lm, 'log(Comadre_Stage$Tc)', level=0.95)
confint(lm, 'log(Compadre_Stage$Tc)', level=0.95)

### 3.3 plots based on OLS ####
### These figures are not presented in the paper ###
#### (1) Tc and sigma #####
myshape = c(4, 21, 21, 21)
all_data$db_source <- factor(all_data$db_source, levels = c("GO_Age","Comadre_Age",
                                                            "Comadre_Stage","Compadre_Stage"))

p_Tc_sigma <- ggplot(all_data, 
                     aes(log(Tc), log(sigma))) +
  geom_point(aes(shape = db_source), alpha = 0.7)+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_shape_manual(values=myshape)+
  theme_bw()+
  xlab(expression(log(Tc)))+ylab(expression(log(sigma)))+
  theme(legend.title=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  facet_wrap(. ~db_sep, nrow = 3)
p_Tc_sigma
ggsave("./plot/Tc and sigma.png", p_Tc_sigma , width = 9, height = 9)

p_Tc_sigma <- ggplot(all_data, 
                     aes(Tc, sigma)) +
  geom_point(aes(shape = db_source), alpha = 0.7)+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_shape_manual(values=myshape)+
  theme_bw()+
  xlab(expression(T[c]))+ylab(expression(S))+
  theme(legend.title=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  facet_wrap(. ~db_sep, nrow = 3, scales="free_y")
p_Tc_sigma
ggsave("./plot/Tc and sigma without log.png", p_Tc_sigma , width = 9, height = 9)

p_Tc_sigma_class1 <- ggplot(Animal, 
                           aes(log(Tc), log(sigma))) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 4) +
  theme_bw()+
  theme(legend.title=element_blank())+
  xlab(expression(log(Tc)))+ylab(expression(log(sigma)))+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  facet_wrap(~Class, ncol = 3)+
  geom_point(shape = 21, alpha = 0.7)
p_Tc_sigma_class1
ggsave("./plot/Tc and sigma (animal), facet by Class.png", p_Tc_sigma_class1, width = 9, height = 10)

p_Tc_sigma_class2 <- ggplot(Plant, 
                       aes(log(Tc), log(sigma))) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 4) +
  theme_bw()+
  theme(legend.title=element_blank())+
  xlab(expression(log(Tc)))+ylab(expression(log(sigma)))+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  facet_wrap(~Class)+
  geom_point(shape = 21, alpha = 0.7)
p_Tc_sigma_class2
ggsave("./plot/Tc and sigma (plant), facet by Class.png", p_Tc_sigma_class2, width = 9, height = 9)

#### (2) Tc and damping #####
p_Tc_damping <- ggplot(all_data, 
                       aes(log(Tc), log(1/damping.cal)))+
geom_point(aes(shape = db_source), alpha = 0.7)+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_shape_manual(values=myshape)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(Tc)))+ylab(expression(log(tau)))+
  facet_wrap(. ~db_sep, nrow = 3)
p_Tc_damping
ggsave("./plot/Tc and damping.png", p_Tc_damping , width = 9, height = 9)

p_Tc_damping <- ggplot(all_data, 
                       aes(Tc, 1/damping.cal))+
  geom_point(aes(shape = db_source), alpha = 0.7)+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_shape_manual(values=myshape)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(T[c]))+ylab(expression(tau))+
  facet_wrap(. ~db_sep, nrow = 3)
p_Tc_damping
ggsave("./plot/Tc and damping without log.png", p_Tc_damping , width = 9, height = 9)

p_Tc_damping_class1 <- ggplot(Animal, 
                        aes(log(Tc), log(1/damping.cal))) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 4) +
  facet_wrap(~Class, ncol = 3)+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  xlab(expression(log(Tc)))+ylab(expression(log(tau)))+
  geom_point(shape = 21, alpha = 0.7)
p_Tc_damping_class1
ggsave("./plot/Tc and damping (animal), facet by Class.png", p_Tc_damping_class1, width = 9, height = 10)

p_Tc_damping_class2 <- ggplot(Plant, 
                        aes(log(Tc), log(1/damping.cal))) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 4) +
  facet_wrap(~Class)+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  xlab(expression(log(Tc)))+ylab(expression(log(tau)))+
  geom_point(shape = 21, alpha = 0.7)
p_Tc_damping_class2
ggsave("./plot/Tc and damping (plant), facet by Class.png", p_Tc_damping_class2 , width = 9, height = 9)

#### (3) Tc and reproductive period #####
p_Tc_repperiod <- ggplot(filter(Animal, 
                                db_source %in% c("GO_Age", "Comadre_Age")), 
                         aes(log(Tc), log(omega-alpha), color = Class)) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 5) +
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        legend.background = element_rect(fill = "white", color = "gray"))+
  xlab(expression(log(omega-alpha)))+ylab(expression(log(Tc)))+
  geom_point(shape = 21, alpha = 0.7)
p_Tc_repperiod
ggsave("./plot/Tc and reproduction period (age-structured animal).png", p_Tc_repperiod , width = 9, height = 9)

#### (4) residual of Tc damping and sigma #####
reg.a <- lm(log(1/damping.cal)~log(Tc), data = Animal)
reg2.a <- lm(reg.a$residuals~ log(Animal$sigma))
summary(reg2.a)
plot(log(reg.a$residuals),log(Animal$sigma))

db.a<-as.data.frame(bind_cols(reg.a$residuals, Animal$sigma, Animal$db_taxa))
names(db.a) <- c("residuals","sigma", "db_taxa")

reg.p <- lm(log(1/damping.cal)~log(Tc), data = Plant)
reg2.p <- lm(reg.p$residuals~ log(Plant$sigma))
summary(reg2.p)
plot(log(reg.p$residuals),log(Plant$sigma))

db.p<-as.data.frame(bind_cols(reg.p$residuals, Plant$sigma, Plant$db_taxa))
names(db.p) <- c("residuals","sigma", "db_taxa")

db<-rbind(db.a, db.p)

p <- ggplot(db, aes(log(sigma),residuals)) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3, size = 5) +
  geom_point(shape = 21, alpha = 0.7)+
  facet_grid(. ~db_taxa)+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.89, 0.15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  ylim(-6,6)+
  ylab(expression(residuals))+xlab(expression(log(sigma)))
p
ggsave("./plot/residuals (of Tc and damping) and sigma.png", p , width = 9, height = 6)

#### (5) Coefficient of Variation #####
p_invariant3 <- ggplot(Animal, 
                       aes(Order, sigma/Tc, color = Class)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 15))+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.89, 0.15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  ylab(expression(sigma/Tc))+
  geom_point(shape = 21, alpha = 0.7)
p_invariant3
ggsave("./plot/coefficient of variation (animal).png", p_invariant3 , width = 9, height = 9)

p_invariant4 <- ggplot(Plant, 
                       aes(Order, sigma/Tc, color = Class)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size = 15))+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.89, 0.15),
        legend.background = element_rect(fill = "white", color = "gray"))+
  ylab(expression(sigma/Tc))+
  geom_point(shape = 21, alpha = 0.7)
p_invariant4
ggsave("./plot/coefficient of variation (plant).png", p_invariant4, width = 9, height = 9)


### (6) compare calculated damping (from PPM) and approximated damping (from analytical approximation)
p_damping <- ggplot(all_data, 
                    aes(log(1/damping.cal), log(1/damping.approx))) +
  geom_point(aes(shape = db_source), alpha = 0.7)+
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_shape_manual(values=myshape)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15))+
  facet_wrap(. ~db_sep, nrow = 3) + 
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(paste(log(tau), " calculated from PPM")))+
  ylab(expression(paste(log(tau), " from analytical approximation")))
p_damping
ggsave("./plot/damping calcuated and approximated.png", p_damping , width = 9, height = 9)

######### (7) main results by R0 #######
p_Tc_sigma_R0 <- ggplot(all_data, 
                        aes(log(Tc), log(sigma), 
                            color = db_source)) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_color_manual(values = c("#F8766D","#C77CFF","#7CAE00","#00BFC4")) +
  geom_point(shape = 21, alpha = 1)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(Tc)))+ylab(expression(log(sigma)))+
  facet_wrap(. ~db_sep+R0_range, nrow = 3)
p_Tc_sigma_R0
ggsave("./plot/Tc and sigma by R0.png", p_Tc_sigma_R0 , width = 9, height = 9)

p_Tc_damping_R0 <- ggplot(all_data, 
                          aes(log(Tc), log(1/damping.cal), 
                              color = db_source)) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_color_manual(values = c("#F8766D","#C77CFF","#7CAE00","#00BFC4")) +
  geom_point(shape = 21, alpha = 1)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  xlab(expression(log(Tc)))+ylab(expression(log(tau)))+
  facet_wrap(. ~db_sep+R0_range, nrow = 3)
p_Tc_damping_R0
ggsave("./plot/Tc and damping by R0.png", p_Tc_damping_R0 , width = 9, height = 9)

p_R0_damping <- ggplot(all_data, 
                       aes(log(R0), log(Tc/sigma), 
                           color = db_source)) +
  stat_poly_eq(formula = y~x,aes(label = paste(..eq.label.., ..rr.label..,..p.value.label.., sep = "~','~")),
               parse = TRUE, color = "blue", rr.digits = 2, coef.digits = 3,size = 5) +
  scale_color_manual(values = c("#F8766D","#C77CFF","#7CAE00","#00BFC4")) +
  geom_point(shape = 21, alpha = 1)+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=15),
        axis.text.x = element_text(color = "grey20", size = 20),
        axis.text.y = element_text(color = "grey20", size = 20), 
        axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20), 
        strip.text.x = element_text(size = 15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  facet_wrap(. ~db_sep, nrow = 2)
p_R0_damping
ggsave("./plot/R0 and damping.png", p_R0_damping , width = 9, height = 9)
