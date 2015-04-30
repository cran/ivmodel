# condPvalue: This is the conditional p-value integral described in Andrews, Moreira, and Stock (2007, Journal of Econometrics)
#             P(m,qT) = 1 - Pr(LR > m |QT = qT)
#          m: test stat value
#         qT:  conditional value of Qt
#          k: number of moment conditions (i.e. number of instruments)
#        eps: epsilon (for numerical stability, default set based on Andrews, Moreira, and Stock's 2007 recommendation (0.02)
# OUTPUT: a p-value
# NOTE: mimics condivreg.ado's new_try function in Mikusheva and Poi (2006)'s STATA function
condPvalue <- function(m,qT,k,eps = 0.02) {
  if(k == 1) { # Eqivalent to AR #
    return(1 - pchisq(m,k))
  }
  K = gamma(k/2) / (sqrt(pi) * gamma((k-1)/2))  
  if(k == 2) {
    return(1 - 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * sin(x)^2/m),k)},lower=0,upper=pi/2)$value) 
  } else if(k == 3) {
    return(1 - 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k)},lower=0,upper=1)$value)
  } else if(k == 4) {
    nonapproxRegion = 2 * K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k) * (1 - x^2)^((k-3)/2)},lower=0,upper=1-eps)$value
    approxRegion =  2 * K * pchisq( (qT + m) / (1 + qT *(1 - eps/2)^2 /m),k) * (1/2 * (pi/2 - asin(1 - eps)) - (1 - eps)/2 * sqrt( 1- (1 - eps)^2)) 
    return(1 - nonapproxRegion - approxRegion)
  } else {
    return(1 - 2*K * integrate(function(x){pchisq( (qT + m)/(1 + qT * x^2/m),k) * (1 - x^2)^((k-3)/2)},lower=0,upper=1)$value)
  }
}


# invTwobyTwoSymMatrix: inverts two-by-two symmetric matrices
# M: symmetric matrix
# OUTPUT: a two-by-two symmetric matrix
invTwobyTwoSymMatrix <- function(M) {
  if(missing(M)) stop("invTwobyTwoSymMatrix: Supply M!") 
  if(nrow(M) != 2 || ncol(M) != 2) stop("invTwobyTwoSymMatrix: The matrix is not 2 by 2")
  if(M[1,2] != M[2,1]) stop("invTwobyTwoSymMatrix: The matrix is not symmetric!")
  detM = M[1,1] * M[2,2] - M[1,2]^2
  Minverse = M
  Minverse[2,2] = M[1,1]
  Minverse[1,1] = M[2,2]
  Minverse[1,2] = Minverse[2,1] = -1 * M[1,2]
  return(1/detM * Minverse)
}

# quadSolver: solves the quadratic equation ax^2 + bx + c = 0
#             Can also be used to solve eigenvalues of two-by-two matrices
#             E.x. if [d,e; f,g] is the two-by-two matrix
#                  eigenvalues are defined as the solution to det( [d,e; f,g] - lambda * I) = 0
#                  Or det([d - lambda, e; f, g - lambda]) = 0
#                  Or (d-lambda) * (g -lambda) - e*f = 0
#                  Or lambda^2 - (d +g) *lambda + d*g - e*f = 0
quadSolver = function(a,b,c) {
  if(missing(a) || missing(b) || missing(c)) stop("Supply the coefficients!")
  if(length(a) != 1 || length(b) != 1 || length(c) != 1) stop("Supply only scalar values to a,b, and c")
  if(a == 0) stop("Something is wrong! Leading coefficient is zero")
  
  discriminant = b^2 - 4 * a * c
  if(discriminant < 0) return(c(NA,NA))
  if(discriminant == 0) return( rep(-b/(2*a),2))
  if(discriminant > 0) {
    if(a > 0) {
	  lowerRoot = -b / (2*a) - sqrt(discriminant) / (2*a)
	  upperRoot = -b / (2*a) + sqrt(discriminant) / (2*a)
	  return(c(lowerRoot,upperRoot))
	} else {
	  lowerRoot = -b / (2*a) + sqrt(discriminant) / (2*a)
	  upperRoot = -b / (2*a) - sqrt(discriminant) /(2*a)
	  return(c(lowerRoot,upperRoot))
	}
  }
}



### kClassEst: Generates point estimates and standard errors 
###                  Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        k, a vector (or scalar) of k values for estimation
###        beta0, a vector (or scalar) of null values for testing and p-values
###        alpha, significance level for confidence intervals
###        heteroSE, use heteroscedastic robust standard errors?
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
kClassEst = function(ivmodel,beta0=0,alpha=0.05,k=c(0,1),heteroSE=FALSE) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("You must supply an ivmodel class. Run ivmodel() and see ivmodel() function for details")
    return(NULL)
  }
  if(missing(k)) {
    print("You must specify a value for k")
    return(NULL)
  }
 
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj =ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR
  degF = ivmodel$n - ivmodel$p - 1
  
  # Compute k-class estimator
  denom = (sum(Dadj^2) - k * sum(qr.resid(ZadjQR,Dadj)^2))
  kPointEst = (sum(Yadj * Dadj) - k * sum(Yadj * qr.resid(ZadjQR,Dadj))) / denom 
  
  # Compute std error, testStat, and confidence intervals
  kVarPointEst = rep(0,length(k))
  kCI = matrix(0,length(k),2)
  for(i in 1:length(k)) {
    if(denom[i] <= 0) { #k exceeds the maximum possible value
	  kVarPointEst[i] = NA; kCI[i,] = c(NA,NA)
	} else {
      if(heteroSE) {
	    kVarPointEst[i] = ivmodel$n/degF * sum( (Yadj - Dadj * kPointEst[i])^2 * (Dadj - k[i]*qr.resid(ZadjQR,Dadj))^2 ) / (denom[i])^2
	  } else {
        kVarPointEst[i] = 1/(degF) *sum((Yadj - Dadj * kPointEst[i])^2) / denom[i]
	  }
	  # Note that R's AER package uses the z-scores instead of the t distribution. They are basically the same as n gets large.
	  kCI[i,] =  c(kPointEst[i] - qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]),
	               kPointEst[i] + qt(1-alpha/2,degF) * sqrt(kVarPointEst[i]))
    }
  }
  
  # Compute test statistics
  # This operation takes a vector of k estimates (pointEst, k*1 dim) and 
  # another vector of null values (beta0,1 * length(beta0))
  # and subtracts pointEst from a vector of null values.
  # The last operation of dividing by sqrt() is done column-wise 
  # E.g. [1,2,3 ; 2,3,4] / (2,3,4) --> [0.5,2/3,0.75;1,1,1] 
  kTestStat = outer(kPointEst,beta0,FUN="-") / sqrt(kVarPointEst)
  
  # Compute p-value
  kPValue = 2*(1 - pt(abs(kTestStat),degF))
  
  # Package output 
  kPointEst = matrix(kPointEst,length(k),1)
  kVarPointEst = matrix(sqrt(kVarPointEst),length(k),1)
  rownames(kPointEst) = rownames(kVarPointEst) = rownames(kCI) = rownames(kPValue) = rownames(kTestStat) = k
  colnames(kPValue) = colnames(kTestStat) = beta0
  colnames(kPointEst) = "Estimate"
  colnames(kVarPointEst) = "Std. Error"
  colnames(kCI) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(point.est = kPointEst,std.err = kVarPointEst,test.stat = kTestStat,p.value = kPValue,ci = kCI))
}

### LIMLEst: Generates point estimates and standard errors 
###                  Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        beta0, a vector (or scalar) of null values for testing and p-values
###        alpha, significance level for confidence intervals
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
LIMLEst = function(ivmodel,beta0=0,alpha=0.05,heteroSE=FALSE) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("LIML: You must supply an ivmodel class. See ivmodel function for details")
    return(NULL)
  }
 
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR
  
  # Value of k for LIML
  LIMLMatrix1 = matrix(0,2,2)
  LIMLMatrix1[1,1] = sum(Yadj^2)
  LIMLMatrix1[1,2] = LIMLMatrix1[2,1] = sum(Dadj * Yadj)
  LIMLMatrix1[2,2] = sum(Dadj^2)
  
  LIMLMatrix2 = matrix(0,2,2); projYadj = qr.resid(ZadjQR,Yadj); projDadj = qr.resid(ZadjQR,Dadj)
  LIMLMatrix2[1,1] = sum(projYadj^2)
  LIMLMatrix2[1,2] = LIMLMatrix2[2,1] = sum(projDadj * projYadj)
  LIMLMatrix2[2,2] = sum(projDadj^2)
  
  kLIML = eigen(LIMLMatrix1 %*% invTwobyTwoSymMatrix(LIMLMatrix2))$values[2]
  output = kClassEst(ivmodel,k=kLIML,beta0=beta0,alpha=alpha,heteroSE=heteroSE)
  
  # Package output
  rownames(output$point.est) = rownames(output$std.err) = rownames(output$test.stat) = rownames(output$ci) = rownames(output$p.value) = NULL
  return(c(output,list(k=kLIML)))
}

### FullerEst: Generates point estimates and standard errors 
###            Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        b, adjustment by Fuller
###        beta0, a vector (or scalar) of null values for testing and p-values
###        alpha, significance level for confidence intervals
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
FullerEst = function(ivmodel,beta0=0,alpha=0.05,b=1,heteroSE=FALSE) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("Fuller: You must supply an ivmodel class. See ivmodel function for details")
	return(NULL)
  } 
  if(b <= 0) {
    print("Fuller: Fuller adjustment cannot be less than 0")
	return(NULL)
  }
 
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; ZadjQR = ivmodel$ZadjQR
  
  # Value of k for LIML
  LIMLMatrix1 = matrix(0,2,2)
  LIMLMatrix1[1,1] = sum(Yadj^2)
  LIMLMatrix1[1,2] = LIMLMatrix1[2,1] = sum(Dadj * Yadj)
  LIMLMatrix1[2,2] = sum(Dadj^2)
  
  LIMLMatrix2 = matrix(0,2,2); projYadj = qr.resid(ZadjQR,Yadj); projDadj = qr.resid(ZadjQR,Dadj)
  LIMLMatrix2[1,1] = sum(projYadj^2)
  LIMLMatrix2[1,2] = LIMLMatrix2[2,1] = sum(projDadj * projYadj)
  LIMLMatrix2[2,2] = sum(projDadj^2)
  
  kLIML = eigen(LIMLMatrix1 %*% invTwobyTwoSymMatrix(LIMLMatrix2))$values[2]
  kFuller = kLIML - b/(ivmodel$n - ivmodel$L - ivmodel$p)
  output = kClassEst(ivmodel,k=kFuller,beta0=beta0,alpha=alpha,heteroSE=heteroSE)
  
  # Package output
  rownames(output$point.est) = rownames(output$std.err) = rownames(output$test.stat) = rownames(output$ci) = rownames(output$p.value) = NULL
  return(c(output,list(k=kFuller)))
}

### CLR.test: Generates the CLR confidnece interval
###         Must run ivmodel before you run this code
### INPUT: ivmodel, an object from ivmodel() function
###        alpha, significance level for confidence intervals
### OUTPUT: a list of point estimate, standard error, test statistic, and p-value
CLR.test = function(ivmodel,beta0=0,alpha=0.05) {
  # Error checking
  if(class(ivmodel) != "ivmodel") {
    print("CLR: You must supply an ivmodel class. See ivmodel function for details")
	return(NULL)
  }
  # Extract objects from ivmodel
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; ZadjQR = ivmodel$ZadjQR
  
  # Qs, Qst, Qst, Qt #
  YadjANDDadj = cbind(Yadj,Dadj)
  PZYadjANDDadj = qr.fitted(ZadjQR,YadjANDDadj); 
  RZYadjANDDadj = qr.resid(ZadjQR,YadjANDDadj); 
  
  sigmaHat = (t(RZYadjANDDadj) %*% RZYadjANDDadj) / (ivmodel$n - ivmodel$p - ivmodel$L)
  sigmaHatInv = invTwobyTwoSymMatrix(sigmaHat) 
  
  a0 = c(beta0,1); b0 = c(1,-beta0)
  denomS = sum((sigmaHat %*% b0)^2); denomT = sum((sigmaHatInv %*% a0)^2)
  numT = PZYadjANDDadj %*% sigmaHatInv 
  QS = sum( (PZYadjANDDadj %*% b0)^2) / denomS
  QT = sum( (numT %*% a0)^2) / denomT
  QTS = sum( b0 * (numT %*% a0)) / (sqrt(denomS) * sqrt(denomT))
  
  LRtest = 1/2 * (QS - QT + sqrt((QS + QT)^2 - 4*(QS *QT - QTS^2)))
  
  test.stat = matrix(LRtest,1,1)
  p.value = matrix(condPvalue(LRtest,QT,ivmodel$L),1,1) 
  
  maxEigen = max(quadSolver(a=1,b=-1*(QS + QT),c=QS *QT - QTS^2)) #of Q matrix
  C = tryCatch({uniroot(function(y){condPvalue(m=maxEigen - y,qT = y,k = ivmodel$L) - alpha},lower=0,upper=maxEigen)$root},
            error=function(e){0}) 
  quadMatrix.CLR = sigmaHatInv %*% t(YadjANDDadj) %*% numT - C * sigmaHatInv
  ci.CLR = quadSolver(a=quadMatrix.CLR[1,1],b=2*quadMatrix.CLR[1,2],c=quadMatrix.CLR[2,2])
  if(quadMatrix.CLR[1,1] > 0) {
    if(is.na(ci.CLR[1])) {
      ci.CLR = matrix(c(-Inf,Inf),1,2)
	  info = c("Whole Real Line")
    } else {
	  if( abs(ci.CLR[1] - ci.CLR[2]) < 10^-6) {
	    ci.CLR = matrix(c(-Inf,Inf),1,2)
		info = c("Whole Real Line")
      }  else {
	    ci.CLR = matrix(c(-Inf,ci.CLR[1],ci.CLR[2],Inf),2,2,byrow=TRUE)
		info = c(paste("(-Infinity, ",ci.CLR[1,2],"] union [",ci.CLR[2,1],", Infinity)",sep=""))
	  }
    }
  } else {
    ci.CLR = matrix(c(ci.CLR[1],ci.CLR[2]),1,2)
	info = paste("[",ci.CLR[1,1],", ",ci.CLR[1,2],"]",sep="")
  }
  
  # Package output 
  colnames(p.value) = colnames(test.stat) = beta0
  colnames(ci.CLR) = c(paste(as.character(round(alpha/2 * 100,1)),"%"),paste(as.character( round((1-alpha/2) * 100,1)),"%"))
  
  return(list(test.stat = test.stat,p.value = p.value,ci = ci.CLR,ci.info = info))
}

#####  AR test and CI

AR.test=function(ivmodel, beta0=0, alpha=0.05){
  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; 
  n=ivmodel$n;  k=ivmodel$p;  l=ivmodel$L

  if(ncol(Yadj)>1){
    print("The outcome variable should be in one dimension!")
	return(NULL)
  }
  if(ncol(Dadj)>1){
    print("The treatment variable should be in one dimension!")
	return(NULL)
  }
  if(l+k>=n){
    print("Too many IVs, AR can't handle!")
    return(NULL)
  }

### test  
  temp=Yadj-beta0*Dadj
  Fstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(Fstat, df1=l, df2=n-k-l)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l)*l/(n-k-l)
  coef.beta0sq=cval*sum(Dadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Dadj)^2)
  coef.beta0=-2*cval*sum(Dadj*Yadj)+2*(cval+1)*sum(Dadj*qr.fitted(ivmodel$ZadjQR, Yadj))
  coef.constant=cval*sum(Yadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Yadj)^2)
  Delta=coef.beta0^2-4*coef.constant*coef.beta0sq

  ci=matrix(NA, ncol=2)
  colnames(ci)<-c("lower", "upper")

  if(coef.beta0sq==0){
    if(coef.beta0>0){
      info=c("[",-coef.constant/coef.beta0,",Infinity)")
      ci[1,]=c(-coef.constant/coef.beta0, Inf)
    }
    if(coef.beta0<0){
      info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
      ci[1,]=c(-Inf, -coef.constant/coef.beta0)
    }
    if(coef.beta0==0){
      if(coef.constant>=0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.constant<0){
        info="Empty Set"
      }
    }
  }
  
  if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",lower.root,",",upper.root,"]")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }  

  return(list(Fstat=Fstat, df=c(l, n-k-l), p.value=p.value,
              ci.info=info, ci=ci))
}



#####  calculate the power of AR test

AR.power=function(n, k, l, beta, gamma, Zadj_sq, 
                  sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
	print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*n*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	    (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  temp = qf(1-alpha, df1=l, df2=n-k-l)
  power = 1-pf(temp, df1=l, df2=n-k-l, ncp=ncp)

  return(power)
}


##### calculate the sample size needed for certain power

AR.size=function(power, k, l, beta, gamma, Zadj_sq, 
                 sigmau, sigmav, rho, alpha=0.05){
  if(length(gamma)!=ncol(as.matrix(Zadj_sq))|length(gamma)!=nrow(as.matrix(Zadj_sq))){
    print("The dimension of Zadj_sq doesn't match gamma")
	stop()
  }
  ncp=beta^2*c(t(as.matrix(gamma))%*%as.matrix(Zadj_sq)%*%as.matrix(gamma))/
	  (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)  

  oldn<-k+l+1
  state<-1
  while(state){
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*oldn)  
    if(temppower < power){
	  oldn <- oldn*2
	}else{
	  state <- 0
	}
  }

  lower <- oldn%/%2
  upper <- oldn
  while((upper-lower)>2){
    new <- (upper+lower)%/%2
    temp = qf(1-alpha, df1=l, df2=oldn-k-l)  
    temppower = 1-pf(temp, df1=l, df2=oldn-k-l, ncp=ncp*new)  
    if(temppower < power){
	  lower <- new
	}else{
	  upper <- new
	}
  }

  return(upper)
}


#####  AR sensitivity test and CI

ARsens.test=function(ivmodel, beta0=0, alpha=0.05, deltarange=NULL){
  if(is.null(deltarange))
    return(NULL)

  Yadj = ivmodel$Yadj; Dadj = ivmodel$Dadj; Zadj = ivmodel$Zadj; 
  n=ivmodel$n;  k=ivmodel$p;  l=ivmodel$L

  if(ncol(Yadj)>1){
    print("The outcome variable should be in one dimension!")
	return(NULL)
  }
  if(ncol(Dadj)>1){
    print("The treatment variable should be in one dimension!")
	return(NULL)
  }
  if(l+k>=n){
    print("Too many IVs, AR can't handle!")
    return(NULL)
  }
  
  if(l!=1){
    print("Please input exact one IV for AR sensitivity analysis")
	return(NULL)
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    ncp=max(deltarange^2)*sum(Zadj^2)
  }else{
    print("Wrong input of the sensitivity range.")
    return(NULL)
  }
  
### test  
  temp=Yadj-beta0*Dadj
  ncFstat=c(sum(qr.fitted(ivmodel$ZadjQR, temp)^2))/c(sum(temp^2)-sum(qr.fitted(ivmodel$ZadjQR, temp)^2))*(n-k-l)/l
  p.value=1-pf(ncFstat, df1=l, df2=n-k-l, ncp=ncp)
  
### confidence interval  
  cval=qf(1-alpha, df1=l, df2=n-k-l, ncp=ncp)*l/(n-k-l)
  coef.beta0sq=cval*sum(Dadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Dadj)^2)
  coef.beta0=-2*cval*sum(Dadj*Yadj)+2*(cval+1)*sum(Dadj*qr.fitted(ivmodel$ZadjQR, Yadj))
  coef.constant=cval*sum(Yadj^2)-(cval+1)*sum(qr.fitted(ivmodel$ZadjQR, Yadj)^2)
  Delta=coef.beta0^2-4*coef.constant*coef.beta0sq

  ci=matrix(NA, ncol=2)
  colnames(ci)<-c("lower", "upper")

  if(coef.beta0sq==0){
    if(coef.beta0>0){
      info=c("[",-coef.constant/coef.beta0,",Infinity)")
      ci[1,]=c(-coef.constant/coef.beta0, Inf)
    }
    if(coef.beta0<0){
      info=c("(-Infinity,",-coef.constant/coef.beta0,"]");
      ci[1,]=c(-Inf, -coef.constant/coef.beta0)
    }
    if(coef.beta0==0){
      if(coef.constant>=0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.constant<0){
        info="Empty Set"
      }
    }
  }
  
  if(coef.beta0sq!=0){
    if(Delta<=0){
      if(coef.beta0sq>0){
        info="Whole Real Line"
        ci[1,]=c(-Inf, Inf)
      }
      if(coef.beta0sq<0){
        info="Empty Set"
      }
    }
    if(Delta>0){
      # Roots of quadratic equation
      root1=(-coef.beta0+sqrt(Delta))/(2*coef.beta0sq)
      root2=(-coef.beta0-sqrt(Delta))/(2*coef.beta0sq)
      upper.root=max(root1,root2)
      lower.root=min(root1,root2)
      if(coef.beta0sq<0){
        info=paste("[",lower.root,",",upper.root,"]")
        ci[1, ]=c(lower.root, upper.root)
      }
      if(coef.beta0sq>0){
        info= paste("(-Infinity,",lower.root,"] union [",upper.root,",Infinity)")
        ci[1, ]=c(-Inf, lower.root)
        ci<-rbind(ci, c(upper.root, Inf))
      }
    }
  }  

  return(list(ncFstat=ncFstat, df=c(l, n-k-l), ncp=ncp, 
              p.value=p.value, ci.info=info, ci=ci, deltarange=deltarange))
}


#####  calculate the power of AR sensitivity analysis

ARsens.power=function(n, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                      rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*n*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     n*Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*n*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  temp = qf(1-alpha, df1=1, df2=n-k-1, ncp=ncp2)
  power = 1-pf(temp, df1=1, df2=n-k-1, ncp=ncp1)

  return(power)
}


##### calculate the sample size needed for certain power of sensitivity analysis

ARsens.size=function(power, k, beta, gamma, Zadj_sq, sigmau, sigmav, 
                     rho, alpha=0.05, deltarange=deltarange, delta=NULL){
  if(!is.numeric(gamma) | length(gamma)!=1){
	print("Wrong input of gamma or gamma is not one dimension")
	stop()
  }
  
  if(!is.numeric(Zadj_sq) | length(Zadj_sq)!=1){
	print("Wrong input of Zadj_sq or Zadj_sq is not one dimension")
	stop()
  }

  if(is.numeric(deltarange) & length(deltarange)==2){
    deltarange<-sort(deltarange)
    ncp2=max(deltarange^2)*Zadj_sq
	if(is.null(delta)){
	  if((deltarange[1] < -gamma*beta/sigmau) & 
	     (deltarange[2] > -gamma*beta/sigmau)){
	    ncp1=0
	  }else{
	    ncp1=min((beta*gamma+deltarange[1]*sigmau)^2,
		         (beta*gamma+deltarange[2]*sigmau)^2)*
		     Zadj_sq/(sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	  }	
	}else{
	  ncp1=(beta*gamma+delta*sigmau)^2*Zadj_sq/
	       (sigmau^2+2*rho*sigmau*sigmav*beta+sigmav^2*beta^2)
	}
  }else{
    print("Wrong input of the sensitivity range.")
    stop()
  }

  if(ncp1<=ncp2){
    print("Sensitivity range too large")
	stop()
  }else{
    oldn<-k+2
    state<-1
    while(state){
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*oldn)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*oldn)  
      if(temppower < power){
	    oldn <- oldn*2
	  }else{
	    state <- 0
	  }
    } 

    lower <- oldn%/%2
    upper <- oldn
    while((upper-lower)>2){
      new <- (upper+lower)%/%2
      temp = qf(1-alpha, df1=1, df2=oldn-k-1, ncp=ncp2*new)  
      temppower = 1-pf(temp, df1=1, df2=oldn-k-1, ncp=ncp1*new)  
      if(temppower < power){
	    lower <- new
	  }else{
	    upper <- new
	  }
    }
  }
  return(upper)
}
