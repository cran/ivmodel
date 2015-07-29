### functions for QR or sparseQR class
qrrank<-function(QR, tol=1e-8){
  if(class(QR)=="sparseQR")
    return(sum(abs(diag(QR@R))>tol))
  else
    return(QR$rank)
}
qrRM<-function(QR){
  if(class(QR)=="sparseQR")
    return(qrR(QR))
  else
    return(qr.R(QR))
}


### ivmodel: Generates an instance of ivmodel
###          In this package, our IV model is
###          Y_i = beta_0 + D_i * beta_1 + X_i^T * gamma + epsilon_i
###          E(epsilon_i | Z_i, X_i) = 0
### INPUT: Y, outcome (n * 1 vector)
###        D, exposure (n * 1 vector)
###        Z, instruments (n * L matrix)
###        X, baseline and exogenous covariates (n * p matrix)
###        intercept, should we include the intercept term?
### OUTPUT: ivmodel object which contains a clean-up version of Y,D,Z,X (if available)

ivmodel <- function(Y,D,Z,X,intercept=TRUE,beta0=0,alpha=0.05,k=c(0,1), heteroSE = FALSE, deltarange=NULL) {
  # Error checking: check to see if necessary inputs are there!
  if(missing(Y)) stop("Y is missing!")
  if(missing(D)) stop("X is missing!")
  if(missing(Z)) stop("Z is missing!")
  
  # Error checking: check Y and D
  if( (!is.vector(Y) && !is.matrix(Y) && !is.data.frame(Y)) || (is.matrix(Y) && ncol(Y) != 1) || (is.data.frame(Y) && ncol(Y) != 1) || (!is.numeric(Y))) stop("Y is not a numeric vector.")
  if( (!is.vector(D) && !is.matrix(D) && !is.data.frame(Y)) || (is.matrix(D) && ncol(D) != 1) || (is.data.frame(D) && ncol(D) != 1) || (!is.numeric(D))) stop("D is not a numeric vector.")
  Y = as.numeric(Y); D = as.numeric(D)
  if(length(Y) != length(D)) stop("Dimension of Y and D are not the same!")
  
  # Error checking: check Z and convert "strings" into factors
  Z = data.frame(Z); stringIndex = sapply(Z,is.character); Z[stringIndex] = lapply(Z[stringIndex],as.factor)
  if(nrow(Z) != length(Y)) stop("Row dimension of Z and Y are not equal!")
  colnames(Z) = paste("Z",colnames(Z),sep="")
  
  # Add intercept as X
  if(intercept && !missing(X)) X = data.frame(X,1)
  if(intercept && missing(X)) X = data.frame(rep(1,length(Y)))
  
  # Error checking: check X and convert "strings" into factors
  if(!missing(X)) {
    X = data.frame(X); stringIndex = sapply(X,is.character); X[stringIndex] = lapply(X[stringIndex],as.factor)
	if(nrow(X) != length(Y)) stop("Row dimension of X and Y are not equal!")
	colnames(X) = paste("X",colnames(X),sep="")	
  }
  
  # Coalesce all data into one data.frame 
  if(!missing(X)) {
    allDataOrig = cbind(Y,D,Z,X)
  } else {
    allDataOrig = cbind(Y,D,Z)
  }
  
  # Fit adjustment model
  ff = terms(Y ~ D + . -1,data=allDataOrig)
  mf = model.frame(ff,allDataOrig)
  
  # Declare Y and D
  Y = as.matrix(mf$Y); colnames(Y) = "Y"
  D = as.matrix(mf$D); colnames(D) = "D"
 
  # Extract Z and X
  allData = sparse.model.matrix(ff,allDataOrig); attr(allData,"assign") = NULL; attr(allData,"contrasts") = NULL
  Zindex = grep(paste("^",colnames(Z),sep="",collapse="|"),colnames(allData))
  Z = allData[,Zindex,drop=FALSE]; colnames(Z) = sub("^Z","",colnames(Z))
  if(!missing(X)) {
    Xindex = grep(paste("^",colnames(X),sep="",collapse="|"),colnames(allData)); X = allData[,Xindex,drop=FALSE]; colnames(X) = sub("^X","",colnames(X))
  }
  # Check to see if there's enough data points after removing missing values
  if(nrow(Y) <= 1) stop("Too much missing data!")

  n = length(Y)
  if(!missing(X)){
    ### clean X and project Y, D, Z
    qrX<-qr(X)
    p = qrrank(qrX)
    if(p==0)
      stop("vector in X are all 0")
    #if(p<ncol(X))
      #X<-qr.Q(qrX)[, 1:p]%*%qrRM(qrX)[1:p, 1:p]
    Yadj = as.matrix(qr.resid(qrX,Y)); Dadj = as.matrix(qr.resid(qrX,D)); Zadj = as.matrix(qr.resid(qrX,Z))
    ### clean Zadj
    ZadjQR = qr(Zadj)
    L = qrrank(ZadjQR)
    if(L==0)
      stop("No useful instrumental variables")
    if(L<ncol(Z)){
      Zadj<-qr.Q(ZadjQR)[, 1:L]%*%qrRM(ZadjQR)[1:L, 1:L]  ### shall we update the Z by doing qr(X, Z)?   
      #qrXZ<-qr(cbind(X, Z))
      #Z<-(qr.Q(qrXZ)[, 1:(p+L)]%*%qrRM(qrXZ)[1:(p+L), 1:(p+L)])[,(p+1):(p+L)]
    }
    ivmodelObject = list(call = match.call(),n=n,L=L,p=p,Y=Y,D=D,Z=Z,X=X,Yadj=Yadj,Dadj=Dadj,Zadj=Zadj, ZadjQR = ZadjQR)

  }else{
    p = 0
    Yadj = Y; Dadj = D; Zadj = as.matrix(Z)
    ### only need to clean Z
    ZadjQR = qr(Zadj)
    L = qrrank(ZadjQR)
    if(L==0)
      stop("No useful instrumental variables")
    if(L<ncol(Z))
      Z<-Zadj<-qr.Q(ZadjQR)[, 1:L]%*%qrRM(ZadjQR)[1:L, 1:L]
	
    ivmodelObject = list(call = match.call(),n=n,L=L,p=p,Y=Y,D=D,Z=Z,X=NA,Yadj=Yadj,Dadj=Dadj,Zadj=Zadj, ZadjQR = ZadjQR)
  }

  class(ivmodelObject) = "ivmodel"

  ivmodelObject$alpha = alpha
  ivmodelObject$beta0 = beta0
  
  ivmodelObject$AR = AR.test(ivmodelObject,beta0=beta0,alpha=alpha)
  ivmodelObject$ARsens = ARsens.test(ivmodelObject,beta0=beta0,alpha=alpha,deltarange=deltarange)
  ivmodelObject$kClass = kClassEst(ivmodelObject,beta0=beta0,alpha=alpha,k=k,heteroSE=heteroSE)
  ivmodelObject$LIML = LIMLEst(ivmodelObject,beta0=beta0,alpha=alpha,heteroSE=heteroSE)  
  ivmodelObject$Fuller = FullerEst(ivmodelObject,beta0=beta0,alpha=alpha,heteroSE=heteroSE)
  ivmodelObject$CLR = CLR.test(ivmodelObject,beta0=beta0,alpha=alpha)

  return(ivmodelObject)
}


coef.ivmodel<-function(object, ...){
  ivmodel<-object
  coefmat <- matrix(NA, ncol=5, nrow=0)
  colnames(coefmat) <- c("k", "Estimate", "Std. Error", "t value", "Pr(>|t|)")
  if(!is.null(ivmodel$kClass)){
    temp<-cbind(as.numeric(rownames(ivmodel$kClass$ci)), ivmodel$kClass$point.est, 
            	ivmodel$kClass$std.err, ivmodel$kClass$test.stat, ivmodel$kClass$p.value)
	rownames(temp) <- rep("k-class", nrow(temp))
    rownames(temp)[temp[,1]==0] <- "OLS"
    rownames(temp)[temp[,1]==1] <- "TSLS"
    coefmat <- rbind(coefmat, temp)
  }
  if(!is.null(ivmodel$LIML)){
    temp<-cbind(ivmodel$LIML$k, ivmodel$LIML$point.est, ivmodel$LIML$std.err,
	            ivmodel$LIML$test.stat, ivmodel$LIML$p.value)
	rownames(temp) <- "LIML"
    coefmat <- rbind(coefmat, temp)
  }
  if(!is.null(ivmodel$Fuller)){
    temp<-cbind(ivmodel$Fuller$k, ivmodel$Fuller$point.est, ivmodel$Fuller$std.err, ivmodel$Fuller$test.stat, ivmodel$Fuller$p.value)
	rownames(temp) <- "Fuller"
    coefmat <- rbind(coefmat, temp)
  }
  
  coefmat<-coefmat[sort(coefmat[,1], index.return=T)$ix,]  
  return(coefmat)
}


summary.ivmodel <- function(object, ...){
  ivmodel<-object
### print formula
  cat("\nCall:\n", paste(deparse(ivmodel$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("sample size: ", ivmodel$n, "\n", sep="")
  
### first stage regression result
  cat(rep("_", 30), "\n")
  cat("\nFirst Stage Regression Result:\n\n")

  SSM=c(sum(qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj)^2))
  SST=sum(ivmodel$Dadj^2)
  SSE=SST-SSM
  DM=ivmodel$L
  DE=ivmodel$n-ivmodel$p-ivmodel$L
  DT=DM+DE
  
  Fstat=SSM/SSE*DE/DM
  pval=1-pf(Fstat, df1=DM, df2=DE)
  RSquare=SSM/SST
  adjRSquare=1-(1-RSquare)*DT/DE
  RMSE=sqrt(SSE/DE)
  
  cat("F=", Fstat, ", df1=", DM, ", df2=", DE, ", p-value is ", format.pval(pval), "\n", sep="")
  cat("R-squared=", RSquare, ",   Adjusted R-squared=", adjRSquare, "\n", sep="")
  cat("Residual standard error: ", RMSE, " on ", DT, " degrees of freedom\n", sep="")
	  
### Sargan test
  if(ivmodel$L>1){

  cat(rep("_", 30), "\n")
  cat("\nSargan Test Result:\n\n")

  TSLS=sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Yadj))/
       sum(ivmodel$Dadj*qr.fitted(ivmodel$ZadjQR, ivmodel$Dadj))
  e=ivmodel$Yadj-ivmodel$Dadj*TSLS
  Sargan=sum(qr.fitted(ivmodel$ZadjQR, e)^2)/(sum(e^2)/length(e))
  pval=1-pchisq(Sargan, df=ivmodel$L-1)
  
  cat("Sargan Test Statistics=", Sargan, ", df=", ivmodel$L-1, ", p-value is ", format.pval(pval), "\n", sep="")
  
  }
  
### print TSLS, kClass, LIML, Fuller
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of k-Class Estimators:\n\n")
  printCoefmat(coef(ivmodel), digits = max(3L, getOption("digits") - 3L))
	  
### print AR, ARsens, CLR
  cat(rep("_", 30), "\n")
  cat("\nAlternative tests for the treatment effect under H_0: beta=", ivmodel$beta0, ".\n",sep = "")
  if(!is.null(ivmodel$AR)){
    cat("\nAnderson-Rubin test:\n")
	cat("F=", ivmodel$AR$Fstat, ", df1=", ivmodel$AR$df[1], ", df2=", 
	    ivmodel$AR$df[2], ", p-value=", format.pval(ivmodel$AR$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$AR$ci.info)
  }
  if(!is.null(ivmodel$ARsens)){
    cat("\n\nSensitivity analysis with deltarange [", ivmodel$ARsens$deltarange[1], 
	    ", ", ivmodel$ARsens$deltarange[2], "]:\n")
	cat("non-central F=", ivmodel$ARsens$ncFstat, ", df1=", ivmodel$ARsens$df[1], 
	    ", df2=", ivmodel$ARsens$df[2], ", ncp=", ivmodel$ARsens$ncp, ", p-value=", format.pval(ivmodel$ARsens$p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$ARsens$ci.info)
  }
  if(!is.null(ivmodel$CLR)){
    cat("\n\nConditional Likelihood Ratio test:\n")
	cat("Test Stat=", ivmodel$CLR$test.stat, ", p-value=", format.pval(ivmodel$CLR$    p.value), "\n", sep="")
	cat(round((1-ivmodel$alpha)*100, digits=1), "percent confidence interval:\n", 
	    ivmodel$CLR$ci.info)
  }
  cat("\n")
}

confint.ivmodel <- function(object, parm,level=NULL, ...){
  ivmodel = object
  alpha=ivmodel$alpha
  if(!is.null(level)){
    if(!is.numeric(level) | level>1 | level<0){
	  print("Wrong input of confidence level!")
	  return()
	}else{
	  alpha=1-level
	}
  }
  temp<-coef(ivmodel)
  result<-matrix(NA, ncol=2, nrow=nrow(temp))
  colnames(result)<-c(paste(round(alpha/2*100, digits=2), "%", sep=""),
                      paste(round((1-alpha/2)*100, digits=2), "%", sep=""))
  rownames(result)<-rownames(temp)
  index<-rownames(temp)=="k-class"
  rownames(result)[index]<-paste("k-class", round(temp[index, 1], digits=2))
  result[,1]<-temp[,2]+temp[,3]*qt(alpha/2, df=ivmodel$n-ivmodel$p-1)
  result[,2]<-temp[,2]+temp[,3]*qt(1-alpha/2, df=ivmodel$n-ivmodel$p-1)

  if(!is.null(ivmodel$AR)){
    if(alpha!=ivmodel$alpha){
	  temp=AR.test(ivmodel, beta0=ivmodel$beta0, alpha=alpha)
	}else{
	  temp=ivmodel$AR
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"AR"
	}else{
	  rownames(temp$ci)<-paste("AR(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  if(!is.null(ivmodel$ARsens)){
    if(alpha!=ivmodel$alpha){
	  temp=ARsens.test(ivmodel, beta0=ivmodel$beta0, alpha=alpha, 
	                   deltarange=ivmodel$ARsens$deltarange)
	}else{
	  temp=ivmodel$ARsens
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"Sensitivity Analysis"
	}else{
	  rownames(temp$ci)<-paste("Sensitivity Analysis(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  if(!is.null(ivmodel$CLR)){
    if(alpha!=ivmodel$alpha){
	  temp=CLR.test(ivmodel, beta0=ivmodel$beta0, alpha=alpha)
	}else{
	  temp=ivmodel$CLR
	}
    if(nrow(temp$ci)==1){
      rownames(temp$ci)<-"CLR"
	}else{
	  rownames(temp$ci)<-paste("CLR(part", 1:nrow(temp$ci), ")", sep="")
	}
    result<-rbind(result, temp$ci)
  }
  
  return(result) 
}



