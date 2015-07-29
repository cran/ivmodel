\name{para}
\alias{para}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameter Estimation from Ivmodel
%%  ~~function to do ... ~~
}
\description{\code{para} computes the estimation of several parameters for the \code{ivmodel} object, which could be used in the power calculation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
para(ivmodel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ivmodel}{\code{ivmodel} object.
%%     ~~Describe \code{ivmodel} here~~
}
}
\details{\code{para} computes the coefficients of 1st and 2nd stage regression (gamma and beta). It also computes the covariance matrix of the error term of 1st and 2nd stage. (sigmau, sigmav, and rho)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \code{para} returns a list containing the following components
    
	\item{gamma}{The coefficient of IV in first stage, calculated by linear regression}
	\item{beta}{The TSLS estimator of the exposure effect}
	\item{sigmau}{Standard deviation of potential outcome under control (structural error for y).}
	\item{sigmav}{Standard deviation of error from regressing treatment on instruments}
	\item{rho}{Correlation between u (potential outcome under control) and v (error from regressing treatment on instrument).}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Yang Jiang, Hyunseung Kang, Dylan Small
%%  ~~who you are~~
}

\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", 
        "reg662", "reg663", "reg664", "reg665", "reg666", "reg667", 
		"reg668", "smsa66")
X=card.data[,Xname]
cardfit=ivmodel(Y=Y, D=D, Z=Z, X=X)
para(cardfit)
}

