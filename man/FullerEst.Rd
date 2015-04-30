\name{FullerEst}
\alias{FullerEst}
\title{Fuller-k Estimator}
\description{\code{FullerEst} computes the Fuller-k (Fuller 1977) estimate for the \code{ivmodel} object.}
\usage{FullerEst(ivmodel, beta0 = 0, alpha = 0.05, b = 1, heteroSE = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ivmodel}{\code{ivmodel} object}
  \item{beta0}{Null value \eqn{\beta_0} for testing null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}. Default is 0}
  \item{alpha}{The significance level for hypothesis testing. Default is 0.05}
  \item{b}{Positive constant \eqn{b} in Fuller-k estimator. Default is 1}
  \item{heteroSE}{Should heteroscedastic-robust standard errors be used? Default is FALSE}
}
\details{
\code{FullerEst} computes the Fuller-k estimate for the instrumental variables model in \code{ivmodel}, specifically for the parameter \eqn{beta}. The computation uses \code{kClassEst} with the value of \eqn{k = k_{LIML} - b/(n - L - p)}. It generates a point estimate, a standard error associated with the point estimate, a test statistic and a p value under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel} along with a \eqn{1-\alpha} confidence interval.}
\value{
    \code{FullerEst} returns a list containing the following components
    
	\item{k}{The k value used when computing the Fuller estimate with the k-Class estimator}
	\item{point.est}{Point estimate of \eqn{\beta}}
    \item{std.err}{Standard error of the estimate}
	\item{test.stat}{The value of the test statistic for testing the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
	\item{p.value}{The p value of the test under the null hypothesis \eqn{H_0: \beta = \beta_0} in \code{ivmodel}}
	\item{ci}{A matrix of one row by two columns specifying the confidence interval associated with the Fuller estimator}
}
\references{
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953.  
}
\author{Yang Jiang, Hyunseung Kang, Dylan Small}
\seealso{
See also \code{\link{ivmodel}} for details on the instrumental variables model. See also \code{\link{kClassEst}} for more information about the k-Class estimator.
}
\examples{
data(card.data)
Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,c("nearc4","nearc2")]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", "reg662", 
        "reg663", "reg664", "reg665", "reg666", "reg667", "reg668", "smsa66")
X=card.data[,Xname]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
FullerEst(card.model2IV,alpha=0.01)
}
\keyword{Fuller-k estimator}
\keyword{Instrumental variables}% __ONLY ONE__ keyword per line
