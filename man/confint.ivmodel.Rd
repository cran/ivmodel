\name{confint.ivmodel}
\alias{confint.ivmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{confint.ivmodel
%%  ~~function to do ... ~~
}
\description{This confint methods returns a matrix of two columns, each row represents a confident interval for different IV approaches, which include k-Class, AR (Anderson and Rubin 1949) and CLR (Moreira 2003) estimations. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{confint}{ivmodel}(object, parm, level = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{ivmodel} object
%%     ~~Describe \code{ivmodel} here~~
}
  \item{parm}{ignored for our code}
  \item{level}{the confidence level
%%     ~~Describe \code{level} here~~
}
  \item{\dots}{additional argument(s) for methods.}
}

\value{
A matrix, each row represents a confidence interval for different IV approaches.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Andrews, D. W. K., Moreira, M. J., and Stock, J. H. (2006). Optimal two-side invariant similar tests for instrumental variables regression. \emph{Econometrica} 74, 715-752.\cr
  Moreira, M. J. (2003). A conditional likelihood ratio test for structural models. \emph{Econometrica} 71, 1027-1048.\cr
  Fuller, W. (1977). Some properties of a modification of the limited information estimator. \emph{Econometrica}, 45, 939-953. \cr
  Anderson, T.W. and Rubin, H. (1949), Estimation of the parameters of a single equation in a complete system of stochastic equations, Annals of Mathematical Statistics, 20, 46-63.
%% ~put references to the literature/web site here ~
}
\author{Yang Jiang, Hyunseung Kang, and Dylan Small
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
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661", "reg662", 
        "reg663", "reg664", "reg665", "reg666", "reg667", "reg668", "smsa66")
X=card.data[,Xname]
foo = ivmodel(Y=Y,D=D,Z=Z,X=X)
confint(foo)
}
