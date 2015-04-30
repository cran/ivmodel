\name{coef.ivmodel}
\alias{coef.ivmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{coef.ivmodel
%%  ~~function to do ... ~~
}
\description{This coef methods returns the point estimation, standard error, test statistic and p value for all specified k-Class estimation from an \code{ivmodel} object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{coef}{ivmodel}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{ivmodel} object
%%     ~~Describe \code{ivmodel} here~~
}
  \item{\dots}{additional coef arguments}
}

\value{
A matrix summarises all k-Class estimations.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Yang Jiang, Hyunseung Kang, and Dylan Small
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also \code{\link{ivmodel}} for details on the instrumental variables model.
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
coef(foo)

}
