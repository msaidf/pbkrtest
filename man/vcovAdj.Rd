\name{vcovAdj}
\alias{vcovAdj}
\alias{vcovAdj.lmerMod}
\alias{vcovAdj_internal}
\alias{vcovAdj0}
\alias{vcovAdj2}

\alias{vcovAdj.mer}


\alias{LMM_Sigma_G}

\alias{get_SigmaG}
\alias{get_SigmaG.lmerMod}

\alias{get_SigmaG.mer}

\title{
	Ajusted covariance matrix for linear mixed models according to Kenward
	and Roger
}
\description{
  Kenward and Roger (1997) describbe an improved small sample approximation to
  the covariance matrix estimate of the fixed parameters in a linear mixed model.
}
\usage{
vcovAdj(object, details=0)
LMM_Sigma_G(object, details=0)
}
	
\arguments{
  \item{object}{An \code{lmer} model}
  \item{details}{If larger than 0 some timing details are printed.}
}

\value{
  \item{phiA}{the estimated covariance matrix, this has attributed P, a
  list of matrices used in \code{KR_adjust} and the estimated matrix  W
  of the variances of the covariance parameters of the random effetcs}
	
  \item{SigmaG}{list: Sigma: the covariance matrix of Y; G: the G
		matrices that sum up to Sigma; n.ggamma: the number (called M in the
		article) of G matrices)
	} 
}

\references{
  Ulrich Halekoh, S�ren H�jsgaard (2014).,
  A Kenward-Roger Approximation and Parametric Bootstrap Methods for
  Tests in Linear Mixed Models - The R Package pbkrtest., 
  Journal of Statistical Software, 58(10), 1-30., \url{http://www.jstatsoft.org/v59/i09/}

  Kenward, M. G. and Roger, J. H. (1997), \emph{Small Sample Inference 
    for Fixed Effects from Restricted Maximum Likelihood}, Biometrics
  53: 983-997. 
}

\note{
	If $N$ is the number of observations, then the \code{vcovAdj()}
	function involves inversion of an $N x N$ matrix, so the
	computations can be relatively slow. 
}

\author{
  Ulrich Halekoh \email{uhalekoh@health.sdu.dk}, Soren Hojsgaard \email{sorenh@math.aau.dk}
}


\seealso{
  \code{\link{getKR}}
  \code{\link{KRmodcomp}}
  \code{\link{lmer}}
  \code{\link{PBmodcomp}}
  \code{\link{vcovAdj}}
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

## Here the adjusted and unadjusted covariance matrices are identical,
## but that is not generally the case
v1 <- vcov(fm1)
v2 <- vcovAdj(fm1,detail=0)
v2 / v1

## For comparison, an alternative estimate of the variance-covariance
## matrix is based on parametric bootstrap (and this is easily
## parallelized): 

\dontrun{
nsim <- 100
sim <- simulate(fm.ml, nsim)
B <- lapply(sim, function(newy) try(fixef(refit(fm.ml, newresp=newy))))
B <- do.call(rbind, B)
v3 <- cov.wt(B)$cov
v2/v1
v3/v1
}


}

\keyword{function}
