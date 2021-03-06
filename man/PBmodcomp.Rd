\name{PBmodcomp}
\alias{PBmodcomp}
\alias{PBmodcomp.lm}

\alias{PBmodcomp.merMod}
\alias{getLRT}
\alias{getLRT.lm}

\alias{getLRT.merMod}


\alias{plot.XXmodcomp}


\alias{PBmodcomp.mer}
\alias{getLRT.mer}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Model comparison using parametric bootstrap methods.
}
\description{
  Model comparison of nested models using parametric bootstrap
  methods. Implemented for some commonly applied model types.
}
\usage{
PBmodcomp(largeModel, smallModel, nsim = 1000, ref = NULL, seed=NULL,
          cl = NULL, details = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{largeModel}{
    A model object. Can be a linear mixed effects model or generalized
    linear mixed effects model (as fitted with \code{lmer()} and \code{glmer()}
    function in the \pkg{lme4} package) or a linear normal model or a
    generalized linear model. The \code{largeModel} must be larger than
    \code{smallModel} (see below).
  }
  \item{smallModel}{
    A model of the same type as \code{largeModel} or a restriction
    matrix.  
  }
  \item{nsim}{
    The number of simulations to form the reference distribution.
  }
  \item{ref}{
    Vector containing samples from the reference distribution. If NULL,
    this vector will be generated using PBrefdist(). }
  \item{seed}{
    A seed that will be passed to the simulation of new datasets.
  }
  \item{cl}{
    A vector identifying a cluster; used for
    calculating the reference distribution using several cores. See
    examples below.
  }
  \item{details}{
    The amount of output produced. Mainly relevant for debugging
    purposes.
  }
}
\details{
	The model \code{object} must be fitted with maximum likelihood (i.e. with
	\code{REML=FALSE}). If the object is fitted with restricted maximum
	likelihood (i.e. with
	\code{REML=TRUE}) then the model is refitted with \code{REML=FALSE}
	before the p-values are calculated. Put differently, the user needs
	not worry about this issue.

	Under the fitted hypothesis (i.e. under the fitted small model)
  \code{nsim} samples of the likelihood ratio test statistic (LRT) are
  generetated.

  Then p-values are calculated as follows:

  LRT: Assuming that LRT has a chi-square distribution.

  PBtest: The fraction of simulated LRT-values that are larger or equal
  to the observed LRT value.

  Bartlett: A Bartlett correction is of LRT is calculated from the mean
  of the simulated LRT-values

  Gamma: The reference distribution of LRT is assumed to be a gamma
  distribution with mean and variance determined as the sample mean and
  sample variance of the simulated LRT-values.

  F: The LRT divided by the number of degrees of freedom is assumed to
  be F-distributed, where the denominator degrees of freedom are
  determined by matching the first moment of the reference
  distribution.

}
%\value{
%  ## To be added...
%}
%\references{
%  ## To be added...
%}

\references{ Ulrich Halekoh, S�ren H�jsgaard (2014)., A Kenward-Roger
  Approximation and Parametric Bootstrap Methods for Tests in Linear
  Mixed Models - The R Package pbkrtest., Journal of Statistical
  Software, 58(10), 1-30., \url{http://www.jstatsoft.org/v59/i09/} }

\author{
  Soren Hojsgaard \email{sorenh@math.aau.dk}
}

\note{It can happen that some values of the LRT statistic in the
  reference distribution are negative. When
  this happens one will see that the number of used samples (those where
  the LRT is positive) are reported (this number is smaller than the
  requested number of samples).

  In theory one can not have a negative value of the LRT statistic but
  in practice on can: We speculate that the reason is as follows: We
  simulate data under the small model and fit both the small and the
  large model to the simulated data. Therefore the large model
  represents - by definition - an overfit; the model has superfluous
  parameters in it. Therefore the fit of the two models will for some
  simulated datasets be very similar resulting in similar values of the
  log-likelihood. There is no guarantee that the the log-likelihood for
  the large model in practice always will be larger than for the
  small (convergence problems and other numerical issues can play a role
  here).

  To look further into the problem, one can use the \code{PBrefdist()}
  function for simulating the reference distribution (this reference
  distribution can be provided as input to \code{PBmodcomp()}). Inspection
  sometimes reveals that while many values are negative, they are
  numerically very small. In this case one may try to replace the
  negative values by a small positive value and then invoke
  \code{PBmodcomp()} to get some idea about how strong influence there
  is on the resulting p-values. (The p-values get smaller this way
  compared to the case when only the originally positive values are
  used).
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{KRmodcomp}}
  \code{\link{PBrefdist}}
}
\examples{
data(beets, package="pbkrtest")
head(beets)

## Linear mixed effects model:
sug   <- lmer(sugpct ~ block + sow + harvest + (1|block:harvest), data=beets, REML=FALSE)
sug.h <- update(sug, .~. -harvest)
sug.s <- update(sug, .~. -sow)

anova(sug, sug.h)
PBmodcomp(sug, sug.h, nsim=50)
anova(sug, sug.h)
PBmodcomp(sug, sug.s, nsim=50)

## Linear normal model:
sug <- lm(sugpct ~ block + sow + harvest, data=beets)
sug.h <- update(sug, .~. -harvest)
sug.s <- update(sug, .~. -sow)

anova(sug, sug.h)
PBmodcomp(sug, sug.h, nsim=50)
anova(sug, sug.s)
PBmodcomp(sug, sug.s, nsim=50)

## Generalized linear model
counts    <- c(18,17,15,20,10,20,25,13,12)
outcome   <- gl(3,1,9)
treatment <- gl(3,3)
d.AD      <- data.frame(treatment, outcome, counts)
head(d.AD)
glm.D93   <- glm(counts ~ outcome + treatment, family = poisson())
glm.D93.o <- update(glm.D93, .~. -outcome)
glm.D93.t <- update(glm.D93, .~. -treatment)

anova(glm.D93, glm.D93.o, test="Chisq")
PBmodcomp(glm.D93, glm.D93.o, nsim=50)
anova(glm.D93, glm.D93.t, test="Chisq")
PBmodcomp(glm.D93, glm.D93.t, nsim=50)

## Generalized linear mixed model (it takes a while to fit these)
\dontrun{
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
(gm2 <- update(gm1, .~.-period))
anova(gm1, gm2)
PBmodcomp(gm1, gm2)
}


\dontrun{
(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Days
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
anova(fmLarge, fmSmall)
PBmodcomp(fmLarge, fmSmall)

## The same test using a restriction matrix
L<-cbind(0,1)
PBmodcomp(fmLarge, L)

## Vanilla
PBmodcomp(beet0, beet_no.harv, nsim=1000)

## Simulate reference distribution separately:
refdist <- PBrefdist(beet0, beet_no.harv, nsim=1000)
PBmodcomp(beet0, beet_no.harv, ref=refdist)

## Do computations with multiple processors:
## Number of cores:
(nc <- detectCores())
## Create clusters
cl <- makeCluster(rep("localhost", nc))

## Then do:
PBmodcomp(beet0, beet_no.harv, cl=cl)

## Or in two steps:
refdist <- PBrefdist(beet0, beet_no.harv, nsim=1000, cl=cl)
PBmodcomp(beet0, beet_no.harv, ref=refdist)

## It is recommended to stop the clusters before quitting R:
stopCluster(cl)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{inference}


%\method{PBmodcomp}{lm}(largeModel, smallModel, nsim = 1000, ref = NULL, cl = NULL, details = 0)
%\method{PBmodcomp}{mer}(largeModel, smallModel, nsim = 1000, ref = NULL, cl = NULL, details = 0)
%\method{PBmodcomp}{lmerMod}(largeModel, smallModel, nsim = 1000, ref = NULL, cl = NULL, details = 0)
