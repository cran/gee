\name{gee}
\alias{gee}
\alias{print.gee}
\alias{summary.gee}
\alias{print.summary.gee}
\title{
  Function to solve a Generalized Estimation Equation Model
}
\description{
  Produces an object of class \code{"gee"} which is a Generalized Estimation 
  Equation fit of the data.
}
\usage{
gee(formula, id,
    data, subset, na.action,
    R = NA, b = NA,
    tol = 0.001, maxiter = 25,
    family = gaussian, corstr = "independence",
    Mv = 1, silent = TRUE, contrasts = NULL,
    scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE)
}
\arguments{
\item{formula}{
  a formula expression as for other regression models, of the form
  \code{response ~ predictors}. See the documentation of
  \code{\link{lm}} and \code{\link{formula}} for details.
}
\item{id}{
  a vector which identifies the clusters.  The length of \code{id} should be
  the same as the number of observations.  Data are assumed to be sorted
  so that observations on a cluster are contiguous rows for all entities
  in the formula.
}
\item{data}{
  an optional data frame in which to interpret the variables occurring
  in the \code{formula}, along with the \code{id} and \code{n} variables.
}
\item{subset}{
  expression saying which subset of the rows of the data should be used
  in the fit.  This can be a logical vector (which is replicated to have
  length equal to the number of observations), or a numeric vector
  indicating which observation numbers are to be included, or a
  character vector of the row names to be included.
  All observations are included by default.
}
\item{na.action}{
  a function to filter missing data.  For \code{gee} only \code{na.omit}
  should be used here.
}
\item{R}{
  a square matrix of dimension maximum cluster size containing the user
  specified correlation.  This is only appropriate if \code{corstr = "fixed"}.
}
\item{b}{
  an initial estimate for the parameters.
}
\item{tol}{
  the tolerance used in the fitting algorithm.
}
\item{maxiter}{
  the maximum number of iterations.
}
\item{family}{
  a \code{family} object: a list of functions and expressions for
  defining link and variance functions.  Families supported
  in \code{gee} are \code{gaussian}, \code{binomial}, \code{poisson},
  \code{Gamma}, and \code{quasi};
  see the \code{\link{glm}} and \code{\link{family}} documentation.
  Some links are not currently available: \code{1/mu^2} and \code{sqrt} have
  not been hard-coded in the cgee engine at present.
  The inverse gaussian variance function is not available.
  All combinations of remaining functions can be obtained
  either by family selection or by the use of \code{quasi}.
}
\item{corstr}{
  a character string specifying the correlation structure.
  The following are permitted:
     \code{"independence"},
     \code{"fixed"},
     \code{"stat_M_dep"},
     \code{"non_stat_M_dep"},
     \code{"exchangeable"},
     \code{"AR-M"} and
     \code{"unstructured"}
}
\item{Mv}{
  When the corstr is \code{"stat_M_dep"}, \code{"non_stat_M_dep"},
  or \code{"AR-M"} then \code{Mv} must be specified.
}
\item{silent}{
  a logical variable controlling whether parameter estimates at each
  iteration are printed.
}
\item{contrasts}{
  a list giving contrasts for some or all of the factors appearing
  in the model formula.  The elements of the list should have the
  same name as the variable and should be either a contrast matrix
  (specifically, any full-rank matrix with as many rows as there are
  levels in the factor), or else a function to compute such a matrix
  given the number of levels.
}
\item{scale.fix}{
  a logical variable; if true, the scale parameter is fixed at
  the value of \code{scale.value}.
}
\item{scale.value}{
  numeric variable giving the value to which the scale parameter
  should be fixed; used only if \code{scale.fix == TRUE}.
}
\item{v4.4compat}{
  logical variable requesting compatibility of correlation
  parameter estimates with previous versions; the current
  version revises to be more faithful to the Liang and Zeger (1986)
  proposals (compatible with the Groemping SAS macro, version 2.03)
}}
\value{
  An object of class \code{"gee"} representing the fit.
}
\section{Side Effects}{
  Offsets must be specified in the model formula, as in \code{\link{glm}}.
}
\details{
  Though input data need not be sorted by the variable
  named \code{"id"}, the program
  will interpret physically contiguous records possessing the
  same value of \code{id} as members of the same cluster.  Thus it
  is possible to use the following vector as an \code{id} vector
  to discriminate 4 clusters of size 4: 
  \code{c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1)}.
}
\note{
  This is version 4.8 of this user documentation file, revised
  98/01/27.  The assistance of Dr B Ripley is gratefully acknowledged.
}
\references{
  Liang, K.Y. and Zeger, S.L. (1986)
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73} 13--22. 

  Zeger, S.L. and Liang, K.Y. (1986)
  Longitudinal data analysis for discrete and continuous outcomes.
  \emph{Biometrics}, \bold{42} 121--130.
}
\seealso{
  \code{\link{glm}}, \code{\link{lm}}, \code{\link{formula}}.
}
\examples{
data(warpbreaks)
## marginal analysis of random effects model for wool
summary(gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="exchangeable"))
## test for serial correlation in blocks
summary(gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="AR-M", Mv=1))

if(require(MASS)) {
data(OME)
## not fully appropriate link for these data.
(fm <- gee(cbind(Correct, Trials-Correct) ~ Loud + Age + OME, id = ID,
           data = OME, family = binomial, corstr = "exchangeable"))
summary(fm)
}}
\keyword{nonlinear}
