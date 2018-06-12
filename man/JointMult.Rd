\name{JointMult}
\alias{JointMult}
\title{
Estimation of multi-dimensional joint models 
}
\description{
This function fits multi-dimensional joint models. Each dimension is
modelled as a multivariate mixed model that handles continuous longitudinal outcomes (Gaussian or non-Gaussian, curvilinear). 
For each dimension, it is assumed that all the outcomes measure the same
underlying latent process defined as their common factor, and each
outcome is related to this latent common factor by a specific
parameterized link function.
The events are modelled with degradation processes, defined as linear
combinations of the different dimensions. One or two events only can be
considered. The model is defined is a competing risk setting, so only
the first entcountered event is considered.
}
\usage{
JointMult(Y, D, data, var.time, RE = "block-diag", BM = "diag", 
    B, posfix, maxiter = 0, nproc = 1, verbose = FALSE, 
    file = "", pred = FALSE, breaks = NULL, delayed = TRUE,  
    eps = c(1e-04, 1e-04, 0.001)) 

}
\arguments{
  \item{Y}{
a list of \code{multlcmm} objects. Each \code{multlcmm} object defines
the outcomes and longitudinal structures of one dimension of the
model. The number of outcomes, covariates or random effects can differ
between dimensions.
}
  \item{D}{
    a list of two-sided formula defining the event part of the
    model. The left side should either \code{longDiag(time0, time, indicator)}
    with \code{time0} the variable indicating entry time (which will only be
    used for delayed entry), \code{time} the repeated measure times,
    \code{indicator} the repeated event indicator
    or \code{survBreak(T0, T, Event)} with \code{T0} the entry time,
    \code{T} the event time and \code{Event} the event indicator.
}
  \item{var.time}{
a character vector indicating the name of the time variable of each
dimension. The scales of these
different time variables should be the same.
}
\item{RE}{
  an indicator of the random effect structure between dimensions. Should
  be either "block-diag" for independant random effects between
  dimensions (the internal structure being difined in the \code{multlcmm}
  object) or "full" for correlated random effects.
}
  \item{BM}{
in the case where Brownian motions are included in the
\code{multlcmm} objects, an indicator of the correlation of the Brownian
motions. Should be "diag" for independance, or "full" for correlated
Brownian motions.
}
  \item{breaks}{
optional vector specifying the break points in the case where the event
time is discretized in time intervals. 
}
\item{delayed}{
  logical vector indicating, for each event, if delayed entry should be
  accounted for.
}
  \item{B}{
optional specification for the initial values for the parameters.
}
\item{posfix}{
  Optional vector specifying the indices in vector B of the parameters that
  should not be estimated. By default, all parameters are estimated.
}
  \item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=?.
  }
  \item{eps}{
    optional thresholds for the convergence criteria. Default to 0.0001
    for the parameters stability, to 0.0001 for the log-likelihood
    stability, and to 0.001 for the criterion based on second derivatives.
}
\item{nproc}{
  optional integer indicating the number of processors to be used for
  parallel computation. Default to 1, the algorithm runs sequentially.
}
  \item{verbose}{
  logical indicating if information about computation should be
  reported. Default to FALSE.
}
\item{file}{
 optional filename in which to report information about computation (if
 \code{verbose=TRUE}).
}
  \item{pred}{
 logical indicating if subject-specific predictions should be
 computed. Default to FALSE.
}
}

%\details{
%% ordre des parametres, var dep du temps sur seuil, factor(t), ...
%}


\value{
  \item{istop}{convergence status. 1 if the model converged properly, 2
  if the maximum number of iterations is reached without convergence, >2
  if an error occured.}
  \item{ni}{number of iterations}
  \item{loglik}{log-likelihood of the model}
  \item{b}{vector of estimated parameters}
  \item{v}{variance matrix of the estimated parameters}
  \item{convcrit}{convergence criteria at stop point}
  \item{time}{estimation time}
  \item{nproc}{number of processors}
  \item{bopt}{vector of estimated an fixed parameters}
  \item{nef}{number of fixed effects paramters}
  \item{ncontr}{number of contrasts parameters}
  \item{nea}{number of random effects variables}
  \item{nvc}{number of random effects parameters}
  \item{idiag}{indicator of correlation of the inner random effects}
  \item{ntr}{number of link function parameters}
  \item{ntrtot}{number of link function parameters}
  \item{ncor}{number of Brownian motion parameters}
  \item{nalea}{number of outcome-specific random effects parameters}
  \item{ny}{number of longitudinal outcomes}
  \item{link}{type of link functions}
  \item{nodes}{nodes for the link functions}
  \item{nRE}{number of between correlation parameters for the random effects}
  \item{nRM}{number of between correlation for the Brownian motions}
  \item{varD}{time independent covariates in the event model}
  \item{vardept}{time dependent covariates in the event model}
  \item{nvarD}{number of time independent covariates in the event model}
  \item{ndept}{number of time dependent covariates in the event model}
  \item{idcause}{indicator of presence of the covariates in the event model}
  \item{call}{the model's call}
  \item{fix}{indicator of fixed parameter}
  \item{chol}{indicator of Cholesky transformation}
  \item{Ynames}{name of the longitudinal outcomes}
  \item{Xnames}{name of the covariates in the longitudinal part}
  \item{ns}{number of subjects}
  \item{nbevt}{number of observed events}
  \item{nbmes}{mean number of measurement}
  \item{entreRetard}{indicator of left truncation}
  \item{discretise}{indicator of dicretization}
  \item{breaks}{list of break points}
  \item{VRE}{variance-covariance matrix of the random effetcs}
  \item{corRE}{correlation matrix of the random effetcs}
  \item{mod}{list of pdated \code{multlcmm} objects}
}


\author{
Cecile Proust-Lima and Viviane Philipps
}


%% \references{

%% }


%% \seealso{
%% }

%% \examples{
%% \dontrun{
%% }
%% }

