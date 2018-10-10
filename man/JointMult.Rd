\name{JointMult}
\alias{JointMult}
\title{
Estimation of a joint model for multivariate latent processes and clinical endpoints using an exact likelihood inference method
}
\description{
This function fits a joint model for multivariate longitudinal markers (possibly summarized into latent processes) and clinical endpoints. The estimation is performed within the Maximum Likelihood Framework benefiting from an exact likelihood formulation. 
More details are given below. 

Each dimension (constituted of a unique longitudinal marker or of several longitudinal markers measuring the same underlying latent process) is modeled according to a mixed model that handles continuous (Gaussian or non-Gaussian, curvilinear) data. The technique for each dimension follows the methodology of lcmm and multlcmm functions (lcmm package). All the dimensions can be correlated through correlated random effects. 
The dimensions are linked with one or two competing clinical endpoints through a joint degradation process model. Specifically, each clinical endpoint is defined as a binary repeated endpoint either measured at various visit times or derived from a time-to-event discretized into time intervals. The model assumes that the clinical endpoint reaches 1 when its underlying degradation process becomes above a threshold to estimate, this degradation process being defined as a linear combination of the longitudinal dimensions. 

The methodology is fully described in :
Proust-Lima C, Philipps V, Dartigues JF (2018). A joint model for multiple dynamic processes and clinical endpoints: application to Alzheimer’s disease. https://arxiv.org/abs/1803.10043
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
the outcomes and longitudinal structure of one dimension of the
model. The number of outcomes, covariates or random effects can differ
between dimensions.
}
  \item{D}{
    a list of two-sided formula defining the event part of the
    model. The left side should be either \code{longDiag(time0, time, indicator)}
    with \code{time0} the variable indicating entry time (which will only be
    used for delayed entry), \code{time} the repeated measurement times,
    \code{indicator} the repeated event indicator at each time
    or \code{survBreak(T0, T, Event)} with \code{T0} the entry time,
    \code{T} the event time and \code{Event} the event indicator.
}
  \item{var.time}{
a character vector indicating the name of the time variable of each
dimension. The scales of these different time variables should be the same.
}
\item{RE}{
  an indicator of the random effect structure between dimensions. Should
  be either "block-diag" for independent random effects between
  dimensions (the internal structure being defined in the \code{multlcmm}
  object) or "full" for correlated random effects between dimensions.
}
  \item{BM}{
in the case where Brownian motions are included in the
\code{multlcmm} objects, an indicator of the correlation of the Brownian
motions. Should be "diag" for independence, or "full" for correlated
Brownian motions between dimensions.
}
  \item{breaks}{
optional vector specifying the break points in the case where the event
time is discretized into time intervals. 
}
\item{delayed}{
  logical vector indicating, for each event, if delayed entry should be
  accounted for.
}
  \item{B}{
optional specification for the initial values for the parameters.
}
\item{posfix}{
  optional vector specifying the indices in vector B of the parameters that
  should not be estimated. By default, all parameters are estimated.
}
  \item{maxiter}{optional maximum number of iterations for the Marquardt iterative algorithm. By default, maxiter=XXX.
  }
  \item{eps}{
    optional thresholds for the convergence criteria. Default is set to 0.0001
    for the parameters stability, to 0.0001 for the log-likelihood
    stability, and to 0.001 for the criterion based on second derivatives.
}
\item{nproc}{
  optional integer indicating the number of processors to be used for
  parallel computation. Default is set to 1 (i.e., the algorithm runs sequentially).
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
 computed. Default is set to FALSE.
}
}

%% \details{
%% }


\value{
  \item{istop}{convergence status: 1 if the model converged properly, 2
  if the maximum number of iterations was reached without convergence, >2
  if an error occurred.}
  \item{ni}{number of iterations}
  \item{loglik}{log-likelihood of the model}
  \item{b}{vector of estimated parameters}
  \item{v}{estimated variance matrix of the estimated parameters}
  \item{convcrit}{convergence criteria at stop point}
  \item{time}{estimation time}
  \item{nproc}{number of processors}
  \item{bopt}{total vector of estimated and fixed parameters}
  \item{nef}{number of fixed effects parameters}
  \item{ncontr}{number of contrasts parameters}
  \item{nea}{number of random effects variables}
  \item{nvc}{number of random effects parameters}
  \item{idiag}{indicator of intra dimension correlation between random effects}
  \item{ntr}{number of link function parameters}
  \item{ntrtot}{number of link function parameters}
  \item{ncor}{number of Brownian motion parameters}
  \item{nalea}{number of outcome-specific random effects parameters}
  \item{ny}{number of longitudinal outcomes}
  \item{link}{type of link functions}
  \item{nodes}{nodes for the link functions}
  \item{nRE}{number of correlation parameters for the random effects between dimensions}
  \item{nRM}{number of correlation for the Brownian motions between dimensions}
  \item{varD}{time independent covariates in the event model}
  \item{vardept}{time dependent covariates in the event model}
  \item{nvarD}{number of time independent covariates in the event model}
  \item{ndept}{number of time dependent covariates in the event model}
  \item{idcause}{indicator of presence of the covariates in the event model}
  \item{call}{the model's call}
  \item{fix}{indicator of fixed parameters}
  \item{chol}{indicator of Cholesky transformation}
  \item{Ynames}{name of the longitudinal outcomes}
  \item{Xnames}{name of the covariates in the longitudinal part}
  \item{ns}{number of subjects}
  \item{nbevt}{number of observed events}
  \item{nbmes}{mean number of measurement}
  \item{entreRetard}{indicator of left truncation}
  \item{discretise}{indicator of dicretization}
  \item{breaks}{list of break points}
  \item{VRE}{variance-covariance matrix of the random effects}
  \item{corRE}{correlation matrix of the random effects}
  \item{mod}{list of updated \code{multlcmm} objects}
}


\author{
Cecile Proust-Lima and Viviane Philipps
}


\references{
Proust-Lima C, Philipps V, Dartigues JF (2018). A joint model for multiple dynamic processes and clinical endpoints: application to Alzheimer’s disease. https://arxiv.org/abs/1803.10043
}


%% \seealso{
%% }

%% \examples{
%% \dontrun{
%% }
%% }

