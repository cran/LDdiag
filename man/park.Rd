\name{park}
\alias{park}
\alias{park.default}
\alias{park.lm}
\alias{park.glm}
\title{Modified Park Test}

\description{
\bold{Modified Park Test} is to identify the potential distribution of the dependent variable. This function first computes tentative parameter
estimates for a generalized linear model based on one's prior belief about the appropriate variance function. Gamma-like regressions
are the most often used option. The linear predictors from the tentative regression are used to get
raw-scale residuals. \cr

The modified Park test is then to regress the squared raw-scale residuals on a constant and the log of linear 
predictors in a generalized linear model with a Gamma distribution and log link. This test also works for linear model. 
For linear model, the dependent variable of squared raw-scale residuals is set to the log scale.\cr

The coefficient on the linear predictor indicates which variance structure is the
most appropriate. \itemize{
   \item[1.] if the coefficient = 0, Gaussian distribution. Variance is constant.
   \item[2.] if the coefficient = 1, Poisson distribution.  Variance is proportional to mean.
   \item[3.] if the coefficient = 2, Gamma distribution. Variance is proportional to square of mean.
   \item[4.] if the coefficient = 3, Inverse Gaussian or Wald distribution. Variance is proportional to cube of mean.
   }

}

\usage{park(object,digits=4)}

\arguments{
    \item{object}{an object created by \code{lm} or \code{glm} function}
    \item{digits}{round to decimals. Default is 4}
      }

\references{
Manning, W.G., Mullahy, J. (2001). Estimating log models: to transform or not to transform? \emph{Journal of Health
Economics}. \bold{20}(4), 461-494

Park, R.E. (1966). Estimation with heteroscedastic error terms. \emph{Econometrica}. \bold{34}, 888}


\examples{

ex <- lm(Days ~ ., data = quine)  #quine is from MASS
park(ex)

quine$Days <- ifelse(quine$Days==0, 1, quine$Days)
ex <- glm(Days ~ ., family = Gamma(link="log"), data = quine)
park(ex)

}

