\name{Copas}
\alias{Copas}
\alias{Copas.lm}
\alias{Copas.glm}
\title{Test for Over-fitting Using Cross Validation}

\description{
\bold{Copas Test} is processed as below \itemize{
\item[1.] Randomly split data into two subsets A and B (eg,50\% versus 50\%, 67\% versus 33\%) for cross validation. This function chooses equal sizes as default.
\item[2.] Fit a linear or generalized linear model \code{Y ~ X*betas} on subset A, retain the coefficients (\code{betas}).
\item[3.] Predict Yhat on subset B using the coefficients (\code{betas}) obtained from step 2.
\item[4.] Fit a linear model of true values (Y) versus predicted values (Yhat) \code{Y ~ alpha + beta*Yhat} on subset B.
\item[5.] Test null hypothesis (\code{H0:alpha=0 and beta=1}) from step 4.
\item[6.] Repeat step 1 to step 5 for multiple times. This function chooses 1000 times as default.
}

The \% of times null hypothesis being rejected is reported. High \% of rejecting times indicates an
 over-fitting problem might exist.\cr

}

\usage{Copas(obj,rate=0.5, ssize=1000)}

\arguments{
    \item{obj}{an object created by \code{lm} or \code{glm} function}
    \item{rate}{proportion of the population, must between 0 and 1. Default is 0.5}
    \item{ssize}{times of cross validation. Default is 1000}
}

\references{
Copas, J.B. (1983). Regression, prediction and shrinkage. \emph{Journal of the Royal Statistical Society: Series B}. \bold{45}, 311-354
}

\examples{
quine$Days <- ifelse(quine$Days==0, 1, quine$Days) #quine from MASS

fit.lm <- lm(Days ~ Age + Sex, data=quine)
fit.glm <- glm(Days ~ Age + Sex, family=Gamma(link="log"), data=quine)

Copas(fit.lm, ssize=1000)
Copas(fit.glm, ssize=1000)

}

