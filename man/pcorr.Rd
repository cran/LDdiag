\name{pcorr}
\alias{pcorr}
\title{Test for Pearson Correlation}

\description{
If a linear or generalized linear model fits the data well, the predicted values and the residuals should be independent.
\bold{pcorr Test} is to test the pearson correlation between the predicted values and the residuals. 
Small correlation coefficient is highly desired. The 95\% confidence
interval of the correlation coefficient is calculated by Fisher transformation.\cr

}

\usage{pcorr(object,digits=4)}

\arguments{
    \item{object}{an object created by \code{lm} or \code{glm} function}
    \item{digits}{round to decimals. Default is 4}
}


\references{
Fisher, R.A.(1915). Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population. \emph{Biometrika}. \bold{10(4)}, 507-521

Fisher, R.A.(1921). On the probable error of a coefficient of correlation deduced from a small sample. \emph{Metron}. \bold{1(4)}, 3-32
}


\examples{
quine$Days <- ifelse(quine$Days==0, 1, quine$Days) #quine from MASS

fit.lm <- lm(Days ~ Age + Sex, data=quine)
fit.glm <- glm(Days ~ Age + Sex, family=Gamma(link="log"), data=quine)

pcorr(fit.lm)
pcorr(fit.glm)

}

