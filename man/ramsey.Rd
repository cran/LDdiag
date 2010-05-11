\name{ramsey}
\alias{ramsey}

\title{Ramsey RESET Test}

\description{
\bold{Ramsey Regression Equation Specification Error Test (RESET)} is a general specification test
for the linear regression model. This test can also be done by \code{resettest} in \bold{lmtest}.\cr

More specifically, function \code{ramsey} tests whether non-linear combinations of the
estimated values help explain the endogenous variable. The intuition behind the test is that,
if non-linear combinations of the explanatory variables have any power in explaining the endogenous variable, then the model is mis-specified.
}

\usage{ramsey(object,power=4)}

\arguments{
    \item{object}{an object created by \code{lm} or \code{glm} function}
    \item{power}{power of the expected values. Minimum is 3. Default is 4}
      }

\references{
Ramsey, J.B.(1969). Tests for Specification Errors in Classical Linear Least Squares Regression Analysis,
 \emph{Journal of the Royal Statistical Society: Series B}, \bold{31(2)}, 350-371.
}

\examples{

ex <- lm(Days ~ ., data = quine) #quine is from MASS
ramsey(ex)

quine$Days <- ifelse(quine$Days==0, 1, quine$Days)
ex <- glm(Days ~ ., family = Gamma(link="log"), data = quine)
ramsey(ex)

}

