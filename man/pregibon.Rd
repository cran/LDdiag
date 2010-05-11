\name{pregibon}
\alias{pregibon}
\alias{pregibon.default}
\alias{pregibon.lm}
\alias{pregibon.glm}
\title{Pregibon Link Diagnostic Test}

\description{
\bold{Pregibon Test for Linearity} was introduced by Daryl Pregibon in 1980. The test is to examine the adequacy of the hypothesized
link used in fitting a linear or generalized linear model.
}

\usage{pregibon(object,digits=4)}

\arguments{
    \item{object}{an object created by \code{lm} or \code{glm} function}
    \item{digits}{round to decimals. Default is 4}
      }

\references{
Pregibon, D.(1980). Goodness of Link Tests for Generalized Linear Models. \emph{Applied Statistics}. \bold{29}, 15-24
}

\examples{

ex <- lm(Days ~ ., data = quine)  #quine is from MASS
pregibon(ex)

quine$Days <- ifelse(quine$Days==0, 1, quine$Days)
ex <- glm(Days ~ ., family = Gamma(link="log"), data = quine)
pregibon(ex)

}

