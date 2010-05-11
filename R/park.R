#  File library/LDdiag/R/Copas.R/park.R
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


park <- function(object,digits=4) UseMethod("park")

park.default <- function(object,digits=4) {
  scip <- options()$scipen
  options(scipen=4)
  if(as.character(object$call)[1]=="lm") park.lm(object)
  else  if(as.character(object$call)[1]=="glm") park.glm(object,digits)
  else cat("The function call doesn't support yet.\n")
  options(scipen=scip)
}

park.lm <- function(object, digits=4) {
  xb <- fitted(object)
  X <- cbind(1,log(xb))
  residsq <- log(object$residuals^2)
  fits <- lm.fit(X,residsq)
  Park <- fits$coefficients[2]
  Park <- as.numeric(Park)
  cat("\n  Park =", Park)
  thresh <- round(Park)
  if (thresh==0) {cat(". Gaussian distribution might be appropriate.\n\n")
  } else if(thresh==1) { cat(". Poisson distribution might be appropriate.\n\n")
  } else if(thresh==2) { cat(". Gamma distribution might be appropriate.\n\n")
  } else if(thresh==3) { cat(". Inverse Gaussian or Wald distribution might be appropriate.\n\n")
  } else {cat(". Unknown distribution.\n\n")}
}

park.glm <- function(object, digits=4){
  testf <- function(object) { #this function modified from summary.glm()
    df.r <- object$df.residual
   	dispersion <-
  	if(object$family$family %in% c("poisson", "binomial"))  1
  	else if(df.r > 0) {
    sum(object$residuals^2)/ df.r
  	}
    p <- object$rank
    p1 <- 1:p
    Qr <- object$qr
    coef.p <- object$coefficients[Qr$pivot[p1]]
    covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
    covmat <- dispersion*covmat.unscaled
    var.cf <- diag(covmat)
      s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    pvalue <- 2*pt(-abs(tvalue), df.r)
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    return(coef.table)
    }
  fm <- object$'family'
  xb <- fitted(object)
  X <- cbind(1,log(xb))
  residsq <- object$residuals^2
  fits <- glm.fit(X,residsq, family=fm)
  ans <- testf(fits)
  Park <- ans[2,1]
  cat("\n  Park =", Park)
  thresh <- round(Park)
  if (thresh==0)       { cat(". Gaussian distribution might be appropriate.\n\n")
  } else if(thresh==1) { cat(". Poisson distribution might be appropriate.\n\n")
  } else if(thresh==2) { cat(". Gamma distribution might be appropriate.\n\n")
  } else if(thresh==3) { cat(". Inverse Gaussian or Wald distribution might be appropriate.\n\n")
  } else {cat(". Unknown distribution.\n\n") }
}


