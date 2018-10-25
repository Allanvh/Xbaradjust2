#' LadjSp
#'
#' Essa funcao  ajusta o L
#'
#'
#'
#' @param x o \code{x} da questao
#'
#' @return L ajustado
#'
#' @examples
#'
#' LadjSp(370,25,5)
#'
#' @export cubature
#' @import cubature
#'

library(cubature)

LadjSp <- function(ARL0nom,m,n) {

  alpha <- 1/ARL0nom

  if (alpha<0.001 || m < 15 || n < 3) {
    print(paste("Please revise your entries according to the following conditions:"))
    print(paste("The nominal in-control ARL must be equal or smaller than 1000"))
    print(paste("The number (m) of Phase I Samples must be equal or larger than 15"))
    print(paste("The size (n) of each Phase I Samples must be equal or larger than 3"))
  }
  else {

    secant <- function(fun, x0, x1, tol=1e-10, niter=100000){
      for ( i in 1:niter ) {
        funx1 <- fun(x1)
        funx0 <- fun(x0)
        x2 <- ( (x0*funx1) - (x1*funx0) )/( funx1 - funx0 )
        funx2 <- fun(x2)
        if (abs(funx2) < tol) {
          return(x2)
        }
        if (funx2 < 0)
          x1 <- x2
        else
          x0 <- x2
      }
      stop("exceeded allowed number of iteractions")
    }

    ARL0XbarSp <- function (alpha1,m,n) {

      CARL <- function (U) {
        a <- 1/(1 - pnorm(((1/sqrt(m))*qnorm(U[1],0,1))+((-1*qnorm(alpha1/2))*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1) + pnorm(((1/sqrt(m))*qnorm(U[1],0,1))-((-1*qnorm(alpha1/2))*sqrt(qchisq(U[2],m*(n-1))/(m*(n-1)))),0,1))
        return(a)
      }
      a <- adaptIntegrate(CARL, lowerLimit = c(0, 0), upperLimit = c(1, 1))$integral
      return (a)
    }

    ARLfunc <- function (alphaf) {
      a <- ARL0XbarSp(alphaf,m,n) - (1/alpha)
      return(a)
    }

    cat("This may take several minutes. Please, wait... ")
    s <- secant(ARLfunc,alpha*0.5,alpha*1.5)

    L <- (-1*qnorm(s/2))

    b <- round(ARL0XbarSp(s,m,n),3)

    print(paste("End of calculations. See results below:"))
    print(paste("L = ", L))
    print(paste("This is the Limit Factor that generates an in-control ARL equal to the specified one for a given m and n for the Xbar chart with Sp estimator"))
    print(paste("In-Control ARL = E[In-Control RL] = E[In-Control CARL] ", b))
    return(L)
  }
}

