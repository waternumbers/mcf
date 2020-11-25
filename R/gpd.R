#' Generalised (Type II) pareto distribution
#'
#' Functions for evaluating the density, cumulative density and quantilies of the GPD along with random number generation.
#' 
#' @param a,b the shape and scale parameters of the GPD
#' @param x vector of samples
#' @param q vector of quantiles
#' @param n number of samples to draw
#' @param log Logical if TRUE returns log of density
#'
#' @details For \eqn{b\ge0} with \eqn{x\ge0} when \eqn{c\le0} or
#' \eqn{0\lex\le\frac{b/a}} when \eqn{c\ge0}, the density function is given by
#' \deqn{ f(x)=\frac{1}{b}(1-\frac{a}{b}x)^{\frac{1}{a}-1} }
#' For \eqn{a} not equal to \eqn{0} the corresponding distribution function is
#' \deqn{ F(x) = 1-(1-\frac{a}{b}x)^{\frac{1}{a}} }
#' while for \eqn{a=0}
#' \deqn{F(x) = 1-\exp(1-\frac{x}{b})}
#'
#' @examples
#' x <- rgpd(3,c(-1,0,1),2)
#' p <- pgpd(x,c(-1,0,1),2)
#' xx <- qgpd(p,c(-1,0,1),2)
#' stopifnot( all(x==xx) )
#' d <- dgpd(x,c(-1,0,1),2)
#' 
#' @name gpd
NULL

#' @rdname gpd
#' @export
dgpd <- function(x,a,b,log=FALSE){
    a <- rep_len(a,length(x))
    b <- rep_len(b,length(x))
    x[x<0 | b<=0] <- NA
    idx <- a!=0
    ## compute log density
    x[!idx] <- -log(b[!idx]) ##case a==0
    x[idx] <- -log(b[idx]) + ((1/a[idx])-1)*log(1-a[idx]*x[idx]/b[idx]) ## case a!=0
    if(!log){
        x <- exp(x)
    }
    return(x)
}

#' @rdname gpd
#' @export
pgpd <- function(x,a,b){
    a <- rep_len(a,length(x))
    b <- rep_len(b,length(x))
    x[x<0 | b<=0] <- NA
    idx <- a!=0
    x[idx] <- 1 - ( (1 - a[idx]*x[idx]/b[idx])^(1/a[idx]))
    x[!idx] <- 1 - exp( 1-(x[!idx]/b[!idx]) )
    return(x)
}


#' @rdname gpd
#' @export
qgpd <- function(q,a,b){
    a <- rep_len(a,length(q))
    b <- rep_len(b,length(q))
    q[q<0 | q>1 | b<=0] <- NA
    idx <- a!=0
    q[idx] <- b[idx]*( 1 - (1-q[idx])^a[idx] )/a[idx]
    q[!idx] <- b[!idx]*(1-log(1-q[!idx]))
    return(q)
}

#' @rdname gpd
#' @export
rgpd <- function(n,a,b){
    return( qgpd(runif(n),a,b) )
}
