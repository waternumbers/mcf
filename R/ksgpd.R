#' Empirical distribution with kernel smoothing for bulk and Generalised (Type II) pareto distribution upper tail
#'
#' Functions for evaluating the density, cumulative density and quantilies of the distribution.
#'
#' @param dst A ksgpd object
#' @param numBreak Maximum number of breaks to consider. The values are the numBreak highest values in the sample
#' @param skipBreak Miss the skipBreak largest values out when testing break points
#' @param mina Minimum value of the GPD shape parameter to try (a<=0 gives infinite tails to the distribution)
#' @param keepEstimate Flag to indicate is estimation information should be kept
#' @param x vector of samples
#' @param q vector of quantiles
#' @param log Logical if TRUE returns log density [FALSE]
#' @param xlb x axes label for plot [expression(paste(paste("Discharge [",m^3),"/s]"))]
#' @param ylb y axes label for plot ['Profile Likelihood']
#' @param ttl Title of the plot [""]
#'
#' @return In the case of ksgpd a list object containing
#'         brk, Pbrk - the selecte break value and the distribution function value at that point
#'         gpd - a list with the shape (a) and scale (b) of the GPD tail of the distribution
#'         ks - a list with the centers (cntr) and bandwidth (bw) of the kernal smooth approximation
#'         estimate - a list of tried break points (brks) and the corresponding log likelihoos (llkl)
#' 
#' @name ksgpd
NULL

#' @rdname ksgpd
ksgpd <- function(x,numBreak = 1000,skipBreak=10,mina=1e-10,keepEstimate){

    ## sort sample for easier looping
    x <- sort(x[is.finite(x)],decreasing=TRUE)

    ## compute cross-validation densities for all points
    bw <- bw.nrd0(x) # this is a standard deviation
    fcv <- function(ii,x,bw){ mean( dnorm(x[ii],x[-ii],bw) ) }
    cvDens <- sapply(1:length(x),fcv,x=x,bw=bw)
    cvDens[cvDens<1e-100] <- 1e-100 # since in some cases might return -Inf when logged

    
    ## define potential breaks
    brks <- unique(x)
    nb <- min(length(brks),numBreak)
    if(nb < 10){ stop('Fewer then 10 breaks') }
    brks <- brks[skipBreak:nb]

    ## initialise loglikelihood record for the breaks
    nb <- length(brks) 
    bllkl <- rep(NA,nb)

    best <- -Inf
    out <- list()
    for(ii in 1:nb){
        c <- brks[ii] ## break point
        H <- mean( pnorm(c,x,bw) ) ## distribution function at c
        h <- mean( dnorm(c,x,bw) ) ## density function at c

        b <- (1-H)/h ## scale parameter of gpd to get matching density

        ## index of points in kernal smoothed bulk
        idx <- x<=c
        
        ## points covered by gpd shifted so start at 0
        z <- (x[!idx] - c)

        maxa <-  b/(x[1]-c) ## maximum value of shape parameter that x is within the support of the gpd
        ## cehck there is a valid range for a
        if(max<mina){
            warning(paste("Cannot try break point at",c,"since no valid shape paremter for gpd - check mina?"))
            next
        }

        ## optimise the likelihood of z to get a
        fl <- function(a,z,b){ sum(dgpd(z,a,b,log=TRUE)) }
        opt <- optimise(fl,c(mina,maxa),z=z,b=b,maximum = TRUE)
        ## on output opt$max is the value of a, opt$obj is the log likelihood of the gpd part

        ## compute the total log-likelihood
        llkl[ii] <- ( log(1-H) + opt$obj ) + sum( log(cvDens[idx]) )

        if(llkl[ii] > best){
            out <- list(brk=c,Pbrk=H,gpd=list(a=opt$max,b=b),ks=list(cntr=x,bw=bw))
            best <- llkl[ii]
        }
    }

    ## finish off making output
    if(keepEstimate){
        out$estimate = list(brks=brks,llkl=llkl)
    }
    return(out)
}

#' @rdname ksgpd
plot_ksgpd <- function(dst,
                       xlb = expression(paste(paste("Discharge [",m^3),"/s]")),
                       ylb = "Profile Likelihood",
                       ttl = ""){
    if(is.null(dst$estimate)){
        stop("No estimation data saved - cannot plot")
    }

    ii <- which(dst$estimate$brks==dst$brk)
    
    plot(dst$estimate$brks,exp(dst$estimate$brks - dst$estimate$brks[ii]),
             xlab=xlb,ylab=ylb,main=ttl,
             pch=20,ylim=c(0,1.1))
    abline(v=dst$brk,lty=2,col='red',lwd=2)
}
    
#' @rdname ksgpd
#' @export
dksgpd <- function(x,dst,log=FALSE){

    ## initialise with the kernalsmooth densisty estimate
    fk <- function(z,m,s){mean( dnorm(z,m,s) )}
    out <- sapply(x,fk,m=dst$ks$cntr,dst$ks$bw)
    if(log){ out <- log(out) }
    
    ## adapt those above break point
    idx <- x > dst$brk
    out[idx] <- dgpd(x[idx]-dst$brk,dst$gpd$a,dst$gpd$b,log)
    if(log){
        out[idx] <- out[idx] + log(1-dst$Pbrk)
    }else{
        out[idx] <- out[idx] * (1-dst$Pbrk)
    }
    return(out)
}

#' @rdname ksgpd
#' @export
pksgpd <- function(x,dst){

    ## initialise with kernal smooth estimate
    fk <- function(z,m,s){mean( pnorm(z,m,s) )}
    out <- sapply(x,fk,m=dst$ks$cntr,s=dst$ks$bw)

    ## apply alternative for points in gpd tail
    idx <- x > dst$brck
    out[idx] <- dst$Pbrk + (1-dst$Pbrk)*pgpd(x[idx]-dst$brk,dst$gpd$a,dst$gpd$b)
    
    return(out)
}

#' @rdname ksgpd
#' @export
qksgpd <- function(q,dst){

    ## initialise output
    out <- q; out[] <- NA

    ## case in tail - evaluate expicitly
    idx <- q >= dst$Pbrk
    out[idx] <- qgpd( (q[idx]-dst$Pbrk)/(1-dst$Pbrk), dst$gpd$a,dst$gpd$b )

    ## case when in bulk of distribution
    xx <- sort(unique(dst$ks$cntr)) ## a sample of x
    pxx <- pksgpd(xx,dst) ## and corresponding percetile - only valid for pxx < dst$Pbrk

    ## if pxx[1] if less then lowest value of q we need to extend the range downwards
    minq <- min(q[q>0])
    if(pxx[1]>minq){
        lw <- xx[1] - dst$ks$bw
        plw <- pksgpd(lw,dst)
        cnt <- 1
        while(plw > minq & cnt <= 100){
            lw <- lw - dst$ks$bw
            plw <- pksgpd(lw,dst)
        }
        if(cnt>100){
            warning(paste("Unable to evaluate for q<",plw))
        }
        xx <- c(lw,xx)
        pxx <- c(plw,pxx)
    }

    ## loop to evaluate
    for(ii in which(!idx)){
        ## handle the case q=0
        if(q[ii]==0){ out[ii] <- -Inf; next }
        ## handle points out of range
        if(q[ii] < pxx[1]){ next }
                
        ## find number of xx below value
        jj <- sum(pxx<=q)
        ## then bracketed by a set of xx
        rg <- c(xx[jj],xx[jj+1])

        ## optimise
        fopt <- function(z,v,d){ pksgpd(z,d) - v }
        opt <- uniroot(fopt,rg,v=q[ii],d=dst)
        out[ii] <- opt$root
    }

    return(out)
}

