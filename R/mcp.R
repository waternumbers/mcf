#' Gaussian model conditional processor from a sample
#'
#' Function produces a Model Conditional Processor, based on relating the various timeseries using a Gaussian Copula. The input time series are expected to be quatiles of the marginal distributions 
#' 
#' @param q An xts object of time series of quantiles. There must be a constant timestep
#' @param series names of series to be used
#' @param maxLag largest lag to apply in the mcp
#' @param method Estimation method to be used, either acf to use auto/cross covariance values or empirical to use empirical estimate derived from arranging the sample as a matrix
#' 
#' @return A list of
#'           mn mean vector
#'           Vr covariance matrix
#'           ts time step in seconds
#' 
#' The mean and variance are labelled by the variable and lag e.g. obs_0, sim_10
#'
#' @export
gmcp <- function(q,series=c("obs","sim"),maxLag,method=c("empirical","acf")){
    ## initialise output
    out <- list()
    
    ## check method
    method = match.arg(method)

    ## check x timestep
    if(!is.xts(q)){ stop("x should be an xts object") }
    dt <- diff(as.numeric(as.POSIXct(index(q))))
    if(!all(dt==dt[1])){ stop("x should have a constant time step") }
    out$ts <- dt[1]
    ## check series are present
    nms <- colnames(x)
    if(!all(series%in%names(x))){ stop('No all variables are present') }
    ## check values lie in 0,1 and convert to gaussian
    x <- q[,series]
    if( any( x<0 | x>1 ) ){ stop("Values of time series so not lie qithin 0,1") }
    x <- qnorm(x)


    ## labels to use for column and row names
    nlbl <- rep(series,each=maxLag+1)
    lbls <- paste(nlbl,rep(0:maxLag,length(nms)),sep="_")
    
    ## empirical method
    if(method=="empirical"){
        Z <- matrix(as.numeric(NA),nrow(x),length(lbl))
        cl <- 1
        for(nm in series){
            for(ii in 0:maxLag){
                Z[,cl] <- filter(x[,nm],c(rep(0,ii),1),method = "convolution",sides=1)
                cl <- cl+1
            }
        }
        
        colnames(Z) <- lbls
        out$mn <- setNames(rep(0,length(lbls)),lbls)
        idx <- which(rowSums(is.finite(Z))==ncol(Z))
        out$vr <- (t(Z[idx,]) %*% Z[idx,])/length(idx)
        out$vr <- cov2cor(out$vr)
    }

    ## acf method
    if(method=="acf"){
        ## estimate mean
        out$mn <- setNames(rep(0,length(lbls)),lbls)

        ## make diagonal elements of Vr
        out$vr <- matrix(NA,length(lbls),length(lbl))
        colnames(vrEst) <- lbls
        rownames(vrEst) <- lbls
        for(ii in series){
            for(jj in series){
                if(ii==jj){
                    tmp <- acf(x[,ii],lag.max=maxLag,
                               type = "covariance",na.action = na.pass,
                               plot=FALSE,demean=FALSE)
                    tmp <- toeplitz(drop(tmp$acf))
                    out$vr[nlbl==ii,nlbl==jj] <- tmp
                }else{
                    tmp <- ccf(x[,ii],x[,jj],lag.max=maxLag,
                               type="covariance",na.action = na.pass,
                               plot=FALSE,demean=FALSE)
                    lwr <- toeplitz(tmp$acf[tmp$lag>=0])
                    
                    uppr <- toeplitz(rev(tmp$acf[tmp$lag<=0]))
                    tmp <- lwr
                    idx <- lower.tri(tmp,diag=FALSE)
                    tmp[idx] <- uppr[idx]
                    rdx <- (ii-1)*(maxLag+1) + (1:(maxLag+1))
                    cdx <- (jj-1)*(maxLag+1) + (1:(maxLag+1))
                    vrEst[rdx,cdx] <- tmp
                    vrEst[cdx,rdx] <- t(tmp)
                }
            }
        }
        vrEst <- cov2cor(vrEst)
    }
    
    out <- list(mn = mnEst,
                Vr =  posdefify(vrEst),
                ts = td)

    return(out)
}
