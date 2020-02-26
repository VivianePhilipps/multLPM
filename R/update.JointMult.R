update.JointMult <- function(object,...)
{
    dots <- list(...)

    K <- length(object$nef)
    
    cl <- sapply(dots,class)
    im <- which(cl=="multlcmm")
        
    if(length(im)!=K) stop("m1, m2, etc?")

    V <- matrix(0,length(object$b),length(object$b))
    V[upper.tri(V,diag=TRUE)] <- object$v

    posfix <- which(object$fix==1)
    if(length(posfix))
    {
        Vtot <- matrix(0,length(object$bopt),length(object$bopt))
        Vtot[setdiff(1:length(object$bopt),posfix),setdiff(1:length(object$bopt),posfix)] <- V
    }
    else
    {
        Vtot <- V
    }

    res <- vector("list",length=K)
    

    if(length(im))
    {
        avt <- 0
        for(k in 1:K)
        {
            mm <- dots[[im[k]]]
            
            vk <- Vtot[avt+1:length(mm$best),avt+1:length(mm$best)]

            if(length(mm$call))
            {
                z <- mm$call
                z$maxiter <- 0
                z$verbose <- FALSE
                z$B <- object$bopt[avt+1:length(mm$best)]

                mm <- eval(z)
                
                mm$V <- vk[upper.tri(vk,diag=TRUE)]
                mm$conv <- object$istop
            }
            else
            {
                mm$best <- object$bopt[avt+1:length(mm$best)]
                
                mm$V <- vk[upper.tri(vk,diag=TRUE)]
                
                mm$conv <- object$istop
                mm$gconv <- NA
                mm$loglik <- NA
                mm$niter <- NA
                mm$pred <- NA
                mm$pprob <- NA
                mm$predRE <- NA
                mm$predRE_Y <- NA
                mm$cholesky <- NA
                mm$AIC <- NA
                mm$BIC <- NA
            }

            avt <- avt+length(mm$best)

            res[[k]] <- mm
        }
        
    }

    return(res)
}
