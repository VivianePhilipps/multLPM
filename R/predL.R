
predL <- function(m,newdata,plot=TRUE,add=FALSE,mfrow=NULL,confint=FALSE,...)
{
    dots <- list(...)

    K <- length(m$nef)

    cl <- sapply(dots,class)
    im <- which(cl=="multlcmm")

    dots.plot <- dots[setdiff(1:length(dots),im)]

    if(length(im)!=K) stop("m1, m2, etc?")

    if(is.null(mfrow)) mfrow <- c(K,1)
    if(isTRUE(plot) & !isTRUE(add)) par(mfrow=mfrow)        
                                        #if(isTRUE(plot) & !isTRUE(add)) layout(matrix(1:K,K,1))

    eltk <- function(v,k)
    {
        if(length(v)>=k) return(v[k])
        else return(v[1])
    }
    
    res <- vector("list",K)

                                        #xlim <- matrix(c(65,97,65,97,65,97),2,3)
                                        #ylim <- matrix(c(-3,8,-1,10,-1,8),2,3)
                                        #xlim <- matrix(c(0,1.2,0,1.2,0,1.2),2,3)
                                        #ylim <- matrix(c(-3,6,-1,10,-1,3),2,3)

    if(length(dots$ylim))
    {
        ylim <- matrix(dots$ylim,2,K)
    }
    else
    {
                                        #ylim <- rbind(rep(-10,K),rep(10,K))
    }
    if(length(dots$xlim))
    {
        xlim <- matrix(dots$xlim,2,K)
    }
    else
    {
                                        #xlim <- rbind(rep(min(newdata[,1]),K),rep(max(newdata[,1]),K))
    }
    
    
    usr <- t(cbind(xlim[1,]-(xlim[2,]-xlim[1,])*4/100,
                   xlim[2,]+(xlim[2,]-xlim[1,])*4/100,
                   ylim[1,]-(ylim[2,]-ylim[1,])*4/100,
                   ylim[2,]+(ylim[2,]-ylim[1,])*4/100))

    avt <- 0
    if(length(im))
    {
        univ <- do.call("update",args=c(list(object=m),dots[im]))
        
        for(k in 1:length(im))
        {
            mm <- univ[[k]] #dots[[im[k]]]
            #mm$best <- m$bopt[avt+1:length(mm$best)]

            dotsk <- lapply(dots.plot,eltk,k=k)

            pr <- predictL(x=mm,newdata=newdata,confint=confint)
            if(isTRUE(plot))
            {
                if(isTRUE(add))
                {
                    jk <- (k-1) %% par("mfg")[3] +1
                    ik <- (k-1) %/% par("mfg")[3]  +1
                    par(mfg=c(ik,jk),usr=usr[,k])
                                        #par(mfg=c(k,1),usr=usr[,k])
                    do.call("plot",c(list(x=pr,add=TRUE,ylim=ylim[,k]),dotsk))
                }
                else
                {
                    do.call("plot",c(list(x=pr,ylim=ylim[,k],xlim=xlim[,k]),dotsk))
                                        #print(par("usr"))
                }
            }
            avt <- avt+length(mm$best)
            res[[k]] <- pr
        }
    }

    return(invisible(res))
}



