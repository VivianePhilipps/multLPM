summary.JointMult <- function(object,latex=FALSE,...)
{
    x <- object
    
    if(isTRUE(latex))
        {
            try(library(xtable),silent=TRUE)
            if("package:xtable" %in% search())
                {
                    latex <- TRUE
                }
            else
                {
                    stop("package xtable must be installed to use latex=TRUE")
                }
        }

    id <- 1:length(x$b)
    indice <- rep(id*(id+1)/2)
    se <- sqrt(x$v[indice])
    wald <- (x$b/se)**2
    z <- abs(qnorm((1 + .95)/2))
    binf <- x$b-1.96*se
    bsup <- x$b+1.96*se
    if(x$istop!=1)
        {
            se <- rep(NA,length(se))
            wald <- rep(NA,length(wald))
            binf <- rep(NA,length(x$b))
            bsup <- rep(NA,length(x$b))
        }

    tmp <- matrix(NA,length(x$bopt),5)
    tmp[which(x$fix==0),] <- cbind(round(x$b,3),round(se,3),round(1 - pchisq(wald, 1),5),round(binf,3),round(bsup,3))
    if(sum(x$fix))
        {
            tmp[which(x$fix==1),1] <- x$bopt[which(x$fix==1)]
        }
    colnames(tmp) <- c("coef","se","pval","binf","bsup")
    rownames(tmp) <- names(x$bopt)
    
    cat(" \n")
    cat("Number of iterations: ", x$ni, "\n")
    cat(" \n")
    cat("Convergence criteria: parameters stability=", round(x$convcrit[1],6), "\n")
    cat("                    : likelihood stability=", round(x$convcrit[2],6), "\n") 
    cat("                    : relative distance to maximum=", round(x$convcrit[3],6), "\n")
    cat(" \n")
    cat("Number of parameters: ", length(x$bopt)-sum(x$fix), "\n")
    cat(" \n")
    cat("Maximum log-likelihood:", round(x$loglik,2)," \n")
    cat(" \n")
    cat("AIC:",round(-2*x$loglik+2*(length(x$bopt)-sum(x$fix)),2),"\n") 
    cat(" \n Convergence criteria:",round(x$convcrit,6))
    cat(" \n")
    if(isTRUE(latex)) cat("\n \\bigskip \n ")    
    cat(" \n")
    
    K <- length(x$nef)
    M <- sum(x$ny)
    
    ifix <- which(x$fix==1)

    err <- matrix(NA,sum(x$ny),1+ifelse(any(x$nalea>0),1,0))
 
    cat("Fixed effects : \n")
    if(isTRUE(latex)) cat("\n")
    n <- 0
    nerr <- 0
    for(k in 1:K)
        {
            cat(" Submodel ",k,":\n")

            if(isTRUE(latex))
                {
                    cat("\n")
                    print(xtable(tmp[n+1:x$nef[k],],digits=c(2,3,3,5,3,3),
                                 align=c("l",rep("r",5))),
                          only.contents=FALSE,include.colnames=TRUE,
                          comment=FALSE,floating=FALSE,
                          hline.after=NULL)
                    cat("\n")
                    if(x$ncontr[k]>0)
                        {
                            cat("Contrasts : \n")
                            cat("\n")
                            print(xtable(tmp[n+x$nef[k]+1:x$ncontr[k],],
                                         digits=c(2,3,3,5,3,3),
                                         align=c("l",rep("r",5))),
                                  only.contents=FALSE,include.colnames=TRUE,
                                  comment=FALSE,floating=FALSE,
                                  hline.after=NULL)
                            cat("\n")
                        }                            
                }
            else
                {
                    print(tmp[n+1:x$nef[k],])
                    if(x$ncontr[k]>0)
                        {
                            cat("Contrasts : \n")
                            print(tmp[n+x$nef[k]+1:x$ncontr[k],])
                            cat("\n")
                        }                  
                }

            err[nerr+1:x$ny[k],1] <- tmp[n+x$nef[k]+x$ncontr[k]+x$nvc[k]+x$ncor[k]+1:x$ny[k]]
            if(x$nalea[k]>0) err[nerr+1:x$ny[k],2] <- tmp[n+x$nef[k]+x$ncontr[k]+x$nvc[k]+x$ncor[k]+x$ny[k]+1:x$nalea[k]]
            

            n <- n+x$nef[k]+x$ncontr[k]+x$nvc[k]+x$ncor[k]+x$ny[k]+x$nalea[k]+x$ntrtot[k]
            nerr <- nerr+x$ny[k]

        }

    npmMM <- sum(c(x$nef+x$ncontr+x$nvc+x$ncor+x$ny+x$nalea+x$ntrtot))
    
    cat(" \n")
    if(any(x$nea>0))
        {
            #colnames(x$VRE) <- c("1V","tV","qV","1M","tM","qM","1F","tF","qF")[1:sum(x$nea)]
            #rownames(x$VRE) <- colnames(x$VRE)
            if(isTRUE(latex)) cat("\n \\bigskip \n \n")
            cat("Matrix of variance of the random effects : \n")
            if(isTRUE(latex))
                {
                    if(isTRUE(latex)) cat("\n \\bigskip \n \n")
                    if(ncol(x$VRE)>9) cat("\\small \n")
                    print(xtable(x$VRE,digits=3,
                                 align=c("l",rep("r",ncol(x$VRE)))),
                          include.rownames=TRUE,include.colnames=TRUE,
                          comment=FALSE,floating=FALSE,hline.after=NULL)
                }
            else
                {
                    print(round(x$VRE,3))
                }
            cat(" \n")
            if(isTRUE(latex)) cat("\n \\bigskip \n \n")
            cat("Correlation and p-values of the random effects : \n")

            ptot <- rep(NA,length(x$bopt))
            ptot[which(x$fix==0)] <- as.vector(tmp[which(x$fix==0),3])

            pcor <- matrix(NA,nrow(x$corRE),ncol(x$corRE))
            imatB <- 0
            iB <- 0
            iRE <- 0
            for(k in 1:K)
                {
                    if(x$nea[k]==0 | x$idiag[k]==1) next


                    if(x$nvc[k]>0)
                        {
                            ii <- upper.tri(x$corRE[imatB+1:x$nea[k],imatB+1:x$nea[k]])
                            jj <- setdiff(1:x$nvc[k],((1:x$nea[k])*(1:x$nea[k]+1))/2-1)
                            pcor[imatB+1:x$nea[k],imatB+1:x$nea[k]][ii] <- ptot[iB+x$nef[k]+jj]
                        }
                            
                    if((imatB>0) & (x$nRE>0))
                        {
                            qavt <- sum(x$nea[1:(k-1)])
                            pcor[1:qavt,imatB+1:x$nea[k]] <- ptot[npmMM+iRE+1:(qavt*x$nea[k])]
                            iRE <- iRE+qavt*x$nea[k]
                        }
                    
                    imatB <- imatB+x$nea[k]
                    iB <- iB+x$nef[k]+x$ncontr[k]+x$nvc[k]+x$ncor[k]+x$ny[k]+x$nalea[k]+x$ntrtot[k]
                }

            corRE <- round(x$corRE,3)
            corRE <- format(corRE,digits=3)
            pcor <- format(pcor,digits=5)

            corRE[upper.tri(corRE)] <- pcor[upper.tri(pcor)]
            corRE <- t(corRE)
            # corRE <- qsub("NA","",corRE)

            colnames(corRE) <- colnames(x$VRE)
            rownames(corRE) <- colnames(x$VRE)

            if(isTRUE(latex))
                {
                    if(isTRUE(latex)) cat("\n \\bigskip \n \n")
                    if(ncol(corRE)>9)
                        {
                            cat("\\small \n")
                            cat("\\hspace*{-1.5cm} \n")
                        }
                    corRE <- gsub("NA","",corRE)
                    print(xtable(corRE,digits=3,
                                 align=c("l",rep("r",ncol(x$VRE)))),
                          include.rownames=TRUE,include.colnames=TRUE,
                          comment=FALSE,floating=FALSE,hline.after=NULL,
                          NA.string="")
                }
            else
                {
                    print(corRE,quote=FALSE,right=TRUE)
                }
        }

    cat(" \n")
    if(sum(x$ncor)>0)
        {# AR pas traite
            if(x$nBM==0)
                {
                    npm <- c(0,cumsum(x$nef+x$ncontr+x$nvc+x$ncor+x$ny+x$nalea+x$ntrtot))
                    
                    bm <- abs(x$bopt[npm[1:K]+x$nef[1:K]+x$ncontr[1:K]+x$nvc[1:K]+x$ncor[1:K]])
                    names(bm) <- paste("BM",which(x$ncor==1),sep="")
                    
                    if(isTRUE(latex))
                        {
                            cat("\n \\bigskip \n \n")
                            cat("Standard errors of the brownian motions: \n")
                            cat("\n \\bigskip \n \n")
                            print(xtable(t(bm),digits=c(2,rep(3,length(bm))),
                                         align=c("l",rep("r",length(bm)))),
                                  comment=FALSE,floating=FALSE,
                                  hline.after=NULL,include.rownames=FALSE)
                        }
                    else
                        {
                            cat("Standard errors of the brownian motions: \n")
                            print(round(bm,3))
                            cat(" \n")  
                        }
                }
            else
                {
                    npm <- c(0,cumsum(x$nef+x$ncontr+x$nvc+x$ncor+x$ny+x$nalea+x$ntrtot))
                    
                    bm <- x$bopt[npm[1:K]+x$nef[1:K]+x$ncontr[1:K]+x$nvc[1:K]+x$ncor[1:K]]
                    vbm <- matrix(NA,length(bm),length(bm))
                    diag(vbm) <- bm^2
                    vbm[upper.tri(vbm)] <- x$bopt[npmMM+x$nRE+1:x$nBM]

          
                    if(isTRUE(latex))
                        {
                            cat("\n \\bigskip \n \n")
                            cat("Variance matrix of the brownian motions: \n")
                            cat("\n \\bigskip \n \n")
                            print(xtable(vbm,digits=c(2,rep(3,ncol(vbm))),
                                         align=c("l",rep("r",ncol(vbm)))),
                                  comment=FALSE,floating=FALSE,
                                  hline.after=NULL)
                        }
                    else
                        {
                            cat("Variance matrix of the brownian motions: \n")
                            print(round(vbm,3),na.print="")
                            cat(" \n")  
                        }
                }
        }
    
    if(isTRUE(latex)) cat("\n \\bigskip \n \n")
    cat(" \n")
    cat("Outcome specific parameters:")
    cat(" \n")
    if(ncol(err)==1) colnames(err) <- c("Residual se")
    if(ncol(err)==2) colnames(err) <- c("Residual se", "Random intercept se")
    rownames(err) <- unlist(x$Ynames)
    ntrmax <- max(x$ntr)
    prmlink <- matrix(NA,M,ntrmax)
    j <- 0
    iB <- 0
    for(k in 1:K)
        {
            iB <- iB+x$nef[k]+x$ncontr[k]+x$nvc[k]+x$ncor[k]+x$ny[k]+x$nalea[k]
            for(m in 1:x$ny[k])
                {
                    j <- j+1
                    prmlink[j,1:x$ntr[j]] <- x$bopt[iB+1:x$ntr[j]]
                    iB <- iB+x$ntr[j]
                }
        }
    err2 <- cbind(err,prmlink)
    rownames(err2) <- rownames(err)
    colnames(err2) <- c(colnames(err),"Link",rep("",ncol(err2)-ncol(err)-1))
    if(isTRUE(latex))
        {
            cat("\n \\bigskip \n \n")
            print(xtable(err2,digits=c(2,rep(3,ncol(err2))),
                         align=c("l",rep("r",ncol(err2)))),
                  comment=FALSE,floating=FALSE,
                  hline.after=NULL)
        }
    else
        {
            print(round(err2,3),na.print="")
            cat(" \n")  
        }

    

    cat(" \n")
    if(isTRUE(latex)) cat("\n \\bigskip \n \n")
    cat("Submodel D: \n")

    avtD <- npmMM+x$nRE+x$nBM   
    tmp2 <- tmp[(avtD+1):nrow(tmp),]
    
    if(isTRUE(latex))
        {
            if(isTRUE(latex)) cat("\n \\bigskip \n \n")
            print(xtable(tmp2,digits=c(2,3,3,5,3,3),align=c("l",rep("r",5))),
                  comment=FALSE,floating=FALSE,
                  hline.after=NULL)
        }
    else
        {                 
            print(tmp2)
            cat(" \n")  
        }

    return(invisible(tmp))
}


