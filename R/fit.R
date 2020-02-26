
fit_ss <- function(m,subm,data,id,y,time,breaks)
    {
        mcall <- match.call()$m
        mcall$maxiter <- 0
        mcall$pred <- TRUE
        mcall$B <- m$bopt
        
        mfit <- eval(mcall)

        K <- length(y) # K designe l'outcome ici (pas la sphere!)

        res1 <- vector("list",K)
        res2 <- vector("list",K)
#browser()
        jpred <- which(names(mfit)=="predcondY")

        for(k in 1:K)
            {
                fitk <- mfit[[jpred]][which(mfit[[jpred]][,2]==k),]
                colnames(fitk) <- c("dimension","outcome",id,paste("H",y[k],sep=""),"pred",y[k])

                obsk <- na.omit(data[,c(id,y[k],time[k])])

                nk <- nrow(fitk)
                if(nrow(obsk)!=nk) stop(paste("nobs != npred pour y",k,sep=""))
                if(sum(fitk[,id]==obsk[,id])!=nk) stop(paste("pas le meme ordre entre pred et obs pour y",k,sep=""))
                if(sum(fitk[,y[k]]==obsk[,y[k]])!=nk) stop(paste("pas les memes obs pour y",k,sep=""))

                predk <- cbind(fitk,obsk[,3])
                colnames(predk) <- c(colnames(fitk),time[k])

                predk <- cbind(predk,itps=cut(predk[,time[k]],breaks=breaks,include.lowest=TRUE))

                res1[[k]] <- predk

                moy_predk <- tapply(predk[,"pred"],predk[,"itps"],mean)
                moy_obsk <- tapply(predk[,paste("H",y[k],sep="")],predk[,"itps"],mean)
                sd_obsk <- tapply(predk[,paste("H",y[k],sep="")],predk[,"itps"],sd)
                nk <- table(cut(predk[,time[k]],breaks=breaks,include.lowest=TRUE))
                inf_obsk <- moy_obsk -1.96*sd_obsk/sqrt(nk)
                sup_obsk <- moy_obsk +1.96*sd_obsk/sqrt(nk)

                res2[[k]] <- cbind(x=(breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,ypred=moy_predk,yobs=moy_obsk,binf=inf_obsk,bsup=sup_obsk)
                               
            }

        res <- list(res1,res2)
        class(res) <- "fit_ss"
        return(res)
    }

