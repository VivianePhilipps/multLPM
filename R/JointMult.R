RloglikJointMult <- function(b0, bfix0, fix0, Y0, X0, Xd0, idD0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, epsY0, nRE0, nBM0, nbevt0, idcause0, entreRetard0, nbTronc0)
    {
        .Call("loglikJointMult", b0, bfix0, fix0, Y0, X0, Xd0, idD0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, epsY0, nRE0, nBM0, nbevt0, idcause0, entreRetard0, nbTronc0)
    }


predcondY <- function(b0, bfix0, fix0, Y0, X0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, nRE0, nBM0, chol0)
    {
        .Call("predcondY",b0, bfix0, fix0, Y0, X0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, nRE0, nBM0, chol0)
    }


## D : liste de 1 ou 2 formules avec a gauche
##      longDiag(age0,age,dem) si on a deja des mesures repetees dans age et dem
##  ou
##      survBreak(T0,T,DC) si on a 1 mesure par sujet dans T et dc (-> a discretiser)



JointMult <- function(Y,D,data,var.time,RE="block-diag",BM="diag",B,posfix,maxiter=0,nproc=1,verbose=FALSE,file="",pred=FALSE,breaks=NULL,delayed=TRUE,eps=c(0.0001,0.0001,0.001))
    {
        ## verif arg
        if(missing(Y)) stop("Y is missing")
        if(missing(D)) stop("D is missing")
        if(missing(data)) stop("data is missing")
        if(missing(var.time)) stop("var.time is missing")

        if(!is.list(Y)) stop("Y should be a list of multlcmm objects")
        if(!(all(sapply(Y,class)=="multlcmm"))) stop("Y should only contain multlcmm objects")
        if(!is.list(D)) stop("D should be a list of 1 or 2 formula")
        if(length(D)>2) stop("D should be of length 1 or 2")
        if(!(all(sapply(D,class)=="formula")))  stop("D should only contain two-sided formula")
        if(!is.data.frame(data)) stop("data should be a data frame")
        if((!nrow(data)) | (!ncol(data))) stop("data is empty")
        if(!(RE %in% c("block-diag","full"))) stop("RE should be either 'block-diag' or 'full'")
        if(!(BM %in% c("diag","full"))) stop("BM should be either 'diag' or 'full'")
        if(!(maxiter>=0)) stop("maxiter should be non negative")
        if(!(nproc>0)) stop("nproc should be positive")
        if(!is.logical(verbose)) stop("verbose should be logical")
        if(!is.character(file)) stop("file should be a character")
        if(!is.logical(pred)) stop("pred should be logical")
        if(!is.null(breaks)){if(!is.vector(breaks)) stop("breaks should be a vector")}
        if(!all(is.logical(delayed))) stop("delayed should be a logical vector")
        if(length(eps)!=3) stop("eps should be of length 3")
        

        ##### donnees Y #####
        K <- length(Y)
        if(length(var.time)==1) var.time <- rep(var.time,K)
        dataY <- NULL
        ny <- rep(NA,K)
        Ynames <- vector("list",K)
        Xnames <- vector("list",K)
        nomsX <- unique(unlist(sapply(Y,function(x) x$Xnames2)))
        
        for(k in 1:K)
            {
                ### modele k
                if(length(Y[[k]]$call))
                    {
                        z <- Y[[k]]$call
                        z$data <- data
                        z$maxiter <- 0
                        z$B <- Y[[k]]$best
                        z$verbose <- FALSE
                        mod <- eval(z)
                    }
                else
                    {
                        mod <- eval(Y[[k]])
                    }
                assign(paste("mod",k,sep=""),mod)

                subject <- mod$call$subject
                if(k>1){if(subject != colnames(dataY)[1]) stop("Subject variable should be the same for all multlcmm models")}
                
                colx <- c(subject,mod$Xnames2)
                
                ny[k] <- length(mod$Ynames)

                Ynames[[k]] <- mod$Ynames
                Xnames[[k]] <- mod$Xnames2

                ## donnees km
                for(m in 1:ny[k])
                    {
                        ## data frame de l'outcome m
                        colx <- c(subject,nomsX,mod$Ynames[m])
                        if(length(mod$na.action[[m]]))
                            {
                                datam <- data[-mod$na.action[[m]],colx,drop=FALSE]
                            }
                        else
                            {
                                datam <- data[,colx,drop=FALSE]
                            }

                        datam <- datam[order(datam[,1],datam[,var.time[k]]),,drop=FALSE]
                        old <- colnames(datam)
                        datam$processK <- k
                        datam$outcomeM <- m
                        #datam <- datam[,c("processK","outcomeM",old)]
                        datam <- datam[,c(old[1],"processK","outcomeM",old[-1])]

                        if((k==1) & (m==1))
                            {
                                dataY <- datam
                                colnames(dataY)[which(colnames(dataY)==var.time[k])] <- "timeT"
                                colnames(dataY)[which(colnames(dataY)==Ynames[[k]][m])] <- "measureY"
                            }
                        else
                            {
                                colnames(datam)[which(colnames(datam)==var.time[k])] <- "timeT"
                                colnames(datam)[which(colnames(datam)==Ynames[[k]][m])] <- "measureY"
                                Xplus <- setdiff(colnames(datam),colnames(dataY))
                                if(length(Xplus))
                                    {
                                        for(l in 1:length(Xplus))
                                            {
                                                old <- colnames(dataY)
                                                dataY <- cbind(dataY,NA)
                                                colnames(dataY) <- c(old,Xplus[l])
                                            }
                                    }

                                Xmqt <- setdiff(colnames(dataY),colnames(datam))
                                if(length(Xmqt))
                                    {
                                        for(l in 1:length(Xmqt))
                                            {
                                                old <- colnames(datam)
                                                datam <- cbind(datam,NA)
                                                colnames(datam) <- c(old,Xmqt[l])
                                            }
                                        datam <- datam[,colnames(dataY),drop=FALSE]
                                    }

                                dataY <- rbind(dataY,datam)
                            }
                    }
                
            }

        
        #### donnees D ####
        nbevt <- length(D)

        ## premier evt
        dem <- NULL
        D1 <- D[[1]]
        tD1 <- terms(D1)        
        formD1 <- gsub("[[:space:]]","",D1[3])
        
        if(!any(c(grep("longDiag",D1),grep("survBreak",D1)))) stop("Left side of formula in D should be either 'longDiag()' or 'survBreak()'")
         
        gauche1 <- all.vars(tD1[[2]]) #t0, t et d
       
        if(length(grep("longDiag",D1[2])))
            {
                discretise <- 0
                
                if(length(gauche1)!=3) stop("Please specified entry time, time and indicator in that order in longDiag")
            
                entry <- gauche1[1]
                var.time <- c(var.time,gauche1[2])
                dem <- gauche1[3]
            }
        else
            {
                if(length(grep("survBreak",D1[2])))
                    {
                        discretise <- 1
                        if(length(gauche1)!=3) stop("Please specified entry time, time and indicator in that order in survBreak")
                        
                        entry <- gauche1[1]
                        var.time <- c(var.time,gauche1[2])
                        dem <- gauche1[3]
                    }
            }
               
        varexplD1 <- all.vars(formula(paste("~",D1[3])))

        
        ## deuxieme evt
        dc <- NULL
        D2 <- NULL
        gauche2 <- NULL
        varexplD2 <- NULL
        if(nbevt==2)
            {
                D2 <- D[[2]]
                tD2 <- terms(D2)
                formD2 <- gsub("[[:space:]]","",D2[3])

                if(!any(c(grep("longDiag",D2),grep("survBreak",D2)))) stop("Left side of formula in D should be either 'longDiag()' or 'survBreak()'")
        
                gauche2 <- all.vars(tD2[[2]]) #t0, t et d
                
                if(length(grep("longDiag",D2[2])))
                    {
                        discretise <- c(discretise,0)

                        if(length(gauche2)!=3) stop("Please specified entry time, time and indicator in that order in longDiag")

                        entry <- c(entry,gauche2[1])
                        var.time <- c(var.time,gauche2[2])
                        dc <- gauche2[3]
                    }
                else
                    {
                        if(length(grep("survBreak",D2[2])))
                            {
                                discretise <- c(discretise,1)

                                if(length(gauche2)!=3) stop("Please specified entry time, time and indicator in that order in survBreak")
                                
                                entry <- c(entry,gauche2[1])
                                var.time <- c(var.time,gauche2[2])
                                dc <- gauche2[3]
                            }
                    }
                
                varexplD2 <- all.vars(formula(paste("~",D2[3])))

                if(all(discretise==1)) stop("At most 1 event time can be discretized")
                if((discretise[1]==1) & (discretise[2]==0)) stop("exchange events!") #!**
            }

        varexplD <- unique(c(varexplD1,varexplD2))
        nvarexplD <- length(varexplD)
        delayed <- rep(delayed,length.out=nbevt)
        
        
        
        dataD1 <- NULL
        dataD2 <- NULL
        ni01 <- NULL
        ni02 <- NULL
        
        ## donnees dem
        if(!is.null(dem))
            {
                if(discretise[1]==0)
                    {
                        dataD1 <- data[,c(subject,entry[1],var.time[K+1],dem,setdiff(unique(c(varexplD1,nomsX)),c(var.time,entry[1])))]
                        dataD1 <- na.omit(dataD1)
                        dataD1 <- dataD1[order(dataD1[,1],dataD1[,3]),]

                        ## ajouter temps T0 en premiere ligne de chaque sujet
                        if(delayed[1]==TRUE)
                            {
                                ntmp <- table(dataD1[,1])
                                dataD1entry <- dataD1[c(1,cumsum(ntmp[-length(ntmp)])+1),]
                                dataD1entry[,3] <- dataD1entry[,2]
                                dataD1entry$measure <- 0
                                dataD1entry[,dem] <- 0
                                
                                dataD1measure <- dataD1
                                dataD1measure$measure <- 1

                                dataD1bis <- merge(dataD1entry,dataD1measure,all=TRUE,sort=FALSE)
                                jinf <- which(!is.finite(dataD1bis[,var.time[K+1]]))
                                if(length(jinf))
                                    {
                                        dataD1bis <- dataD1bis[-jinf,]
                                    }
                                dataD1 <- dataD1bis[order(dataD1bis[,subject],dataD1bis[,"measure"],dataD1bis[,var.time[K+1]]),]
                            }
                        
                        ## garder jusqu au premier 1
                        nd <- as.vector(table(dataD1[,1]))
                        nok <- tapply(dataD1[,dem],dataD1[,1],function(x){return(1:ifelse(any(x==1),which(x==1)[1],length(x)))})
                        garde <- unlist(sapply(1:length(nok),function(k,nd,nok) sum(nd[1:k])-nd[k]+nok[[k]],nd=nd,nok=nok))
                        
                        dataD1 <- dataD1[garde,,drop=FALSE]
                        
                        ## temps max pr chaque sujet
                        max1 <- tapply(dataD1[,3],dataD1[,1],max)
                        t1max <- data.frame(id=names(max1),t1max=as.vector(max1),stringsAsFactors=FALSE)
                        maxd1 <- as.numeric(tapply(dataD1[,4],dataD1[,1],max))
                        t1max <- data.frame(t1max,maxd1)
                        
                        colnames(t1max) <- c(colnames(dataY)[1],"t1max","d1")
                        
                        ## nb de mesures
                        ni01 <- table(dataD1[,1]) # mesure T0 dedans! revoir !**

                        ## entree retardee
                        if(delayed[1]==TRUE)
                            {
                                nbTronc <- matrix(1,length(ni01),ncol=1)

                                dataD1T0 <- unique(dataD1[,c(subject,entry[1])])
                                nbTronc[which(!is.finite(dataD1T0[,2])),1] <- 0
                            }
                        else
                            {
                                nbTronc <- matrix(0,length(ni01),ncol=1)
                            }

                        ## enlever prevalents
                        ## ipreval <- which((ni01==1) & (t1max[,3]==1))
                        ## if(length(ipreval))
                        ##     {
                        ##         #stop("Prevalent subjects for event 1") #!**
                        ##         dataD1 <- dataD1[-which(dataD1[,1] %in% t1max[ipreval,1]),,drop=FALSE]
                        ##         ni01 <- ni01[-ipreval]
                        ##         t1max <- t1max[-ipreval,,drop=FALSE]
                        ##     }

                        ## enlever les Y post-diag
                        ## dataYbis <- merge(dataY,t1max[,1:2],all.x=TRUE)
                        ## dataYbis <- dataYbis[order(dataYbis[,2],dataYbis[,3],dataYbis[,1],dataYbis[,"timeT"]),]
                        ## jpost <- which(dataYbis$timeT>dataYbis$t1max)
                        ## if(length(jpost))
                        ##     {
                        ##         dataYbis <- dataYbis[-jpost,,drop=FALSE]
                        ##     }
                        ## dataY <- dataYbis[,-which(colnames(dataYbis)=="t1max")]
                    }
                else
                    {
                        ## evt absorbant et en temps continu
                        
                        dataD1 <- data[,c(subject,entry[1],var.time[K+1],dem,setdiff(unique(c(varexplD1,nomsX)),c(var.time,entry[1])))]
                        dataD1 <- na.omit(dataD1)
                        dataD1 <- dataD1[order(dataD1[,1],dataD1[,3]),]
                        
                        ## garder la premiere/derniere ligne de chaque sujet de dataD1
                        ntmp <- table(dataD1[,1])
                        #dataD1 <- dataD1[c(1,cumsum(ntmp[-length(ntmp)])+1),,drop=FALSE]# premiere
                        dataD1 <- dataD1[cumsum(ntmp),,drop=FALSE]# derniere
                        
                        ## discretiser
                        breaks <- sort(breaks)
                        tbr <- (breaks[-length(breaks)]+breaks[-1])/2
                        nbr <- length(breaks)

                        t1max <- dataD1[,c(1,3,4,2)]
                        colnames(t1max) <- c(colnames(dataY)[1],"t1max","d1","t1min")

                        t1max <- data.frame(t1max,intervalle=as.numeric(cut(x=t1max$t1max,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))
                        t1max <- data.frame(t1max,intervT0=as.numeric(cut(x=t1max$t1min,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))

                        ## si derniere mesure=0, s arreter a l intervalle  precedent
                        t1max$intervalle[which(t1max$d1==0)] <- t1max$intervalle[which(t1max$d1==0)]-1

                        ## entree retardee
                        nbTronc <- matrix(rep(as.numeric(delayed[1]==TRUE),nrow(t1max)),ncol=1)
                        nbTronc[which(t1max$intervT0==1),1] <- 0

                        ## recreer les donnees dc
                        dataD1bis <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))
                        ni01 <- NULL
                        for(i in 1:nrow(t1max))
                            {
                                if(delayed[1]==TRUE)
                                    {
                                        t1 <- tbr[(t1max[i,6]-1):t1max[i,5]]
                                        t01 <- rep(-Inf,length(t1))
                                        if(t1max$intervT0[i]>1)
                                            {
                                                t01 <- rep(tbr[t1max[i,6]-1],length(t1))
                                            }
                                    }
                                else
                                    {
                                        t1 <- tbr[t1max[i,6]:t1max[i,5]]
                                        t01 <- rep(-Inf,length(t1))
                                    }

                                di1 <- data.frame(rep(t1max[i,1],length(t1)),t01,t1,rep(0,length(t1)))
                                if(t1max[i,3]==1) di1[length(t1),4] <- 1
                                if(ncol(dataD1)>4)
                                    {
                                        ## boucle pour garder les differents types (character, numeric)
                                        for(k in 5:ncol(dataD1))
                                            {
                                                di1 <- data.frame(di1,rep(dataD1[which(dataD1[,1]==t1max[i,1])[1],k],length(t1)))
                                            }
                                    }
                                        
                                colnames(di1) <- colnames(dataD1)
                                dataD1bis <- rbind(dataD1bis,di1)

                                ni01 <- c(ni01,length(which(dataD1bis[,1]==t1max[i,1])))
                            }
                        
                        dataD1 <- dataD1bis
                    }

                D <- t1max[,c(1,3)]
            }
        
        ## donnees dc
        if(!is.null(dc))
            {                
                if(discretise[2]==0)
                    {
                        ## 2e evt de type diagnostic
                        dataD2 <- data[,c(subject,entry[2],var.time[K+2],dc,setdiff(unique(c(varexplD2,nomsX)),c(var.time,entry[2])))]
                        dataD2 <- na.omit(dataD2)
                        dataD2 <- dataD2[order(dataD2[,1],dataD2[,3]),]

                        ## ajouter temps T0 en premiere ligne de chaque sujet
                        if(delayed[2]==TRUE)
                            {
                                ntmp <- table(dataD2[,1])
                                dataD2entry <- dataD2[c(1,cumsum(ntmp[-length(ntmp)])+1),]
                                dataD2entry[,3] <- dataD2entry[,2]
                                dataD2entry$measure <- 0
                                dataD2entry[,dc] <- 0
                                
                                dataD2measure <- dataD2
                                dataD2measure$measure <- 1

                                dataD2bis <- merge(dataD2entry,dataD2measure,all=TRUE,sort=FALSE)
                                
                                jinf <- which(!is.finite(dataD2bis[,var.time[K+1]]))
                                if(length(jinf))
                                    {
                                        dataD2bis <- dataD2bis[-jinf,]
                                    }
                                dataD2 <- dataD2bis[order(dataD2bis[,subject],dataD2bis[,var.time[K+2]]),]
                            }
                        
                        ## garder jusqu au premier 1
                        nd <- as.vector(table(dataD2[,1]))
                        nok <- tapply(dataD2[,dc],dataD2[,1],function(x){return(1:ifelse(any(x==1),which(x==1)[1],length(x)))})
                        garde <- unlist(sapply(1:length(nok),function(k,nd,nok) sum(nd[1:k])-nd[k]+nok[[k]],nd=nd,nok=nok))
                        
                        dataD2 <- dataD2[garde,,drop=FALSE]

                        ## temps max pr chaque sujet
                        max2 <- tapply(dataD2[,3],dataD2[,1],max)
                        t2max <- data.frame(id=names(max2),t2max=as.vector(max2),stringsAsFactors=FALSE)
                        maxd2 <- as.numeric(tapply(dataD2[,4],dataD2[,1],max))
                        t2max <- data.frame(t2max,maxd2)
                        colnames(t2max) <- c(colnames(dataY)[1],"t2max","d2")

                        ## enlever prevalents
                        ## ni02 <- table(dataD2[,1])
                        ## ipreval <- which((ni02==1) & (t2max[,3]==1))
                        ## if(length(ipreval))
                        ##     {
                        ##         #stop("Prevalents subjects for event 2") #!**
                        ##         dataD2 <- dataD2[-which(dataD2[,1] %in% t2max[ipreval,1]),,drop=FALSE]
                        ##         ni02 <- ni02[-ipreval]
                        ##         t2max <- t2max[-ipreval,,drop=FALSE]
                        ##     }
                                     
                        ## ## enlever les Y post-diag           
                        ## dataYbis <- merge(dataY,t2max[,1:2,drop=FALSE],all.x=TRUE)
                        ## dataYbis <- dataYbis[order(dataYbis[,2],dataYbis[,3],dataYbis[,1]),]
                        ## jpost <- which(dataYbis$timeT>dataYbis$t2max)
                        ## if(length(jpost))
                        ##     {
                        ##         dataYbis <- dataYbis[-jpost,,drop=FALSE]
                        ##     }
                        ## dataY <- dataYbis[,-which(colnames(dataYbis)=="t2max")]

                        ## matrice sujet indicateur evt
                        D <- merge(t1max,t2max,all=FALSE,sort=FALSE) #selection ici
                        evt <- (D[,"t1max"]<=D[,"t2max"])*D[,"d1"]+2*(D[,"t1max"]>=D[,"t2max"])*D[,"d2"]
                        evt[which(evt==3)] <- 1 # ? !** si ex aequo, on prend D1
                        Tevt <- apply(D[,c("t1max","t2max")],1,min) # ? a voir si ok  si evt=0!**
                        D <- data.frame(D[,1],evt,Tevt)
                        colnames(D) <- c(colnames(dataY)[1],"evt","Tevt")
                        ##D <- D[order(D[,1]),] #!** remettre si sort=FALSE enleve

                        ## dataD1 et dataD2 sans les temps post-diag
                        dataD1bis <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))
                        dataD2bis <- as.data.frame(matrix(NA,0,ncol(dataD2),dimnames=list(NULL,colnames(dataD2))))
                        ni01 <- NULL
                        ni02 <- NULL
                        for(i in 1:nrow(D))
                        {
                            if(D[i,2]==0)
                            {
                                dataD1bis <- rbind(dataD1bis,dataD1[which(dataD1[,1]==D[i,1]),])
                                dataD2bis <- rbind(dataD2bis,dataD2[which(dataD2[,1]==D[i,1]),])
                            }
                            
                            if(D[i,2]==1)
                            {
                                dataD1bis <- rbind(dataD1bis,dataD1[which(dataD1[,1]==D[i,1]),])
                                
                                di2 <- dataD2[which((dataD2[,1]==D[i,1]) & (dataD2[,3]<=D[i,3])),] # ? a voir si < ou <= !**
                                if(nrow(di2)==0) stop(paste("Subject",D[i,1],"has no observation for event 2 \n"))
                                di2[,4] <- 0
                                dataD2bis <- rbind(dataD2bis,di2)
                            }
                            if(D[i,2]==2)
                            {
                                dataD2bis <- rbind(dataD2bis,dataD2[which(dataD2[,1]==D[i,1]),])
                                
                                di1 <- dataD1[which((dataD1[,1]==D[i,1]) & (dataD1[,3]<=D[i,3])),]
                                if(nrow(di1)==0) stop(paste("Subject",D[i,1],"has no observation for event 1 \n"))
                                di1[,4] <- 0
                                dataD1bis <- rbind(dataD1bis,di1)
                            }
                            
                            ni01 <- c(ni01,length(which(dataD1bis[,1]==D[i,1])))
                            ni02 <- c(ni02,length(which(dataD2bis[,1]==D[i,1])))
                        }
                        
                        dataD1 <- dataD1bis
                        dataD2 <- dataD2bis
                        D <- D[,1:2]                        

                        ## entree retardee
                        nbTronc <- matrix(0,nrow(D),2)
                        if(delayed[1]==TRUE)
                            {
                                nbTronc[,1] <- 1

                                dataD1T0 <- unique(dataD1[,c(subject,entry[1])])
                                nbTronc[which(!is.finite(dataD1T0[,2])),1] <- 0
                            }
                        
                        if(delayed[2]==TRUE)
                            {
                                nbTronc[,2] <- 1

                                dataD2T0 <- unique(dataD2[,c(subject,entry[2])])
                                nbTronc[which(!is.finite(dataD2T0[,2])),2] <- 0
                            }
                    }
                else
                    {
                        ## evt absorbant et en temps continu
                        
                        dataD2 <- data[,c(subject,entry[2],var.time[K+2],dc,setdiff(unique(c(varexplD2,nomsX)),c(var.time,entry[2])))]
                        dataD2 <- na.omit(dataD2)
                        dataD2 <- dataD2[order(dataD2[,1],dataD2[,3]),]
                        
                        ## garder la premiere ligne de chaque sujet de dataD2
                        ntmp <- table(dataD2[,1])
                        #dataD2 <- dataD2[c(1,cumsum(ntmp[-length(ntmp)])+1),,drop=FALSE]#premiere
                        dataD2 <- dataD2[cumsum(ntmp),,drop=FALSE] # derniere
                        
                        ## discretiser
                        breaks <- sort(breaks)
                        tbr <- (breaks[-length(breaks)]+breaks[-1])/2
                        nbr <- length(breaks)

                        t2max <- dataD2[,c(1,3,4,2)]
                        colnames(t2max) <- c(colnames(dataY)[1],"t2max","d2","t2min")

                        t2max <- data.frame(t2max,intervalle=as.numeric(cut(x=t2max$t2max,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))
                        t2max <- data.frame(t2max,intervT0=as.numeric(cut(x=t2max$t2min,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))

  
                        tmax <- merge(t2max,t1max,all=FALSE,sort=FALSE) #selection ici
                        tmax <- data.frame(tmax,intervDem=as.numeric(cut(x=tmax$t1max[which(t1max[,1] %in% t2max[,1])],breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))
                        
                        ## si dem=1, mettre dc=0
                        tmax$d2[which(tmax$d1==1)] <- 0
                        ## si dc=0, s arreter avt l intervalle de dem
                        tmax$intervalle[which(tmax$d2==0)] <- tmax$intervDem[which(tmax$d2==0)]-1
                        
                        ## indicateur evt
                        D <- matrix(tmax[,1],ncol=1)
                        evt <- rep(0,nrow(tmax))
                        evt[which(tmax$d1==1)] <- 1
                        evt[which(tmax$d2==1)] <- 2
                        D <- data.frame(D,evt,stringsAsFactors=FALSE)
                        colnames(D) <- c(colnames(dataY)[1],"evt")
                        D <- D[order(D[,1]),]
                        
                        ## recreer les donnees dc
                        dataD1bis <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))
                        dataD2bis <- as.data.frame(matrix(NA,0,ncol(dataD2),dimnames=list(NULL,colnames(dataD2))))
                        ni01 <- NULL
                        ni02 <- NULL
                        for(i in 1:length(evt))
                            {
                                dataD1bis <- rbind(dataD1bis,dataD1[which(dataD1[,1]==D[i,1]),])
                                if(delayed[nbevt]==TRUE)
                                    {
                                        t2 <- tbr[(tmax[i,6]-1):tmax[i,5]]
                                        t02 <- rep(-Inf,length(t2))
                                        if(tmax$intervT0[i]>1)
                                            {
                                                t02 <- rep(tbr[tmax[i,6]-1],length(t2))
                                            }
                                    }
                                else
                                    {
                                        t2 <- tbr[tmax[i,6]:tmax[i,5]]
                                        t02 <- rep(-Inf,length(t2))
                                    }
                                
                                di2 <- data.frame(rep(D[i,1],length(t2)),t02,t2,rep(0,length(t2)))
                                if(D[i,2]==2) di2[length(t2),4] <- 1
                                if(ncol(dataD2)>4)
                                    {
                                        ## boucle pour garder les differents types (character, numeric)
                                        for(k in 5:ncol(dataD2))
                                            {
                                                di2 <- data.frame(di2,rep(dataD2[which(dataD2[,1]==D[i,1])[1],k],length(t2)))
                                            }
                                        #di2X <- matrix(rep(dataD2[which(dataD2[,1]==D[i,1])[1],4:ncol(dataD2)],each=length(t2)),nrow=length(t2))
                                        #di2 <- data.frame(di2,di2X)
                                    }
                                colnames(di2) <- colnames(dataD2)
                                dataD2bis <- rbind(dataD2bis,di2)
                                
                                ni01 <- c(ni01,length(which(dataD1bis[,1]==D[i,1])))
                                ni02 <- c(ni02,length(which(dataD2bis[,1]==D[i,1])))
                                
                            }
                        
                        dataD1 <- dataD1bis
                        dataD2 <- dataD2bis

                        ## entree retardee
                        nbTronc <- matrix(0,nrow(D),2)
                        if(delayed[1]==TRUE)
                            {
                                nbTronc[,1] <- 1

                                dataD1T0 <- unique(dataD1[,c(subject,entry[1])])
                                nbTronc[which(!is.finite(dataD1T0[,2])),1] <- 0
                            }
                        if(delayed[2]==TRUE)
                            {
                                nbTronc[,2] <- 1
                                nbTronc[which(tmax$intervT0==1),2] <- 0
                            }
                    }
            }
        
       

        ## union des temps de mesure evts
        if(nbevt==1)
            {
                idD <- rep(1,sum(ni01))
                dataD <- dataD1
                dataD$tdemdc <- dataD[,var.time[K+1]]
            }
        if(nbevt==2)
            {
                #colnames(dataD2)[which(colnames(dataD2)==var.time[K+2])] <- colnames(dataD1)[which(colnames(dataD1)==var.time[K+1])]
                dataD1$tdemdc <- dataD1[,var.time[K+1]]
                dataD2$tdemdc <- dataD2[,var.time[K+2]]
                dataD <- merge(dataD1,dataD2,all=TRUE,sort=FALSE)
                #dataD <- dataD[order(dataD[,1],dataD[,var.time[K+1]]),]
                dataD <- dataD[order(dataD[,1],dataD[,"tdemdc"]),]
                idD <- c(as.numeric(!is.na(dataD[,dem])),as.numeric(!is.na(dataD[,dc])))
            }
        

        ## nb de sujets (select au moins 1 mesure dc, non prevalents diag)       
        ni0 <- table(dataD[,1])
        ns <- length(ni0)

        ## enlever les sujets exclus de dataY
        ntmp <- table(dataY[,1])
        ienlev <- setdiff(names(ntmp),names(ni0))
        if(length(ienlev))
            {
                dataY <- dataY[-which(dataY[,1] %in% ienlev),]
            }
        

        ## variables sur D
        nvarD <- 0
        ndept <- 0
        nvarD1 <- 0
        ndept1 <- 0
        nvarD2 <- 0
        ndept2 <- 0
        varD <- NULL
        varD1 <- NULL
        varD2 <- NULL
        vardept <- NULL
        vardept1 <- NULL
        vardept2 <- NULL
        idcause <- matrix(0,nbevt,nvarexplD)
        modmat1 <- NULL
        modmat2 <- NULL
        Xseuil <- NULL

        namesevt1 <- NULL
        if(D1[[3]]!=1)
            {
                form1 <- formula(paste("~",subject,"+",D1[3]))
                modmat1 <- model.matrix(form1,data=dataD[which(!is.na(dataD[,dem])),])[,-1,drop=FALSE]
                vars1 <- labels(terms(form1))
                placea <- NULL
                placeb <- NULL
                for(l in 2:ncol(modmat1))
                    {
                        tmp <- unique(na.omit(modmat1[,c(1,l)]))
                        if(nrow(tmp)==ns) # variable pas dep du tps
                            {
                                varD1 <- c(varD1,vars1[attr(modmat1,"assign")[l]])
                                nvarD1 <- nvarD1+1
                                placea <- c(placea,l-1)
                            }
                        else # variable dep du tps
                            {
                                vardept1 <- c(vardept1,vars1[attr(modmat1,"assign")[l]])
                                ndept1 <- ndept1+1
                                placeb <- c(placeb,l-1)
                            }
                    }
                namesevt1 <- colnames(modmat1)[-1]
                namesevt1 <- c(namesevt1[placea],namesevt1[placeb])
                
                modmat1 <- modmat1[,c(1,placea+1,placeb+1)]
            }

        namesevt2 <- NULL
        if(nbevt==2)
            {
                if(D2[[3]]!=1)
                    {
                        form2 <- formula(paste("~",subject,"+",D2[3]))
                        modmat2 <- model.matrix(form2,data=dataD[which(!is.na(dataD[,dc])),])[,-1,drop=FALSE]
                        vars2 <- labels(terms(form2))
                        placea <- NULL
                        placeb <- NULL
                        for(l in 2:ncol(modmat2))
                            {
                                tmp <- unique(na.omit(modmat2[,c(1,l)]))
                                if(nrow(tmp)==ns) # variable pas dep du tps
                                    {
                                        varD2 <- c(varD2,vars2[attr(modmat2,"assign")[l]])
                                        nvarD2 <- nvarD2+1
                                        placea <- c(placea,l-1)
                                    }
                                else # variable dep du tps
                                    {
                                        vardept2 <- c(vardept2,vars2[attr(modmat2,"assign")[l]])
                                        ndept2 <- ndept2+1
                                        placeb <- c(placeb,l-1)
                                    }
                            }
                        namesevt2 <- colnames(modmat2)[-1]
                        namesevt2 <- c(namesevt2[placea],namesevt2[placeb])

                        modmat2 <- modmat2[,c(1,placea+1,placeb+1)]
                    }
            }

        if(nbevt==1)
            {
                if(nvarD1>0)
                    {
                        D <- cbind(D,unique(modmat1[,c(1:(nvarD1+1))])[,-1])
                    }

                if(ndept1>0)
                    {
                        Xseuil <- modmat1[,1+nvarD1+1:ndept1,drop=FALSE]
                    }

                idcause <- matrix(1,nrow=1,ncol=nvarD1+ndept1)
            }

        if(nbevt==2)
            {
                if((nvarD1+nvarD2)>0)
                    {
                        tmp1 <- unique(modmat1[,1:(nvarD1+1),drop=FALSE])
                        tmp2 <- unique(modmat2[,1:(nvarD2+1),drop=FALSE])
                        modmat <- cbind(tmp1,tmp2[,-1,drop=FALSE])
                        modmat <- modmat[,c(1,1+which(!duplicated(colnames(modmat)[-1])))]

                        D <- cbind(D,modmat[,-1])

                        idcause <- matrix(0,nrow=nbevt,ncol=ncol(modmat)-1)
                        idcause[1,1:nvarD1] <- 1
                        idcause[2,which(colnames(modmat1)[-1] %in% colnames(modmat2)[-1])] <- 1
                        if(ncol(idcause)>nvarD1)
                            {
                                idcause[2,(nvarD1+1):ncol(idcause)] <- 1
                            }
                    }

                if((ndept1+ndept2)>0)
                    {
                        Xseuil1 <- NULL
                        Xseuil2 <- NULL

                        if(ndept1>0)
                            {
                                Xseuil1 <- modmat1[,1+nvarD1+1:ndept1,drop=FALSE]
                            }
                    
                        if(ndept2>0)
                            {
                                Xseuil2 <- modmat2[,1+nvarD2+1:ndept2,drop=FALSE]
                            }

                        Xseuil1 <- cbind(Xseuil1,un=rep(1,sum(ni01)))
                        Xseuil2 <- cbind(Xseuil2,deux=rep(2,sum(ni02)))
                        Xseuil <- merge(Xseuil1,Xseuil2,by.x=c("un"),by.y=c("deux"),all=TRUE,sort=FALSE)
                        Xseuil <- Xseuil[,setdiff(colnames(Xseuil),c("un","deux")),drop=FALSE]
                        if(any(is.na(Xseuil)))
                            {
                                Xseuil <- apply(Xseuil,2,function(x) { x[which(is.na(x))] <- -3000; return(x) })
                            }

                        idcause <- cbind(idcause,
                                         rbind(c(rep(1,ndept1),rep(0,ndept2)),
                                               c(rep(0,ndept1),rep(1,ndept2))))
                    }
            }
        

        
###############################
        ## if(nvarexplD>0)
        ##     {
        ##         for (l in 1:nvarexplD)
        ##             {
        ##                 tmp <- unique(na.omit(dataD[,c(subject,varexplD[l])]))
        ##                 if(nrow(tmp)==ns) # variable pas dep du tps
        ##                     {
        ##                         varD <- c(varD,varexplD[l])
        ##                         nvarD <- nvarD+1
        ##                     }
        ##                 else # variable dep du tps
        ##                     {
        ##                         vardept <- c(vardept,varexplD[l])
        ##                         ndept <- ndept+1
        ##                     }
        ##             }
        ##         varexplD <- c(varD,vardept) # les prm seront dans cet ordre
                
        ##         if(nbevt==1)
        ##             {
        ##                 nvarD1 <- nvarD
        ##                 nvarD2 <- 0
        ##                 ndept1 <- ndept
        ##                 ndept2 <- 0
                        
        ##                 idcause[1,] <- 1
        ##             }
        ##         else
        ##             {
        ##                 ## variables non dep du temps sur D1 et D2
        ##                 varD1 <- intersect(varD,varexplD1)
        ##                 nvarD1 <- length(varD1)
        ##                 varD2 <- intersect(varD,varexplD2)
        ##                 nvarD2 <- length(varD2)
                        
        ##                 ## variables dep du temps sur (les seuils de) D1 et D2
        ##                 vardept1 <- intersect(vardept,varexplD1)
        ##                 ndept1 <- length(vardept1)
        ##                 vardept2 <- intersect(vardept,varexplD2)
        ##                 ndept2 <- length(vardept2)

        ##                 ## quelle variable joue sur quel evt
        ##                 idcause[1,] <- as.numeric(varexplD %in% c(varD1,vardept1))
        ##                 idcause[2,] <- as.numeric(varexplD %in% c(varD2,vardept2))
        ##             }
        ##     }

        ## ## ajouter varD a D et faire Xseuil
        ## if(nbevt==1)
        ##     {
        ##         namesevt1 <- NULL
        ##         if(nvarD>0)
        ##             {
        ##                 f <- sapply(paste("factor(",varD,")",sep=""),function(x) grepl(x,formD1,fixed=TRUE))
                        
        ##                 varDform <- varD
        ##                 if(any(f))
        ##                     {
        ##                         varDform[which(f==TRUE)] <- paste("factor(",varD[which(f==TRUE)],")",sep="")
        ##                     }
                        
        ##                 DvarD <- model.matrix(formula(paste("~",paste(varDform,collapse="+"))),data=unique(dataD[which(!is.na(dataD[,dem])),c(subject,varD)]))

        ##                 nvarD <- ncol(DvarD)-1
        ##                 nvarD1 <- nvarD

        ##                 namesevt1 <- c(namesevt1,colnames(DvarD)[-1])
                        
        ##                 D <- cbind(D[,2],DvarD[,-1,drop=FALSE])
        ##             }
        ##         else
        ##             {
        ##                 D <- matrix(D[,2],ncol=1)
        ##             }

        ##          if(ndept>0)
        ##              {
        ##                  f <- sapply(paste("factor(",vardept,")",sep=""),function(x) grepl(x,formD1,fixed=TRUE))
                        
        ##                 vardeptform <- vardept
        ##                 if(any(f))
        ##                     {
        ##                         vardeptform[which(f==TRUE)] <- paste("factor(",vardept[which(f==TRUE)],")",sep="")
        ##                     }
                      
        ##                  Xseuil <- model.matrix(formula(paste("~",paste(vardeptform,collapse="+"))),data=dataD)
        ##                  Xseuil <- Xseuil[,-1,drop=FALSE]

        ##                  ndept <- ncol(Xseuil)
        ##                  ndept1 <- ndept

        ##                  namesevt1 <- c(namesevt1,colnames(Xseuil))
        ##              }
        ##         else
        ##             {
        ##                 Xseuil <- NULL
        ##             }

        ##         ## attention : nvarD+ndept != nvarexplD a cause des factors
        ##         idcause <- matrix(1,nrow=1,ncol=nvarD+ndept) 
        ##     }
        ## else
        ##     {
        ##         if(nvarD>0)
        ##             {
        ##                 DvarD1 <- NULL
        ##                 DvarD2 <- NULL
        ##                 namesevt1 <- NULL
        ##                 namesevt2 <- NULL
                        
        ##                 if(nvarD1>0)
        ##                     {
        ##                         f <- sapply(paste("factor(",varD1,")",sep=""),function(x) grepl(x,formD1,fixed=TRUE))
                        
        ##                         varD1form <- varD1
        ##                         if(any(f))
        ##                             {
        ##                                 varD1form[which(f==TRUE)] <- paste("factor(",varD1[which(f==TRUE)],")",sep="")
        ##                             }
                                
        ##                         DvarD1 <- model.matrix(formula(paste("~",paste(varD1form,collapse="+"))),data=unique(dataD[which(!is.na(dataD[,dem])),c(subject,varD1)]))
        ##                         DvarD1 <- DvarD1[,-1,drop=FALSE]

        ##                         nvarD1 <- ncol(DvarD1)
        ##                         namesevt1 <- c(namesevt1,colnames(DvarD1))

        ##                         idcause <- matrix(0,2,nvarD1)
        ##                         idcause[1,] <- 1
        ##                     }

        ##                 if(nvarD2>0)
        ##                     {
        ##                         ## sert juste a faire nvarD2
        ##                         f <- sapply(paste("factor(",varD2,")",sep=""),function(x) grepl(x,formD2,fixed=TRUE))
                        
        ##                         varD2form <- varD2
        ##                         if(any(f))
        ##                             {
        ##                                 varD2form[which(f==TRUE)] <- paste("factor(",varD2[which(f==TRUE)],")",sep="")
        ##                             }
                                
        ##                         DvarD2 <- model.matrix(formula(paste("~",paste(varD2form,collapse="+"))),data=unique(dataD[which(!is.na(dataD[,dc])),c(subject,varD2)]))
        ##                         DvarD2 <- DvarD2[,-1,drop=FALSE]

        ##                         nvarD2 <- ncol(DvarD2)
        ##                         namesevt2 <- c(namesevt2,colnames(DvarD2))

        ##                         idcause[2,which(colnames(DvarD1 %in% DvarD2))] <- 1
        ##                     }

        ##                 DvarD2 <- NULL
        ##                 if(length(setdiff(varD2,varD1))>0)
        ##                     {
        ##                         f <- sapply(paste("factor(",setdiff(varD2,varD1),")",sep=""),function(x) grepl(x,formD2,fixed=TRUE))
                        
        ##                         varD2form <- setdiff(varD2,varD1)
        ##                         if(any(f))
        ##                             {
        ##                                 varD2form[which(f==TRUE)] <- paste("factor(",varD2form[which(f==TRUE)],")",sep="")
        ##                             }
                                
        ##                         DvarD2 <- model.matrix(formula(paste("~",paste(varD2form,collapse="+"))),data=unique(dataD[which(!is.na(dataD[,dc])),c(subject,varD2)]))
        ##                         DvarD2 <- DvarD2[,-1,drop=FALSE]

        ##                         idcause <- cbind(idcause,rbind(rep(0,ncol(DvarD2)),rep(1,ncol(DvarD2))))
        ##                     }

        ##                 DvarD <- cbind(DvarD1,DvarD2)
        ##                 nvarD <- ncol(DvarD)
                        
        ##                 D <- cbind(D[,2],DvarD)
        ##             }
        ##         else
        ##             {
        ##                 D <- matrix(D[,2],ncol=1)
        ##             }

        ##         if(ndept>0)
        ##              {
        ##                  Xseuil1 <- NULL
        ##                  Xseuil2 <- NULL
        ##                  if(ndept1>0)
        ##                      {
        ##                          f <- sapply(paste("factor(",vardept1,")",sep=""),function(x) grepl(x,formD1,fixed=TRUE))
                        
        ##                          vardept1form <- vardept1
        ##                          if(any(f))
        ##                              {
        ##                                  vardept1form[which(f==TRUE)] <- paste("factor(",vardept1[which(f==TRUE)],")",sep="")
        ##                              }
                                 
        ##                          Xseuil1 <- model.matrix(formula(paste("~",paste(vardept1form,collapse="+"))),data=dataD[which(!is.na(dataD[,dem])),])
        ##                          Xseuil1 <- Xseuil1[,-1,drop=FALSE]

        ##                          ndept1 <- ncol(Xseuil1)
        ##                          namesevt1 <- c(namesevt1,colnames(Xseuil1))

        ##                          idcause <- cbind(idcause,rbind(rep(1,ncol(Xseuil1)),rep(0,ncol(Xseuil1))))
        ##                      }
        ##                  if(ndept2>0)
        ##                      {
        ##                          f <- sapply(paste("factor(",vardept2,")",sep=""),function(x) grepl(x,formD2,fixed=TRUE))

        ##                          #colnames(dataD)[2] <- var.time[K+2]
                        
        ##                          vardept2form <- vardept2
        ##                          if(any(f))
        ##                              {
        ##                                  vardept2form[which(f==TRUE)] <- paste("factor(",vardept2form[which(f==TRUE)],")",sep="")
        ##                              }
                                 
        ##                          Xseuil2 <- model.matrix(formula(paste("~",paste(vardept2form,collapse="+"))),data=dataD[which(!is.na(dataD[,dc])),])
        ##                          Xseuil2 <- Xseuil2[,-1,drop=FALSE]

        ##                          colnames(dataD)[2] <- var.time[K+1]
                                 
        ##                          ndept2 <- ncol(Xseuil2)
        ##                          namesevt2 <- c(namesevt2,colnames(Xseuil2))

        ##                          idcause <- cbind(idcause,rbind(rep(0,ncol(Xseuil2)),rep(1,ncol(Xseuil2))))
        ##                      }
                         
        ##                  Xseuil1 <- cbind(Xseuil1,un=rep(1,sum(ni01)))
        ##                  Xseuil2 <- cbind(Xseuil2,deux=rep(2,sum(ni02)))
        ##                  Xseuil <- merge(Xseuil1,Xseuil2,by.x=c("un"),by.y=c("deux"),all=TRUE)
        ##                  Xseuil <- Xseuil[,setdiff(colnames(Xseuil),c("un","deux")),drop=FALSE]
        ##                  if(any(is.na(Xseuil)))
        ##                      {
        ##                          Xseuil <- apply(Xseuil,2,function(x) { x[which(is.na(x))] <- -3000; return(x) })
        ##                      }

        ##                  ndept <- ncol(Xseuil)
        ##              }
        ##         else
        ##             {
        ##                 Xseuil <- NULL
        ##             }
                
        ##     }
########################################
        
        ## extraire les infos de chaque sous-modele
        link <- NULL
        nodes <- NULL
        nmes <- NULL
        btot <- NULL
        idx <- NULL
        fix <- NULL
        epsY <- NULL

        nv <- rep(0,K)
        nobs <- rep(0,K+nbevt)
        idiag <- rep(0,K)
        npmtot <- rep(0,K)
        p <- rep(0,K)
        q <- rep(0,K)
        ctr <- rep(0,K)
        ncor <- rep(0,K)
        ntrtot <- rep(0,K)
        ntr <- NULL
        nalea <- rep(0,K)

        X <- vector("list",K)
        Xd <- vector("list",K)
        Y <- NULL

        namesmod <- NULL

        for(k in 1:K)
            {   
                mod <- get(paste("mod",k,sep=""))
                colx <- c(mod$call$subject,mod$Xnames2,var.time[K+1])
                namesmod <- c(namesmod,names(mod$best))


                ## formule k       
                formk <- paste(paste(mod$call$fixed[3],collapse="+"),paste(mod$call$random[2],collapse="+"),"1",sep="+")
                if(!is.null(mod$call$cor))
                    {
                        formk <- paste(formk,as.character(mod$call$cor)[2],sep="+")
                    }
                
                ## les donnees X et Y
                yk <- NULL
                xk <- NULL
                id <- NULL
                outcome <- NULL
                levid <- unique(dataY[,1])
                for(m in 1:ny[k])
                    {
                        datam <- dataY[which((dataY$processK==k) & (dataY$outcomeM==m)),]
                        colnames(datam)[which(colnames(datam)=="timeT")] <- var.time[k]
                        colnames(datam)[which(colnames(datam)=="measureY")] <- mod$Ynames[m]
                        yk <- c(yk,datam[,mod$Ynames[m]])
                        xk <- rbind(xk,model.matrix(formula(paste("~",formk)),data=datam))
                        id <- c(id,as.numeric(factor(datam[,1],levels=levid))) # ? va poser pb si pas les memes sujets pr chaque outcome !**
                        outcome <- c(outcome,rep(m,nrow(datam)))

                        ## nb mesures
                        #nmesm <- table(datam[,1])
                        #nmesm <- data.frame(id=names(nmesm),nm=as.vector(nmesm))
                        nmesm <- rle(datam[,1])
                        nmesm <- data.frame(id=nmesm[[2]],nm=nmesm[[1]])
                        if((k==1) & (m==1))
                            {
                                nmes <- nmesm
                                colnames(nmes) <- c("id","k1m1")
                            }
                        else
                            {
                                old <- colnames(nmes)
                                nmes <- merge(nmes,nmesm,by="id",all=TRUE)#,sort=FALSE)
                                colnames(nmes) <- c(old,paste("k",k,"m",m,sep=""))
                            }
                        
                        ## ntr
                        if(mod$linktype[m]==0) ntr <- c(ntr,2)
                        if(mod$linktype[m]==1) ntr <- c(ntr,4)
                        if(mod$linktype[m]==2) ntr <- c(ntr,mod$nbnodes[m]+2)
                    }

                ## mettre dans l ordre
                yx <- cbind(id,outcome,yk,xk)
                yxord <- yx[order(yx[,1],yx[,2]),]
                Y <- c(Y,as.vector(yxord[,3]))
                X[[k]] <- yxord[,-c(1,2,3),drop=FALSE]

                ## Xd
                ty <- paste("\\b",var.time[k],"\\b",sep="")
                form0 <- gsub(ty,"tdemdc",formk)
                Xd[[k]] <- model.matrix(formula(paste("~",form0)),data=dataD)

                ## infos du modele k
                nv[k] <- ncol(xk)
                nobs[k] <- nrow(xk) #mod$N[9] # plus pareil si prevalents enleves
                idiag[k] <- mod$idiag
                npmtot[k] <- length(mod$best)
                ncor[k] <- mod$N[7]
                p[k] <- sum(mod$idg)
                ctr[k] <- sum(mod$idcontr)
                q[k] <- sum(mod$idea)
                ntrtot[k] <- npmtot[k]-sum(mod$N[c(3,4,5,6,7,8)])
                nalea[k] <- mod$N[6]

                idx <- cbind(idx,rbind(mod$idg,mod$idcontr,mod$idea,mod$idcor))

                link <- c(link,mod$linktype)
                nodes <- c(nodes,as.vector(mod$linknodes))
                epsY <- c(epsY,mod$epsY)
                
                ## parametres (estimes ou fixes)
                btotk <- as.vector(mod$best)
                if(mod$N[4]>0)
                    {
                        if(!mod$idiag)
                            {
                                vark <- matrix(0,q[k],q[k])
                                vark[upper.tri(vark,diag=TRUE)] <- c(1,btotk[mod$N[3]+1:mod$N[4]])
                                vark <- t(vark)
                                vark[upper.tri(vark,diag=TRUE)] <- c(1,btotk[mod$N[3]+1:mod$N[4]])
                                
                                sigk <- sqrt(diag(vark))
                                sigcorrk <- sweep(vark,1,sigk,"/")
                                sigcorrk <- sweep(sigcorrk,2,sigk,"/")
                                diag(sigcorrk) <- sigk
                                rho <- sigcorrk[upper.tri(sigcorrk)]
                                sigcorrk[upper.tri(sigcorrk)] <- log((1+rho)/(1-rho))
                                btotk[mod$N[3]+1:mod$N[4]] <- sigcorrk[upper.tri(sigcorrk,diag=TRUE)][-1]
                            }
                        else
                            {
                                btotk[mod$N[3]+1:mod$N[4]] <- sqrt(btotk[mod$N[3]+1:mod$N[4]])
                            }

                        #if(any(!is.finite(btotk))) stop(paste("Infinite parameters for multlcmm model ",k,". Check the correlations.")) # a remettre !**
                    }
                btot <- c(btot,btotk)

                fix <- c(fix,rep(0,length(mod$best)))
                if(length(mod$call$posfix)) fix[length(fix)-length(mod$best)+eval(mod$call$posfix)] <- 1
            } # fin boucle sur k
                        
        ## pas de link beta
        if(any(ntr==4)) stop("Use only linear or splines links")
        if(any(ncor==2)) stop("Autoregressive correlations are not allowed")
        
        ## tous les parametres du modele
        if(K==1 | RE=="block-diag") nRE <- 0
        if(K>1 & RE!="block-diag") nRE <- sum(sapply(1:(length(q)-1), function(k,x) sum(x[k]*x[(k+1):length(x)]),x=q))
        bRE <- rep(0,nRE)
        BM <- "diag"
        if((K==1) | (BM=="diag")) nBM <- 0
        if((K>1) & (BM=="full"))
            {
                ncortot <- sum(ncor==1)
                nBM <- ncortot*(ncortot-1)/2
            }
        bBM <- rep(0,nBM)
        eta1 <- rep(0,1+nvarD1+ndept1)
        gamma1 <- rep(0,K*(1+nvarD1))
        btot <- c(btot,bRE,bBM,eta1,gamma1)
        if(nbevt==2)
            {
                eta2 <- rep(0,1+nvarD2+ndept2)
                gamma2 <- rep(0,K*(1+nvarD2))
                btot <- c(btot,eta2,gamma2)
            }
        fix <- c(fix,rep(0,nRE+nBM+nbevt+nvarD1+nvarD2+ndept1+ndept2+K*(nbevt+nvarD1+nvarD2))) #par defaut, fixer les prm qui ont ete fixes dans les sous-modeles et ne pas fixer les autres

        ## parametres fixes
        if(!missing(posfix))
            {
                if(missing(B)) stop("B must be specify to use posfix")

                if(!all(posfix %in% 1:length(B))) stop("posfix must contain indexes between 1 and length(B)")

                fix <- rep(0,length(B))
                fix[posfix] <- 1
            }

        
        ## rassembler nb mesures
        uniqueid <- nmes[,1]
        nmes <- cbind(nmes[,-1],ni01=as.vector(ni01))
        if(nbevt==2) nmes <- cbind(nmes,ni02=as.vector(ni02))
        nmes <- cbind(nmes,ni0=as.vector(ni0))
        nmes <- as.matrix(nmes)
        if(any(is.na(nmes))) nmes[which(is.na(nmes))] <- 0 # peut arriver pr Y
        #if(any(nmes==0)) stop("nmes=0")
        nobs[K+1] <- sum(ni01)
        if(nbevt==2) nobs[K+2] <- sum(ni02)
        niparK <- matrix(0,ns,K)
        m0 <- 0
        for(k in 1:K)
            {
                for(m in m0+1:ny[k])
                    {
                        niparK[,k] <- niparK[,k]+nmes[,m]
                    }
                m0 <- m0+ny[k]
            }

        M <- sum(ny)
        
        ## max 20 mesures
        #if(any(nmes[,M+nbevt+1]>20)) stop("The number of repeated measures for D is limited to 20")
        #if(nbevt==2) if(any(nmes[,M+1]+nmes[,M+2]>20)) stop("The number of repeated measures for D is limited to 20")
        nmesbis <- nmes[,M+1:nbevt]-nbTronc
        if(nbevt==1) if(any(nmesbis[,1]>20)) stop("The number of repeated measures for D is limited to 20")
        if(nbevt==2) if(any(nmesbis[,1]+nmesbis[,2]>20)) stop("The number of repeated measures for D is limited to 20")

        ## si valeurs initiales specifiees
        if(!missing(B))
            {
                if(length(B)!=length(btot)) stop(paste(c("B should be of length",length(btot),"and it is of length",length(B)),collapse=" "))
                
                btot <- B

                ## remplacer varcov par cholesky/corr
                iB <- 0
                qcurr <- 0
                iRE <- 0
                kBM <- 0
                if(nRE>0) vcK <- matrix(0,sum(q),sum(q))
                if(nBM>0) sigBM <- matrix(0,ncortot,ncortot)
                for(k in 1:K)
                    {
                        if(nRE==0)
                            {
                                if(q[k]>1 & idiag[k]==0)
                                    {
                                        nvc <- q[k]*(q[k]+1)/2-1
                                        vc <- matrix(0,q[k],q[k])
                                        if(nvc==0)
                                            {
                                                vc[1,1] <- 1
                                            }
                                        else
                                            {
                                                vc[upper.tri(vc,diag=TRUE)] <- c(1,B[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc])
                                                vc <- t(vc)
                                                vc[upper.tri(vc,diag=TRUE)] <-c(1,B[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc])
                                            }


                                        sigk <- sqrt(diag(vc))
                                        sigcorrk <- sweep(vc,1,sigk,"/")
                                        sigcorrk <- sweep(sigcorrk,2,sigk,"/")
                                        diag(sigcorrk) <- sigk
                                        
                                        rho <- sigcorrk[upper.tri(sigcorrk)]
                                        sigcorrk[upper.tri(sigcorrk)] <- log((1+rho)/(1-rho))
                                        
                                        btot[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc] <- sigcorrk[upper.tri(sigcorrk,diag=TRUE)][-1]
                                    }
                                
                                if(q[k]>1 & idiag[k]==1)
                                    {
                                        nvc <- q[k]
                                        btot[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:q[k]] <- sqrt(B[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:q[k]])
                                        
                                    }
                            }
                        
                        if(nRE>0)
                            {
                                nvc <- q[k]*(q[k]+1)/2-1
                                vc <- matrix(0,q[k],q[k])
                                if(nvc==0)
                                    {
                                        vc[1,1] <- 1
                                    }
                                else
                                    {
                                        vc[upper.tri(vc,diag=TRUE)] <- c(1,B[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc])
                                        vc <- t(vc)
                                        vc[upper.tri(vc,diag=TRUE)] <- c(1,B[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc])
                                    }
                                
                                vcK[qcurr+1:q[k],qcurr+1:q[k]] <- vc

                                if(k>1)
                                    {
                                        nbr <- sum(q[1:(k-1)])
                                        nbc <- q[k]
                                        vcK[1:nbr,qcurr+1:nbc] <- B[sum(npmtot)+iRE+1:(nbr*nbc)]
                                        vcK <- t(vcK)
                                        vcK[1:nbr,qcurr+1:nbc] <- B[sum(npmtot)+iRE+1:(nbr*nbc)]
                                        iRE <- iRE+nbr*nbc
                                    }
                            }

                        if((nBM>0) & (ncor[k]==1))
                            {
                                kBM <- kBM+1
                                sigBM[kBM,kBM] <- B[iB+p[k]-1+ctr[k]*(ny[k]-1)+nvc+1]
                            }
                        
                        iB <- iB+npmtot[k]
                        qcurr <- qcurr+q[k]
                    }

                
                if(nRE>0)
                    {
                        #cat("varcov : \n ")
                        #print(vcK)

                        sigk <- sqrt(diag(vcK))
                        sigcorrk <- sweep(vcK,1,sigk,"/")
                        sigcorrk <- sweep(sigcorrk,2,sigk,"/")
                        diag(sigcorrk) <- sigk
                        
                                        #cat("sigma et  corr : \n ")
                                        #print(sigcorrk)

                        rho <- sigcorrk[upper.tri(sigcorrk)]
                        sigcorrk[upper.tri(sigcorrk)] <- log((1+rho)/(1-rho))
                        
                                #cat("sigma et prm corr non contraints: \n ")
                                #print(sigcorrk)
                                
                        iB <- 0
                        qcurr <- 0
                        iRE <- 0
                        for(k in 1:K)
                            {
                                if(q[k]>0)
                                    {
                                        isigcorrk <- qcurr+1:q[k]
                                        if(q[k]>1)
                                            {
                                                btot[iB+p[k]-1+ctr[k]*(ny[k]-1)+1:nvc] <- sigcorrk[isigcorrk,isigcorrk][upper.tri(sigcorrk[isigcorrk,isigcorrk],diag=TRUE)][-1]
                                            }
                                        
                                        iB <- iB+npmtot[k]
                                    }
                                
                                        # pas possible de faire idiag=TRUE avec cov interdim

                                if(k>1)
                                    {
                                        nbr <- sum(q[1:(k-1)])
                                        nbc <- q[k]
                                        btot[sum(npmtot)+iRE+1:(nbr*nbc)] <- sigcorrk[1:nbr,nbr+1:nbc]
                                        iRE <- iRE+nbr*nbc
                                    }
                                
                                qcurr <- qcurr+q[k]
                            }
                        
                    }
                
                if(nBM>0)
                    {
                        sigBM[upper.tri(sigBM)] <- B[sum(npmtot)+nRE+1:nBM]
                        corBM <- sweep(sigBM,1,diag(sigBM),"/")
                        corBM <- sweep(corBM,2,diag(sigBM),"/")
                        corBM <- corBM[upper.tri(corBM)]
                        prmBM <- log((1+corBM)/(1-corBM))
                        btot[sum(npmtot)+nRE+1:nBM] <- prmBM
                    }
            }
        
        ## differencier prm estimes et prm fixes
        b <- btot[fix==0]
        bfix <- 0
        if(length(fix)) bfix <- btot[fix==1]
        
        nef <- p-1
        ncontr <- ctr*(ny-1)
        nea <- q

        ch <- 0
        D <- as.matrix(D[,-1,drop=FALSE])
        if(is.null(Xseuil)) Xseuil <- matrix(0,0,0)
        Xseuil <- as.matrix(Xseuil)
        idcause <- as.vector(t(idcause))
        entreRetard <- as.numeric(delayed)
        algo <- 0
        cond <- 0
        decesExact <- 0

        #browser()
        #return(list(dataD1=dataD1,D=D,Xseuil=Xseuil,Xd=Xd[[1]],X=X[[1]],idcause=idcause))
        ## evaluer une fois la vraisemblance si maxiter=0
        if(maxiter==0)
            {
               # vrais <- vraisR(b,bfix,fix,Y,X,Xd,idD,D,Xseuil,nmes,nv,idx,idiag,ncor,ny,nalea,ntr,link,nodes,epsY,nRE,nBM,ch,nbevt,idcause,decesExact,algo,cond,entreRetard,discretise,nbTronc)
                vrais <- RloglikJointMult(b,bfix,fix,Y,X,Xd,idD,D,Xseuil,nmes,nv,idx,idiag,ncor,ny,nalea,ntr,link,nodes,epsY,nRE,nBM,nbevt,idcause,entreRetard,nbTronc)

                res <- list(istop=2,ni=0,loglik=vrais,b=b,v=rep(NA,length(b)),convcrit=rep(NA,3),time=0,nproc=1)
            }
        else
            {
                res <- marqLevAlg(b=b,m=FALSE,fn=RloglikJointMult,gr=NULL,hess=NULL,
                                  epsa=eps[1],epsb=eps[2],epsd=eps[3],
                                  digits=8,print.info=verbose,blinding=FALSE,
                                  multipleTry=25,file=file,
                                  nproc=nproc,maxiter=maxiter,minimize=FALSE, 
                                  bfix0=bfix,fix0=fix,Y0=Y,X0=X,Xd0=Xd,D0=D,
                                  idD0=idD,nbevt0=nbevt,idcause0=idcause,
                                  Xseuil0=Xseuil,nmes0=nmes,nv0=nv,idx0=idx,
                                  idiag0=idiag,ncor0=ncor,ny0=ny,
                                  nalea0=nalea,ntr0=ntr,link0=link,
                                  nodes0=nodes,epsY0=epsY,nRE0=nRE,nBM0=nBM,
                                  entreRetard0=entreRetard,nbTronc0=nbTronc)
                
                res <- list(istop=res$istop,ni=res$ni,loglik=res$fn.value,
                            b=res$b,v=res$v,convcrit=c(res$ca,res$cb,res$rdm),
                            time=res$time,nproc=nproc)
            }

        bopt <- rep(0,length(btot))
        bopt[fix==0] <- res$b
        bopt[fix==1] <- bfix
        
        res$bopt <- bopt
        npmMM <- sum(npmtot)
        names(res$bopt)[1:npmMM] <- namesmod #noms Y
        if(nRE>0) names(res$bopt)[npmMM+1:nRE] <- paste("cov",1:nRE,sep="")
        if(nBM>0) names(res$bopt)[npmMM+nRE+1:nBM] <- paste("covBM",1:nBM,sep="")
        if(nbevt==1)
            {
                idcause <- matrix(idcause,nrow=nbevt,byrow=TRUE)
                if(nvarD1==0)
                    {
                        if(ndept1==0)
                            {
                                names(res$bopt)[npmMM+nRE+nBM+1] <- "Threshold"
                                names(res$bopt)[npmMM+nRE+nBM+1+1:K] <- paste("Gamma",1:K,sep="")
                            }
                        else
                            {
                                names(res$bopt)[npmMM+nRE+nBM+1:(ndept1+1)] <- paste("Threshold",c("",namesevt1))
                                names(res$bopt)[npmMM+nRE+nBM+1+ndept1+1:K] <- paste("Gamma",1:K,sep="")
                            }
                    }
                else
                    {
                        if(ndept1==0)
                            {
                                names(res$bopt)[npmMM+nRE+nBM+1:(nvarD1+1)] <- paste("Threshold",c("",namesevt1),sep=" ")
                                names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+1:(K*(1+nvarD1))] <- paste(rep(paste("Gamma",1:K,sep=""),each=nvarD1+1),c("",namesevt1),sep=" ")
                            }
                        else
                            {  
                                names(res$bopt)[npmMM+nRE+nBM+1:(1+nvarD1+ndept1)] <- paste("Threshold",c("",namesevt1),sep=" ")
                                names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+1:(K*(1+nvarD1))] <- paste(rep(paste("Gamma",1:K,sep=""),each=nvarD1+1),c("",namesevt1[1:nvarD1]),sep=" ")
                            }
                    }
            }
        if(nbevt==2)
            {
                idcause <- matrix(idcause,nrow=nbevt,byrow=TRUE)

                ## noms dem
                names(res$bopt)[npmMM+nRE+nBM+1] <- paste("Threshold",dem)
                names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+seq(1,K*(nvarD1+1),by=nvarD1+1)] <- paste("Gamma",dem,1:K,sep=" ")
                
                if(nvarD1!=0)
                    {
                        names(res$bopt)[npmMM+nRE+nBM+1+1:nvarD1] <- paste("Threshold",dem,namesevt1[1:nvarD1])
                        names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+1:(K*(nvarD1+1))] <- paste("Gamma",dem,paste(rep(1:K,each=nvarD1+1),c("",namesevt1[1:nvarD1]),sep=" "),sep=" ")
                    }
                      
                if(ndept1!=0) names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+1:ndept1] <- paste("Threshold",dem,namesevt1[nvarD1+1:ndept1])

                ## noms dc
                names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+K*(nvarD1+1)+1] <- paste("Threshold",dc)
                names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+K*(nvarD1+1)+1+nvarD2+ndept2+seq(1,K*(nvarD2+1),by=nvarD2+1)] <- paste("Gamma",dc,1:K,sep=" ")

                if(nvarD2!=0)
                    {
                        names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+K*(nvarD1+1)+1+1:nvarD2] <- paste("Threshold",dc,namesevt2[1:nvarD2])
                        names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+K*(nvarD1+1)+1+nvarD2+ndept2+1:(K*(nvarD2+1))] <- paste("Gamma",dc,paste(rep(1:K,each=nvarD2+1),c("",namesevt2[1:nvarD2]),sep=" "),sep=" ")
                    }
                      
                if(ndept2!=0) names(res$bopt)[npmMM+nRE+nBM+1+nvarD1+ndept1+K*(nvarD1+1)+1+nvarD2+1:ndept2] <- paste("Threshold",dc,namesevt2[nvarD2+1:ndept2])
            }

        res$nef <- nef  # nef ne contient pas ncontr
        res$ncontr <- ncontr
        res$nea <- nea
        res$nvc <- nea-1
        res$idiag <- idiag
        res$nvc[which(res$idiag==0)] <- nea[which(res$idiag==0)]*(nea[which(res$idiag==0)]+1)/2-1
        res$ntr <- ntr
        res$ntrtot <- ntrtot
        res$ncor <- ncor
        res$nalea <- nalea
        res$ny <- ny
        res$link <- link
        res$nodes <- nodes
        res$nRE <- nRE
        res$nBM <- nBM
        res$varexplD <- varexplD
        res$nvarD <- c(nvarD1,nvarD2)
        res$ndept <- c(ndept1,ndept2)
        res$idcause <- idcause
        res$call <- match.call()
        res$fix <- fix
        res$Ynames <- Ynames
        res$Xnames <- Xnames
        res$Dnames <- c(dem,dc)
        res$ns <- ns
        res$events <- table(D[,1])
        ## 1 mesure en trop si nbTronc=1
        if(any(nbTronc>0))
            {
                if(any(nbTronc[,1]>0))
                    {
                        nmes[,M+1] <- nmes[,M+1] - nbTronc[,1]
                    }
                if((nbevt==2) & (any(nbTronc[,nbevt]>0)))
                    {
                        nmes[,M+nbevt] <- nmes[,M+nbevt] - nbTronc[,nbevt]
                    }
            }
        res$nbmes <- apply(nmes[,-ncol(nmes)],2,mean)
        res$entreRetard <- entreRetard
        res$breaks <- breaks
        
        ## remplacer corr par varcov
        nvc <- res$nvc
        U <- matrix(0,sum(nea),sum(nea))
        iB <- 0
        imatB <- 0
        iRE <- 0
                                        #browser()
        for(k in 1:K)
            {
                if(nvc[k]==0)
                    {
                        U[imatB+1,imatB+1] <- 1
                    }
                else
                    {
                        if(nvc[k]>0 & idiag[k]==1)
                            {
                                diag(U[imatB+1:nea[k],imatB+1:nea[k]]) <- c(1,bopt[iB+nef[k]+ncontr[k]+1:nvc[k]])
                            }
                        else
                            {
                                U[imatB+1:nea[k],imatB+1:nea[k]][upper.tri(U[imatB+1:nea[k],imatB+1:nea[k]],diag=TRUE)] <- c(1,bopt[iB+nef[k]+ncontr[k]+1:nvc[k]])
                                
                            }
                    }
              
                                
                if((imatB>0) & (nRE>0))
                    {
                        qavt <- sum(nea[1:(k-1)])
                        U[1:qavt,imatB+1:nea[k]] <- bopt[npmMM+iRE+1:(qavt*nea[k])]
                        iRE <- iRE+qavt*nea[k]
                    }
                
                imatB <- imatB+nea[k]
                iB <- iB+npmtot[k]
            }
        
        sig <- abs(diag(U))
        U <- -1+2*(exp(U)/(1+exp(U)))
        corr <- t(U)
        corr[upper.tri(U)] <- U[upper.tri(U)]
        diag(corr) <- 1
        
        vc <- sweep(corr,1,sig,"*")
        vc <- sweep(vc,2,sig,"*")
        
        res$VRE <- vc
        
        res$corRE <- corr
        

        iB <- 0
        imatB <- 0
        qavt <- 0
        covRE <- NULL
        ll <- vector("list",K)
        for(k in 1:K)
            {
                if(res$nea[k]>0)
                    {
                        if(res$nvc[k]>0)
                            {
                                vtmp <- res$VRE[imatB+1:res$nea[k],imatB+1:res$nea[k]]
                                if(res$idiag[k]==0)
                                    {
                                        res$bopt[iB+res$nef[k]+res$ncontr[k]+1:res$nvc[k]] <- vtmp[upper.tri(vtmp,diag=TRUE)][-1]
                                    }
                                else
                                    {
                                        res$bopt[iB+res$nef[k]+res$ncontr[k]+1:res$nvc[k]] <- diag(vtmp)[-1]
                                    }
                            }
                
                        if((imatB>0) & (nRE>0))
                            {
                                qavt <- sum(res$nea[1:(k-1)])
                                covRE <- c(covRE,as.vector(res$VRE[1:qavt,imatB+1:nea[k]]))
                            }
                        
                        imatB <- imatB+res$nea[k]
                    }
                
                llk <- get(paste("mod",k,sep=""))
                llk$binit <- llk$best
                llk$best <- res$bopt[iB+1:(res$nef[k]+res$ncontr[k]+res$nvc[k]+res$ncor[k]+res$ny[k]+res$nalea[k]+res$ntrtot[k])]
                ll[[k]] <- llk
                iB <- iB+res$nef[k]+res$ncontr[k]+res$nvc[k]+res$ncor[k]+res$ny[k]+res$nalea[k]+res$ntrtot[k]
            }
        
        res$mod <- ll
        
        if(nRE>0)
            {
                res$bopt[iB+1:nRE] <- covRE
            }
        
        class(res) <- "JointMult"

        if((verbose==TRUE) & (file!=""))
            {
                if(length(grep(".txt",file)))
                    {
                        ressscall <- res
                        ressscall$call <- NULL
                        ressscall$mod <- NULL
                        fileres <- sub(".txt","_last.txt",file)
                        dput(ressscall,file=fileres)
                    }
            }



        if(pred==TRUE)
            {
                ## indicateur de dimension et d'outcome
                dimension <- rep(1:K,nobs[1:K])
                outcome <- NULL
                m0 <- 0
                for(k in 1:K)
                    {
                        outcome <- c(outcome,unlist(apply(nmes[,m0+1:ny[k],drop=FALSE],1,function(n) rep(m0+1:ny[k],n))))
                        m0 <- m0+ny[k]
                    }
                
                
                ## transformer les obs par H
                HY <- NULL
                iB <- 0
                m0 <- 0
                for(k in 1:K)
                    {
                        jntr <- 0
                        for(m in m0+1:ny[k])
                            {
                                HY[which(outcome==m)] <- transfo_spl(Y[which(outcome==m)], get(paste("mod",k,sep=""))$linknodes[,m-m0],res$bopt[iB+res$nef[k]+res$ncontr[k]+res$nvc[k]+res$ncor[k]+res$ny[k]+res$nalea[k]+jntr+1:res$ntr[m]])
                                jntr <- jntr+res$ntr[m]
                            }
                        iB <- iB+res$nef[k]+res$ncontr[k]+res$nvc[k]+res$ncor[k]+res$ny[k]+res$nalea[k]+res$ntrtot[k]
                        m0 <- m0+ny[k]
                    }
               
                ## predictions conditionnelles (sachant Y slt)
                nmes <- nmes[,1:sum(ny),drop=FALSE]
                pred_condY <- predcondY(b0=b,bfix0=bfix,fix0=fix,
                                        Y0=Y,X0=X,
                                        nmes0=nmes,nv0=nv,
                                        idx0=idx,idiag0=idiag,
                                        ncor0=ncor,ny0=ny,nalea0=nalea,
                                        ntr0=ntr,link0=link,nodes0=nodes,
                                        nRE0=nRE,nBM0=nBM,chol0=ch)
                sujet <- as.vector(apply(niparK,2,function(n) rep(uniqueid,n)))
                #browser()
                predHY <- data.frame(dimension,outcome,as.vector(unlist(sujet)),HY,pred_condY,Y)
                colnames(predHY) <- c("dimension","outcome","subject","obs","pred","Y")
                rownames(predHY) <- NULL
                
                res <- as.list(res)
                res$predcondY <- predHY
            }
        
        
        
        return(res)
    }




########################

## dtest <- data.frame(numero=c(1,2,2,3,3,4,4,4,5,5,5,6,6,6,6),
##                     temps=c(0.6,0.5,1,0.5,1,0.2,0.4,0.6,0.3,0.5,0.7,0.5,0.6,0.8,0.9),
##                     ty1=c(0.6,0.5,1,0.5,1,0.2,0.4,0.6,0.3,0.5,0.7,0.5,0.6,0.8,NA),
##                     ty2=c(0.6,0.5,1,0.5,1,0.2,0.4,0.6,0.3,0.5,0.7,0.5,0.6,0.8,NA),
##                     ty3=c(0.6,0.5,1,0.5,1,0.2,0.4,0.6,0.3,0.5,0.7,0.5,0.6,0.8,NA),
##                     tdiag1=c(0.6,0.5,1,0.5,1,0.2,0.4,0.5,0.4,0.6,0.8,0.5,0.6,0.8,NA),
##                     tdiag2=c(0.6,0.5,1,0.5,1,0.2,0.4,0.5,0.2,0.4,0.7,0.5,0.6,0.8,NA),
##                     tdeces=c(NA,0.5,1,0.5,1,0.2,0.4,0.6,0.2,0.4,0.6,0.5,0.6,0.8,0.9),
##                     y1=c(1.1,12,12,11,11,10,9,8,7,7,7,13,13,12,NA),
##                     y2=c(2.2,11,11,12,12,10,10,9,6,7,6,3,4,9,NA),
##                     y3=c(3.3,13,13,13,13,10,NA,NA,12,10,12,3,NA,9,NA),
##                     diag1=c(1,0,0,0,1,0,0,0,0,1,1,0,0,0,NA),
##                     diag2=c(0,0,0,0,1,0,0,1,0,0,0,0,0,0,NA),
##                     deces=c(NA,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
##                     X=c(0,0,0,1,1,0,0,0,1,1,1,0,0,0,0),
##                     Xt=c(0,0,1,0,1,0,0,0,0,0,1,0,0,1,1))


## m1 <- multlcmm(y1+y2~temps,random=~1,subject="numero",data=dtest)
## m2 <- multlcmm(y3~temps+X,random=~1+temps,subject="numero",data=dtest)

## JointMult(Y=list(m1,m2),D=diag1~Xt,var.time="temps",data=dtest)

## JointMult(Y=list(m1,m2),D=deces~X,var.time="temps",data=dtest,breaks=c(0,0.3,0.5,0.8,1),B=c(-0.3,0.1,1.5,9,4,8,2,0.8,0,-0.5,1,0.1,5,9,1,0.2,0.1,0,0.1,0.05))

