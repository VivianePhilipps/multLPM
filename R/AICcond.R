
RloglikUACV <- function(b0, bfix0, fix0, Y0, X0, Xd0, idD0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, epsY0, nRE0, nBM0, nbevt0, idcause0, entreRetard0, nbTronc0)
{
  .Call("loglikUACV", b0, bfix0, fix0, Y0, X0, Xd0, idD0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, epsY0, nRE0, nBM0, nbevt0, idcause0, entreRetard0, nbTronc0)
}


## AICcond pour calculer le fit de la partie survie ##
#' Computation of the Conditionnal Akaike Information Criterion (AICcond) for a joint model estimated by JointMult function
#'
#' @param model a JointMult model
#' @param Y a list of \code{multlcmm} objects
#' @param D a list of two-sided formula defining the event part of the model
#' @param data data.frame containing the observations and variables
#' @param var.time a character vector indicating the name of the time variables
#' @param RE an indicator of the random effect structure between dimensions
#' @param BM an indicator of the correlation of the Brownian motions
#' @param B vector containing initial values for the parameters
#' @param posfix optional vector specifying the indices in vector B of the parameters that are not estimated
#' @param breaks optional vector specifying the break points in the case where the event time is discretized
#' @param delayed logical vector indicating if delayed entry should be accounted for
#' @return the value of the Universal Approximate Cross Validation criterion
#' @author Cecile Proust-Lima and Viviane Philipps
#' @references
#' 
#' @export

AICcond <- function(model,Y,D,data,var.time,RE="block-diag",BM="diag",B,posfix,breaks=NULL,delayed=TRUE)
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
  if(!is.null(breaks)){if(!is.vector(breaks)) stop("breaks should be a vector")}
  if(!all(is.logical(delayed))) stop("delayed should be a logical vector")
  
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
        
        dataD1bis <- merge(dataD1entry,dataD1measure,all=TRUE)
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
      dataD2 <- data[,c(subject,entry,var.time[K+2],dc,setdiff(unique(c(varexplD2,nomsX)),c(var.time,entry[2])))]
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
        
        dataD2bis <- merge(dataD2entry,dataD2measure,all=TRUE)
        
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
      
      ## matrice sujet indicateur evt
      D <- merge(t1max,t2max,all=FALSE) #selection ici
      evt <- (D[,"t1max"]<=D[,"t2max"])*D[,"d1"]+2*(D[,"t1max"]>=D[,"t2max"])*D[,"d2"]
      evt[which(evt==3)] <- 1 # ? !** si ex aequo, on prend D1
      Tevt <- apply(D[,c("t1max","t2max")],1,min) # ? a voir si ok  si evt=0!**
      D <- data.frame(D[,1],evt,Tevt)
      colnames(D) <- c(colnames(dataY)[1],"evt","Tevt")
      D <- D[order(D[,1]),]
      
      ## dataD1 et dataD2 sans les temps post-diag
      dataD1bis <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))
      dataD2bis <- as.data.frame(matrix(NA,0,ncol(dataD2),dimnames=list(NULL,colnames(dataD2))))
      ni01 <- NULL
      ni02 <- NULL
      for(i in 1:nrow(D))
      {
        if(D[i,2]==1)
        {
          dataD1bis <- rbind(dataD1bis,dataD1[which(dataD1[,1]==D[i,1]),])
          
          di2 <- dataD2[which((dataD2[,1]==D[i,1]) & (dataD2[,2]<=D[i,3])),] # ? a voir si < ou <= !**
          di2[,3] <- 0
          dataD2bis <- rbind(dataD2bis,di2)
        }
        if(D[i,2]==2)
        {
          dataD2bis <- rbind(dataD2bis,dataD2[which(dataD2[,1]==D[i,1]),])
          
          di1 <- dataD1[which((dataD1[,1]==D[i,1]) & (dataD1[,2]<=D[i,3])),]
          di1[,3] <- 0
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
      
      
      tmax <- merge(t2max,t1max,all=FALSE) #selection ici
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
    dataD1$tdemdc <- dataD1[,var.time[K+1]]
    dataD2$tdemdc <- dataD2[,var.time[K+2]]
    dataD <- merge(dataD1,dataD2,all=TRUE)
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
    modmat1 <- model.matrix(form1,data=dataD)[,-1,drop=FALSE]
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
      modmat2 <- model.matrix(form2,data=dataD)[,-1,drop=FALSE]
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
      Xseuil <- merge(Xseuil1,Xseuil2,by.x=c("un"),by.y=c("deux"),all=TRUE)
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
    for(m in 1:ny[k])
    {
      datam <- dataY[which((dataY$processK==k) & (dataY$outcomeM==m)),]
      colnames(datam)[which(colnames(datam)=="timeT")] <- var.time[k]
      colnames(datam)[which(colnames(datam)=="measureY")] <- mod$Ynames[m]
      yk <- c(yk,datam[,mod$Ynames[m]])
      xk <- rbind(xk,model.matrix(formula(paste("~",formk)),data=datam))
      id <- c(id,as.numeric(factor(datam[,1]))) # ? va poser pb si pas les memes sujets pr chaque outcome !**
      outcome <- c(outcome,rep(m,nrow(datam)))
      
      ## nb mesures
      nmesm <- table(datam[,1])
      nmesm <- data.frame(id=names(nmesm),nm=as.vector(nmesm))
      if((k==1) & (m==1))
      {
        nmes <- nmesm
        colnames(nmes) <- c("id","k1m1")
      }
      else
      {
        old <- colnames(nmes)
        nmes <- merge(nmes,nmesm,by="id",all=TRUE,sort=FALSE)
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
  fix <- model$fix
  
  
  ## a voir : demander 1 mesure a chaque Y ?
  
  ## rassembler nb mesures
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
  nmesbis <- nmes[,M+1:nbevt]-nbTronc
  if(nbevt==1) if(any(nmesbis[,1]>20)) stop("The number of repeated measures for D is limited to 20")
  if(nbevt==2) if(any(nmesbis[,1]+nmesbis[,2]>20)) stop("The number of repeated measures for D is limited to 20")
  
  nef <- p-1
  ncontr <- ctr*(ny-1)
  nea <- q
  
  ch <- 0
  D <- as.matrix(D[,-1,drop=FALSE])
  if(is.null(Xseuil)) Xseuil <- matrix(0,0,0)
  Xseuil <- as.matrix(Xseuil)
  idcause <- as.vector(t(idcause))
  entreRetard <- as.numeric(delayed)
  
  
  ## prm estimes et prm fixes
  b <- model$b
  bfix <- 0
  if(any(model$fix==1)) bfix <- model$bopt[which(model$fix==1)]
  
  
  ## vraisemblances conditionnelles ###
  
  ## vrais cond (P(dem/Y)) a l'optimum
  vrais2i <- RloglikUACV(b0=b, bfix0=bfix, fix0=fix, Y0=Y, X0=X, Xd0=Xd, idD0=idD, D0=D, Xseuil0=Xseuil, nmes0=nmes, nv0=nv, idx0=idx, idiag0=idiag, ncor0=ncor, ny0=ny, nalea0=nalea, ntr0=ntr, link0=link, nodes0=nodes, epsY0=epsY, nRE0=nRE, nBM0=nBM, nbevt0=nbevt, idcause0=idcause, entreRetard0=entreRetard, nbTronc0=nbTronc) # vrais conditionnelles puis vrais totales (condi1,condi2,.., toti1, toti2...)

  
  
  ### AICcond ###
  npm <- length(b) # nbre de prms du modele complet
  npmMM <- sum(npmtot) # nbre de prms des modeles longitudinaux
  AICcond <- -2*sum(vrais2i[1:(length(vrais2i)/2)])+2*(npm-npmMM)
  
  
  res <- list(AICcond=AICcond,vrais2i=vrais2i,npm=npm,npmtot=npmtot)
  
  
  
  return(res)
}

