RprobaEC <- function(Bs0, Bst0, b0, Y0, X0, Xd0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, nRE0, nBM0)
{
  .Call("probaEC", Bs0, Bst0, b0, Y0, X0, Xd0, D0, Xseuil0, nmes0, nv0, idx0, idiag0, ncor0, ny0, nalea0, ntr0, link0, nodes0, nRE0, nBM0)
} # appel de fonction probaEC codee en C++



## calcul de probabilite d'evenement ds intervalle de prediction pr un sujet specifique, cas continu ##
#' Probability of continuous event in the prediction interval ]s,s+t] for a specific subject
#'
#' @param num subject ID for whom we want to predict
#' @param s landmark s, lower bound of the prediction interval
#' @param t horizon t of prediction defining the prediction interval
#' @param model a JointMult model
#' @param Y a list of \code{multlcmm} objects. Each multlcmm object defines the outcomes and longitudinal structure of one dimension of the model. The number of outcomes, covariates or random effects can differ between dimensions.
#' @param D a list of two-sided formula defining the event part of the model. The left side should be survBreak(T0,T,Event) with T0 the entry time, T the event time and Event the event indicator
#' @param data data.frame containing the observations and variables
#' @param ID name of the subjects'ID variable
#' @param var.time a character vector indicating the name of the time variable of each dimension. The scales of these different time variables should be the same.
#' @param RE an indicator of the random effect structure between dimensions. Should be either "block-diag" for independent random effects between dimensions (the internal structure being defined in the multlcmm object) or "full" for correlated random effects between dimensions. Default is "block-diag".
#' @param BM in the case where Brownian motions are included in the multlcmm objects, an indicator of the correlation of the Brownian motions. Should be "diag" for independence, or "full" for correlated Brownian motions between dimensions. Default is "diag".
#' @param breaks a vector specifying the break points for event time's discretization
#' @param B optional vector containing values for the parameters
#' @return the probability of continuous event in the prediction interval ]s,s+t] for a specific subject
#' @author Cécile Proust-Lima, Viviane Philipps and Tiphaine Saulnier
#' @references
#' Proust-Lima C, Philipps V, Dartigues JF.
#' A joint model for multiple dynamic processes and clinical endpoints: application to Alzheimer's disease.
#' arXiv preprint arXiv:1803.10043, 2018.
#' 
#' @export




PredEC <- function(num,s,t,model,Y,D,data,ID,var.time,RE="block-diag",BM="diag",breaks,B)
{
  #------------------------------#
  #- verification des arguments -#
  #------------------------------#
  
  if(missing(num)) stop("num is missing")
  if(missing(s)) stop("s is missing")
  if(missing(t)) stop("st is missing")
  if(missing(model)) stop("model is missing")
  if(missing(Y)) stop("Y is missing")
  if(missing(D)) stop("D is missing")
  if(missing(data)) stop("data is missing")
  if(missing(ID)) stop("ID is missing")
  if(missing(var.time)) stop("var.time is missing")
  if(missing(breaks)) stop("breaks is missing")
  
  if(!is.list(Y)) stop("Y should be a list of multlcmm objects")
  if(!(all(sapply(Y,class)=="multlcmm"))) stop("Y should only contain multlcmm objects")
  if(!is.list(D)) stop("D should be a list of 1 formula")
  if(length(D)>1) stop("D should be of length 1")
  if(!(all(sapply(D,class)=="formula")))  stop("D should only contain two-sided formula")
  if(!length(grep("survBreak",D[[1]][2]))) stop("Left side of formula in D should be 'survBreak()'")
  if(!is.data.frame(data)) stop("data should be a data frame")
  if((!nrow(data)) | (!ncol(data))) stop("data is empty")
  if(!(RE %in% c("block-diag","full"))) stop("RE should be either 'block-diag' or 'full'")
  if(!(BM %in% c("diag","full"))) stop("BM should be either 'diag' or 'full'")
  if(!is.null(breaks)){if(!is.vector(breaks)) stop("breaks should be a vector")}
  
  
  if(!(num %in% data[,ID])) stop("num should be a subject in data")
  if(!(s %in% breaks) || (s==breaks[1])) stop("landmark s should be a bound of discretization intervals (except first bound), meaning among :",paste(breaks[-1],collapse =", ") ) 
  st <- s+t
  if(!(st %in% breaks)) stop("landmark+horizon s+t should be a bound of discretization intervals, meaning among :",paste(breaks[-1],collapse =", ") )
  if(!(s < st)) stop("s should be strictly inferior to s+t")
  # pblm limite a breaks, on ne va pas plus loin
  
  
  
  #####################
  ##### donnees Y #####
  #####################
  
  
  K <- length(Y) #nbre de dimensions
  if(length(var.time)==1) var.time <- rep(var.time,K) #variable de tps de chq dim
  dataY <- NULL
  ny <- rep(NA,K) #initialisation vect de taille K
  Ynames <- vector("list",K) #initialisation K listes
  Xnames <- vector("list",K)
  nomsX <- unique(unlist(sapply(Y,function(x) x$Xnames2))) # noms des VarExplLong EF ttes dims confondues et sans repet
  
  for(k in 1:K) #pr chq dimension
  {
    ### modele k
    if(length(Y[[k]]$call)) # si la dimension a ete estimee, ie. length != 0
    {
      z <- Y[[k]]$call  # z = multlcmm ou
      z$data <- data        # data = data
      z$maxiter <- 0        # maxiter = 0
      z$B <- Y[[k]]$best    # B = prms de la dimension
      z$verbose <- FALSE    # verbose = FALSE
      mod <- eval(z)        # descriptif de la dimension : nbre obs/sujets/prms, type de link functions, loglik, AIC
    }
    else  # si la dimension n a pas ete estimee
    {
      mod <- eval(Y[[k]])    
    }
    assign(paste("mod",k,sep=""),mod) # cree nv element ds RData, nomme modk dt valeur = mod
    
    subject <- mod$call$subject # nom de variable indicatrice de l identifiant des sujets
    if(k>1){if(subject != colnames(dataY)[1]) stop("Subject variable should be the same for all multlcmm models")}
    
    colx <- c(subject,mod$Xnames2) # nom de variable identifiant + noms des VarExpLong EF
    
    ny[k] <- length(mod$Ynames) # nbre de marqueurs dans la dimension
    
    Ynames[[k]] <- mod$Ynames # noms des marqueurs de la dimension
    Xnames[[k]] <- mod$Xnames2 # noms des VarExpLong EF
    
    ## donnees km
    for(m in 1:ny[k])  # pr chq marqueur de la dimension
    {
      ## data frame de l outcome m
      colx <- c(subject,nomsX,mod$Ynames[m]) # vect = noms des variables identifiant, VarExpLong EF (ttes dims) et du marqueur
      if(length(mod$na.action[[m]])) # si on supprime les donnees manquantes
      {
        datam <- data[-mod$na.action[[m]],colx,drop=FALSE]
      }
      else  # si on garde les donnees manquantes
      {
        datam <- data[,colx,drop=FALSE]  # datam = donnees id sujets et varexpl EF     #drop = reste un dataframe m eme si 1 colonne
      }
      
      datam <- datam[order(datam[,1],datam[,var.time[k]]),,drop=FALSE]  # donnees ordonnees par sujet et par temps
      old <- colnames(datam) 
      datam$processK <- k # dim
      datam$outcomeM <- m # marqueur
      datam <- datam[,c(old[1],"processK","outcomeM",old[-1])] # ordre des colonnes : sujet, dim, marqueur, varexplEF
      
      
      if((k==1) & (m==1)) # si 1ere dimension, 1er marqueur
      {
        dataY <- datam
        colnames(dataY)[which(colnames(dataY)==var.time[k])] <- "timeT"   # variable de temps nommee timeT
        colnames(dataY)[which(colnames(dataY)==Ynames[[k]][m])] <- "measureY"  # variable des donnees du marqueur nommee measureY
      }
      else # si pas 1er marqueur de 1ere dimension
      {
        colnames(datam)[which(colnames(datam)==var.time[k])] <- "timeT"
        colnames(datam)[which(colnames(datam)==Ynames[[k]][m])] <- "measureY"
        Xplus <- setdiff(colnames(datam),colnames(dataY))  # renvoie les noms des colonnes de datam qu il n y a pas dans dataY
        if(length(Xplus)) # si il y en a, alors ces colonnes vont être ajoutees a dataY
        {
          for(l in 1:length(Xplus))
          {
            old <- colnames(dataY)
            dataY <- cbind(dataY,NA)
            colnames(dataY) <- c(old,Xplus[l])
          }
        }
        
        Xmqt <- setdiff(colnames(dataY),colnames(datam))  # renvoie les noms des colonnes de dataY qu il n y a pas dans datam
        if(length(Xmqt)) # si il y en a, alors ces colonnes vont être ajoutees a datam
        {
          for(l in 1:length(Xmqt))
          {
            old <- colnames(datam)
            datam <- cbind(datam,NA)
            colnames(datam) <- c(old,Xmqt[l])
          }
          datam <- datam[,colnames(dataY),drop=FALSE] # reordonner les colonnes dans le m eme ordre que dataY
        }
        
        dataY <- rbind(dataY,datam) # ajout des donnees de la dimension
      }
    } #fin boucle marqueur
    
  } #fin boucle dim
  
  
  #----------------------#
  #- donnees Y du sujet -#
  #----------------------#
  
  ## Y_i = vecteur des obs des marqueurs du sujet collectees avant temps s
  dataY_i <- dataY[which(dataY[,1]==num),] # selection des donnees Y du sujet
  dataY_i <- dataY_i[which(dataY_i$timeT <= s),] # selection des donnees collectees avant le temps s
  Y_i <- dataY_i$measureY
  
  
  
  
  
  
  ###################
  #### donnees D ####
  ###################
  
  ## evt absorbant et en temps continu -> tps a discretiser (fct surBreaks)
  
  dc <- NULL
  D1 <- D[[1]]  # recupere def de evenement
  tD1 <- terms(D1)  # descriptif      
  
  gauche1 <- all.vars(tD1[[2]]) # recupere noms des variables T0:temps initial, T:temps d event et d:indicateur event
  if(length(gauche1)!=3) stop("Please specified entry time, time and indicator in that order in survBreak") 
  
  entry <- gauche1[1] #T0:temps initial
  var.time <- c(var.time,gauche1[2]) #T:temps d event
  dc <- gauche1[3] #d:indicateur event
  
  varexplD <- all.vars(formula(paste("~",D1[3])))  # noms VarExplSurv
  nvarexplD <- length(varexplD) # nbre de VarExplSurv
  
  
  dataD1 <- NULL

  
  dataD1 <- data[,c(subject,entry,var.time[K+1],dc,setdiff(unique(c(varexplD,nomsX)),c(var.time,entry)))]
  dataD1 <- na.omit(dataD1) # suppression des valeurs manquantes
  dataD1 <- dataD1[order(dataD1[,1],dataD1[,3]),] # ordonne par sujet puis par tps d event
  
  ## garder la derniere ligne de chaque sujet de dataD1
  ntmp <- table(dataD1[,1])
  dataD1 <- dataD1[cumsum(ntmp),,drop=FALSE]# derniere   #cumsum indique pour chq sujet a quelle ligne dans tte le base correspond sa derniere visite
  
  ## verif si le sujet a deja eu l event avant s
  dataD1_i <- dataD1[which(dataD1[,1]==num),]
  # alerte si sujet a deja eu l'evenement 1 avant s
  if(dataD1_i[,4]==1 & dataD1_i[,3]<=s)  stop("subject already had the event at time  ",dataD1_i[,3],"  (before s = ",s, ")")
  
  
  #------------------#
  #- discretisation -#
  #------------------#

  breaks <- sort(breaks) # ordonner par ordre croissant
  tbr <- (breaks[-length(breaks)]+breaks[-1])/2 # calcul des milieux des intervalles
  nbr <- length(breaks)
  
  t1max <- dataD1[,c(1,3,4,2)] # sujet, tps event, indicateur event, t0  _  1 par patient
  colnames(t1max) <- c(colnames(dataY)[1],"t1max","d1","t1min")
  #ajout colonne intervalle qui indique dans quel intervalle tombe le temps d event max observe de chq patient
  t1max <- data.frame(t1max,intervalle=as.numeric(cut(x=t1max$t1max,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))
  #ajout colonne intervT0 qui indique dans quel intervalle tombe le temps initial de chq patient
  t1max <- data.frame(t1max,intervT0=as.numeric(cut(x=t1max$t1min,breaks=breaks,labels=1:(nbr-1),include.lowest=TRUE,right=FALSE)))
  ## si la patient n a pas l event a la derniere mesure, on s arr ete a la fin de l intervalle precedent   # UTILE ???
  t1max$intervalle[which(t1max$d1==0)] <- t1max$intervalle[which(t1max$d1==0)]-1
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~#
  #~ tps d event du sujet ~#
  #~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # t1max du sujet
  t1max_i <- t1max[which(t1max[,1]==num),]
  
  mip <- tbr[which(tbr < st)] # milieux des intervalles avant s+t
  mip <- mip[which(mip >= tbr[t1max_i$intervT0])] # milieux des intervalles entre tps initial (tbr[t1max_i$intervT0]) du sujet et s+t
  Bs <- length(mip[which(mip < s)])-1 # identifiant de l intervalle juste avant s _ sachant 1er intervalle identifie 0
  Bst <- length(mip[which(mip < st)])-1 # identifiant de l intervalle juste avant s+t
  
  
  #----------------------------------------------------------------------#
  #- reconstruction des donnees dc + construction du dataframe du sujet -#
  #----------------------------------------------------------------------#

  dataD1bis <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))
  dataD1bis_i <- as.data.frame(matrix(NA,0,ncol(dataD1),dimnames=list(NULL,colnames(dataD1))))

  for(i in 1:nrow(t1max)) # pr chq patient
  {
    t1 <- tbr[t1max[i,6]:t1max[i,5]]  # milieux des intervalles entre tps initial et tps maxi d event obs
    t01 <- rep(-Inf,length(t1))
    
    di1 <- data.frame(rep(t1max[i,1],length(t1)),t01,t1,rep(0,length(t1))) # creation d un dataframe avec chq milieu d intervalles du patient
    if(t1max[i,3]==1) di1[length(t1),4] <- 1  #si le patient a eu l event, alors on associe 1 au dernier intervalle, 0 aux autres
    if(ncol(dataD1)>4) #ie. si il y a des varexpl
    {
      for(k in 5:ncol(dataD1))
      {
        di1 <- data.frame(di1,rep(dataD1[which(dataD1[,1]==t1max[i,1])[1],k],length(t1))) # di1 = di1 + varExpl
      }
    }
    
    colnames(di1) <- colnames(dataD1)
    dataD1bis <- rbind(dataD1bis,di1)
    
    
    if(t1max[i,1]==num){
      t1_i <- mip
      t01_i <- rep(-Inf,length(t1_i))
      
      di1_i <- data.frame(rep(t1max[i,1],length(t1_i)),t01_i,t1_i,rep(0,length(t1_i)))
      if(ncol(dataD1)>4) #ie. si il y a des varexpl
      {
        for(k in 5:ncol(dataD1))
        {
          di1_i <- data.frame(di1_i,rep(dataD1[which(dataD1[,1]==t1max[i,1])[1],k],length(t1_i))) # di1 = di1 + varExpl
        }
      }
      
      colnames(di1_i) <- colnames(dataD1)
      dataD1bis_i <- rbind(dataD1bis_i,di1_i)
      
    } 
  }
  
  dataD1 <- dataD1bis
  dataD1_i <- dataD1bis_i
  
  D <- t1max[,c(1,3)] #sujet, indicateur event
  
  
  ## union des temps de mesure evts
  dataD <- dataD1
  dataD$tdemdc <- dataD[,var.time[K+1]]
  
  dataD_i <- dataD1_i
  dataD_i$tdemdc <- dataD_i[,var.time[K+1]]
  
  
  ## nb de sujets (select au moins 1 mesure dc, non prevalents diag)       
  ni0 <- table(dataD[,1]) #nbre de mesures d event par patient
  ns <- length(ni0) #nbre de patients
  
  
  ## enlever les sujets exclus de dataY
  ntmp <- table(dataY[,1])
  ienlev <- setdiff(names(ntmp),names(ni0))
  if(length(ienlev))
  {
    dataY <- dataY[-which(dataY[,1] %in% ienlev),]
  }
  if(!(num %in% dataY[,1])) stop("no survival information for subject choosen")
  
  
  
  #--------------#
  #- VarExpSurv -#
  #--------------#
  
  nvarD1 <- 0
  ndept1 <- 0
  varD1 <- NULL
  vardept1 <- NULL
  modmat1 <- NULL
  Xseuil_i <- NULL
  
  
  if(D1[[3]]!=1) # si il y a des VarExplSurv
  {
    form1 <- formula(paste("~",subject,"+",D1[3])) # ID VarExpSurv
    modmat1 <- model.matrix(form1,data=dataD[which(!is.na(dataD[,dc])),])[,-1,drop=FALSE] #construit automatiquement les variables des interactions et des contrastes
    modmat1_i <- model.matrix(form1,data=dataD_i[which(!is.na(dataD_i[,dc])),])[,-1,drop=FALSE]
    vars1 <- labels(terms(form1))
    placea <- NULL
    placeb <- NULL
    for(l in 2:ncol(modmat1)) #pr chq colonne de modmat (sauf la colonne intercept)
    {
      tmp <- unique(na.omit(modmat1[,c(1,l)])) #valeurs prises par la variable
      if(nrow(tmp)==ns) # si la variable n est pas dependante du temps
      {
        varD1 <- c(varD1,vars1[l])
        nvarD1 <- nvarD1+1
        placea <- c(placea,l-1)
      }
      else # si la variable est dependante du tps
      {
        vardept1 <- c(vardept1,vars1[l])
        ndept1 <- ndept1+1
        placeb <- c(placeb,l-1)
      }
    }
    
    modmat1 <- modmat1[,c(1,placea+1,placeb+1)] #ordre = ID, VarExplSurv non-dep, VarExplSurv dep du temps
    modmat1_i <- modmat1_i[,c(1,placea+1,placeb+1)]
  }
  
  
  if(nvarD1>0) #si il y a des variables non-dep du tps
  {
    D <- cbind(D,unique(modmat1[,c(1:(nvarD1+1))])[,-1]) # D = ID, ind event, VarExpSurv ind du tps
  }
  
  if(ndept1>0) #si il y a des variables dep du tps
  {
    Xseuil_i <- modmat1_i[,1+nvarD1+1:ndept1,drop=FALSE]
  }
  
  
  ##########################
  # Infos des sous-modeles #
  ##########################
  
  link <- NULL
  nodes <- NULL
  nmes_i <- NULL
  idx <- NULL
  ntr <- NULL
  
  nv <- rep(0,K)
  idiag <- rep(0,K)
  q <- rep(0,K)
  ncor <- rep(0,K)
  nalea <- rep(0,K)
  
  X_i <- vector("list",K)
  Xd_i <- vector("list",K)
  
  
  for(k in 1:K) # pr chq dimension
  {   
    mod <- get(paste("mod",k,sep=""))  # recupere le mod/descriptif de la dimension k

    ## formule k       
    formk <- paste(paste(mod$call$fixed[3],collapse="+"),paste(mod$call$random[2],collapse="+"),"1",sep="+") # VarExplEF1 + VarExplEF2 + ... + VarExplEA1 + VarExplEA2 + ... + 1
    if(!is.null(mod$call$cor)) # si la dimension comporte un processus d autocorrelation
    {
      formk <- paste(formk,as.character(mod$call$cor)[2],sep="+") #VarExplEF1 + VarExplEF2 + ... + VarExplEA1 + VarExplEA2 + ... + 1 + BM
    }
    
    ## les donnees X et Y
    yk_i <- NULL
    xk_i <- NULL
    outcome_i <- NULL
    
    for(m in 1:ny[k])  # pr chaque marqueur de la dimension
    {
      datam_i <- dataY_i[which((dataY_i$processK==k) & (dataY_i$outcomeM==m)),]
      colnames(datam_i)[which(colnames(datam_i)=="timeT")] <- var.time[k]
      colnames(datam_i)[which(colnames(datam_i)=="measureY")] <- mod$Ynames[m]
      yk_i <- c(yk_i,datam_i[,mod$Ynames[m]])  # obs des marqueurs
      xk_i <- rbind(xk_i,model.matrix(formula(paste("~",formk)),data=datam_i))  # obs des variables des EF, EA, BM
      outcome_i <- c(outcome_i,rep(m,nrow(datam_i)))
      
      
      #-----------------------#
      #- construction nmes_i -#
      #-----------------------#
      
      nmesm_i <- table(datam_i[,1])  #nbre de mesures des marqueurs
      nmesm_i <- data.frame(id=names(nmesm_i),nm=as.vector(nmesm_i))  # ces infos dans un data.frame
      if((k==1) & (m==1)) # si 1er marqueur de la 1ere dimension -> initialisation
      {
        nmes_i <- nmesm_i
        colnames(nmes_i) <- c("id","k1m1")
      }
      else # si pas 1er marqueur de la 1ere dimension
      {
        old_i <- colnames(nmes_i)
        nmes_i <- merge(nmes_i,nmesm_i,by="id",all=TRUE,sort=FALSE)
        colnames(nmes_i) <- c(old_i,paste("k",k,"m",m,sep=""))
      }
      
      
      
      ## ntr # nbre de prms de la fonction de lien
      if(mod$linktype[m]==0) ntr <- c(ntr,2) # linear
      if(mod$linktype[m]==1) ntr <- c(ntr,4) # beta
      if(mod$linktype[m]==2) ntr <- c(ntr,mod$nbnodes[m]+2) # splines
    }
    
    
    
    #--------------------#
    #- construction X_i -#
    #--------------------#
    
    yx_i <- cbind(outcome_i,yk_i,xk_i)
    yxord_i <- yx_i[order(yx_i[,1]),] # ordre des lignes = sujet, outcome  
    X_i[[k]] <- yxord_i[,-c(1,2),drop=FALSE] # valeurs des VarExp du sujet avant tps s
    
    
    #---------------------#
    #- construction Xd_i -#
    #---------------------#
    
    ty <- paste("\\b",var.time[k],"\\b",sep="") # "\\bTPS\\b" ! sans espace
    form0 <- gsub(ty,"tdemdc",formk) # remplace ty par tdemdc dans formk
    Xd_i[[k]] <- model.matrix(formula(paste("~",form0)),data=dataD_i) # valeurs des VarExp du sujet aux temps d event
    
    
    #---------------------#
    #- infos du modele k -#
    #---------------------#
    
    nv[k] <- ncol(xk_i) #nbre de colonnes dans chacun des X pr chq dim
    idiag[k] <- mod$idiag #vecteur indiquant si les EA sont correles pr chq dim
    ncor[k] <- mod$N[7] # 1 si BM, 0 sinon
    q[k] <- sum(mod$idea) # nbre EA
    nalea[k] <- mod$N[6] # nbre EA marker-spe
    idx <- cbind(idx,rbind(mod$idg,mod$idcontr,mod$idea,mod$idcor)) # EA, contrastes, EA, BM
    link <- c(link,mod$linktype)
    nodes <- c(nodes,as.vector(mod$linknodes))
    
  } # fin boucle sur k
  
  
  
  
  
  
  ########################
  ## ARGUMENTS DE PROBA ##
  ########################
  
  ## identifiant de l intervalle juste avant s  _ sachant que le 1er intervalle est identifie 0
  Bs
  
  ## identifiant de l intervalle juste avant s+t  _ sachant que le 1er intervalle est identifie 0
  Bst
  
  ## vecteur des parametres du modele conjoint
  b <- model$bopt
  if(!missing(B)){
    if(length(B)!=length(b))
      stop("argument B must be of length ",length(b)," but it's of length ",length(B))
    b <- B
  }
  
  ## vecteur des obs des marqueurs collectees avant s
  Y_i
  
  ## liste des matrices des VarExplLong de chq dim aux tps de visites avant s
  X_i
  
  ## liste des matrices des VarExplLong de chq dim aux tps d event
  Xd_i
  
  ## valeurs des VarExpSurv independantes du temps pr sujet i
  D_i <- NULL # si pas de VarExpSurv independante du temps -> pose pblm pr proba ??
  D_i <- D[which(D[,1]==num),-2] # on ne garde que les VarExplSurv indep. du tps, on supprime son statut d event
  D_i <- unlist(D_i)
  
  ## valeurs des VarExpSurv dependantes du temps pr sujet i
  if(is.null(Xseuil_i)) Xseuil_i <- matrix(0,0,0)
  Xseuil_i <- as.matrix(Xseuil_i)
  
  ## nbre d observations pr chq marqueur pr sujet i
  nmes_i <- as.vector(nmes_i[,-1]) # 1ere colonne = id a enlever # nbre obs par marqueur pr individu
  if(any(is.na(nmes_i))) nmes_i[which(is.na(nmes_i))] <- 0 # ie. pas d obs pr 1 marqueur
  nmes_i <- unlist(nmes_i)
  
  ## nbre de colonnes de chq matrice de la liste X
  nv
  
  ## nature des VarExplLong (EF,contrast,EA,BM)
  idx
  
  ## indicateur de la correlation des EAs pr chq dim
  idiag
  
  ## indicateur de presence et nature du processus d autocorrelation de chq dim
  ncor
  if(any(ncor==2)) stop("Autoregressive correlations are not allowed") # pas AR
  
  ## nbre de marqueur(s) ds chq dim
  ny
  
  ## nbre d EAs marker-specific ds chq dim
  nalea
  
  ## nbre de prms des fcts de liens des marqueurs
  ntr 
  if(any(ntr==4)) stop("Use only linear or splines links") # pas beta
  
  ## types des fcts de lien des marqueurs
  link
  
  ## noeuds ou min/max des fcts de lien des marqueurs
  nodes
  
  ## nbre de correlations entre les EAs des dims 
  if(K==1 | RE=="block-diag") nRE <- 0  # nbre de prms de corr entre les dims # 1ere dim ou EA indep
  if(K>1 & RE!="block-diag") nRE <- sum(sapply(1:(length(q)-1), function(k,x) sum(x[k]*x[(k+1):length(x)]),x=q)) # pas 1ere dim & EA indep
  
  ## nbre de correlations entre les BMs des dims 
  BM <- "diag"  #on ne traite pas le cas ou les BM sont correles ??
  if((K==1) | (BM=="diag")) nBM <- 0  # 1ere dim ou BM indep
  if((K>1) & (BM=="full"))    # pas 1ere dim & EA corr
  {
    ncortot <- sum(ncor==1) # nbre de dim dt les BM sont correles
    nBM <- ncortot*(ncortot-1)/2
  }
  
  
  #########################################
  ## AFFICHAGE DES ARGUMENTS DE LA PROBA ##
  #########################################
  
  # cat("Bs = ",Bs,"\n")
  # cat("Bst = ",Bst,"\n")
  # cat("b = ",b,"\n")
  # cat("Y_i = ",Y_i,"\n")
  # cat("X_i = \n")
  # print(X_i)
  # cat("Xd_i = \n")
  # print(Xd_i)
  # cat("D_i = ",D_i,"\n")
  # cat("Xseuil_i = \n")
  # print(Xseuil_i)
  # cat("nmes_i = ",nmes_i,"\n")
  # cat("nv = ",nv,"\n")
  # cat("idx = \n")
  # print(idx)
  # cat("idiag = ",idiag,"\n")
  # cat("ncor = ",ncor,"\n")
  # cat("ny = ",ny,"\n")
  # cat("nalea = ",nalea,"\n")
  # cat("ntr = ",ntr,"\n")
  # cat("link = ",link,"\n")
  # cat("nodes = ",nodes,"\n")
  # cat("nRE = ",nRE,"\n")
  # cat("nBM = ",nBM, "\n")

  cat("num = ",num,"\n")
  
  #####################
  ## CALCUL DE PROBA ##
  #####################
  
  res <- RprobaEC(Bs0=Bs, Bst0=Bst, b0=b, Y0=Y_i, X0=X_i, Xd0=Xd_i, D0=D_i, Xseuil0=Xseuil_i, nmes0=nmes_i, nv0=nv, idx0=idx, idiag0=idiag, ncor0=ncor, ny0=ny, nalea0=nalea, ntr0=ntr, link0=link, nodes0=nodes, nRE0=nRE, nBM0=nBM)
  
  return(res)
}





############################################################################################################################
############################################################################################################################
############################################################################################################################