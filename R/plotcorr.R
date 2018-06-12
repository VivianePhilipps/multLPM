


#### graph des correlations des effets aléatoires, avec les p-values,
#### pour un modele a K processus

plotcorr <- function(m)
    {
        npm <- length(m$b)
        indice <- c((1:npm)*(1:npm+1)/2)
        se <-sqrt(m$v[indice])
        wald <- (m$b/se)**2
        pval <- 1 - pchisq(wald, 1)
        binf <- m$b-1.96*se
        bsup <- m$b+1.96*se

        K <- length(m$nef)
        npmMM <- sum(c(m$nef+m$ncontr+m$nvc+m$ncor+m$ny+m$nalea+m$ntrtot))
        q <- sum(m$nea)
        
        ##p 
        ptot <- rep(NA,length(m$bopt))
        ptot[which(m$fix==0)] <- pval

        pcor <- matrix(NA,nrow(m$corRE),ncol(m$corRE))
        imatB <- 0
        iB <- 0
        iRE <- 0
        for(k in 1:K)
            {
                if(m$nea[k]==0 | m$idiag[k]==1) next
                
                if(m$nvc[k]>0)
                    {
                        ii <- upper.tri(m$corRE[imatB+1:m$nea[k],imatB+1:m$nea[k]])
                        jj <- setdiff(1:m$nvc[k],((1:m$nea[k])*(1:m$nea[k]+1))/2-1)
                        pcor[imatB+1:m$nea[k],imatB+1:m$nea[k]][ii] <- ptot[iB+m$nef[k]+jj]
                    }
                
                if((imatB>0) & (m$nRE>0))
                    {
                        qavt <- sum(m$nea[1:(k-1)])
                        pcor[1:qavt,imatB+1:m$nea[k]] <- ptot[npmMM+iRE+1:(qavt*m$nea[k])]
                        iRE <- iRE+qavt*m$nea[k]
                    }
                
                imatB <- imatB+m$nea[k]
                iB <- iB+m$nef[k]+m$ncontr[k]+m$nvc[k]+m$ncor[k]+m$ny[k]+m$nalea[k]+m$ntrtot[k]
            }
        diag(pcor) <- NA
        pcor <- t(pcor)
        
        pm2 <- t(pcor)
        pm2[which(pm2==0)] <- 0.00001
        pm2[lower.tri(pm2)] <- 0
        pm2[which(m$corRE<0)] <- -pm2[which(m$corRE<0)]
        pm2[which(pm2>0.95)] <- 0.9
        pm2[which(is.na(pm2))] <- 0.9
        diag(pm2) <- 1


#browser()
        
        colorp <- c(rep("white",19),"seagreen1","steelblue1",rep("white",18),"black")
        gris <- paste("grey",c(50,rep(70,18)),sep="")


        corrplot(pm2,method="color",number.digits=3,type="upper",tl.col="white",na.label=" ",col=colorp,is.corr=FALSE,cl.pos="n",addgrid.col=NULL)


        corrplot(m$corRE,method="number",number.digits=3,type="upper",tl.col="white",na.label=" ",col="black",add=TRUE,bg="transparent",addgrid.col=NULL,cl.pos="n")

        try(corrplot(corr=pcor,p.mat=pcor,method="number",number.digits=3,type="lower",tl.pos="n",add=TRUE,col=c(rev(gris),"black","black",gris),diag=FALSE,legend=NULL,cl.length=-1,na.label=" ",sig.level=0.05,insig="n",addgrid.col="grey",cl.pos="n"),silent=TRUE)

        segments(0.5,q+0.5,q+0.5,0.5,col="white",lwd=2,xpd=TRUE)

        text(1.2,q+0.3,"corr",col="white",xpd=TRUE)

        text(0.9,q-0.3,"p value",col="white",xpd=TRUE)


    }
