## splines lineaires

linspl <- function(time,nodes)
    {
        nodes <- sort(nodes)
        
        nb <- length(nodes)-1

        if(any(na.omit(time)<nodes[1])) stop("time < tmin")
        if(any(na.omit(time)>nodes[nb+1])) stop("time > tmax")

        
        spl <- matrix(0,length(time),nb)

        for(k in 1:nb)
            {
                ik <- which((time>=nodes[k]) & (time<nodes[k+1]))
                if(length(ik))
                   {
                       spl[ik,k] <- time[ik]-nodes[k]
                   }

                ikplus <- which(time>=nodes[k+1])
                if(length(ikplus))
                    {
                        spl[ikplus,k] <- nodes[k+1]-nodes[k]
                    }
            }

        ina <- which(is.na(time))
        if(length(ina)) spl[ina,] <- NA

        colnames(spl) <- paste("s",1:nb,sep="")
        return(spl)
    }


