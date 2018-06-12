print.JointMult <- function(x,...)
{
    cat(" \n")
    cat("Values:", "\n")
    cat(" \n")
    id <- 1:length(x$b)
    indice <- rep(id*(id+1)/2)
    se <-sqrt(x$v[indice])
    wald <- (x$b/se)**2
    z <- abs(qnorm((1 + .95)/2))
    binf <- x$b-1.96*se
    bsup <- x$b+1.96*se
    tmp <- data.frame("coef"=format(round(x$b,3)),"SE coef"=format(round(se,3)),"Wald"=format(wald,4),"P-value"=round(1 - pchisq(wald, 1),5),"binf"=round(binf,3),"bsup"=round(bsup,3))
    print(tmp,row.names=FALSE)

}
