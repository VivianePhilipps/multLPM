
plot.fit_ss <- function(x,which="residuals",outcome=1,...)
    {
        if(which=="residuals")
            {
                res <- x[[1]][[outcome]]
                plot(x=res[,5],y=res[,4]-res[,5],...)
            }
        else
            {
                res <- x[[2]][[outcome]]
                matplot(x=res[,1],y=res[,2:5],...)
            }
        
    }

