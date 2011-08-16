dotplot<-function(fdr,th=.1,add=0,col=1,pch="-",iso=1,ii,sw=0,...){
# This essentially rehashes Brad's function qfd1.

    if(ncol(fdr)<=3){u=fdr}
    else{
        if(!missing(ii))fdr=fdr[ii,]

        # U=.zer(fdr);
        U=matrix(nrow=nrow(fdr),ncol=ncol(fdr),data=0);
        U[fdr<=th]=1
        if(iso>0){
            s=diff(U)
            ss=diff(s)
            sss=ifelse(ss==-2,1,0)
            s=1-.rb(0,sss,0)
            U=s*U
            if(sw==3)return(s)
        }
        u=.index(U,1)
        if(!missing(ii))u[,1]=u[,1]+min(ii)-1
    }



    if(add==1)  points(u,pch=pch,col=col)
    else  plot(u,pch=pch,col=col,xlab="marker",ylab="subject",...)
    if(sw==1)return(u)

}

