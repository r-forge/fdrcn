fdrcn <-
function(x,win=11,scaleSamples=TRUE,nulltype=0,J=1,sw=1,B=100,L=NA,alpha=.05,plots=TRUE,th=0.01){
    
    cat("Computing median smoothed matrix ...")
    nsamples=ncol(x)
    nsnps = nrow(x)
    z = matrix(nrow=nsnps, ncol=nsamples)
    for(i in 1:nsamples){
        z[,i] = runmed(x[,i],k=win)
    }
    # robust standardize each column 
    diffrec=rep(0,nsamples)
    for(i in 1:nsamples){
        med = median(z[,i])
        z[,i] = z[,i]-med
        if(scaleSamples){
            diff=quantile(z[,i],.841)-quantile(z[,i], .158)
            diffrec[i] = diff
            z[,i] = 2*z[,i]/diff
        }
    }
    cat(" done.\n")


    cat("Computing FDRs ...")
    res1=qfdrk(z, nulltype=nulltype,sgn=1,th=-1, J=J, sw=3)
    res2=qfdrk(z, nulltype=nulltype,sgn=2,th=-1, J=J, sw=3)
    K1=.colsum(res1$Tdr,1) 
    K2=.colsum(res2$Tdr,2) 
    cat(" done.\n")


    if(B>0){
        cat("Bootstrapping K ")
        if(is.na(L)){
            L = round(nrow(z)/50)
        }
        Tdr1 = res1$Tdr
        Tdr2 = res2$Tdr
        Kstar1 = matrix(nrow=nrow(z),ncol=B) 
        Kstar2 = matrix(nrow=nrow(z),ncol=B)
        for(b in 1:B){
#            cat("--- Bootstrap iteration ",b," ---\n")
#            ptm <- proc.time()
            cat(".")
            X = blockbootstrap(Tdr1,L)
            Kstar1[,b] = .colsum(X,1) 
            X = blockbootstrap(Tdr2,L)
            Kstar2[,b] = .colsum(X,1)        
#            cat(proc.time() - ptm)
        }
        Kmax1 = apply(Kstar1,2,max) 
        Kmax2 = apply(Kstar2,2,max)
        cat(" done.\n")
    } else {
        Kmax1=rep(0,0)
        Kmax2=rep(0,0)
    }
    
    if(plots){
        par(mfrow=c(2,1))
        dotplot(fdr=1-res1$Tdr,th=th, col="red", main=paste("Aberrations at FDR Threshold",th))
        dotplot(fdr=1-res2$Tdr,th=th, add=TRUE, col="blue")
        legend(x="topright",pch=c("-","-"),col=c("red","blue"),legend=c("Losses","Gains"))
        Kplot(K2,K1,Kmax2,Kmax1)
    }
    
    return(list(tdrGains=res2$Tdr,tdrLosses=res1$Tdr, KGains = K2, KLosses = K1, maxKGainsBoot =Kmax2, maxKLossesBoot = Kmax1))

}
