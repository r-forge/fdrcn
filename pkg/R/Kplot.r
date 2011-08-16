Kplot <-
function(KGains, KLosses, maxKGainsBoot, maxKLossesBoot, alpha=0.95){

    ymax = max(max(KGains),max(KLosses))
    plot(KGains,col="blue",type="l",ylim=c(-1.2*ymax,1.2*ymax),main="K Profile")
    lines(-KLosses,col="red")
    qKmaxGains = quantile(maxKGainsBoot,alpha)
    qKmaxLosses = quantile(maxKLossesBoot,alpha)
    abline(qKmaxGains,0,lty=2,col="blue")
    abline(-qKmaxLosses,0,lty=2,col="red")
}
