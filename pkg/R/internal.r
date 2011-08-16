`qfdrk` <-
function(z,th=c(.1,.01),sgn=0,nulltype=0,bre=seq(-6.1,6.1,.2),J=5,sw=1,plotlocfdr=FALSE){
    ## fdrk of matrix z   1/1/10
    ## sgn=0,1,2  does both,left,right side;  ##th<=0 avoids plotting
    
    N=nrow(z); n=ncol(z)
    
    vl=locfdr(z,bre,nulltype=nulltype,plot=plotlocfdr)[-1]
    ma=vl$mat
    x=ma[,1]; fdr=ma[,2]
    
    if(sgn==1)fdr[x>0]=1
    if(sgn==2)fdr[x<0]=1
    
    Tdr=(1-matrix(.approx(x,fdr,z),N))
    if(sw==3)return(.list(x,fdr,Tdr))
    
    K=.colsum(Tdr,1) 
    Kbar=mean(K)
    S1=K/Kbar-1
    rat=1+(Tdr)*S1
    fdrk=(1-Tdr)/rat
    
    if(J>1){ok=qfd2(Tdr,J,sw=1)
    K=ok$K; fdrk=ok$fdrk}
    
    if(th[1]>0)for(k in 1:length(th)){qfd1(fdrk,th[k],add=min(k,2)-1,col=k)}
    
    if(sw==2) return(.list(fdrk,K,x,fdr,Tdr,N,n))
    if(sw==1)return(.list(K,fdrk))
    
    fdrk
}

`qfd1` <-
function(fdrk,th=.1,add=0,col=1,pch="-",iso=1,ii,sw=0,...){
#plotting program for qfdrk output fdrk 1/2/10
# can enter fdrk=u;  ii takes those rows of fdrk

    if(ncol(fdrk)<=3){u=fdrk}
    else{
        if(!missing(ii))fdrk=fdrk[ii,]
        
        # U=.zer(fdrk); 
        U=matrix(nrow=nrow(fdrk),ncol=ncol(fdrk),data=0); 
        U[fdrk<=th]=1
        if(iso>0){
        s=diff(U)
        ss=diff(s)
        sss=ifelse(ss==-2,1,0)
        s=1-.rb(0,sss,0)
        U=s*U
        if(sw==3)return(s)}
        u=.index(U,1)
        if(!missing(ii))u[,1]=u[,1]+min(ii)-1
    }
    if(add==1)points(u,pch=pch,col=col)
        else  plot(u,pch=pch,col=col,xlab="marker",ylab="subject", ...)
    if(sw==1)return(u)
}

blockbootstrap<-function(X,L){
    N=nrow(X)
    n=ncol(X)
    
    nbl = floor(N/L)
    Lrem = N-L*nbl
    
    stmat=matrix(ncol=n,nrow=nbl,data=sample(N-L+1, nbl*n,replace=TRUE))
    
    Lmat = matrix(data=rep(c(0:(L-1)),nbl),byrow=TRUE,nrow=nbl)
    Xstar = matrix(data=0,nrow=N,ncol=n)
    
    strem = sample(N-Lrem+1,n,replace=TRUE)
    
    for(i in 1:n){
        st = stmat[,i]
        st = matrix(nrow=nbl,data=rep(st,L),byrow=FALSE)
        inds= st+Lmat
        Xi = X[,i]   
        Xistar = Xi[t(inds)]
        Xistar = c(Xistar,Xi[strem[i]+c(0:(Lrem-1))])
        Xstar[,i] = Xistar
    }
    Xstar
}

compute.var<-function (y) 
{
    N = dim(y)[2]
    T = dim(y)[1]
    y.diff <- y[1:(T - 1), ] - y[2:T, ]
    y.var <- apply(y.diff^2, 2, mean)/2
    y.var
}

`qfd2` <-
function(Tdr,J=5,simp=1,sw=0){
#iterated fdrk and K 1/16/10
# Tdr from qfdrk(z,sw=2); simp=0 uses simpler formula

    n=ncol(Tdr); 
    N=nrow(Tdr)
    pi1=mean(Tdr); pi0=1-pi1
    Fdr=1-Tdr
    tdrk=Tdr
    #    KK=.zer(N,J)
    KK = matrix(nrow=N, ncol=J, data=0)
    
    for(j in 1:J){
        #   K=.colsum(tdrk,1)
        K = apply(tdrk,1,sum)
        
        KK[,j]=K
        if(simp==1){
            R=(n*pi1/K)/(n*pi0/(n-K))-1
        } else {
            R=mean(K)/K -1
        }
        tdrk=Tdr/(1+Fdr*R)
    }
    
    if(sw==1){
        fdrk=1-tdrk; return(list(KK=KK,fdrk=fdrk))
    }
    KK
}

`.ma` <-
function(x,m,edge=m,sides=2,sw=0){
    #moving averages using "filter"  12/19/09
    #sides as in "filter": 2 for centered, 1 for past
    # edge=0 avoids end corrections
    #if x is matrix, does ma columnwise
    
    
    if(is.vector(x))x=matrix(x)
    nc=ncol(x);nr=nrow(x)
    
    if(edge>0){o1=x[m:1,]; o2=x[nr:(nr-m+1),]; x=.rb(o1,x,o2)}
    
    for(j in 1:nc){
    x[,j]=filter(x[,j],.one(m)/m, sides=sides)
    }
    x[(m+1):(nr+m),]
}

".approx" <-
function(x, y, xout)
{
if(is.matrix(x)){xout=y;y=x[,2];x=x[,1]}
	approx(x, y, xout, rule = 2.,ties=mean)$y
}

`.colsum`<-function(x, tr = 0.)
{
        #colsums of available elements
        if(tr == 0.) x <- t(x)
        nc = ncol(x)
        xav = is.na(x)
        xav = as.vector(xav %*% .one(nc))
        x[is.na(x)] = 0
        xs = as.vector(x %*% .one(nc))
        xs[xav == nc] = NA
        xs
}

".one" <-
function(m, x=1,n, o) {
  ##matrix(x,m,n) or a one's vector
  if(length(c(m)) > 1.) {
    if(is.matrix(m)) {
      n <- dim(m)[2.]
      m <- dim(m)[1.]
      return(matrix(x, m, n))
    }
    else {
      m <- length(c(m))
    }
  }
  if(missing(n))
    return(rep(x, m))
  if(!missing(o))
    return(array(x, c(m, n, o)))
  matrix(x, m, n)
}

`.rb`<-function(..., deparse.level = 1.) {
  rbind(..., deparse.level = deparse.level)
}


`.index`<-function(x, val, rel = 1.)
{
        #rel: 1= ; 2 <=  ;  3 < ;  4 >=  ; 5 > ;[Or enter "<" etc]
        #missing val==> coords min(x),max(x)
        #coords matrix or vector x having relation rel with val
        if(missing(val)) {
                v <- (x == min(x) | x == max(x))
        }
        else {
                if(is.character(rel)) {
                        rell <- get(rel)
                        v <- rell(x, val)
                }
                else {
                        switch(rel,
                                v <- (x == val),
                                v <- (x <= val),
                                v <- (x < val),
                                v <- (x >= val),
                                v <- (x > val))
                }
        }
        if(is.matrix(x)) {
                rowx <- .v(row(x))
                colx <- .v(col(x))
                rc <- .b(rowx, colx)[v,  ]
                xx <- x[rc]
                rc <- .b(rc, xx)
                return(rc)
        }
        rx <- seq(x)
        coord <- rx[v]
        xx <- x[coord]
        .b(coord, xx)
}

`.list`<-function(...) {
 call.list <- as.character(sys.call())
 result <- list(...)
 names(result) <- call.list[-1]
 result
}

`.v`<-function (x = stop("Argument `x' is missing")) 
{
    as.vector(x)
}

`.b`<-function(..., deparse.level=1.) {
  cbind(..., deparse.level = deparse.level)
}


