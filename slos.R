library(NonpModelCheck)
slos <- function(X,y,K=1,D=NULL,lambda=0,gamma=0,intercept=T,wrap=FALSE,beta.basis=NULL,cutoff=1e-4,
                 max.iter=1000,tol=1e-12,a=3.7,domain=NULL,tuning=NULL)
{
    # D     : the design points
    # X     : predictor(s), could be a matrix (fd object) or list of matrix (fd objects)
    #         if matrix, then it is organized in the way that each column is an independent observation
    # y     : a column vector, the response
    
    ################ sanity check #####################
    
    if(!require('fda'))
        stop("The slos function depends on the fda package which is missing")
    if(is.null(D))
    {
        if(is.matrix(X))
            stop("The design should be provided when X is matrix")
    }
    
    if(is.matrix(X) || is.fd(X)) # there is only one predictor, convert it into a fd object
    {
        X <- list(x=X)
    }
    
    if(!is.list(X))
        stop("X must be a matrix, fd object or a list of matrices/fd objects")
    
    ################### prepare data in the following way #####################
    # K     :   the number of predictors
    # xfd   :   list of K matrices / fd objects
    # design:   the design points for each objects
    
    xfd <- slos.preprocess.x(X,D)
    
    ######################### prepare basis for coefficient functions #############
    
    if(is.null(beta.basis))
    {
        beta.basis <- list()
        for(k in 1:K)
        {
            x <- xfd[[k]]
            rng <- getbasisrange(x$basis)
            beta.basis[[k]] <- create.bspline.basis(rangeval=rng,nbasis=73)
        }
    }
    
    Jmat <- list()
    for(k in 1:K)   Jmat[[k]] <- inprod(xfd[[k]]$basis,beta.basis[[k]])
    
    if(length(gamma)>1 || length(lambda)>1)
    {
        tuned <- slos.tune(xfd,y,lambda,gamma,intercept,wrap,beta.basis,cutoff,
                           max.iter,tol,a,Jmat,tuning)
        result <- slos.alg(xfd,y,tuned$lambda,tuned$gamma,intercept,wrap,beta.basis,cutoff,
                           max.iter,tol,a,Jmat)
        result$lambda <- tuned$lambda
        result$gamma <- tuned$gamma
        result$score <- tuned$score
        result
    }
    else
    {
        slos.alg(xfd,y,lambda,gamma,intercept,wrap,beta.basis,cutoff,
                 max.iter,tol,a,Jmat)
    }
    
}

slos.tune <- function(xfd,y,lambda,gamma,intercept,wrap,beta.basis,cutoff,
                      max.iter,tol,a,Jmat,tuning,Vs=NULL,Ws=NULL)
{
    K <- length(xfd)
    # compute weight and roughness penalty matrices if they are NULL
    if(is.null(Vs))
    {
        Vs <- list()
        for(k in 1:K)
        {
            Vs[[k]] <- eval.penalty(beta.basis[[k]],int2Lfd(2))
        }
    }
    
    if(is.null(Ws))
    {
        Ws <- list()
        for(k in 1:K)
        {
            Ws[[k]] <- slos.compute.weights(beta.basis[[k]])
        }
    }
    
    dm <- dim(xfd[[1]]$coef)
    n <- dm[2]
    
    if(tuning == 'VALIDATION')
    {
        val.xfd <- slos.preprocess.x(attr(tuning,'X'))
        val.y <- attr(tuning,'y')
    }
    
    err.score <- matrix(0,length(lambda),length(gamma))
    for(i in 1:length(lambda))
    {
        for(j in 1:length(gamma))
        {
            if(tuning != 'CV' && tuning != 'VALIDATION')
            {
                fit <- slos.alg(xfd,y,lambda[i],gamma[j],intercept,wrap,beta.basis,cutoff,
                                max.iter,tol,a,Jmat,Vs,Ws,use.index=NULL)
                resid <- y-fit$fitted.values
                rss <- sum(resid^2)
                df <- sum(diag(fit$projMat))
                
                if(tuning == 'RIC')
                {
                    sig2 <- rss/(n-df)
                    err <- (n-df)*log(sig2)+df*(log(n)-1)+4/(n-df-2)
                }
                else if(tuning == 'AIC')
                {
                    err = n*log(rss/n)+2*df
                }
                else if(tuning == 'BIC')
                {
                    err = n*log(rss/n) + log(n)*df
                }
                else stop("Unrecognized tuning criterion")
            }
            else if(tuning == 'VALIDATION')
            {
                fit <- slos.alg(xfd,y,lambda[i],gamma[j],intercept,wrap,beta.basis,cutoff,
                                max.iter,tol,a,Jmat,Vs,Ws,use.index=NULL)
                yhat <- predict(fit,val.xfd)
                resid <- val.y - yhat
                err <- sum(resid^2)
            }
            else if(tuning == 'CV')
            {
                cv.type <- attr(tuning,'type')
                if (is.null(cv.type)) cv.type <- "k-fold"
                if (cv.type == 'leave-one-out') Kfold <- n
                else Kfold <- attr(tuning,'K')
                if (is.null(K)) Kfold <- 10
                idx <- sample(1:n,n) # shuffle data
                
                s <- ceiling(n/Kfold)
                
                err.cv <- matrix(0,1,Kfold)
                for(k in 1:Kfold)
                {
                    test.index <- idx[(s*(k-1)+1):min(s*k,n)]
                    train.index <- idx[-(s*(k-1)+1):min(s*k,n)]
                    fit <- slos.alg(xfd,y,lambda[i],gamma[j],intercept,wrap,beta.basis,cutoff,
                                    max.iter,tol,a,Jmat,Vs,Ws,use.index=train.index)
                    yhat <- predict(fit,xfd,test.index)
                    resid <- y[test.index] - yhat
                    err.cv[k] <- sum(resid^2)
                }
                
                err <- sum(err.cv)
            }
            
            err.score[i,j] <- err
        }
    }
    
    h <- which(err.score==min(err.score),arr.ind=T)
    list(score=err.score,lambda=lambda[h[1,1]],gamma=gamma[h[1,2]],min.score=min(err.score))
}

slos.mse <- function(y,fitted.values)
{
    mean((y-fitted.values)^2)
}

slos.preprocess.x <- function(X,D)
{
    if(is.matrix(X) || is.fd(X))
    {
        Xtmp <- list()
        Xtmp[[1]] <- X
        X <- Xtmp
        K <- 1
    }
    else K <- length(X)
    
    xfd <- list()
    for(k in 1:K)
    {
        x <- X[[k]]
        if(!is.fd(x))
        {
            if(is.matrix(D)||is.vector(D))
            {
                x <- Data2fd(argvals=D,y=x)
            }
            else
            {
                x <- Data2fd(argvals=D[[k]],y=x)
            }
        }
        
        if(!is.fd(x))
            stop("X must be a matrix, fd object or a list of matrices/fd objects")
        xfd[[k]] <- x
    }
    
    xfd
    
}

slos.compute.weights <- function(basis)
{
    L <- basis$nbasis
    rng <- getbasisrange(basis)
    breaks <- c(rng[1],basis$params,rng[2])
    M <- length(breaks) - 1
    norder <- L-M+1
    W <- array(0,dim=c(norder,norder,M))
    for (j in 1:M)
    {
        temp <- inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
        W[,,j] <- temp[j:j+norder-1,j:j+norder-1]
    }
    W
}

slos.alg <- function(xfd,y,lambda,gamma,intercept,wrap,beta.basis,cutoff,
                     max.iter,tol,a,Jmat,Vs=NULL,Ws=NULL,use.index=NULL)
{
    # get some constants
    K <- length(xfd)
    Ls <- rep(0,K)
    Ms <- rep(0,K)
    ds <- rep(0,K)
    L2NNer <- rep(0,K)
    
    for (k in 1:K) Ls[k] <- beta.basis[[k]]$nbasis
    for (k in 1:K) Ms[k] <- length(beta.basis[[k]]$params) + 1
    for (k in 1:K) ds[k] <- Ls[k] - Ms[k]
    for (k in 1:K) L2NNer[k] <- sqrt(Ms[k]/diff(getbasisrange(xfd[[k]]$basis)))
    P <- c(0,cumsum(Ls))
    
    # compute weight and roughness penalty matrices if they are NULL
    if(is.null(Vs))
    {
        Vs <- list()
        for(k in 1:K)
        {
            Vs[[k]] <- eval.penalty(beta.basis[[k]],int2Lfd(2))
        }
    }
    
    if(is.null(Ws))
    {
        Ws <- list()
        for(k in 1:K)
        {
            Ws[[k]] <- slos.compute.weights(beta.basis[[k]])
        }
    }
    
    # compute design matrix
    U <- NULL
    VV <- NULL
    
    if(!is.null(use.index)) y <- y[use.index]
    
    n <- length(y)
    
    for (k in 1:K)
    {
        xcoef <- xfd[[k]]$coefs # M by n
        if(!is.null(use.index))
        {
            U <- cbind(U,t(xcoef[,use.index])%*%Jmat[[k]])
            VV <- bdiag(VV,n*gamma*Vs[[k]])
        }
        else
        {
            U <- cbind(U,t(xcoef)%*%Jmat[[k]])
            if(is.null(VV))  VV <- n*gamma*Vs[[k]]
            else VV <- bdiag(VV,n*gamma*Vs[[k]])
        }
    }
    
    if (wrap)
    {
        oldU <- U
        oldV <- VV
        dimU <- dim(U)
        dimV <- dim(VV)
        U <- cbind(U[,1]+U[,dimU[2]],U[,2:(dimU[2]-1)])
        VV[,1] <- VV[,1] + VV[,dimV[2]]
        VV[1,] <- VV[1,] + VV[dimV[1],]
        VV <- VV[-dimV[1],-dimV[2]]
    }
    
    if(intercept)
    {
        U <- cbind(matrix(1,n,1),U)
        VV <- bdiag(0,VV)
    }
    
    # compute the initial estimate
    bHat <- solve(t(U)%*%U+VV,t(U)%*%y)
    
    if(intercept)
    {
        muhat <- bHat[1]
        bHat <- bHat[-1]
    }
    
    # now, potential intercept has been removed from bHat
    # if periodic, patch the bHat in the end
    if(wrap)
    {
        bHat <- c(bHat,bHat[1])
    }
    
    bTilde <- bHat
    
    # if sparse, perform SLoS, otherwise, bTilde is the solution to \beta
    if(lambda > 0)
    {
        changeThres <- tol
        bTilde <- slos.lqa(U,y,bHat,intercept,wrap,Ws,Ms,Ls,L2NNer,ds,
                           lambda,a,VV,changeThres,max.iter)
        bZero <- (abs(bTilde) < cutoff)
        bTilde[bZero] <- 0
        
        if(intercept)
        {
            bTilde <- rbind(muhat,bTilde)
            bNonZero <- c(TRUE,!bZero)
        }
        else bNonZero <- !bZero
        
        if(wrap)
        {
            bNonZeroPrime <- bNonZero[-length(bNonZero)]
            U1 <- U[,bNonZeroPrime]
            V1 <- VV[bNonZeroPrime,bNonZeroPrime]
            bb <- solve(t(U1)%*%U1+V1,t(U1)%*%y)
            bTilde <- matrix(0,dim(U)[2],1)
            bTilde[bNonZeroPrime,1] <- matrix(bb,length(bb),1)
            if(intercept)
            {
                bTilde <- c(bTilde,bTilde[2])
            }
            else
            {
                bTilde <- c(bTilde,bbTilde[1])
            }
        }
        else
        {
            U1 <- U[,bNonZero]
            V1 <- VV[bNonZero,bNonZero]
            bb <- solve(t(U1)%*%U1+V1,t(U1)%*%y)
            if(intercept)
            {
                bTilde <- matrix(0,1+sum(Ls),1)
                bTilde[bNonZero,1] <- matrix(bb,length(bb),1)
            }
            else
            {
                bTilde <- matrix(0,sum(Ls),1)
                bTilde[bNonZero,1] <- matrix(bb,length(bb),1)
            }
        }
    }
    
    bNonZero <- as.vector((bTilde != 0))
    if(wrap)
    {
        U1 <- U[,bNonZero[-length(bNonZero)]]
        V1 <- VV[bNonZero[-length(bNonZero)],bNonZero[-length(bNonZero)]]
    }
    else
    {
        U1 <- U[,bNonZero]
        V1 <- VV[bNonZero,bNonZero]
    }
    
    projMat <- U1 %*% solve(t(U1)%*%U1+V1,t(U1))
    result <- list(beta=NULL,projMat=projMat,intercept=0,fitted.values=projMat%*%y)
    if(intercept)
    {
        result$intercept <- bTilde[1]
        bTilde <- bTilde[-1]
    }
    betaobj <- list()
    for(k in 1:K)
    {
        bfd <- fd(coef=bTilde[(P[k]+1):(P[k+1])],basisobj=beta.basis[[k]])
        betaobj[[k]] <- bfd
    }
    result$beta <- betaobj
    class(result) <- 'slos'
    result
    
}

predict.slos <- function(slosobj,newX,use.index=NULL,domain=NULL)
{
    xfd <- slos.preprocess.x(newX,domain)
    K <- length(xfd)
    y <- 0
    for(k in 1:K)
    {
        x <- xfd[[k]]
        b <- slosobj$beta[[k]]
        G <- inprod(x$basis,b$basis)
        if(!is.null(use.index))
            y <- y + t(x$coef[,use.index]) %*% G %*% b$coef
        else
            y <- y + t(x$coef) %*% G %*% b$coef
    }
    y <- y + slosobj$intercept
    y
    
}

slos.vec.norm <- function(v)
{
    sqrt(sum(v^2))
}

slos.lqa <- function(U,y,bHat,intercept,wrap,Ws,Ms,Ls,L2NNer,ds,lambda,a,V,changeThres,max.iter)
{
    Mmax <- max(Ms)
    Lmax <- max(Ls)
    K <- length(Ms)
    P <- c(0,cumsum(Ls))
    
    n <- length(y)
    
    betaNormj <- matrix(0,Mmax,K)
    bZeroMat <- matrix(FALSE,1,sum(Ls))
    
    betaNorm <- rep(Inf,K)
    it <- 1
    
    while (it <= max.iter)
    {
        betaNormOld <- betaNorm
        for (k in 1:K)
        {
            betaNorm[k] <- slos.vec.norm(bHat[(P[k]+1):P[k+1]])
        }
        change <- max((betaNormOld-betaNorm)^2)
        if(change < changeThres) break
        
        lqaW <- NULL
        for (k in 1:K)
        {
            W <- Ws[[k]]
            lqaWk <- matrix(0,Ls[k],Ls[k])
            for(j in 1:Ms[k])
            {
                index <- c(j:(j+ds[k]))
                bkj <- bHat[P[k]+c(j:(j+ds[k]))]
                betaNormj[j,k] <- sqrt(t(bkj) %*% W[,,j] %*% bkj)
                cjk <- Dpfunc(betaNormj[j,k]*L2NNer[k],lambda,a)
                if(cjk != 0)
                {
                    if(betaNormj[j,k] < changeThres) bZeroMat[P[k]+index] <- TRUE
                    else lqaWk[index,index] <- lqaWk[index,index] + cjk*(L2NNer[k]/betaNormj[j,k])*W[,,j]
                }
            }
            
            if(is.null(lqaW)) lqaW <- lqaWk
            else lqaW <- bdiag(lqaW,lqaWk)
        }
        
        bZeroVec <- bZeroMat
        bNonZeroVec <- !bZeroVec
        
        if(wrap)
        {
            if(bNonZeroVec[length(bNonZeroVec)]==FALSE || 
              bNonZeroVec[1+intercept]==FALSE)
            {
                bNonZeroVec[1+intercept] <- FALSE
                bNonZeroVec[length(bNonZeroVec)] <- FALSE
            }
            
            dimW <- dim(lqaW)
            lqaW[,1] <- lqaW[,1] + lqaW[,dimW[2]]
            lqaW[1,] <- lqaW[1,] + lqaW[dimW[1],]
            lqaW <- lqaW[-dimW[1],-dimW[2]]
        }
        
        if(intercept)
        {
            lqaW <- bdiag(0,lqaW)
        }
        
        lqaW <- lqaW / 2
        
        
        if(intercept)  bNonZeroVec <- c(TRUE,bNonZeroVec)
        
        if(wrap)
        {
            bNonZeroVecPrime <- bNonZeroVec[-length(bNonZeroVec)]
            UtU <- t(U[,bNonZeroVecPrime]) %*% U[,bNonZeroVecPrime]
            Ut <- t(U[,bNonZeroVecPrime])
            Vp <- V[bNonZeroVecPrime,bNonZeroVecPrime]
            
            theta <- solve(UtU+Vp+n*lqaW[bNonZeroVecPrime,bNonZeroVecPrime],Ut %*% y)
            bHat <- matrix(0,length(bNonZeroVecPrime),1)
            bHat[bNonZeroVecPrime,1] <- matrix(theta,length(theta),1)
        }
        else
        {
            UtU <- t(U[,bNonZeroVec]) %*% U[,bNonZeroVec]
            Ut <- t(U[,bNonZeroVec])
            Vp <- V[bNonZeroVec,bNonZeroVec]
            
            theta <- solve(UtU+Vp+n*lqaW[bNonZeroVec,bNonZeroVec],Ut %*% y)
            bHat <- matrix(0,length(bNonZeroVec),1)
            bHat[bNonZeroVec,1] <- matrix(theta,length(theta),1)
        }
        
        
        
        if(intercept) bHat <- bHat[-1]
        if(wrap) bHat <- c(bHat,bHat[1])
        
        it <- it + 1
        
    }
    
    bHat
}

Dpfunc <- function(u,lambda,a)
{
    if(u<=lambda) Dpval <- lambda
    else if(u<a*lambda) Dpval <- -(u-a*lambda)/(a-1)
    else Dpval <- 0
    Dpval
}
