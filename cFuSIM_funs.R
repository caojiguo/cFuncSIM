library(ggplot2)
library(fda)
library(dplyr)
# library(NonpModelCheck)
library(DiceKriging)
library(reshape2)

wtemp <- function (xbasis0) 
{
  xbasis0 %>% knots(, interior = FALSE) %>% unique %>% data.frame(knot = .) %>% 
    mutate(knotlead = lead(knot)) %>% dplyr::filter(!is.na(knotlead)) %>% 
    rowwise() %>% do(temp = eval.penalty(xbasis0, int2Lfd(0), 
                                         rng = c(.$knot, .$knotlead)))
}

# cv to select the optimal parameters given y scalar responses and xfds the functional predictor


cv_cFuSIM = function(y,xfds, log10lambda=c(),gamma=c(),bandwidth=c(),spline_basis=xfds$basis, k=3, ...){
  parameters = expand.grid( log10lambda =  log10lambda, gamma=gamma, bandwidth=bandwidth)
  cv_res = c()
  cv_folds = caret::createFolds(1:length(y),  k, list = F)
  for (j in 1:nrow(parameters)){
    log10lambda = parameters$log10lambda[j]
    gamma = parameters$gamma[j]
    bandwidth = parameters$bandwidth[j]
    cv_errors = c()
    for (i in 1:max(cv_folds)){
      test = which(cv_folds==i)
      train = which(cv_folds!=i)
      ytrain=y[train]
      xfdtrain = xfds[train]
      ytest=y[test]
      xfdtest = xfds[test]
      cv_errors = c(cv_errors,predict_cFuSIM(ytrain, xfdtrain, ytest, xfdtest,spline_basis,
                                              log10lambda=log10lambda,
                                              gamma =gamma,
                                              bandwidth = bandwidth,
                                              ...))
    }
    cv_res = c(cv_res,mean(cv_errors))
  }
  parameters$cv_errors=cv_res
  parameters = parameters%>%arrange(cv_errors)
  return(parameters)
}

# W0 <- function (beta_j, lambda, bspi, W_m, normthres = 10^-2) 
# {
#   W = inprod(bspi, bspi, 1, 1)
#   R = inprod(bspi, bspi)
#   r = inprod(bspi, bspi, 2, 2)
#   range_knots = fda:::knots.basisfd(bspi, interior = FALSE) %>% 
#     range
#   wi = function(i) {
#     zero = NULL
#     beta_m_norm = t(as.matrix(beta_j)) %*% W_m[[i]] %*% as.matrix(beta_j) %>% 
#       as.numeric %>% sqrt
#     #print(beta_m_norm)
#     if (beta_m_norm < normthres) {
#       zero = c(i:(i + 3))
#       res = matrix(0, nrow = bspi$nbasis, ncol = bspi$nbasis)
#     }
#     else {
#       temp1 = SCAD.derivative(sqrt((length(W_m))/((bspi$rangeval[2] - 
#                                                      bspi$rangeval[1]))) * beta_m_norm, lambda = lambda)
#       temp2 = beta_m_norm * sqrt(((bspi$rangeval[2] - bspi$rangeval[1]))/(length(W_m)))
#       res = temp1/temp2 * W_m[[i]]
#     }
#     return(list(res = res, zero = zero))
#   }
#   res_list = lapply(1:length(W_m), wi)
#   res = lapply(res_list, function(x) x[[1]])
#   zero = do.call(c, lapply(res_list, function(x) x[[2]])) %>% 
#     unique
#   W_j0 = Reduce("+", res) * 0.5
#   return(list(W = W_j0, zero = zero))
# }



cFuSIM_index <- function(y, xfd, spline_basis, threshold=1e-5, maxit=150,log10lambda=6, gamma=1e5, bandwidth=NULL, verbose=FALSE){
  bandwidth = ifelse(is.null(bandwidth),"CV", bandwidth)
  lambda = 10^log10lambda
  (betac= rep(1,spline_basis$nbasis))
  #@betac = runif(spline_basis$nbasis)
  #betac = 20+rnorm(spline_basis$nbasis)
  B = inprod(spline_basis,spline_basis,0,0)
  betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
  GAMMA = inprod(spline_basis,spline_basis,2,2)
  W_m = wtemp(spline_basis )$temp
  thes = 1
  i = 1
  thresb= 3
  while(thes>threshold){
    i=i+1
    if(i>maxit) break
    betac0 = betac
    pmatf = pmatf_locploy(betac,y,xfd,spline_basis,bandwidth=bandwidth)
    yweighted = as.numeric(pmatf$yweighted)
    nrowx = nrow(coef(xfd))
    xfds_weighted = fd(rep(pmatf$weights, each=nrowx)%>%matrix(nrow=nrowx)*coef(xfd), xfd$basis)
    beta_basis = list()
    beta_basis[[1]] = spline_basis
    fit <- try(slos(xfds_weighted,yweighted%>%as.numeric,D=NULL,lambda=lambda,gamma=gamma,intercept=F,beta.basis=beta_basis,cutoff=1e-4,
                          max.iter=1000,tol=1e-12,a=3.7,domain=NULL, tuning = "BIC"), silent=TRUE)
    if(class(fit)=="try-error") 
    {
      print('error')
      betac = rep(0, spline_basis$nbasis) 
      thes = 0
    } else {
      bandwidth = pmatf$bandwidth
      betac = fit$beta[[1]]%>%coef%>%as.numeric
      pc_fit = fit$beta[[1]]
      betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
      pc_fit = fd(betac,spline_basis)
      score_fit = inprod(xfd,pc_fit)%>%as.numeric
      if(as.numeric(inprod(pc_fit))<0) {
        betac= -betac
        pc_fit = fd(betac,spline_basis)
      }
      
      thes  = abs(betac-betac0)%>%max
      thresb =thes
      if(verbose){
      plot(pc_fit)
      print(thresb)
      }
    }
  }
  
  converaged=(thes< threshold)
  if(!all(betac==0)){
    betac = betac/sqrt(as.numeric(t(betac)%*%B%*%betac))
  }
  pc_fit = fd(betac,spline_basis)
  score_fit = inprod(xfd,pc_fit)%>%as.numeric
  list(coefBeta=betac, basisBeta = spline_basis, score_fit = score_fit, Converaged=converaged, threshold=threshold, maxit=maxit, bandwidth = bandwidth)
}


predict_cFuSIM = function(ytrain, xfdtrain, ytest=ytrain, xfdtest=xfdtrain, ... ){
  res_cFuSIM = cFuSIM_index(ytrain,xfdtrain,...)
  if(all(res_cFuSIM$coefBeta==0)) {
    print('all coefs are zero!')
    return(mean((ytest-mean(ytrain))^2))
  } else {
    beta_fd = fd(res_cFuSIM$coefBeta, res_cFuSIM$basisBeta)
    plot(beta_fd)
    score_fit_test = inprod(xfdtest,beta_fd)%>%as.numeric
    score_fit_train = inprod(xfdtrain,beta_fd)%>%as.numeric
    score_fit = score_fit_train
    sd_train = sd(score_fit)
    mean_train = mean(score_fit)
    # score_fit = scale(score_fit)
    score_fit_train = scale(score_fit_train, center = mean_train, scale = sd_train )[,1]
    score_fit_test = scale(score_fit_test, center = mean_train, scale = sd_train )[,1]
    y = ytrain
    y = y[order(score_fit_train)]
    score_fit_train = score_fit_train[order(score_fit_train)]
    # plot the index function 
    pred_y = localpoly.reg(score_fit_train, y, points=score_fit_test,degree.pol = 1, kernel.type = "gaussian",bandwidth = res_cFuSIM$bandwidth,deriv=0)
    # gbasis= create.bspline.basis(rangeval=range(c(score_fit_test,score_fit_train)),nbasis=10,norder=3)
    # gfd = Data2fd(argvals=score_fit_train, y=pred_y$predicted, basisobj = gbasis)
    # plot(gfd)
    # points(x=score_fit_train, y=pred_y$predicted,col=2)
    # y_pred_test = as.numeric(eval.fd(gfd,score_fit_test))
    # return(mean((ytest-y_pred_test )^2))
    return(mean((ytest-pred_y$predicted )^2))
  }
}


pmatf_locploy <- function (betac, ytrain = ytrain, xtrain = xtrain, spline_basis, bandwidth="CV")
{
  B = inprod(spline_basis, spline_basis, 0, 0)
  pc_fit = fd(betac, spline_basis)
  Z = inprod(xtrain, spline_basis)
  score_fit = inprod(xtrain, pc_fit) %>% as.numeric
  score_fit = score_fit%>%scale
  bandwidth0 = bandwidth
  loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1,
                          kernel.type = "gaussian", bandwidth = bandwidth0 , deriv = 0)
  select_bandwith = loc_fit$bandwidth
  fity = loc_fit$predicted
  loc_fit = localpoly.reg(score_fit, ytrain, degree.pol = 1,
                          kernel.type = "gaussian", bandwidth = loc_fit$bandwidth, deriv = 1)
  derfity = loc_fit$predicted
  Gmat = list()
  yweighted = ytrain - fity + derfity*(Z%*%betac)
  list( bandwidth = select_bandwith, weights = -derfity, yweighted = yweighted, Z = Z)
}



yunlongtheme <- function(fontscale=1) {
  # library(fonts)
  theme_bw() + 
    theme(legend.position="none", 
          # legend.text = element_text(size = 22*fontscale, face = "bold"), 
          # legend.title=element_text(color = "#DB4437", size=22*fontscale,face = "bold"),
          axis.title=element_text(color = "#DB4437",size=22*fontscale,face="bold"), 
          axis.text = element_text(color = "#4285F4",face="bold",size=20*fontscale),
          plot.title = element_text(color="#0F9D58",size=25,face="bold",hjust=0.5),
          strip.text.x = element_text(color = "#0F9D58", size=25*fontscale,  angle = 0, face="bold"),
          strip.background =element_rect(fill="grey95")
    )
}

theme_set(yunlongtheme(fontscale=1.1))


fdagg  = function(fdobj, ylab, xlab,addzero=TRUE, fontscale=10){
  range = fdobj$basis$rangeval
  t = seq(range[1], range[2], len = 1000)
  names_y = fdobj$fdnames$reps
  if(names_y%>%is.null) names_y = paste0("rep",1:(eval.fd(t, fdobj) %>%ncol))
  data = data.frame(Time = t, eval.fd(t, fdobj))
  names(data) = c("Time", names_y)
  plotdata = melt(data, id = 1)
  plot = ggplot(plotdata, aes(x = Time, y = value, group = factor(variable), 
                              color = factor(variable))) + geom_line(size=3, aes(linetype=factor(variable))) +  ylab(ylab) + xlab(xlab)+xlim(c(1,24))
  plot
  if (addzero) {
    plot = plot + geom_hline(yintercept = 0, linetype = 2, 
                             col = 1)+ theme(legend.position = "none")
  }
  plot
}
