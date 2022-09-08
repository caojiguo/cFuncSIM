# install the following package 
# install.packages('fda')
# install.packages('dplyr')
# install.packages('ggplot2')
# install.packages('DiceKriging')
# install.packages('NonpModelCheck') # this package is not available in cran anymore, but you can install it using 
# install.packages('NonpModelCheck_3.0.tar.gz') # you can find NonpModelCheck_3.0.tar.gz in the folder
# install.packages('dr')
# install.packages('caret')
# $version.string
# [1] "R version 3.6.1 (2019-07-05)"

## Please change the work directory to the folder containing these R files
#setwd('/Users/jiguo/Dropbox/Research/Yunlong2/single_index_compact_support_draft/EnvironmetricsRevision2022Aug12/Rcode2022Sep07')


source('slos.R')
source('cFuSIM_funs.R')

## load the real data 
source('bike_dataset.R')

# define the time points on which the functional predictor is observed. 
timepts = bike$timepts

# define the bspline basis 
norder= 8 ## cubic B-spline
nbasis=norder+length(timepts)-2; 
# 
knots = fda:::knots.basisfd
spline_basis=create.bspline.basis(rangeval=c(1,24),nbasis,norder,timepts)
plot(spline_basis)
# convert the functional predictor into a fda object
wull = bike$temp
xfds=  Data2fd(y=wull%>%t, argvals=bike$timepts, basis=spline_basis)

# define the response value 
y = bike$y

# use cv to obtain the optimal tunning parameter settings 
#res_cv = cv_cFuSIM(y,xfds,log10lambda = c(2,4,6,8), gamma = c(1e1,1e3,1e5),bandwidth = c(0.2,0.5,0.8))
#para_cv = res_cv%>%dplyr::arrange(cv_errors)%>%head(1)
#cFuSIM_fit = cFuSIM_index(y, xfds, spline_basis, log10lambda=para_cv$log10lambda,  
#                           gamma = para_cv$gamma, bandwidth=para_cv$bandwidth)
# we skip the cv part for now 

# estimate the index function and link function
cFuSIM_fit = cFuSIM_index(y, xfds, spline_basis, log10lambda=6, gamma = 1e1, bandwidth=0.5) order 12 

# index function estimate
beta_fd = fd(cFuSIM_fit$coefBeta, cFuSIM_fit$basisBeta)

plot = fdagg(beta_fd ,ylab="", xlab='hour',addzero=TRUE)+ scale_x_continuous(breaks=seq(0,24, by=2))+
  ggtitle(expression("Estimated Index Function"~hat(beta)(.)))
plot

# link function estimate
score_fit = inprod(xfds,beta_fd12)%>%as.numeric%>%scale%>%as.numeric
pred_y = localpoly.reg(score_fit, y, points=seq(-2,2,length.out = 50),degree.pol = 1, kernel.type = "gaussian",
                       bandwidth = 0.5,deriv=0)

pred_g = data.frame(x=seq(-2,2,length.out = 50),y=pred_y$predicted)

gplot = qplot(x=score_fit, y, geom = "point")+geom_point(size=3)+geom_line(data=pred_g, aes(x=x,y=y), col=4, size=2)+
  xlab(expression(bold(X)^T~hat(beta)))+ylab('')+ggtitle(expression("Estimated "~hat(g)(.)~"with order 8 B-spline"))
gplot

