# set the folder containing the funs_SOAP.R file 
setwd('.')
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(locfit))
suppressPackageStartupMessages(library(fdapace))
suppressPackageStartupMessages(library(lsei))

# install these packages if any of them is missing 
# sapply(c("dplyr","fda","parallel","locfit","fdapace","lsei"), install.packages)
source('functions_soap.R')

##  generate the eigen functions used to generate the observations
tempdat = daily$tempav
all_matrix=tempdat;all_matrix = sweep(all_matrix,1,rowMeans(all_matrix))
timepts=0:364;
norder=4 ## cubic B-spline
nbasis=norder+length(timepts)-2;
spline.basis=create.bspline.basis(rangeval=c(0,364),nbasis=10,norder=4)
D2Lfd <- int2Lfd(m=2)
D2fdPar <- fdPar(spline.basis, D2Lfd, 1e1)
all.fd =  Data2fd(y=all_matrix, argvals=timepts)
fpca = pca.fd(all.fd,5,D2fdPar)
seed_set=1010010
set.seed(seed_set)
nsim = 300
nobs = sample(x=2:5,size=nsim, replace=TRUE)
sdsn = c(sqrt(800), sqrt(160), sqrt(32), sqrt(6.4), sqrt(1.6))
scores_mat = lapply(sdsn, function(x) {
	temp = rnorm(nsim, mean=1,sd=x)
	return(temp)
})%>%do.call(cbind, .)
yfds = fd(t(scores_mat%*%(fpca$harmonics%>%coef%>%t)), fpca$harmonics$basis)
timepoints = lapply(1:length(nobs), function(x) runif(nobs[x],0,364)%>%sort)
true_y = lapply(1:length(nobs), function(x) eval.fd(timepoints[[x]], yfds[x])%>%as.numeric)
sigma = 0.05
observed = lapply(1:length(nobs), function(x) true_y[[x]]+rnorm(length(true_y[[x]]), sd=sigma ))

########################################
################PACE###################
########################################

res_pace<- FPCA(observed, timepoints,list(dataType='Sparse',error=TRUE, kernel='epan', 
	verbose=FALSE,methodBwCov="GCV",methodBwMu="GCV"))
select_k = SelectK(res_pace, criterion = 'AIC')$K
pred_pace = predict(res_pace, observed, timepoints,K=select_k )
error_pace = c()
for(i in 1:length(observed)){
ytest_fit = rowSums(matrix(rep(pred_pace[i,],each=nrow(res_pace$phi[,1:select_k])),nrow=nrow(res_pace$phi[,1:select_k]))*res_pace$phi[,1:select_k])+res_pace$mu
basist = create.bspline.basis(rangeval=c(0,364),nbasis=length(res_pace$workGrid)-4)
yfit = Data2fd(argvals=res_pace$workGrid,y=ytest_fit,basisobj = basist)
error_pace= c(error_pace,inprod(yfit-yfds[i],yfit-yfds[i]))
}
coef_mat0 = coef(Data2fd(argvals = res_pace$workGrid,y=res_pace$phi,spline.basis))
mean(error_pace)
summary(error_pace)


########################################
################SOAP###################
########################################
spline_basis=create.bspline.basis(rangeval=c(0,364),nbasis=10,norder=4)
observed%>%do.call(c,.)%>%mean
pc1s = first_FPC(coef_mat0[,1],observed=observed, timepoints=timepoints,minit=6,gamma=1e1,threshold=1e-5)
previous_beta = list()
previous_beta[[1]] = pc1s$beta
pc2s = third_FPC_conditional(coef_mat0[,2], observed=observed, timepoints=timepoints, pc_index=2, gamma=3e4,betalist =previous_beta,threshold=1e-4)
previous_beta[[2]] = pc2s$beta
pc3s = third_FPC_conditional(coef_mat0[,3], observed=observed, timepoints=timepoints, pc_index=3, gamma=2e2,betalist =previous_beta,threshold=1e-3)
previous_beta[[3]] = pc3s$beta
pc4s = third_FPC_conditional(coef_mat0[,4], observed=observed, timepoints=timepoints, pc_index=4, gamma=1e4,betalist =previous_beta,threshold=1e-3)
previous_beta[[4]] = pc4s$beta


#########
###AIC###
#########
previous_beta0 = previous_beta
observed2 = observed[which(sapply(observed,length)>1)]
timepoints2 = timepoints[which(sapply(observed,length)>1)]
tempy = observed2
sd_score	 = c()
beta_samples= list()
sigma_est = c()
for (i in 1:length(previous_beta0)){
	print(i)
	res = pred_SOAP_step(previous_beta0[i],	tempy, timepoints2, spline_basis,sigma=0,nminus=0)
	tempy  =res$residuals
	sigma_est = c(sigma_est ,mean((res$residuals%>%do.call(c,.))^2))
	beta_samples[[i]]=as.numeric(res$sfit)
	sd_score = c(sd_score, res$sfit%>%apply(.,2,sd))
}

N = sapply(observed2,length)%>%sum
n = length(observed2)
AIC = N*log(sigma_est ) + N  + 2*n*c(1:length(previous_beta0))
(k_selet = as.numeric(which.min(AIC)))
observed2 = observed[which(sapply(observed,length)==5)]
timepoints2 = timepoints[which(sapply(observed,length)==5)]
tempy = observed2
res = pred_SOAP_step(previous_beta0[1:min(4,k_selet)],	tempy, timepoints2, spline_basis,sigma=0,nminus=0)
(sig = sqrt(mean((res$residuals%>%do.call(c,.))^2)))
betalist=previous_beta0

###############################
### estimating the individual scores###
###############################

scores_est = function(x) {
library(locfit)
(yi=observed[[x]])
(timei=timepoints[[x]])
xmat = lapply(1:k_selet,function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
betai = c()
llike = function(beta10) {
liklog = dnorm(yi, as.numeric(xmat%*%beta10), sd=sig,log=TRUE)
sum(liklog)+ sapply(1:k_selet, function(x) {
	log(density.lf(beta_samples[[x]],ev=beta10[x])$y)
})%>%sum
}
beta10 = sapply(beta_samples[1:k_selet],mean)
betap = lapply(1:k_selet, function(x) {c()})
for (i in 1:100){
for (j in 1:k_selet){
# betap[[j]]=c()
beta11 = beta10
# beta11[j] = sample(size=1,x = beta_samples[[j]])
beta11[j] = rnorm(1,beta10[j],sd= 0.1*sd(beta_samples[[j]]))
if (runif(1)<exp(llike(beta11)-llike(beta10))){
	beta10=beta11
} 
betap[[j]]  =c(betap[[j]],beta10[j])
}
}
sapply(betap, function(x) {x%>%tail(500)%>%mean})%>%return
}

## estimated scores 
score_pred= mclapply(1:length(observed),scores_est)%>%do.call(rbind,.)
coefs_fd = t(score_pred%*%do.call(rbind,betalist[1:k_selet]))
## estimated trajectories 
yfitsfd = fd(coefs_fd, spline_basis)
## compute the errors 
error_soap = sapply(1:nrow(score_pred), function(x){
	inprod(yfds[x] - (yfitsfd)[x],yfds[x] - (yfitsfd)[x])
})
error_soap
error_soap%>%summary%>%print
error_pace%>%summary%>%print