#!/usr/bin/Rscript --vanilla

###########################################################
##                                                       ##
##   bayes-reg.R                                         ##
##                                                       ##
##                Author: Tony Fischetti                 ##
##                        tony.fischetti@gmail.com       ##
##                                                       ##
###########################################################

# workspace cleanup
rm(list=ls())

# options
options(echo=TRUE)
options(stringsAsFactors=FALSE)
options(datatable.fread.datatable=FALSE)
options(mc.cores = parallel::detectCores())

# cli args
args <- commandArgs(trailingOnly=TRUE)

# libraries
library(glmnet)
library(magrittr)
library(dplyr)
library(rethinking)
library(pbapply)
library(cvTools)
library(ggplot2)
library(tidyr)
library(boot)





#--------------------------------------------------#
#--------------------------------------------------#
#                    FUNCTIONS                     #
#--------------------------------------------------#
#--------------------------------------------------#
make.map.formula <- function(startofalist, numofcoefs, stddev){
  # doesn't count the intercept
  gencoeffs <- paste0("b", 1:numofcoefs)
  this <- list()
  for(i in 1:length(gencoeffs)){
    this <- c(this, list(bquote(.(as.name(gencoeffs[i])) ~ dnorm(0, .(stddev)))))
  }
  return(c(startofalist, this))
}

get.folds <- function(df, k=5){
  tmp <- cvFolds(nrow(df), K=k)
  folds <- lapply(1:length(unique(tmp$which)),
                  function(x) tmp$subsets[which(tmp$which!=x),])
  return(folds)
}

pevaluate <- function(X, y, design.matrix, basealist, prior.stddev,
                      startlist, k=5, parallel=TRUE, show.result=TRUE){
  if(parallel){
    mapcarfn <- function(...){ unlist(mclapply(...)) }
  } else{
    mapcarfn <- sapply
  }
  df <- cbind(X, y)
  form <- make.map.formula(basealist, ncol(X), prior.stddev)
  folds <- get.folds(df, k=k)
  mserrors <- mapcarfn(1:length(folds),
                       function(x){
                         rows <- folds[[x]]
                         fit <- map(form, data=df[rows,], start=startlist)
                         best.bayes.coefs <- coef(fit)[-length(coef(fit))] %>% matrix
                         pred.bayes <- design.matrix %*% best.bayes.coefs
                         the.errors <- pred.bayes - y[,1]
                         serrors <- the.errors[-rows]^2
                         return(mean(serrors))
                       })
  erfit <- map(alist(mserrors ~ dnorm(mu, stddev),
                     mu ~ dnorm(0,100),
                     stddev ~ dunif(0,100)),
               data=data.frame(mserrors=mserrors),
               start=list(mu=8, stddev=10))
  emean <- mean(extract.samples(erfit)$mu)
  hdi95 <- HPDI(extract.samples(erfit)$mu, prob=.95)
  fullfit <- map(form, data=df, start=startlist)
  best.bayes.coefs <- coef(fullfit)[-length(coef(fullfit))] %>% matrix
  retlist <- list(coefs=best.bayes.coefs, cv.mse=emean, hdi=hdi95)
  if(show.result){
    cat(sprintf("\n%s\n  LOOCV MSE ESTIMATE FOR PRIOR WIDTH of %.3f: %.3f\n%s\n",
                paste0(rep("*", 55), collapse=""),
                prior.stddev,
                emean,
                paste0(rep("*", 55), collapse="")))
  }
  return(retlist)
}



#==================================================#
#==================================================#
#=============         MTCARS         =============#
#=============         MTCARS         =============#
#=============         MTCARS         =============#
#==================================================#
#==================================================#

# first let's use ridge regression and get info
# on the highest performing model fit

mtstd <- mtcars %>% lapply(scale) %>% data.frame
row.names(mtstd) <- row.names(mtcars)
mtstd$mpg <- mtcars$mpg

design.matrix <- model.matrix(mpg ~ ., data=mtstd)

X <- design.matrix[,-1]
y <- mtstd[, 1, drop=FALSE]

cvfits <- cv.glmnet(X, y[,1], alpha=0, nfolds=10)
# plot(cvfits)


loc <- which(cvfits$lambda==cvfits$lambda.min)
best.coefs <- coef(cvfits, s="lambda.min")
mse <- cvfits$cvm[loc]                           # 7.336113





####### messing around
basealist <- alist(mpg ~ dnorm(mu, sigma),
                   mu <- b0 + b1*cyl + b2*disp + b3*hp + b4*drat + b5*wt + b6*qsec + b7*vs + b8*am + b9*gear + b10*carb,
                   sigma ~ dunif(0,100),
                   b0 ~ dnorm(20, 100))
form <- make.map.formula(basealist, 10, 5)

startlist <- list(b0=20, b1=0, b2=0, b3=0, b4=0,
                  b5=0, b6=0, b7=0, b8=0, b9=0,
                  b10=0, sigma=30)


this <- map(form, data=mtstd, start=startlist)
best.bayes.coefs <- coef(this)[-length(coef(this))] %>% matrix
preds.bayes <- design.matrix %*% best.bayes.coefs
serrors <- (preds.bayes - mtstd$mpg)^2
mse <- mean((preds.bayes - mtstd$mpg)^2)         # 4.622117
mse

tmp <- map(alist(serrors ~ dnorm(mu, stddev),
                 mu ~ dnorm(5, 10),
                 stddev ~ dunif(0, 20)),
           data=data.frame(serrors=serrors),
           start=list(mu=4, stddev=6))


# example
pevaluate(X, y, design.matrix, basealist, 0.575, startlist, k=10)




#--------------------------------------------------#
# 10-fold CV with mtcars predicting `mpg` from     #
# all other variables and gaussian priors of       #
# increasing precision                             #
#--------------------------------------------------#

prior.widths <- seq(0.05, 5, by=0.025)

results <- pblapply(prior.widths, function(x){
                  tryCatch({
                    pevaluate(X, y, design.matrix, basealist, x, startlist, k=10)
                  }, error = function(e){
                    return(list(coeffs=rep(NA, 11), cv.mse=NA, hdi=c(NA, NA)))
                  })
               })
MTCARSK10RESULTS <- results

tresults.cv <- lapply(results, function(x) x$cv.mse) %>% unlist
tresults.lower <- lapply(results, function(x) x$hdi[1]) %>% unlist
tresults.upper <- lapply(results, function(x) x$hdi[2]) %>% unlist
tresults.width <- prior.widths

tresults.df <- data.frame(mse=tresults.cv, width=tresults.width,
                          lower=tresults.lower,
                          upper=tresults.upper)

plot(mse ~ width, data=tresults.df, type="l", ylim=c(0, 22),
     main="10K CV MSE as a function of prior width of coefficients")
lines(upper ~ width, data=tresults.df, type="l", col="red")
lines(tresults.df$width, tresults.df$lower, type="l", col="red")
abline(v=0.775, lty=2)




#--------------------------------------------------#
# LOOCV with mtcars predicting `mpg` from          #
# all other variables and gaussian priors of       #
# increasing precision                             #
#--------------------------------------------------#

prior.widths <- seq(0.05, 5, by=0.025)

results <- pblapply(prior.widths, function(x){
                  tryCatch({
                    pevaluate(X, y, design.matrix, basealist, x, startlist, k=31)
                  }, error = function(e){
                    return(list(coeffs=rep(NA, 11), cv.mse=NA, hdi=c(NA, NA)))
                  })
               })
## up to here :)
MTCARSLOORESULTS <- results

tresults.cv <- lapply(results, function(x) x$cv.mse) %>% unlist %>% .[5:length(prior.widths)]
tresults.lower <- lapply(results, function(x) x$hdi[1]) %>% unlist %>% .[5:length(prior.widths)]
tresults.upper <- lapply(results, function(x) x$hdi[2]) %>% unlist %>% .[5:length(prior.widths)]
tresults.width <- prior.widths %>% .[5:length(prior.widths)]


tresults.df <- data.frame(mse=tresults.cv, width=tresults.width,
                          lower=tresults.lower,
                          upper=tresults.upper)

minxintercept <- tresults.width[which.min(tresults.cv)]

plot(mse ~ width, data=tresults.df, type="l", ylim=c(0, 22),
     main="LOOCV MSE as a function of prior width of coefficients")
lines(upper ~ width, data=tresults.df, type="l", col="red")
lines(tresults.df$width, tresults.df$lower, type="l", col="red")
abline(v=minxintercept, lty=2)


ggplot(tresults.df, aes(width, mse)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, alpha=0.5), fill = "grey70",
              show.legend=FALSE) +
  geom_line() +
  geom_vline(xintercept=minxintercept, linetype=2) +
  ylab("LOOCV MSE") +
  xlab("standard deviation of gaussian coefficient priors") +
  ggtitle('Bayesian "ridge regression" LOOCV MSE with different "penalties" (mtcars)') +
  ggsave("./plots/mtcars-loocv-mse.png") +
  ggsave("./plots/mtcars-loocv-mse.pdf")



tmp <- lapply(1:length(results),
         function(x){ df <- data.frame(pred=names(mtcars)[-1],
                                 cvalue=results[[x]]$coefs[-1,1])
           df %>% spread(pred, cvalue) %>% cbind(data.frame(width=prior.widths[x]), .)
         })

do.call(rbind, tmp) -> coefdf

gcoefdf <- coefdf %>% gather(width, cvalue)
names(gcoefdf)[2] <- "predictor"

ggplot(gcoefdf, aes(x=width, y=cvalue, color=predictor)) +
  geom_line() +
  ylab("coefficient value") +
  xlab("standard deviation of gaussian coefficient priors") +
  geom_vline(xintercept=minxintercept, linetype=2, color="grey") +
  ggtitle('Coefficient shrinkage in bayesian "ridge regression" (mtcars)') +
  ggsave("./plots/mtcars-coef-shrinkage.png") +
  ggsave("./plots/mtcars-coef-shrinkage.pdf")



netbest.coefs <- best.coefs %>% as.matrix %>% .[-1,] %>% data.frame

netbest.coefs$predictor <- row.names(best.coefs)[-1]
row.names(netbest.coefs) <- NULL
names(netbest.coefs)[1] <- "ncvalue"

agcoefdf <- gcoefdf %>% left_join(netbest.coefs)

gagcoefdf <- agcoefdf %>% gather(width, predictor)
names(gagcoefdf)[3] <- "bayes_or_net"
gagcoefdf$bayes_or_net <- ifelse(gagcoefdf$bayes_or_net=="cvalue",
                                 "bayes", "elastic net")
names(gagcoefdf)[4] <- "value"

ggplot(gagcoefdf, aes(x=width, y=value, color=predictor, linetype=bayes_or_net)) +
  geom_line() +
  ylab("coefficient value") +
  xlab("standard deviation of gaussian coefficient priors") +
  geom_vline(xintercept=minxintercept, linetype=2, color="grey") +
  ggtitle('Coefficient shrinkage in bayesian "ridge regression" (mtcars)') +
  ggsave("./plots/mtcars-coef-shrinkage-net-overlay.png") +
  ggsave("./plots/mtcars-coef-shrinkage-net-overlay.pdf")


