###### Li, Gao, and Du 2022 ######
###### Comparison of random design with hetGP package (in Binois et al. 2019) ######

library(DiceKriging)
library(hetGP)
library(foreach)

### Objective function, a bimodal function in [0,1]^2
ftest <- function(x, std = 0.1, center1 = c(0.2,0.8), center2 = c(0.8,0.2)) {
	m1 <- dnorm(x[1], mean = center1[1], sd = std) * dnorm(x[2], mean = center1[2], sd = std)
	m2 <- dnorm(x[1], mean = center2[1], sd = std) * dnorm(x[2], mean = center2[2], sd = std)
	return(m1 - m2)
}

### plot the contour of the objective function ftest()
ngrid <- 51
xgrid <- seq(0,1, length.out = ngrid)
f.true <- apply(expand.grid(xgrid, xgrid), 1, ftest)
contour(x = xgrid, y = xgrid, z = matrix(f.true, ngrid),
			main = "True function", nlevels = 20)	


###### Comparison 1: When the covariate distribution P_X is Uniform[0,1]^2

n <- 200
err.std <- 0.1
x.test.grid <- seq(0,1,by=0.01)		# 1e4 test grid for calculating IMSE
test.grid <- data.frame(expand.grid(x.test.grid,x.test.grid))
colnames(test.grid) <- c("x1","x2")

### Random design: randomly draw 100 points from Uniform[0,1]^2, each points with 2 replicates
### repeated for M = 100 times
set.seed(123)
M <- 100
ran.IMSE <- vector()
test.unif <- data.frame(matrix(runif(2e5), ncol=2))
colnames(test.unif) <- c("x1","x2")
y.true <- apply(test.unif, 1, ftest)

ptm1 <- proc.time()
for (k in 1:M) {
	x <- data.frame(do.call("rbind",rep(list(matrix(runif(n),ncol=2)),2)))
	colnames(x) <- c("x1","x2")
	y <- apply(x, 1, ftest) +  rnorm(n, sd = err.std)
	model.fit <- km(~1, design=x, response=data.frame(y=y), 
					nugget.estim = TRUE,, covtype="gauss", multistart = 5)	
	# covtype="gauss", covtype="matern5_2", covtype=matern3_2"
	model.pred <- predict.km(model.fit, newdata = test.unif, type = "SK")$mean
	ran.IMSE[k] <- mean((model.pred - y.true)^2)		# IMSE of random design
}
ptm2 <- proc.time()

summary(ran.IMSE)	# summary of the IMSE from random design
sd(ran.IMSE)
(ptm2 - ptm1)/M		# average run time of random design including prediction


### Sequential design using Binois et al. (2019)
### Draw 50 initial points with 2 replicates on a 5x5 grid
### And then draw the other 150 points sequentially (may replicate existing points)
n.init <- 25
xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n.init))
designs <- as.matrix(expand.grid(xgrid0, xgrid0))
X <- designs[rep(1:n.init, 2),]
Z <- apply(X, 1, ftest) + rnorm(2*n.init, sd=err.std)
model <- mleHomGP(X, Z, lower = rep(0.01, 2), upper = rep(100, 2), covtype="Gaussian")
# covtype="Gaussian", covtype="Matern5_2", covtype="Matern3_2"
nsteps <- 150 			# nsteps = n - n.init


ncores <- 1
set.seed(4231)
ptm1 <- proc.time()
for(i in 1:nsteps){
	res <- IMSPE_optim(model, h = 5, control = list(multi.start = 100, maxit = 50),
			ncores = ncores)
	# If a replicate is selected
	if(!res$path[[1]]$new) cat("Add replicate \n ")
	newX <- res$par
	newZ <- ftest(newX) + rnorm(1, sd=err.std)
	model <- update(object = model, Xnew = newX, Znew = newZ)
	cat(newX, "\n")
	
	## Plots
	# ngrid <- 51
	# xgrid <- seq(0,1, length.out = ngrid)
	# Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
	# preds <- predict(x = Xgrid, object = model)
	# contour(x = xgrid, y = xgrid, z = matrix(preds$mean, ngrid),
			# main = "Sequential design: Predicted mean", nlevels = 20)
	# points(model$X0, col = 'blue', pch = 20, cex = 1.5)
	# points(newX, col = "red", pch = 20, cex = 1.5)
}
ptm2 <- proc.time()
ptm2 - ptm1

test.unif <- data.frame(matrix(runif(2e5), ncol=2))
colnames(test.unif) <- c("x1","x2")
y.true <- apply(test.unif, 1, ftest)
mean((predict(x = as.matrix(test.unif), object = model)$mean - y.true)^2)		# IMSE of sequential design



###### Comparison 2: When the covariate distribution P_X is bivariate truncated normal 
###### with mean (0.5,0.5), covariance matrix given by Sigma below

### Random design: randomly draw 100 points from bivariate truncated normal, each points with 2 replicates
set.seed(1234)
n <- 200		# total size = 100 x 2 = 200
err.std <- 0.1
Sigma <- matrix(c(1,0.7,0.7,1),2,2) * 0.016
L.Sigma <- chol(Sigma)
mean.x <- c(0.5,0.5)

M <- 100
ran.IMSE <- vector()
test.normal <- data.frame(matrix(rnorm(2e5),ncol=2) %*% L.Sigma + mean.x)	
colnames(test.normal) <- c("x1","x2")	

ptm1 <- proc.time()	
for (k in 1:M) {
	x <- data.frame(do.call("rbind", rep(list(matrix(rnorm(n),ncol=2) %*% L.Sigma + mean.x),2)))
	colnames(x) <- c("x1","x2")
	y <- apply(x, 1, ftest) + rnorm(n, sd = err.std)
	model.fit <- km(~1, design=x, response=data.frame(y=y), nugget.estim = TRUE,
                    covtype="matern3_2",  multistart = 5) 
	# covtype="gauss", covtype=""matern5_2", covtype=""matern3_2" , multistart = 5
	model.pred <- predict.km(model.fit, newdata = test.normal, type = "SK")$mean
	y.true <- apply(test.normal, 1, ftest)
	ran.IMSE[k] <- mean((model.pred - y.true)^2)		# IMSE of random design
}
ptm2 <- proc.time()
ptm2 - ptm1
summary(ran.IMSE)	# summary of the IMSE from random design
sd(ran.IMSE)
(ptm2 - ptm1)/M		# average run time of random design including prediction


### Sequential design using Binois et al. (2019)
### Draw 50 initial points with 2 replicates on a 5x5 grid
### And then draw the other 150 points sequentially (may replicate existing points)
set.seed(1234)
n.init <- 25
xgrid0 <- seq(0.1, 0.9, length.out = sqrt(n.init))
designs <- as.matrix(expand.grid(xgrid0, xgrid0))
X <- designs[rep(1:n.init, 2),]
Z <- apply(X, 1, ftest) + rnorm(2*n.init, sd=err.std)
model <- mleHomGP(X, Z, lower = rep(0.01, 2), upper = rep(100, 2), covtype="Gaussian")
# covtype="Gaussian", covtype="Matern5_2", covtype="Matern3_2" 
nsteps <- 150 		# nstep = n - n.init

ncores <- 1
ptm1 <- proc.time()
for(i in 1:nsteps){
	res <- IMSPE_optim(model, h = 5, control = list(multi.start = 100, maxit = 50),
			ncores = ncores)
	# If a replicate is selected
	if(!res$path[[1]]$new) cat("Add replicate \n ")
	newX <- res$par
	newZ <- ftest(newX) + rnorm(1, sd=err.std)
	model <- update(object = model, Xnew = newX, Znew = newZ)
	cat(newX, "\n")
	
	## Plots
	# ngrid <- 51
	# xgrid <- seq(0,1, length.out = ngrid)
	# Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
	# preds <- predict(x = Xgrid, object = model)
	# contour(x = xgrid, y = xgrid, z = matrix(preds$mean, ngrid),
			# main = "Sequential design: Predicted mean", nlevels = 20)
	# points(model$X0, col = 'blue', pch = 20, cex = 1.5)
	# points(newX, col = "red", pch = 20, cex = 1.5)
}
ptm2 <- proc.time()
ptm2 - ptm1

test.normal <- data.frame(matrix(rnorm(2*1e5),1e5,2) %*% L.Sigma + mean.x)
y.true <- apply(test.normal, 1, ftest)
mean((predict(x = as.matrix(test.normal), object = model)$mean - y.true)^2)		# IMSE of sequential design


### Draw plots
pdf("binois.pdf", width=18, height=6)
par(mfrow=c(1,3))

# Plot the objective function and the covariate distribution
set.seed(123)
ngrid <- 101
xgrid <- seq(0,1, length.out = ngrid)
f.true <- apply(expand.grid(xgrid, xgrid), 1, ftest)
test.normal <- data.frame(matrix(rnorm(2*10000),10000,2) %*% L.Sigma + mean.x)	
contour(x = xgrid, y = xgrid, z = matrix(f.true, ngrid),
			main = "Contour of true f(x) and truncated normal", nlevels = 20,
			cex.main=2, cex.axis=2)
points(test.normal, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2), pch = 20, cex = 1)

# Plot the surface fitted by a random design
set.seed(4231)
x <- data.frame(do.call("rbind", rep(list(matrix(rnorm(n),n/2,2) %*% L.Sigma + mean.x),2)))
colnames(x) <- c("x1","x2")
y <- apply(x, 1, ftest) + rnorm(n, sd = err.std)
model.fit <- km(~1, design=x, response=data.frame(y=y), nugget.estim = TRUE,
                covtype="gauss",  multistart = 5) 
test.grid <- data.frame(expand.grid(xgrid,xgrid))		
colnames(test.grid) <- c("x1","x2")	
model.pred.grid <- predict.km(model.fit, newdata = test.grid, type = "SK")$mean
model.pred.grid

contour(x = xgrid, y = xgrid, z = matrix(model.pred.grid, length(xgrid)),
	    main = "Random design: Predicted function and data cloud", nlevels = 20,
	    cex.main=2, cex.axis=2)	
points(x, col = "blue", pch = 20, cex = 1.5)

# Plot the surface fitted by sequential design
preds <- predict(x = as.matrix(test.grid), object = model)
contour(x = xgrid, y = xgrid, z = matrix(preds$mean, ngrid),
		main = "Sequential design: Predicted function and data points", nlevels = 20,
		cex.main=2, cex.axis=2)
points(model$X0, col = 'blue', pch = 20, cex = 1.5)

dev.off()

