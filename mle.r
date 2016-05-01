#This set of functions works to produce an MLE for
#a Mixture of 2 Gammas, single gamma, and mix of 2 normals
#Primary coders: Connor Melson, Jake Martinez
#support coders: Magda Moses

#function for the log liklihood of gamma which is to
#be optimized over with optim.
mixGamma <- function(theta, x)
{
	#sets passed parameters to internals that
	#make them more easily identified
	A1 = theta[1];
	B1 = theta[2];
	A2 = theta[3];
	B2 = theta[4];
	p = theta[5];

	#floor the values of alphas and betas
	#as gamma distributions for alpha or beta
	#values less than 0 do not exist
	if(A1 <= 0 || B1 <= 0 || A2 <= 0 || B2 <= 0)
	{
		return (999999999);
	}

	#constrains omega which can only be between 0 and 1
	if(p < 0 || p > 1)
	{
		return (999999999);
	}
	 
	#returns the negative log sum of the densities 
	#which is used by optim for optimization
	logl <- sum(log(p * dgamma(x, shape = A1,  scale = B1) 
		+ (1 - p) * dgamma(x, shape = A2,  scale = B2)))
	return(-logl)
}

#Function for producing log likelihood of singular gamma
#This funcions similarly to mix gamma, but with less params
singleGamma <- function(theta, x)
{
	#sets params to locals
	A <- theta[1];
	B <- theta[2];
	  
	#floors values for alpha and beta
	if(A <= 0 || B <= 0)
	{
		return (999999999);
	}

	#returns -logsum result which optim uses to optimize
	logl <- sum(log(dgamma(x, shape=A, scale=B)))
	return(-logl)
}

#function for mix normal that is optimized
mixNormal <- function(theta,x)
{
	mu1 = theta[1];
	sd1 = theta[2];
	mu2 = theta[3];
	sd2 = theta[4];
	p = theta[5];
	 
	#floor the values of omega
	if(p < 0 || p > 1)
	{
		return (999999999);
	}
	#returns the negative log sum result for optim to optimize
	logl <- sum(log(p * dnorm(x, mean = mu1, sd = sd1) + 
		(1 - p)* dnorm(x, mean = mu2, sd = sd2)))
	return(-logl)
}

#primary function called by the user
#param 1 is alpha1
#param 2 is beta1
#param 3 is alpha2
#param 4 is beta2
#param 5 is omega
#param 6 is number of trials
#param 7 is samples per trial
#param 8 is distribution type

MLE <-function(p1,p2,p3,p4,p5,p6,p7,p8)
{
	#sets internal values to passed parameters
	a1 = p1;
	b1 = p2;
	a2 = p3
	b2 = p4
	omega = p5
	numTests = p6
	isNormal = identical("normal",p8)
	n = p7
	#arrays for our results to go to
	param1=c(0)
	param2=c(0)
	param3=c(0)
	param4=c(0)
	param5=c(0)

	#defaulted sample size per test

	for(test in 1:numTests)
	{
		#gets a sample of value which will decide
		#what distribution the observation will come from
		sample.dist <- sample(1:100, n, replace = T);

		#a variable for the number of occurrences in the given side
		left.num = 0;
		right.num = 0;
		#iterates through the generated sample 
		#and adds 1 to the appropriate side
		for(idx in 1:n)
		{
			if (sample.dist[idx]/100 <= omega)
			{
				left.num = left.num +1;
			}
			else
			{
				right.num = right.num+1;
			}  
		}


		if(isNormal) {
			#generates normal1 whose contribution coefficient is omega
			leftSet = rnorm(left.num, mean = a1, sd = b1)
			#generates Normal2 whose contribution coefficient is 1-omega
			rightSet = rnorm(right.num, mean = a2, sd = b2)
			isNormal = TRUE;
		}
		else {
			#generates data based on the number of observations
			#found above and uses them for the gamma distribution functions

			#generates Gamma1 whose contribution coefficient is omega
			leftSet = rgamma(left.num, shape = a1, scale = b1)
			#generates Gamma2 whose contribution coefficient is 1-omega
			rightSet = rgamma(right.num, shape = a2, scale = b2)
		}

		#generates an array for the pooled data of Gamma1 and Gamma2
		result = vector (length = n);
		for(x in 1:right.num)
		{
			result[x] = rightSet[x];
		}

		if(omega != 0)
		{
			for (y in 1:left.num)
			{
				result[y+right.num] = leftSet[y];
			}
		}

		#gets a rough estimate of omega for the distribution, pOmega
		pOmega = left.num/n;

		#obtains the mean and variance of gamma1
		lVar = var(leftSet);
		lMean = mean(leftSet);

		#obtains the mean and variance of gamma2
		rVar = var(rightSet);
		rMean = mean(rightSet);

		#estimates the alphas and betas of the two functions, 
		#based on the method of moments.
		#Since these values will be optimized over, 
		#they will be overwritten
		pBeta1 = lVar/lMean;
		pAlpha1 = lMean/pBeta1;

		pBeta2 = rVar/rMean;
		pAlpha2 = rMean/pBeta2;

		#checks to see if the distribution type is
		#a normal, and if so, optimized using those mixNormal
		if(isNormal){
			runTest = optim(theta <- c(lMean,sqrt(lVar),rMean,sqrt(rVar),pOmega)
				,mixNormal, x = result)$par
		}
		#if not a normal distribution, checks is omega is 1 or 0
		#thus telling us it is a single gamma, and to use
		#the function for a single gamma
		else if(pOmega == 1)
		{
			runTest = optim(theta <- c(pAlpha1,pBeta1),singleGamma, x = result)$par
			runTest[5] = 1
			runTest[3] = 0
			runTest[4] = 0
		}
		else if(pOmega == 0)
		{
			runTest = optim(theta <- c(pAlpha2,pBeta2,pOmega),singleGamma, x = result)$par
			runTest[5] = 0
			runTest[3] = runTest[1]
			runTest[4] = runTest[2]
			runTest[1] = 0
			runTest[2] = 0
		}
		#if not normal or single gamma, runs the mix gamma
		#optimization
		else 
		{
			runTest = optim(theta <- c(pAlpha1,pBeta1,pAlpha2,pBeta2,pOmega),mixGamma, x = result)$par
		}  #runTest = (summary)$par
		param1[test] = runTest[1]
		param2[test] = runTest[2]
		param3[test] = runTest[3]
		param4[test] = runTest[4]
		param5[test] = runTest[5]
	}
	#returns out results as a table of values
	return(table <-list(alpha1=param1, beta1=param2, alpha2=param3, beta2=param4, omega=param5))
}


#function that graphs an estimate against the actual, given the parameters for each
graph <- function(omegaHat, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, omega, alpha1, beta1, alpha2, beta2,imgName)
{
	x <- seq(0, 30, 0.01)
	fit <- omegaHat * dgamma(x, shape = alpha1Hat,  scale = beta1Hat) + (1 - omegaHat) * dgamma(x, shape = alpha2Hat,  scale = beta2Hat)
	actual <- omega * dgamma(x, shape = alpha1,  scale = beta1) + (1 - omega) * dgamma(x, shape = alpha2,  scale = beta2)

	fstring = paste(imgName,".png",sep="")
	png(filename=fstring)
	plot(x,fit, type="l",col="blue")
	lines(x,actual, col="red ")
	legend('topright', c('estimated','true'), 
	   lty=1, col=c('red', 'blue', 'green',' brown'), bty='n', cex=.75)

	dev.off()
}

#graphs the mean of an estimated parameter against N
graphc <- function(N,estimate,r,imgName)
{
	fstring = paste(imgName,".png",sep="")
	png(filename=fstring)
	plot(N,estimate,col="blue")
	lines(N,r, col="red")
	lines(lowess(N,estimate), col="green")
	legend('topright', c('actual','estimated','fit'), 
	   lty=1, col=c('red', 'blue', 'green',' brown'), bty='n', cex=.75)

	dev.off()
}

#graphs the variance of an estimated parameter against N
graphv <- function(N,variance,imgName)
{
	fstring = paste(imgName,".png",sep="")
	png(filename=fstring)
	plot(N,variance,col="blue")
	lines(lowess(N,variance), col="red")
	legend('topright', c('fit','variance'), 
	   lty=1, col=c('red', 'blue', 'green',' brown'), bty='n', cex=.75)

	dev.off()
}

########################
#looking for convergence
########################

#initialize a bunch of vectors to store calculated values
alpha1s = vector (length = 100);
a1v = vector (length = 100);
beta1s = vector (length = 100);
b1v = vector (length = 100);
alpha2s = vector (length = 100);
a2v = vector (length = 100);
beta2s = vector (length = 100);
b2v = vector (length = 100);
omegas = vector (length = 100);
ov = vector (length = 100);
for (N in seq(1, 100, 1))
{
	summary = MLE(2,6,4,3,0.6,N,100,FALSE)
	alpha1Hat = mean(summary[,1])
	alpha1Var = var(summary[,1])

	beta1Hat = mean(summary[,2])
	beta1Var = var(summary[,3])

	alpha2Hat = mean(summary[,3])
	alpha2Var = var(summary[,3])

	beta2Hat = mean(summary[,4])
	beta2Var = var(summary[,4])

	omegaHat = mean(summary[,5])
	omegaVar = var(summary[,5])

	#store the values for mean and variance of each est. parameter
	alpha1s[N] = alpha1Hat
	a1v[N] = alpha1Var

	beta1s[N] = beta1Hat
	b1v[N] = beta1Var

	alpha2s[N] = alpha2Hat
	a2v[N] = alpha2Var

	beta2s[N] = beta2Hat
	b2v[N] = beta2Var

	omegas[N] = omegaHat
	ov[N] = omegaVar
	print(paste('iteration ',N))
}

#graph the mean of each parameter as N increases to look for convergence
graphc(seq(1, 100, 1), alpha1s, rep(2, 100),"a1_conv01")
graphc(seq(1, 100, 1), beta1s, rep(6, 100),"b1_conv01")
graphc(seq(1, 100, 1), alpha2s, rep(4, 100),"a2_conv01")
graphc(seq(1, 100, 1), beta2s, rep(3, 100),"b2_conv01")
graphc(seq(1, 100, 1), omegas, rep(0.6, 100),"omega_conv01")

#graph the variance of each parameter as N increases
graphv(seq(1, 100, 1), a1v,"a1_var01")
graphv(seq(1, 100, 1), b1v,"b1_var01")
graphv(seq(1, 100, 1), a2v,"a2_var01")
graphv(seq(1, 100, 1), b2v,"b2_var01")
graphv(seq(1, 100, 1), ov,"omega_var01")


#re-initialize vectors for new data
alpha1s = vector (length = 100);
a1v = vector (length = 100);
beta1s = vector (length = 100);
b1v = vector (length = 100);
alpha2s = vector (length = 100);
a2v = vector (length = 100);
beta2s = vector (length = 100);
b2v = vector (length = 100);
omegas = vector (length = 100);
ov = vector (length = 100);
for (N in seq(1, 100, 1))
{
	summary = MLE(1,3,2,4,0.6,N,100,FALSE)
	alpha1Hat = mean(summary[,1])
	alpha1Var = var(summary[,1])

	beta1Hat = mean(summary[,2])
	beta1Var = var(summary[,3])

	alpha2Hat = mean(summary[,3])
	alpha2Var = var(summary[,3])

	beta2Hat = mean(summary[,4])
	beta2Var = var(summary[,4])

	omegaHat = mean(summary[,5])
	omegaVar = var(summary[,5])

	#store the values for mean and variance of each est. parameter
	alpha1s[N] = alpha1Hat
	a1v[N] = alpha1Var

	beta1s[N] = beta1Hat
	b1v[N] = beta1Var

	alpha2s[N] = alpha2Hat
	a2v[N] = alpha2Var

	beta2s[N] = beta2Hat
	b2v[N] = beta2Var

	omegas[N] = omegaHat
	ov[N] = omegaVar
	print(paste('iteration ',N))
}

#graph the mean of each parameter as N increases to look for convergence
graphc(seq(1, 100, 1), alpha1s, rep(1, 100),"a1_conv02")
graphc(seq(1, 100, 1), beta1s, rep(3, 100),"b1_conv02")
graphc(seq(1, 100, 1), alpha2s, rep(2, 100),"a2_conv02")
graphc(seq(1, 100, 1), beta2s, rep(4, 100),"b2_conv02")
graphc(seq(1, 100, 1), omegas, rep(0.6, 100),"omega_conv02")

#graph the variance of each parameter as N increases
graphv(seq(1, 100, 1), a1v,"a1_var02")
graphv(seq(1, 100, 1), b1v,"b1_var02")
graphv(seq(1, 100, 1), a2v,"a2_var02")
graphv(seq(1, 100, 1), b2v,"b2_var02")
graphv(seq(1, 100, 1), ov,"omega_var02")


###########################
# any anomalies?
###########################

#compares difference between two numbers against a threshold
cmpval <- function(x,y,threshold=1) {
	return (abs(x - y) < threshold)
}

#loops through a large number of test values for each parameter
#    note: this take a very long time to run
for (a1 in seq(1, 10, 0.1))
{
	for (a2 in seq(1, 10, 0.1))
	{
		for (b1 in seq(1, 10, 0.1))
		{
			for (b2 in seq(1, 10, 0.1))
			{
				for (p in seq(0, 1, 0.1))
				{
					summary = MLE(a1,b1,a2,b2,p,15,100,FALSE)
					alpha1Hat = mean(summary[,1])
					beta1Hat = mean(summary[,2])
					alpha2Hat = mean(summary[,3])
					beta2Hat = mean(summary[,4])
					omegaHat = mean(summary[,5])

					# checks if any parameter is off by a certain threshold
					if (!(cmpval(alpha1Hat,a1) && cmpval(alpha2Hat,a2) 
						&& cmpval(beta1Hat,b1) && cmpval(beta2Hat,b2) 
						&& cmpval(omegaHat,p)))
					{
						#outputs information about the estimate and the actual
						# value if the threshold is exceeded
						print(paste('a1    est:',alpha1Hat,' real:',a1))
						print(paste('b1    est:',beta1Hat,' real:',b1))
						print(paste('a2    est:',alpha2Hat,' real:',a2))
						print(paste('b2    est:',beta2Hat,' real:',b2))
						print(paste('omega est:',omegaHat,' real:',p))
						print('')
					}
				}
			}
		}
	}
}





#############################
# testing for a specific case
#############################

x <- seq(0, 30, 0.01)

#weird case that isn't handled well
summary = MLE(2,5,3,3,0.6,10,100,FALSE)
alpha1Hat = mean(summary[,1])
beta1Hat = mean(summary[,2])
alpha2Hat = mean(summary[,3])
beta2Hat = mean(summary[,4])
omegaHat = mean(summary[,5])

print(alpha1Hat)
print(beta1Hat)
print(alpha2Hat)
print(beta2Hat)
print(omegaHat)
print("")


graph(omegaHat, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, 0.6, 2, 5, 3, 3,"MLE2-5-3-3-10")

#stores first estimated function
fit1 <- omegaHat * dgamma(x, shape = alpha1Hat,  scale = beta1Hat) + (1 - omegaHat) * dgamma(x, shape = alpha2Hat,  scale = beta2Hat)


summary = MLE(2,5,3,3,0.6,100,100,FALSE)
alpha1Hat = mean(summary[,1])
beta1Hat = mean(summary[,2])
alpha2Hat = mean(summary[,3])
beta2Hat = mean(summary[,4])
omegaHat = mean(summary[,5])

print(alpha1Hat)
print(beta1Hat)
print(alpha2Hat)
print(beta2Hat)
print(omegaHat)
print("")


graph(omegaHat, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, 0.6, 2, 5, 3, 3,"MLE2-5-3-3-100")

#stores second estimated function
fit2 <- omegaHat * dgamma(x, shape = alpha1Hat,  scale = beta1Hat) + (1 - omegaHat) * dgamma(x, shape = alpha2Hat,  scale = beta2Hat)


summary = MLE(2,5,3,3,0.6,1000,100,FALSE)
alpha1Hat = mean(summary[,1])
beta1Hat = mean(summary[,2])
alpha2Hat = mean(summary[,3])
beta2Hat = mean(summary[,4])
omegaHat = mean(summary[,5])

print(alpha1Hat)
print(beta1Hat)
print(alpha2Hat)
print(beta2Hat)
print(omegaHat)
print("")


graph(omegaHat, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, 0.6, 2, 5, 3, 3,"MLE2-5-3-3-1000")

#stores third estimated function
fit3 <- omegaHat * dgamma(x, shape = alpha1Hat,  scale = beta1Hat) + (1 - omegaHat) * dgamma(x, shape = alpha2Hat,  scale = beta2Hat)
actual <- 0.6 * dgamma(x, shape = 2,  scale = 5) + (1 - 0.6) * dgamma(x, shape = 3,  scale = 3)


#output plot comparing the 3 stored estimated functions against the actual
fstring = "comparison2-5-3-3.png"
png(filename=fstring)
plot(x,actual, type="l",col="red")
lines(x,fit1, type="l",col="blue")
lines(x,fit2, type="l",col="green")
lines(x,fit3, type="l", col="brown")
legend('topright', c('actual','MLE 10 iter.','MLE 100 iter.','MLE 1000 iter.'), 
   lty=1, col=c('red', 'blue', 'green',' brown'), bty='n', cex=.75)

dev.off()
