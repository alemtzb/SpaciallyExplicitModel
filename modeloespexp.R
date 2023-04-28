setwd("~/Desktop/SpExModel")

require(cubature)
require(rgl)

datdisp=read.csv("disp.csv") #dispersal kernel parameters
datheta=read.csv("Thetas.csv") #moodel dispersal
DD=read.csv("DDs.csv") #mortality parameter
BBs1=read.csv("BBs1.csv") #parameter a
BBs2=read.csv("BBs2.csv") #parameter b
BBs3=read.csv("BBs3.csv") #parametro c
alfas1=read.csv("alfas1.csv") #parameter g 
alfas2=read.csv("alfas2.csv") #parameter h
betas=read.csv("betas.csv")#facilitation parameter


#convert all data to matrices
BB1=as.matrix(BBs1[,2:ncol(BBs1)])
rownames(BB1)=BBs1[,1]
BB2=as.matrix(BBs2[,2:ncol(BBs2)])
rownames(BB2)=BBs2[,1]
BB3=as.matrix(BBs3[,2:ncol(BBs3)])
rownames(BB3)=BBs3[,1]
alphas1=as.matrix(alfas1[,2:ncol(alfas1)])
rownames(alphas1)=alfas1[,1]
alphas2=as.matrix(alfas2[,2:ncol(alfas2)])
rownames(alphas2)=alfas2[,1]
bbetas=as.matrix(betas[,2:ncol(betas)])
rownames(bbetas)=betas[,1]

#We take into account only the birth rate parameters of the last seven years
BB1=BB1[,8:14]
BB2=BB2[,8:14]
BB3=BB3[,8:14]
#Detransform parameters
BB2=1/(1+exp(-BB2))-.5
BB3=-1/(1+exp(-BB3))*5/1000
-999->alphas1[which(is.na(alphas1[,])=="TRUE")]
0->alphas2[which(is.na(alphas2[,])=="TRUE")]
bbetas=1/(1+exp(-bbetas))
0->bbetas[which(is.na(bbetas[,])=="TRUE")]
theta=exp(datheta[,2])

#To obtain total abundance of the species  for which we did not calculate pairwise interactions following Bayesian Variable Selection
spnum=ncol(alphas1) #to obtain the number of species in our study
tx0=matrix(ncol=ncol(alphas1),nrow=nrow(alphas1))
0->tx0[which(is.na(alphas1[,])=="FALSE")]
1->tx0[which(is.na(alphas1[,])=="TRUE")]
tx0=tx0[,-37]

#separate the parameters of the species for which we have abundance data
alpabu1=alphas1[1:33,1:37]
alpabu2=alphas2[1:33,1:37]
betabu=bbetas[1:33,1:37]
DDabu=DD[1:33,]
BB1abu=BB1[1:33,]
BB2abu=BB2[1:33,]
BB3abu=BB3[1:33,]
#separate the parameters of the species for which we have presence/absence data
alppa1=alphas1[34:36,1:37]
alppa2=alphas2[34:36,1:37]
betpa=bbetas[34:36,1:37]
DDpa=DD[34:36,]
BB1pa=BB1[34:36,]
BB2pa=BB2[34:36,]
BB3pa=BB3[34:36,]

#Function for species with abundance data to simulate the abundance of species in the next year
lam=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[1:33]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth)
	alphaspre=alphas[,1:36]%*%tx 
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
	txmas=log(tx+1)
	txmasother=log(txothers+1)
	betaspre=bbetas[,1:36]%*%txmas
	betasothers=bbetas[,37]*txmasother
	betassum=betaspre+betasothers
	fac=exp(betassum)
	new=BB*tx[1:33]/(1+alphasum)*fac 
	cbind(surv,new)
	#published function returns t2, a sum of surv and new instead of cbind object
	#there was also a difference between the columns that are taken into account to calculate alphas but I changed them
	#to match the published code and checked that we were taking the right columns and rows into account
}


#Function for species with presence absence daata to simulate if the species is present of absent in the next yearv
lampa=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[34:36]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth) 
	alphaspre=alphas[,1:36]%*%tx 
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
	txmas=log(tx+1)
	txmasother=log(txothers+1)
	betaspre=bbetas[,1:36]%*%txmas
	betasothers=bbetas[,37]*txmasother
	betassum=betaspre+betasothers
	fac=exp(betassum)
	new=BB*tx[34:36]/(1+alphasum)*fac
	cbind(surv,new)
	#published function returns t2, a sum of surv and new instead of cbind object
	#there was also a difference between the columns that are taken into account to calculate alphas but I changed them
	#to match the published code and checked that we were taking the right columns and rows into account
}

#here initial species abundances are assignes in published code

#Function that integrates both species with abuncance data and species with presence absence data into the same simulation
lamx=function(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,year,depth){
	txothers=as.vector(tx0%*%as.matrix(tx))
	t2abu=lam(DDabu,BB1abu,BB2abu,BB3abu,alpabu1,alpabu2,betabu,tx,txothers[1:33],year,depth)
	t2pa=lampa(DDpa,BB1pa,BB2pa,BB3pa,alppa1,alppa2,betpa,tx,txothers[34:36],year,depth)
	list(t2abu,t2pa)
	#here a list of data for presence/absence species and abundance speciesi is returned instead of an object that includes them both
	}





#End of population model****************************

#Computing of dispersal kernel
#probability density of the dispersal kernel centered at 0, x is a two dimensional vector with coordinates in dmc
	kerExp <- function(x,alfa,cc){
        #Define a radius
        r <- sqrt((x[1]**2)+(x[2]**2))
        # A normalizing term 
        nrm <- 2*pi*(alfa**2)*gamma(2/cc)/cc
		# Returns the value of f(r), where f is the exponential kernel with parameters alfa and cc
        return(exp(-((r/alfa)**cc))/nrm) # Returns the value of f(r), where f is the exponential kernel with parameters alfa and cc
        #return(kerExp(r,alfa,cc)) #***** Check this!!!!!!!!!!!!
} 

#**********Pregunta**************
#Returns the probability of dispersal of a seed to any 1 dm2 square in a dims * dims latice where the seed is produced in the central square
kerBidim=function(alfa,cc,dims){
	if(dims/2==floor(dims/2)) stop("dims must be odd") #why?
	lims=seq(-dims*0.5,dims*0.5,1)
	ncat=length(lims)-1
	mat=matrix(ncol=ncat,nrow=ncat)
	for(i in 1:ncat){
		for(j in 1:ncat){
			mat[i,j]=hcubature(kerExp,c(lims[i],lims[j]),c(lims[i+1],lims[j+1]),alfa=alfa,cc=cc)$integral
		}
	}
	mat
}

#makes a dims * dims lattice for every species and creates an array containing all latices in the 3rd dimension. 
#The array is rearranged so that the seeds are produced in the first cell [1,1] of each lattice assuming a toroidal configuration in space.
allker=function(alfa,cc,dims){
	nsp=length(alfa)
	arr=array(dim=c(dims,dims,nsp))
	for(i in 1:nsp) arr[,,i]=kerBidim(alfa[i],cc[i],dims)
	for(i in 1:nsp) arr[,,i]=arr[,,i]/sum(arr[,,i])
	arr[c((ceiling(dims/2)):dims,1:(ceiling(dims/2)-1)),c((ceiling(dims/2)):dims,1:(ceiling(dims/2)-1)),]
}



rearrange=function(mat,i,j){
	n=dim(mat)[1]
	if(i>1) mat=mat[c((n-i+2):n,1:(n-i+1)),,]
	if(j>1) mat=mat[,c((n-j+2):n,1:(n-j+1)),]
	mat
}


#End of dispersal kernel***********************

#simulates 2D population growth using dispersal kernel 
simu=function(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,dispker,theta,yr,depth){
	dims=dim(txini)[1] #txini defined bellow
	nsp=dim(txini)[3]
	txseed=array(0,dim=c(dims,dims,nsp))
	txabu=array(0,dim=c(dims,dims,nrow(alpabu1)))
	txpa=array(0,dim=c(dims,dims,nrow(alppa1)))
	t2=txseed
	for(i in 1:dims){
		for(j in 1:dims){
			dkij=rearrange(dispker,i,j) #dispker is an object created with the allker function (the one that makes the toroidal object). It uses the dispersal parameters.
		pret2=lamx(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini[i,j,],tx0,yr,depth) 
			seed=c(pret2[[1]][,2],pret2[[2]][,2]) #Here we use the seed production obtained using lamx (abundance data [[1]] and presence/absencce data [[2]])
			for(k in 1:nsp) txseed[,,k]=seed[k]*dkij[,,k]+txseed[,,k] #las semillas que se producieron se multiplican por el kernel de dispersión y se suma ¿el número de semillas del tiempo anterior?
			txabu[i,j,]=pret2[[1]][,1] #txabu es ahora el número de individuos que había en el tiempo anterior más las que llegan en preparación por el siguiente tiempo en la simulación
			txpa[i,j,]=pret2[[2]][,1] #mismo que arriba pero ahora para presencia/ausencia
		}
	}
	t2[,,1:nrow(alpabu1)]=txseed[,,1:nrow(alpabu1)]+txabu #suma de viejos individuos más nuevos para el tiempo 2? Abundance species
	t2[,,(nrow(alpabu1)+1):nsp]=1-(1-txpa)*exp(-txseed[,,(nrow(alpabu1)+1):nsp]) #suma de viejos individuos más nuevos para el tiempo 2? Presence/absence species
	for(k in 1:nrow(alpabu1)) t2[,,k]=matrix(rnbinom(dims^2,size=theta[k],mu=t2[,,k]),ncol=dims) #something related to model dispersion? 
	for(k in (nrow(alpabu1)+1):nsp) t2[,,k]=matrix(rbinom(dims^2,size=1,prob=t2[,,k]),ncol=dims)
	t2
}

#*******I still need to go over simut*******
simut=function(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,dispker,theta,depth,burn,iter){
	dims=dim(txini)[1]
	nsp=dim(txini)[3]
	nyr=dim(BB1abu)[2]
	sal=array(dim=c(dims,dims,nsp,iter))
	for(t in 1:burn){
		yr=ceiling(runif(1)*nyr)
		txini=simu(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,aa,theta,yr,depth)
	}
	for(t in 1:iter){
		yr=ceiling(runif(1)*nyr)
txini=simu(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,aa,theta,yr,depth)
sal[,,,t]=txini
	}
	sal
}


densobs=c(0.11,0.03,0.06,0.03,0.17,0.45,0.29,0.02,0.02,0.02,0.1,0.02,0.02,0.23,1.02,0.05,0.05,0.03,0.02,0.42,0.06,1.05,0.02,0.18,0.05,0.02,0.15,0.02,0.65,0.52,0.59,0.12,0.09,0.09,0.18,0.17,0.02)

densini=function(dims,theta,densobs){
	txini=txini=array(0,dim=c(dims,dims,36))
	densobs[which(densobs<0.1)]=0.5
	for(i in 1:33) txini[,,i]=rnbinom(dims^2,size=theta[i],mu=densobs[i])
	for(i in 34:36) txini[,,i]=rbinom(dims^2,size=1,prob=densobs[i])
	txini
}

densini2=function(dims,theta,densobs){
	txini=txini=array(0,dim=c(dims,dims,36))
	densobs[which(densobs<0.1)]=0.5
	for(i in 1:33) txini[,,i]=rpois(dims^2,2)
	for(i in 34:36) txini[,,i]=rbinom(dims^2,size=1,prob=0.5)
	txini
}
aa=allker(datdisp[,2],datdisp[,3],33)
txini=densini2(33,theta,densobs) 

bb=simu(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,aa,theta,3,20)
bb=simut(alpabu1,alpabu2,betabu,DDabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,txini,tx0,dispker,theta,20,1,100)

