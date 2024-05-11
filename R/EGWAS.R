
#' @title Ensemble-like Genome Wide Associated Analysis Strategy
#' @description a post GWAS strategy aims to control false positive rate inspired by ensemble learning method.
#'
#'
#' @param y vector of phenotype, NA is NOT allowed                               
#' @param geno numeric genotype, row for individual and column for SNP           
#' @param CV additional fixed effects or covariates                              
#' @param map map data, 3 columns for SNPs name, Chromosome, Position            
#' @param GWAS.results p values of GWAS methods, ncol equal to nmethod           
#' @param Bayesian.method column for Bayesian method in GWAS.results             
#' @param Bin.size size of bin in union-intersection step                        
#' @param permu.num number of permutation test                                   
#' @param p.threshold significant level, Bonferroni/PPI/Permutation              
#' @param Kinship genetic relationship matrix                                    
#' @param LD.threshold correlation threshold for cluster                         
#' @param ncpus number of cpus used for parallel
#' @param verbose whether to print detail.
#'
#' @return a list of pseudo QTNs detected by base methods and EGWAS
#' @export
#'
#' @examples
#' \donttest{
#' myEGWAS <- EGWAS(y=phe, geno=geno, map=map,CV=NULL,
#'   GWAS.results=GWAS,Bayesian.method=c(1),
#'   Bin.size=0, permu.num=1000,p.threshold=C(.05,.05,.05),
#'   Kinship=G,LD.threshold=0.7,ncpus=10,verbose=TRUE)
#' }
EGWAS <- function(y=NULL, geno=NULL, map=NULL,CV=NULL, GWAS.results=NULL,Bayesian.method=c(),
				Bin.size=0, permu.num=1000,p.threshold=c(0.05,0.05,0.05),Kinship=NULL,
				LD.threshold=0.7,ncpus=1,verbose=TRUE){
###Date: 2024/04/30
#-------------------------------------------------------------------------#
#--------- Ensemble-like Genome Wide Associated Analysis Strategy --------#
#-------------------------------------------------------------------------#
# Input:                                                                  #
# y: vector of phenotype, NA is NOT allowed                               #
# geno: numeric genotype, row for individual and column for SNP           #
# CV: additional fixed effects or covariates                              #
# map: map data, 3 columns for SNPs name, Chromosome, Position            #
# GWAS.results: p values of GWAS methods, ncol equal to nmethod           #
# Bayesian.method: column for Bayesian method in GWAS.results             #
# Bin.size: size of bin in union-intersection step                        #
# permu.num: number of permutation test                                   #
# p.threshold: significant level, Bonferroni/PPI/Permutation              #
# Kinship: genetic relationship matrix                                    #
# LD.threshold: correlation threshold for cluster                         #
# ncpus: number of cpus used for parallel                                 #
#-------------------------------------------------------------------------#

na.num <- length(which(is.na(y)))
if(na.num != 0) stop("NA is not allowed in phenotype!")
n <- nrow(geno); M<- ncol(geno)
if( (n != length(y)) | (M!= nrow(map)) ) stop("input not match!")
if(!is.null(CV)){
	CV = as.matrix(CV)
	cvna.num <- length(which(is.na(CV)))
	if(cvna.num != 0) stop("NA is not allowed in CV!")
	if( nrow(CV) != length(y)) stop("CV not match!")
}
if(  (M!= nrow(GWAS.results)) ) stop("GWAS results not match geno!")

time1 = as.numeric(Sys.time())
if(verbose){
	cat(c("==========================================================================","\n"))
	cat(c("=========","Ensemble-like Genome Wide Association Anlysis Strategy","=========","\n"))
	cat(c("==========================================================================","\n"))
	cat(c("Number of individual:",length(y),"\n"))
	cat(c("Number of Markers:",nrow(map),"\n"))
	cat(c("Number of base methods:",ncol(GWAS.results),"\n"))
	cat(c("Number of permutation:",permu.num,"\n"))
	cat(c("Significant level:",p.threshold,"\n"))
	cat(c("Number of CPUs:",ncpus,"\n"))
	cat(c("==========================================================================","\n"))
}
if(is.null(Kinship)){
	G=rMVP::MVP.K.VanRaden(bigmemory::as.big.matrix(t(as.matrix(geno))),verbose = FALSE)
}else{
	G=Kinship
}

bonferroni <- p.threshold[1]/nrow(map)
ppi <- p.threshold[2]
p.threshold<- p.threshold[3]
signal=list()
#####################
###union-intersect###
#####################
if(ncol(GWAS.results)==1){
	i=1
	if(i %in% Bayesian.method){
		sig.index = which((1-GWAS.results[,1])<=ppi)
	}else{
		sig.index = which(GWAS.results[,1]<=bonferroni)
	}
	if(length(sig.index)==0){
		stop("no signal exist!")
	}
	if(verbose){
	cat("only 1 method existed","\n")
	cat(length(sig.index),"top markers:","\n")
	cat( sig.index,"\n")
	}

	signal[[1]]<-sig.index
	un.in <- sig.index ###candidate snp for permutation
	specific <-c()
}else{
	for(i in 1:ncol(GWAS.results)){
		if(i %in% Bayesian.method){
			sig.i = which((1-GWAS.results[,i])<=ppi)
		}else{
			sig.i = which(GWAS.results[,i]<=bonferroni)
		}

		if(verbose) cat(length(sig.i),"top markers in method",colnames(GWAS.results)[i],"\n")
		signal[[i]]<-sig.i
	}
	###LD remove in each method
	if("linux"%in%unlist(strsplit(R.Version()$platform,"-"))){
		ncpus=ncpus
	}else{
		ncpus=1
	}
	signal.ld <- parallel::mclapply(1:length(signal), function(x){
		ld.index <- signal[[x]]
		if(length(ld.index)>=2){
			X1<-as.matrix(geno[, c(ld.index)],ncol=length(ld.index))
			colnames(X1)<-ld.index
			sigma = stats::cov(X1)
			sigma.distance = stats::as.dist(1 - abs(stats::cov2cor(sigma)))
			fit = stats::hclust(sigma.distance, method="single")
			clusters = stats::cutree(fit, h=1-LD.threshold)
			if(x %in% Bayesian.method){
				pvals.vector<- 1-GWAS.results[signal[[x]],x]
			}else{
				pvals.vector<- GWAS.results[signal[[x]],x]
			}

			names(pvals.vector)<- colnames(X1)
			top.selected = sapply(1:max(clusters), function(c) {
				cluster_elements = clusters==c
				top_within = which.min(pvals.vector[cluster_elements])
				if( length(top_within)==0 ) top_within = 1
				which(cluster_elements)[top_within]
			})
			ld.index<-as.numeric(names(top.selected))
			rm(sigma);rm(sigma.distance);rm(fit);rm(clusters);rm(pvals.vector);rm(top.selected)
		}
		return(ld.index)
	}, mc.cores=ncpus, mc.preschedule=length(signal)>ncpus,mc.allow.recursive=TRUE)
	un.in <- c()
	for(i in 1:(length(signal)-1)){
	   for(j in (i+1):length(signal)){
		   GWASi=signal.ld[[i]]
		   GWASj=signal.ld[[j]]
		   if((length(GWASi)>0)&(length(GWASj)>0)){
			   for(m in GWASi){
				   for(a in GWASj){
					  if((map[m,2]==map[a,2])&(abs(map[m,3]-map[a,3])<=Bin.size)) un.in=c(un.in,m,a)

					}
				}
			}
		}
	}
	un.in <- unique(un.in)
	whole <- unique(unlist(signal))
	if(verbose) cat(ncol(GWAS.results),"methods detected",length(whole),"different markers","\n")
	if(verbose) cat(length(un.in),"markers overlap","\n")
	if(verbose) cat(un.in,"\n")
	specific <- whole[-which(whole %in% un.in)]
	if(length(un.in)==0) {
		warning("no overlap result among methods, filter the whole signal")
		un.in <- whole
	}
	if(length(un.in)==0) stop("no method has signal!")
}
###############
###LD remove###
###############

if(length(un.in)<2){
	ld.index <- un.in
	if(verbose) cat("LD remove is not needed","\n")
}else{
	if(verbose) cat("LD remove is working...")
	ld.index <- un.in
	X1<-as.matrix(geno[, c(ld.index)],ncol=length(ld.index))
	colnames(X1)<-ld.index
	sigma = stats::cov(X1)
	sigma.distance = stats::as.dist(1 - abs(stats::cov2cor(sigma)))
	fit = stats::hclust(sigma.distance, method="single")
	clusters = stats::cutree(fit, h=1-LD.threshold)

	x<-as.matrix(cbind(X1,rep(1, length(y)),CV), nrow = length(y))
	fit.mv<-gaston::lmm.aireml(y, X=x, K=list(G), EMsteps_fail= 10L, verbose= FALSE)
	V <- G*fit.mv$tau[1]+diag(fit.mv$sigma2,dim(G)[1])
	Vi <- try(solve(V+diag(1,ncol(V))*(1e-10)),silent=T)
    if(inherits(Vi, "try-error")){
        warning("Singular matrix V!")
        Vi <- MASS::ginv(V)
    }
	XVX=t(x)%*%Vi%*%x
    XVXi <- try(solve(XVX+diag(1,ncol(XVX))*(1e-10)),silent=T)
    if(inherits(XVXi, "try-error")){
        warning("Singular matrix XVX!")
        XVXi <- MASS::ginv(XVX)
    }
    se=sqrt(diag(XVXi))
	SNPX_ad<-fit.mv$BLUP_beta/se
	P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(X1)-1,lower.tail=FALSE)
	#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
	P_X<-P_X[1:ncol(X1)]
	pvals.vector<- P_X
	names(pvals.vector)<- colnames(X1)
	top.selected = sapply(1:max(clusters), function(c) {
		cluster_elements = clusters==c
		top_within = which.min(pvals.vector[cluster_elements])
		if( length(top_within)==0 ) top_within = 1
		which(cluster_elements)[top_within]
	})
	ld.index<-as.numeric(names(top.selected))
	rm(sigma);rm(sigma.distance);rm(fit);rm(clusters);rm(pvals.vector);rm(top.selected)
	if(verbose) cat(c("(DONE)","\n"))
	if(verbose) cat(length(ld.index),"markers retained after LD filter","\n")
	if(verbose) cat(ld.index,"\n")
}

###true p values###
X1<-as.matrix(geno[, c(ld.index)],ncol=length(ld.index))
x<-as.matrix(cbind(X1,rep(1, length(y)),CV), nrow = length(y))
fit.mv<-gaston::lmm.aireml(y, X=x, K=list(G), EMsteps_fail= 10L, verbose= FALSE)
V <- G*fit.mv$tau[1]+diag(fit.mv$sigma2,dim(G)[1])
Vi <- try(solve(V+diag(1,ncol(V))*(1e-10)),silent=T)
if(inherits(Vi, "try-error")){
    warning("Singular matrix V!")
    Vi <- MASS::ginv(V)
}
XVX=t(x)%*%Vi%*%x
XVXi <- try(solve(XVX+diag(1,ncol(XVX))*(1e-10)),silent=T)
if(inherits(XVXi, "try-error")){
    warning("Singular matrix XVX!")
    XVXi <- MASS::ginv(XVX)
}
se=sqrt(diag(XVXi))
SNPX_ad<-fit.mv$BLUP_beta/se
P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(X1)-1,lower.tail=FALSE)
#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
P_true<-P_X[1:ncol(X1)]
gc()
######################
###permutation test###
######################
if(verbose) cat("permutation test (",permu.num,"repeats) is working...","\n")
wh.index <- c(ld.index,specific)
X2<-as.matrix(geno[, c(wh.index)],ncol=length(wh.index))
xzgeno<-t(X2)
method=c("GLM")
if(ncol(X2)==1){
	method=c('MLM')}
genotype<-bigmemory::as.big.matrix(xzgeno)
mapn<-map[wh.index,]
if(method=='GLM') K=NULL
if(method=='MLM') K=G
if("linux"%in%unlist(strsplit(R.Version()$platform,"-"))){
	ncpus=ncpus
}else{
	ncpus=1
}
if(verbose & ncpus>1) cat(c("Fork process in Linux","\n"))
p_matrix <- parallel::mclapply(1:permu.num, function(x){
		index = 1:length(y)
		index.shuffle = sample(index,length(index),replace=FALSE)
		myY.shuffle = data.frame(ID=paste0("ID",1:length(y)),phe=y)
		myY.shuffle[,2] = myY.shuffle[index.shuffle,2]

		imMVP <- rMVP::MVP(
			phe=myY.shuffle,
			geno=genotype,
			map=mapn,
			K=K,
			CV.GLM=CV,
			CV.MLM=CV,
			#CV.FarmCPU=Covariates,
			#nPC.GLM=5,
			nPC.MLM=NULL,
			#nPC.FarmCPU=3,
			#perc=1,
			priority="speed",
			#vc.method="EMMA",
			#maxLoop=10,
			#method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
			file.output=FALSE,
			ncpus =1,
			verbose=FALSE,
			method=method
		)
		unit.pern=floor(permu.num/10)
		if(verbose & x%%unit.pern==0) cat(c("=",""))
		if(method=='GLM') return( imMVP$glm.results[,3])
		if(method=='MLM') return( imMVP$mlm.results[,3])

	}, mc.cores=ncpus, mc.preschedule=permu.num>ncpus,mc.allow.recursive=TRUE)
gc()

p_matrix<-unlist(p_matrix)
p_matrix<-matrix(p_matrix,nrow=permu.num,byrow=TRUE)
colnames(p_matrix) <- wh.index

if(length(P_true)>1){
	p_matrix1<-as.matrix(rbind(P_true, p_matrix[,c(1:length(ld.index))]),ncol=length(P_true))
}else{
	p_matrix1<-as.matrix(c(P_true, p_matrix[,c(1:length(ld.index))]),ncol=length(P_true))
}
per.p<-apply(p_matrix1,2,function(x){( length(which(x<=x[1] ) ) -1 )/permu.num})
per.index<-ld.index[which(per.p<=p.threshold)]

# minp<-apply(p_matrix[,c(1:length(ld.index))],1,min)
# pthre<-quantile(minp,p.threshold)
# per.index<-ld.index[which(P_true<=pthre)]

if(verbose) cat(c("(DONE)","\n"))
if(verbose) cat(length(per.index),"markers retained after permutation","\n")
if(verbose) cat(per.index,"\n")

if(length(per.index)==0) {
	stop("Permutation remove all overlap markers!")
}

if(length(specific)==0){
	egwas.index<-per.index
	egwas.index<-sort(egwas.index)
	names(egwas.index)=paste0("Union",1:length(egwas.index))
}else{
	if(verbose) cat(c("==========================================================================","\n"))
	if(verbose) cat(length(specific),"specific markers exist, start filtering...","\n")

	X3<-as.matrix(geno[, c(specific)],ncol=length(specific))
	per.x <- as.matrix(geno[, c(per.index)],ncol=length(per.index))
	CVmlm<-as.matrix(cbind(per.x,CV), nrow = length(y))
	x<-as.matrix(cbind(rep(1, length(y)),CVmlm), nrow = length(y))
	fit.mv<-gaston::lmm.aireml(y, X=x, K=list(G), EMsteps_fail= 10L, verbose= FALSE)
	V <- G*fit.mv$tau[1]+diag(fit.mv$sigma2,dim(G)[1])
	Vi <- try(solve(V+diag(1,ncol(V))*(1e-10)),silent=T)
	if(inherits(Vi, "try-error")){
		warning("Singular matrix V!")
		Vi <- MASS::ginv(V)
	}
	x<-as.matrix(cbind(X3,rep(1, length(y)),CVmlm), nrow = length(y))
	XVX=t(x)%*%Vi%*%x
	XVXi <- try(solve(XVX+diag(1,ncol(XVX))*(1e-10)),silent=T)
	if(inherits(XVXi, "try-error")){
		warning("Singular matrix XVX!")
		XVXi <- MASS::ginv(XVX)
	}
	se=sqrt(diag(XVXi))
	SNPX_ad<-fit.mv$BLUP_beta/se
	P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(X1)-1,lower.tail=FALSE)
	#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
	spe.true<-P_X[1:ncol(X3)]

	gc()
	if(length(specific)>1){
		if(verbose) cat("LD remove is working...")
		ld.index <- specific
		X1<-as.matrix(geno[, c(ld.index)],ncol=length(ld.index))
		colnames(X1)<-ld.index
		sigma = stats::cov(X1)
		sigma.distance = stats::as.dist(1 - abs(stats::cov2cor(sigma)))
		fit = stats::hclust(sigma.distance, method="single")
		clusters = stats::cutree(fit, h=1-LD.threshold)
		pvals.vector<- spe.true
		names(pvals.vector)<- colnames(X1)
		top.selected = sapply(1:max(clusters), function(c) {
			cluster_elements = clusters==c
			top_within = which.min(pvals.vector[cluster_elements])
			if( length(top_within)==0 ) top_within = 1
			which(cluster_elements)[top_within]
		})
		speld.index<-as.numeric(names(top.selected))
		rm(sigma);rm(sigma.distance);rm(fit);rm(clusters);rm(pvals.vector);rm(top.selected)
		if(verbose) cat(c("(DONE)","\n"))
		if(verbose) cat(length(speld.index),"specific markers retained after LD filter","\n")
		if(verbose) cat(speld.index,"\n")

		ld.index.s<-specific[which(specific %in% speld.index)]
		###true p values###
		X1<-as.matrix(geno[, c(ld.index.s)],ncol=length(ld.index.s))
		x<-as.matrix(cbind(X1,rep(1, length(y)),CVmlm), nrow = length(y))
		fit.mv<-gaston::lmm.aireml(y, X=x, K=list(G), EMsteps_fail= 10L, verbose= FALSE)
		V <- G*fit.mv$tau[1]+diag(fit.mv$sigma2,dim(G)[1])
		Vi <- try(solve(V+diag(1,ncol(V))*(1e-10)),silent=T)
		if(inherits(Vi, "try-error")){
			warning("Singular matrix V!")
			Vi <- MASS::ginv(V)
		}
		XVX=t(x)%*%Vi%*%x
		XVXi <- try(solve(XVX+diag(1,ncol(XVX))*(1e-10)),silent=T)
		if(inherits(XVXi, "try-error")){
			warning("Singular matrix XVX!")
			XVXi <- MASS::ginv(XVX)
		}
		se=sqrt(diag(XVXi))
		SNPX_ad<-fit.mv$BLUP_beta/se
		P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(X1)-1,lower.tail=FALSE)
		#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
		spe_true<-P_X[1:ncol(X1)]
	}else{
		spe_true=spe.true
		speld.index=specific
		ld.index.s=specific[which(specific %in% speld.index)]
	}
	if(length(spe_true)>1){
		p_matrix2<-as.matrix(rbind(spe_true,
					p_matrix[,which(wh.index %in% speld.index)]),ncol=length(spe_true))
	}else{
		p_matrix2<-as.matrix(c(spe_true,
					p_matrix[,which(wh.index %in% speld.index)]),ncol=length(spe_true))
	}
	spe.p<-apply(p_matrix2,2,function(x){( length(which(x<=x[1] ) ) -1 )/permu.num})
	spe.index<-ld.index.s[which(spe.p<=p.threshold)]

	# minp<-apply(p_matrix[,which(wh.index %in% speld.index)],1,min)
	# ppthre<-quantile(minp,p.threshold)
	# spe.index1<-specific[which(specific %in% speld.index)]
	# spe.true1<-spe.true[which(specific %in% speld.index)]
	# spe.index<-spe.index1[which(spe.true1<=ppthre)]

	if(verbose) cat(length(spe.index),"specific markers retained after permutation with overlap fixed","\n")
	if(verbose) cat(spe.index,"\n")

	egwas.index<-c(per.index,spe.index)
	seqQTN<-egwas.index
	if (verbose) cat("Total",length(seqQTN),"markers was detected", "\n")
	if(length(seqQTN)>=2){

		if (verbose) cat("Stepwise variable selected is working....")

		stepdata <- data.frame(Y=y, geno[,seqQTN])
		colnames(stepdata) <- c("Y", as.character(paste0("M",seqQTN)))
        ##diff version stepwise
		reg <- StepReg::stepwise(formula=Y~.,
								data=stepdata,
								selection="bidirection",
								select="BIC",
								sle = p.threshold,
								sls = p.threshold)
		gc()
		name <- as.data.frame(reg[[5]][[1]][-1,])
		name<-subset(name,name[,5]<=p.threshold)[,1]
		#name<-subset(name,name[,5]<=p.threshold/nrow(name))[,1]
		egwas.index<-as.numeric(gsub("M","",name))
		if(length(egwas.index)==0) {
			stop("StepReg remove all pseudo QTNs!")
		}
		egwas.index<-sort(egwas.index)
		ui=0;si=0
		for(i in 1:length(egwas.index)){
			if(egwas.index[[i]]%in%per.index){
				ui=ui+1
				names(egwas.index)[i]=paste0("Union",ui)
			}
			if(egwas.index[[i]]%in%spe.index){
				si=si+1
				names(egwas.index)[i]=paste0("Specific",si)
			}
		}
		if(verbose) cat(c("(DONE)","\n"))
		if(verbose) cat(length(egwas.index),"markers retained after StepReg","\n")
		if(verbose) cat(egwas.index,"\n")
	}
}
signal[[ (length(signal)+1) ]]<-egwas.index
names(signal)=c(colnames(GWAS.results),"EGWAS")

time2 = as.numeric(Sys.time())
time.cal <- round(time2-time1)
times <- function(x){
			h <- x %/% 3600
			m <- (x %% 3600) %/% 60
			s <- ((x %% 3600) %% 60)
			index <- which(c(h, m, s) != 0)
			num <- c(h, m, s)[index]
			char <- c("h", "m", "s")[index]
			return(paste(round(num), char, sep="", collapse=""))
}
if(verbose){
	cat(c("=========================","Analysis  Accomplished","=========================","\n"))
	cat(c("E-GWAS is DONE within total run time: ", times(time.cal), "\n"))
}
return(signal)
}

#' Title
#' @description MLM scan based on pseudo QTNs (as fix effects) from EGWAS
#'
#'
#' @param y vector of phenotype, NA is NOT allowed
#' @param geno numeric genotype, row for individual and column for SNP
#' @param map map data, 3 columns for SNPs name, Chromosome, Position
#' @param CV additional fixed effects or covariates
#' @param EGWAS.QTN index of pseudo QTNs detected by E-GWAS
#' @param Kinship genetic relationship matrix
#' @param ncpus number of cpus used for parallel
#'
#' @return a data.frame contain map, beta, se, pvalues
#' @export
#'
#' @examples
#' \donttest{
#' myscan <- EGWAS.SCAN(y=phe, geno=geno, map=map,CV=NULL,
#'   EGWAS.QTN=myEGWAS$EGWAS,Kinship=G,ncpus=10)
#' }
EGWAS.SCAN <- function(y=NULL, geno=NULL, map=NULL,CV=NULL,EGWAS.QTN=NULL,Kinship=NULL,
						ncpus=1){

###Date: 2024/04/30
#-------------------------------------------------------------------------#
#--------- Ensemble-like Genome Wide Associated Anlysis Strategy ---------#
#-------------------------------------------------------------------------#
# Input:                                                                  #
# y: vector of phenotype, NA is NOT allowed                               #
# geno: numeric genotype, row for individual and column for SNP           #
# CV: additional fixed effects or covariates                              #
# map: map data, 3 columns for SNPs name, Chromosome, Position            #
# EGWAS.QTN: index of pseudo QTNs detected by E-GWAS                      #
# Kinship: genetic relationship matrix                                    #
# ncpus: number of cpus used for parallel                                 #
#-------------------------------------------------------------------------#

na.num <- length(which(is.na(y)))
if(na.num != 0) stop("NA is not allowed in phenotype!")
n <- nrow(geno); M<- ncol(geno)
if( (n != length(y)) | (M!= nrow(map)) ) stop("input not match!")
if(!is.null(CV)){
	CV = as.matrix(CV)
	cvna.num <- length(which(is.na(CV)))
	if(cvna.num != 0) stop("NA is not allowed in CV!")
	if( nrow(CV) != length(y)) stop("CV not match!")
}

if(is.null(Kinship)){
	G=rMVP::MVP.K.VanRaden(bigmemory::as.big.matrix(t(as.matrix(geno))),verbose = FALSE)
}else{
	G=Kinship
}

####pseudo QTN significant level
QTN <- EGWAS.QTN
X1<-as.matrix(geno[, c(QTN)],ncol=length(QTN))
x<-as.matrix(cbind(X1,rep(1, length(y)),CV), nrow = length(y))
fit.mv<-gaston::lmm.aireml(y, X=x, K=list(G), EMsteps_fail= 10L, verbose= FALSE)
V <- G*fit.mv$tau[1]+diag(fit.mv$sigma2,dim(G)[1])
Vi <- try(solve(V+diag(1,ncol(V))*(1e-10)),silent=T)
if(inherits(Vi, "try-error")){
	warning("Singular matrix V!")
	Vi <- MASS::ginv(V)
}
XVX=t(x)%*%Vi%*%x
XVXi <- try(solve(XVX+diag(1,ncol(XVX))*(1e-10)),silent=T)
if(inherits(XVXi, "try-error")){
	warning("Singular matrix XVX!")
	XVXi <- MASS::ginv(XVX)
}
se=sqrt(diag(XVXi))
SNPX_ad<-fit.mv$BLUP_beta/se
P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(X1)-1,lower.tail=FALSE)
QTN.P<-P_X[1:ncol(X1)]
Qmap <- data.frame(map[QTN,],Beta=fit.mv$BLUP_beta[1:ncol(X1)],SE=se[1:ncol(X1)],P.value=QTN.P)

CVmlm<-as.matrix(cbind(X1,CV), nrow = length(y))
X2<-as.matrix(geno[, -QTN])
xzgeno<-t(X2)
genotype<-bigmemory::as.big.matrix(xzgeno)
mapn<-map[-QTN,]
myY<-data.frame(ID=paste0("ID",1:length(y)),trait=y)

imMVP <- rMVP::MVP(
			phe=myY,
			geno=genotype,
			map=mapn,
			K=G,
			CV.MLM=CVmlm,
			#CV.FarmCPU=Covariates,
			#nPC.GLM=5,
			nPC.MLM=NULL,
			#nPC.FarmCPU=3,
			#perc=1,
			priority="speed",
			#vc.method="EMMA",
			#maxLoop=10,
			#method.bin="FaST-LMM",#"FaST-LMM","EMMA", "static"
			file.output=FALSE,
			ncpus =ncpus,
			verbose=FALSE,
			method="MLM"
		)
mvp.result <- data.frame(mapn,imMVP$mlm.results)
colnames(mvp.result) <- colnames(Qmap)
result <- rbind(mvp.result,Qmap)

result <- result[order(result[,2],result[,3]),]
return(result)
}
