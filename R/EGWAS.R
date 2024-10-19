#' @title Ensemble-like Genome Wide Associated Analysis Strategy
#' @description a post GWAS strategy aims to control false positive rate inspired by ensemble learning method.
#'
#'
#' @param y vector of phenotype, NA is NOT allowed                               
#' @param geno numeric genotype, row for MARKERS and column for IND           
#' @param CV additional fixed effects or covariates                              
#' @param map map data, 3 columns for SNPs name, Chromosome, Position            
#' @param pvalue p values of GWAS methods, ncol equal to nmethod           
#' @param pip posterior inclusive probability (PIP) derived from Bayes method             
#' @param bin.size size of bin in union-intersection step                        
#' @param permu.num number of permutation test                                   
#' @param pv.threshold significant level, Bonferroni   
#' @param pip.threshold significant level, PIP   
#' @param permu.threshold significant level, Permutation              
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
EGWAS <- function(y=NULL, geno=NULL, map=NULL,CV=NULL, pvalue=NULL, pip=NULL,
				  bin.size=0, permu.num=1000,pv.threshold=0.05,pip.threshold=0.05,
				  permu.threshold=0.05,Kinship=NULL,LD.threshold=0.7,
				  ncpus=1,verbose=TRUE){
###Date: 2024/04/30

EGWAS.version <-  function(verbose = TRUE) {
  if (verbose) {
    cat(paste("+",
              paste(rep("=", 18), collapse = ""),
              "< Welcome to E-GWAS >",
              paste(rep("=", 18), collapse = ""),
              "+", sep = ""), "\n")
    cat("|            A Ensemble-like Strategy for GWAS            |\n")
	cat("|         _ _ _  _ _ _ _    __    __  _ _      _ _        |\n")
	cat("|        / _ _/ / _ _ _/   / /   / / /   |    / _ /       |\n")
	cat("|       / /_   / /  __    / /   / / / /| |   / /          |\n")
	cat("|      / _ _/ / /  /_ \\  / /_  / / / /_| |   \\ \\          |\n")
	cat("|     / /_ _  \\ \\_ _/ / / /  |/ / /  __  |  _ / /         |\n")
	cat("|    /_ _ _/   \\_ _ _/ /_ _/| _/ /_/   |_| /_ _/          |\n")   
    cat("|                                                         |\n")
    cat("|    E-mail: fjxu@webmail.hzau.edu.cn                     |\n")
    cat(paste("+", paste(rep("=", 57), collapse = ""), "+", sep = ""), "\n")
  }
}
EGWAS.version(verbose = verbose)

time.trans <- function(x) {
  h <- x %/% 3600
  m <- (x %% 3600) %/% 60
  s <- ((x %% 3600) %% 60)
  index <- which(c(h, m, s) != 0)
  num <- c(h, m, s)[index]
  char <- c("h", "m", "s")[index]
  return(paste(round(num), char, sep = "", collapse = ""))
}

CovMatrix <- function(CV = NULL, n = NULL) {
  Cov <- matrix(1, n)
  if (!is.null(CV)) {
    if (sum(is.na(CV)) != 0) stop("'NA' isn't allowed in CV!")
    for (i in 1:ncol(CV)){
      if (is.numeric(CV[, i])) {
        #quantitative
        Cov <- cbind(Cov, CV[, i])
      }else {
        #discrete
        if (length(unique(CV[, i])) == 1) {
          stop("No groups in column ", i, " of CV!")
        }
        ci <- nnet::class.ind(CV[, i])
        Cov <- cbind(Cov, ci[, -1])
      }
    }
  }
  return(Cov)
}

LD.remove <- function(index = NULL, geno = NULL, value = NULL,
                      LD.threshold = 0.7, verbose = TRUE) {
  if (length(index) < 2) {
    return(index)
  }else {
    ldr <- function(index, value) {
      X1 <- geno[index, ]
      X1<-sweep(X1[, ], 1, rowMeans(X1))
      rownames(X1) <- index
      sigma <- tcrossprod(X1) / (dim(geno)[2] - 1)
      sigma.distance <- stats::as.dist(1 - abs(stats::cov2cor(sigma)))
      fit <- stats::hclust(sigma.distance, method = "single")
      clusters <- stats::cutree(fit, h = 1 - LD.threshold)
      names(value) <- rownames(X1)
      top.selected <- sapply(1:max(clusters), function(i) {
	cluster_elements <- clusters == i
	top_within <- which.min(value[cluster_elements])
	if (length(top_within) == 0) top_within <- 1
	return(which(cluster_elements)[top_within])
      })
      index<-as.numeric(names(top.selected))
      rm(sigma, sigma.distance, fit, clusters, top.selected)
      return(index)
    }
    if (length(index) <= 10000) {
      if (verbose) {
        cat("LD removing...\n")
        pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      index <- ldr(index = index, value = value)
      if (verbose) {
        pbapply::setTimerProgressBar(pb, 1)
      }
    }else {
      nbreak <- ceiling(length(index) / 10000)
      cutind <- cut(1:length(index), breaks = nbreak, labels = FALSE)
      if (verbose) {
        cat("LD removing...\n")
        pb <- pbapply::timerProgressBar(max = nbreak + 1, width = 30, char = "-", style = 3)
        on.exit(close(pb))
      }
      indi <- lapply(1:nbreak, function(i) {
        idi <- index[which(cutind == i, arr.ind = TRUE)]
        vi <- value[which(cutind == i, arr.ind = TRUE)]
	if (verbose) {
          pbapply::setTimerProgressBar(pb, i)
        }
        return(ldr(index = idi, value = vi))
      })
      indi <- unlist(indi)
      vali <- value[which(index %in% indi)]
      index <- ldr(index = indi, value = vali)
      if (verbose) {
        pbapply::setTimerProgressBar(pb, nbreak + 1)
      }
    }	
    return(index)
  }
}

Kinship.Van <- function(geno = NULL, weight = NULL, ncut = 1, ncpus = 1) {

  #---------------------------------------------------------#
  #Calculate the genomic Kinship by VanRaden's method
  if (!bigmemory::is.big.matrix(geno) && !is.matrix(geno)) {
    stop("Wrong format of geno inputed! matrix/big.matrix is allowed.")
  }
  if (sum(is.na(geno[, ])) != 0) stop("'NA' isn't allowed in geno!")
  if (!is.null(weight)) {
    if (sum(is.na(weight)) != 0) stop("'NA' isn't allowed in weight!")
  }
  if (ncut < 1 || ncut %% 1 != 0) stop("ncut should be a positive integer!")

  n <- dim(geno)[2] #NUM OF INDIVIDUALS
  m <- dim(geno)[1] #NUM OF MARKERS
  if (m == 1) {
    maf <- mean(geno[, ]) / 2
  }else {
    maf <- rowMeans(geno[, ]) / 2
  }
  Z <- bigmemory::big.matrix(m, n, type = "double",
                             init = NULL, shared = FALSE)
  if (!is.null(weight)) {
    Z[, ] <- (geno[, ] - 2 * maf) * sqrt(weight)
  }else {
    Z[, ] <- geno[, ] - 2 * maf
  }
  rm(geno, weight)
  K <- matrix(NA, n, n)
  if (ncut == 1 || ncpus == 1) {
    ##recommend##
    if (m == 1) {
      K[, ] <- tcrossprod(Z[, ]) / (2 * sum(maf * (1 - maf)))
    }else {
      K[, ] <- crossprod(Z[, ]) / (2 * sum(maf * (1 - maf)))
    }
    rm(Z)
  }else {
    cut.lab <- cut(1:n, breaks = ncut, labels = FALSE)
    cut.list <- list()
    for (i in 1:ncut) {
      for (j in i:ncut) {
        cuti <- which(cut.lab == i)
        cutj <- which(cut.lab == j)
        cut.list <- c(cut.list, list(list(cuti, cutj)))
      }
    }
    rm(cut.lab, cuti, cutj)
    block.prod <- function(i) {
      row <- cut.list[[i]][[1]]
      col <- cut.list[[i]][[2]]
      block <- crossprod(Z[, row], Z[, col])
      return(block)
    }
    if (Sys.info()[["sysname"]] == "Linux") {
      bloc <- parallel::mclapply(1:length(cut.list), FUN = block.prod,
                                 mc.cores = ncpus,
                                 mc.preschedule = FALSE)
    }else {
      suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
      Z <- Z[, ]
      suppressMessages(snowfall::sfExport("Z", "cut.list"))
      bloc <- snowfall::sfLapply(1:length(cut.list), fun = block.prod)
      suppressMessages(snowfall::sfStop())
    }
    rm(Z)
    for (i in 1:length(cut.list)) {
      row <- cut.list[[i]][[1]]
      col <- cut.list[[i]][[2]]
      K[row, col] <- bloc[[i]]
      if (!identical(row, col)) {
        K[col, row] <- t(bloc[[i]])
      }
    }
    K[, ] <- K[, ] / (2 * sum(maf * (1 - maf)))
  }
  return(K)
}

if (sum(is.na(y)) != 0) stop("'NA' is not allowed in y!")
n <- ncol(geno); m<- nrow(geno)
if (n != length(y)) stop("Number of individuals does not match between y and geno")
if (m != nrow(map)) stop("Number of markers does not match between map and geno")
if (is.null(pvalue) & is.null(pip)) stop("P value should be provided")
nm <- ifelse(is.null(pvalue),0,ncol(pvalue)) + ifelse(is.null(pip),0,ncol(pip))
if (nm < 2) stop("At least 2 method results needed")

time1 = as.numeric(Sys.time())
if(verbose){
	cat(c("Number of individual:",n,"\n"))
	cat(c("Number of Markers:",m,"\n"))
	cat(c("Number of GWAS methods:",nm,"\n"))
	cat(c("Number of permutation:",permu.num,"\n"))
	cat(c("Number of covariates:",ifelse(is.null(CV),0,ncol(CV)),"\n"))
	cat(c("Number of threads:",ncpus,"\n"))
}

#####################
###union-intersect###
#####################
if (verbose) {
    cat(paste(paste(rep("-", 19), collapse = ""),
              "Union-Intersect  Step",
              paste(rep("-", 19), collapse = ""), sep = ""), "\n")
}
###pick up signal from gwas results
bonfer <- pv.threshold / m
GWAS <- cbind(pvalue, pip)
nma <- ifelse(is.null(pvalue),0,ncol(pvalue))
signal <- lapply(1:ncol(GWAS), function(i){
	if (i <= nma) {
		sig.i <- which(GWAS[,i] <= bonfer)
	}else {
		sig.i <- which((1 - GWAS[,i]) <= pip.threshold)
	}
	return(sig.i)
})
if (length(unlist(signal)) == 0) stop("No method has significant markers")
if (verbose) {
	cat(sapply(signal, length), "significant markers in the gwas methods", "\n")
}
###union-intersect
un.in <- c()
for(i in 1:(length(signal)-1)){
	for(j in (i+1):length(signal)){
		GWASi<-signal[[i]]
		GWASj<-signal[[j]]
		if((length(GWASi)>0)&(length(GWASj)>0)){
			for(m in GWASi){
				for(a in GWASj){
					if((map[m,2]==map[a,2])&(abs(map[m,3]-map[a,3])<=bin.size)) un.in<-c(un.in,m,a)
				}
			}
		}
	}
}
un.in <- unique(un.in)
whole <- unique(unlist(signal))
if(length(un.in)==0) {
	warning("No overlaps among gwas results, filter the whole signal")
	un.in <- whole
	specific <- NULL
}else {
	if(verbose) cat(length(un.in),"overlaps in",nm,"method results","\n")
	specific <- whole[-which(whole %in% un.in)]
}
###mlm
if (verbose) {
    cat("Re-calculating p values...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
}
if(is.null(Kinship)){
	G <- Kinship.Van(geno = geno)
}else{
	G <- Kinship
}
Cov <- CovMatrix(CV , n)
fit.mlm <- gaston::lmm.aireml(y,X=Cov,K=G,EMsteps_fail = 10L,eps = 1e-02,verbose = FALSE)
vg <-fit.mlm$tau[1]
ve <- fit.mlm$sigma2
V <- G*vg + diag(ve, dim(G)[1])
Vi <- try(solve(V + diag(1,ncol(V))*(1e-10)),silent=T)
if(inherits(Vi, "try-error")){
    warning("Singular matrix V!")
    Vi <- MASS::ginv(V)
}
mlm <- function(i){
	indexi <- un.in[i]
	x <- geno[indexi, ]
	xi <- cbind(x,Cov)
    XVX <- crossprod(xi, crossprod(Vi, xi))
    XVXi <- try(solve(XVX + diag(1, ncol(XVX))*(1e-10)), silent=T)
    if(inherits(XVXi, "try-error")){
        warning("Singular matrix XVX!")
        XVXi <- MASS::ginv(XVX)
    }
    Beta <- tcrossprod(tcrossprod(XVXi, xi), Vi)%*%matrix(y,ncol=1)
    #se=XVXi
    #tvalue <- beta/se
    se <- sqrt(XVXi[1,1])
    tvalue <- Beta[1,1] / se
    pvalue <- 2 * stats::pt(abs(tvalue), df=length(y)-ncol(xi),lower.tail = FALSE)
    
    result <- c(Beta,se,pvalue)
    return(result)
}

if (length(un.in) == 1 || ncpus == 1) {
    revl <- lapply(1:length(un.in), FUN = mlm)
}else {
    if (Sys.info()[["sysname"]] == "Linux") {
        revl <- parallel::mclapply(1:length(un.in), FUN = mlm, mc.cores = ncpus)
    }else {
        suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
        suppressMessages(snowfall::sfExport("geno", "y", "un.in", "Cov","Vi"))
        suppressMessages(snowfall::sfLibrary(stats))
        suppressMessages(snowfall::sfLibrary(MASS))
        revl <- snowfall::sfLapply(1:length(un.in), fun = mlm)
        suppressMessages(snowfall::sfStop())
    }
}
revl <- do.call(rbind, revl)
if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
}
ld.index <- LD.remove(index = un.in, geno = geno, value = revl[,3],
                      LD.threshold = LD.threshold, verbose = verbose)
if(verbose) cat(length(ld.index),"markers retained after LD removing","\n")

######################
###permutation test###
######################
if (verbose) {
    cat(paste(paste(rep("-", 19), collapse = ""),
              "Permutation Test Step",
              paste(rep("-", 19), collapse = ""), sep = ""), "\n")
}
if (verbose) {
    cat("Calculating initial p values...\n")
    pb <- pbapply::timerProgressBar(width = 30, char = "-", style = 3)
    on.exit(close(pb))
}
###initial p values###
X1<-t(as.matrix(geno[c(ld.index), ], ncol=length(ld.index)))
x<-cbind(X1,Cov)
fit.mv <- gaston::lmm.aireml(y, X=x, K=G, 
							 EMsteps_fail= 10L,eps = 1e-02,
							 get.P =TRUE,verbose= FALSE)
Beta <- fit.mv$BLUP_beta
if (!is.null(fit.mv$varbeta)) {
	varbeta <- as.vector(diag(fit.mv$varbeta))
}else {
	vg <- fit.mv$tau[1]
	ve <- fit.mv$sigma2
	V <- G*vg + diag(ve, dim(G)[1])
	P <- fit.mv$P
    V_VPV <- V - tcrossprod(crossprod(V, P), V )
	X <- x
    XXiX <- try(tcrossprod(solve(crossprod(X)), X), silent = TRUE)
    if (inherits(XXiX, "try-error")) {
        warning("The CV matrix is singular!\nUsing general inverse insteaded")
        XXiX <- tcrossprod(MASS::ginv(crossprod(X)), X)
    }
    rm(X)
    varbeta <- tcrossprod(tcrossprod(XXiX, V_VPV), XXiX)
    rm(XXiX, V_VPV)
}
se <- sqrt(diag(varbeta))
SNPX_ad <- Beta / se
P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(x),lower.tail=FALSE)
#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
P_true<-P_X[1:ncol(X1)]
 if (verbose) {
    pbapply::setTimerProgressBar(pb, 1)
    cat("\n")
}

###permutation
if (permu.num > ncpus && ncpus > 1) {
    pbseq <- seq(permu.num, 1, -ceiling(permu.num / ncpus))[1:ncpus]
}else {
    pbseq <- 1:permu.num
}
if (verbose) {
    cat("Permutation repeating", permu.num, "times..."," \n")
    maxpb <- length(pbseq) + 1
    pb <- pbapply::timerProgressBar(max = maxpb, width = 30,
                                    char = "-", style = 3)
    on.exit(close(pb))
}
y <- matrix(y, ncol = 1)
wh.index <- c(ld.index,specific)
X2<-as.matrix(geno[c(wh.index),],nrow=length(wh.index))
permu.glm <- function(j){
	index <- 1:length(y)
	index.shuffle = sample(index,length(index),replace=FALSE)
	
	myy = y[index.shuffle,1]
	myy <- matrix(myy, ncol = 1)
	w <- Cov
    DF <- length(myy) - ncol(w) - 1
    wwi <- try(solve(crossprod(w)), silent = TRUE)
    if (inherits( wwi, "try-error")) {
        warning("The cor-matrix of CV is singular!")
        wwi <- MASS::ginv(crossprod(w))
    }
    My <- myy - tcrossprod(tcrossprod(w, wwi), w) %*% myy
    yy <- sum(My ^ 2)
	OL <- function(i) {
		x <- matrix(X2[i, ], ncol = 1)
		xx <- sum(x ^ 2)
		xw <- crossprod(x, w)
		xy <- crossprod(x, My)
		xhxh <- xx - tcrossprod(tcrossprod(xw, wwi), xw)
		
		b <- xy / xhxh
		RSS <- yy - b * xy	
		se <- sqrt(RSS / (xhxh * DF))
		
		pvalue <- 2 * stats::pt(abs(b / se), DF, lower.tail = FALSE)
		pvalue[which(is.na(pvalue))] <- 1
		revl <- c(b, se, pvalue)
		
		return(revl)
    }
  
	revl0 <- lapply(1:length(wh.index), FUN = OL)
	revl0 <- do.call(rbind, revl0)
	colnames(revl0) <- c("Effect", "SE", "Pvalue")
	if (verbose && j %in% pbseq) {
		pbapply::setTimerProgressBar(pb, which(sort(pbseq) %in% j))
	}
	return(revl0[, 3])
}
set.seed(6666)
if (ncpus == 1) {
    p_matrix <- lapply(1:permu.num, FUN = permu.glm)
}else {
    if (Sys.info()[["sysname"]] == "Linux") {
        p_matrix <- parallel::mclapply(1:permu.num, FUN = permu.glm,
                                 mc.cores = ncpus,
                                 mc.preschedule = TRUE)
    }else {
        suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
        suppressMessages(snowfall::sfExport("X2", "y", "Cov", "wh.index", 
		                                    "pbseq", "pb", "verbose"))
        suppressMessages(snowfall::sfLibrary(pbapply))
        suppressMessages(snowfall::sfLibrary(stats))
		suppressMessages(snowfall::sfLibrary(MASS))
        p_matrix <- snowfall::sfLapply(1:permu.num, fun = permu.glm)
        suppressMessages(snowfall::sfStop())
    }
}
p_matrix <- do.call(rbind, p_matrix)
colnames(p_matrix) <- wh.index

if(length(P_true)>1){
	p_matrix1<-as.matrix(rbind(P_true, p_matrix[,c(1:length(ld.index))]),ncol=length(P_true))
}else{
	p_matrix1<-as.matrix(c(P_true, p_matrix[,c(1:length(ld.index))]),ncol=length(P_true))
}
per.p<-apply(p_matrix1,2,function(x){( length(which(x<=x[1] ) ) -1 )/permu.num})
per.index<-ld.index[which(per.p<=permu.threshold)]

# minp<-apply(p_matrix[,c(1:length(ld.index))],1,min)
# pthre<-quantile(minp,p.threshold)
# per.index<-ld.index[which(P_true<=pthre)]
if (verbose) {
    pbapply::setTimerProgressBar(pb, maxpb)
	cat("\n")
}
if(length(per.index)==0) {
	stop("Permutation test remove all candidates!")
}
if(verbose) cat(length(per.index),"markers retained after permutation","\n")

######################
###specific signals###
######################
egwas.index<-per.index
egwas.index<-sort(egwas.index)
if(is.null(specific) || length(specific)==0){
	result <- list(Base.signal = signal, EGWAS = egwas.index, Specific = NULL)
}else{
	if (verbose) {
		cat(paste(paste(rep("-", 18), collapse = ""),
              "Specific Markers Filter",
              paste(rep("-", 18), collapse = ""), sep = ""), "\n")
		cat(length(specific),"specific markers exist","\n")
	}
	
	un.in <- specific
	if (length(un.in) == 1 || ncpus == 1) {
		revl <- lapply(1:length(un.in), FUN = mlm)
	}else {
		if (Sys.info()[["sysname"]] == "Linux") {
			revl <- parallel::mclapply(1:length(un.in), FUN = mlm, mc.cores = ncpus)
		}else {
			suppressMessages(snowfall::sfInit(parallel = TRUE, cpus = ncpus))
			suppressMessages(snowfall::sfExport("geno", "y", "un.in", "Cov","Vi"))
			suppressMessages(snowfall::sfLibrary(stats))
			suppressMessages(snowfall::sfLibrary(MASS))
			revl <- snowfall::sfLapply(1:length(un.in), fun = mlm)
			suppressMessages(snowfall::sfStop())
		}
	}
	revl <- do.call(rbind, revl)
	if (verbose) {
		pbapply::setTimerProgressBar(pb, 1)
		cat("\n")
	}
	ld.index <- LD.remove(index = un.in, geno = geno, value = revl[,3],
						  LD.threshold = LD.threshold, verbose = verbose)
	if(verbose) cat(length(ld.index),"specific markers retained after LD removing","\n")

	X1<-t(as.matrix(geno[c(ld.index), ], ncol=length(ld.index)))
	Cov <- cbind(Cov,t(geno[egwas.index,]))
	x<-cbind(X1,Cov)
	fit.mv <- gaston::lmm.aireml(y, X=x, K=G, 
								 EMsteps_fail= 10L,eps = 1e-02,
								 get.P =TRUE,verbose= FALSE)
	Beta <- fit.mv$BLUP_beta
	if (!is.null(fit.mv$varbeta)) {
		varbeta <- as.vector(diag(fit.mv$varbeta))
	}else {
		vg <- fit.mv$tau[1]
		ve <- fit.mv$sigma2
		V <- G*vg + diag(ve, dim(G)[1])
		P <- fit.mv$P
		V_VPV <- V - tcrossprod(crossprod(V, P), V )
		X <- x
		XXiX <- try(tcrossprod(solve(crossprod(X)), X), silent = TRUE)
		if (inherits(XXiX, "try-error")) {
			warning("The CV matrix is singular!\nUsing general inverse insteaded")
			XXiX <- tcrossprod(MASS::ginv(crossprod(X)), X)
		}
		rm(X)
		varbeta <- tcrossprod(tcrossprod(XXiX, V_VPV), XXiX)
		rm(XXiX, V_VPV)
	}
	se <- sqrt(diag(varbeta))
	SNPX_ad <- Beta / se
	P_X<-2*stats::pt(abs(SNPX_ad),length(y)-ncol(x),lower.tail=FALSE)
	#P_X<-2*pnorm(abs(SNPX_ad),lower=FALSE)
	spe_true <-P_X[1:ncol(X1)]
	speld.index<-ld.index
	ld.index.s<-specific[which(specific %in% speld.index)]
	
	if(length(spe_true)>1){
		p_matrix2<-as.matrix(rbind(spe_true,
					p_matrix[,which(wh.index %in% speld.index)]),ncol=length(spe_true))
	}else{
		p_matrix2<-as.matrix(c(spe_true,
					p_matrix[,which(wh.index %in% speld.index)]),ncol=length(spe_true))
	}
	spe.p<-apply(p_matrix2,2,function(x){( length(which(x<=x[1] ) ) -1 )/permu.num})
	spe.index<-ld.index.s[which(spe.p<=permu.threshold)]

	# minp<-apply(p_matrix[,which(wh.index %in% speld.index)],1,min)
	# ppthre<-quantile(minp,p.threshold)
	# spe.index1<-specific[which(specific %in% speld.index)]
	# spe.true1<-spe.true[which(specific %in% speld.index)]
	# spe.index<-spe.index1[which(spe.true1<=ppthre)]

	if(verbose) cat(length(spe.index),"specific markers retained after permutation with overlap fixed","\n")
	result <- list(Base.signal = signal, EGWAS = egwas.index, Specific = spe.index)
	# egwas.index<-c(per.index,spe.index)
	# seqQTN<-egwas.index
	# if (verbose) cat("Total",length(seqQTN),"markers was detected", "\n")
	# if(length(seqQTN)>=2){

		# if (verbose) cat("Stepwise variable selected is working....")

		# stepdata <- data.frame(Y=y, geno[,seqQTN])
		# colnames(stepdata) <- c("Y", as.character(paste0("M",seqQTN)))
        # ##diff version stepwise
		# reg <- StepReg::stepwise(formula=Y~.,
								# data=stepdata,
								# selection="bidirection",
								# select="BIC",
								# sle = p.threshold,
								# sls = p.threshold)
		# gc()
		# name <- as.data.frame(reg[[5]][[1]][-1,])
		# name<-subset(name,name[,5]<=p.threshold)[,1]
		# #name<-subset(name,name[,5]<=p.threshold/nrow(name))[,1]
		# egwas.index<-as.numeric(gsub("M","",name))
		# if(length(egwas.index)==0) {
			# stop("StepReg remove all pseudo QTNs!")
		# }
		# egwas.index<-sort(egwas.index)
		# ui=0;si=0
		# for(i in 1:length(egwas.index)){
			# if(egwas.index[[i]]%in%per.index){
				# ui=ui+1
				# names(egwas.index)[i]=paste0("Union",ui)
			# }
			# if(egwas.index[[i]]%in%spe.index){
				# si=si+1
				# names(egwas.index)[i]=paste0("Specific",si)
			# }
		# }
		# if(verbose) cat(c("(DONE)","\n"))
		# if(verbose) cat(length(egwas.index),"markers retained after StepReg","\n")
		# if(verbose) cat(egwas.index,"\n")
	# }
# }
# signal[[ (length(signal)+1) ]]<-egwas.index
# names(signal)=c(colnames(GWAS.results),"EGWAS")

}
time2 <- as.numeric(Sys.time())
time.cal <- round(time2-time1)

if(verbose){
	cat(paste(paste(rep("=", 20), collapse = ""),
              "EGWAS  ACCOMPLISHED",
              paste(rep("=", 20), collapse = ""), sep = ""), "\n")
	cat("Time used:", time.trans(time.cal), "\n")
}
return(result)
}

