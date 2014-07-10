##############################################################
##
##         PoDA script - bare bones computation
##  
##############################################################
#
# podaDS() is the main function to compute PoDA distinction 
# score for a set of pathways.
#
# Usage: 
#   podaDS(genotypes, status, pathGeneList, geneSnpList, NsamplePerms, resamplePaths)
# 
# For more details see:
#   podaDS() comments below, 
#   PoDA-example.R script, 
#   Paper:
#     R Braun & K Buetow, "Pathways of Distinction Analysis: A 
#     new technique for multi-SNP analysis of GWAS data", PLoS 
#     Genetics, 2011.
#
# Rosemary Braun <braunr@mail.nih.gov>
##############################################################


#-------------------------------------------------------------
# ANCILLARY FUNCTIONS
#-------------------------------------------------------------


.chi2 <- function(x,y){
	# simple but fast computation of chi2 statistic (no continuity correction)
	x <- table(x,y)
	E <- outer(rowSums(x),colSums(x),"*")/sum(x)
	return( sum( (x-E)^2/E ))
}

.podaD <- function(genotypes, boolstat, foldCV="LOO"){
	# compute Di (Eq. 1) for all SNPs given in genotypes
	# for all samples either with leave-one-out or coarser
	# CV folds.
	
	maxF <- floor(min(sum(boolstat,na.rm=TRUE),sum(!boolstat,na.rm=TRUE))/2)
	if (is.numeric(foldCV) & foldCV>=2 & foldCV<=maxF){
		# use foldCV number of folds
		foldmask <- boolstat
		foldmask[boolstat] <- sample(rep(1:foldCV,length=sum(boolstat)))
		foldmask[!boolstat] <- sample(rep(1:foldCV,length=sum(!boolstat)))
		out <- genotypes
		for (fold in 1:foldCV) {
			gi <- colMeans(genotypes[boolstat&(foldmask!=fold),,drop=FALSE],na.rm=TRUE)
			fi <- colMeans(genotypes[(!boolstat)&(foldmask!=fold),,drop=FALSE],na.rm=TRUE)
			out[foldmask==fold,] <- t(
				apply( genotypes[foldmask==fold,,drop=FALSE],1,function(yi){
					abs(yi-fi) - abs(yi-gi)
				})/2
			)
		}
	} else {
		# use a LOO approach
		out <- t(sapply(1:length(boolstat),function(Y){
			boolstat[Y] <- NA # remove the proband Y
			yi <- genotypes[Y,]  # individual Y's MAFs
			gi <- colMeans(genotypes[boolstat,,drop=F],na.rm=TRUE)  # case MAFs
			fi <- colMeans(genotypes[!boolstat,,drop=F],na.rm=TRUE) # control MAFs
			abs(yi-fi) - abs(yi-gi)
		})/2)
		if (ncol(genotypes)==1) {out<-t(out)}
		dimnames(out)<-dimnames(genotypes)
	}
	return(out)
}


.podaW <- function(genotypes, status, pathGene2snp,withS=FALSE) {  
	# compute S (Eq. 2) and W (Eq. 3) for a pathway

	getSnpForGene <- function(geneSnps){
		geneSnps <- geneSnps[geneSnps%in%colnames(genotypes)]
		chis<-apply(genotypes[,geneSnps,drop=FALSE], 2, .chi2, status)
		return(names(chis)[which.max(chis)])
	}
	keptSnps <- unlist(sapply(pathGene2snp,getSnpForGene))
	if (length(keptSnps)<3) {
		# don't bother id there are fewer than 3 snps
		W <- NA
		S <- NA
	} else {
		D <- .podaD(genotypes[,keptSnps], status, foldCV=100)
		S <- rowMeans(D,na.rm=TRUE)/apply(D,1,sd,na.rm=TRUE)
		Nca <- sum(status)
		W <- sum(rank(S)[status])-Nca*(Nca+1)/2
	}
	if (withS) {
		return(list(W=W,S=S))
	} else {
		return(W)
	}
}


.podaDS.truePaths <- function(genotypes, status, pathGeneList, geneSnpList, NsamplePerms){
	# compute DS (Eq. 4) for all pathways
		
	DS <- numeric(length(pathGeneList)); names(DS)<-names(pathGeneList)
	S <- matrix(nrow=length(status),ncol=length(pathGeneList),dimnames=list(names(status),names(pathGeneList)))
	
	for (path in names(pathGeneList)){
		pathGenes <- pathGeneList[[path]]
		pathGenes <- pathGenes[pathGenes%in%names(geneSnpList)] # only keep genes with SNPs
		pathGene2snp <- geneSnpList[pathGenes] #only keep the pathway's genes
		if (length(pathGene2snp)<3){
			warning(paste("pathway",path,"has fewer than three genes probed; returning NA"))
			DS[path] <- NA
			S[,path] <- NA
		} else {
			WS <- .podaW(genotypes,status,pathGene2snp,withS=TRUE)
			Wstar <- sapply(1:NsamplePerms,function(i){
				.podaW(genotypes,sample(status),pathGene2snp)
			})
			DS[path] <- (WS$W - mean(Wstar,na.rm=TRUE))/sd(Wstar,na.rm=TRUE)
			S[,path] <- WS$S
		}
	}
	return(list(S=S,DS=DS))
}


.podaDS.resampPaths <- function(genotypes, status, pathGeneList, geneSnpList, NsamplePerms){
	# compute DS (Eq. 4) for random pathways (for significance testing)

	DS <- rep(NA,length(pathGeneList))
	names(DS)<-paste("RandomSetWithLengthOf:",names(pathGeneList),sep="  ")
	pathLengths <- sapply(pathGeneList,function(pathGenes){sum(pathGenes%in%names(geneSnpList))})
	resampLengths <- unique(pathLengths[pathLengths>=3])
	for (pathLength in resampLengths){
		pathGenes <- sample(names(geneSnpList),pathLength)
		pathGene2snp <- geneSnpList[pathGenes] #only keep the pathway's genes
		WS <- .podaW(genotypes,status,pathGene2snp)
		Wstar <- sapply(1:NsamplePerms,function(i){
			.podaW(genotypes,sample(status),pathGene2snp)
		})
		DS[pathLengths==pathLength] <- (WS - mean(Wstar,na.rm=TRUE))/sd(Wstar,na.rm=TRUE)
	}
	return(DS)
}



#-------------------------------------------------------------
# MAIN PoDA FUNCTION
#-------------------------------------------------------------


podaDS <- function(genotypes, status, pathGeneList, geneSnpList, NsamplePerms=100, resamplePaths=FALSE){
	
	######################################################################
	# podaDS(): function to compute the PoDA distinction score DS for 
	#           all pathways in pathGeneList (cite: R Braun & K Buetow, 
	#           "Pathways of Distinction Analysis: A new technique for 
	#           multi-SNP analysis of GWAS data", PLoS Genetics, 2011.)
	# usage: 
	#   podaDS(genotypes, status, pathGeneList, geneSnpList, NsamplePerms, resamplePaths)
	#
	# ---- INPUTS ------------------------------------------------------
	# genotypes = named matrix of 0/1/2 genotypes, with SNPs in 
	#             columns and samples in rows
	# status = boolean vector indicating case/control status for the 
	#          samples (TRUE=case)
	# pathGeneList = named list mapping pathways to genes
	# geneSnpList = named list mapping genes to SNPs; expected that 
	#               all SNPs in geneSnpList are in genotypes matrix
	# NsamplePerms = number of case/control permutations for the 
	#                computation of DS
	# resamplePaths = boolean flag indicating whether or not pathways 
	#                 are to be randomly generated corresponding to 
	#                 pathway lengths in pathGeneList (for p-values)
	#
	# ---- OUPUT ------------------------------------------------------
	# If resamplePaths=FALSE, a list with elements:
	# 	S = a matrix of PoDA S values for each pathway (col) in each 
	#       sample (rows)
	#	DS = a vector with the distinction score for each pathway
	#
	# If resamplePaths=TRUE, a vector with distinction scores for 
	#   random gene sets of the same length as the true pathways in 
	#   pathGeneList
	######################################################################
	
	# check input
	if (!is.list(pathGeneList)) {stop("expected a list mapping pathways to genes")}
	if (!is.list(geneSnpList)) {stop("expected a list mapping genes to SNPs")}
	if (length(status)!=nrow(genotypes)) {stop("number of samples in genotypes and status do not match")}

	# remove NA-status individuals
	genotypes <- genotypes[!is.na(status),]
	status <- status[!is.na(status)]
	
	# compute podaDS
	if (resamplePaths){
		out <- .podaDS.resampPaths(genotypes, status, pathGeneList, geneSnpList, NsamplePerms)
	} else {
		out <- .podaDS.truePaths(genotypes, status, pathGeneList, geneSnpList, NsamplePerms)
	}
	return(out)
}

