load("nash-dataForPoDA.RData")
source("PoDA.R")

# we need a boolean status variable where "cases" are TRUE;
# here is an example where nashdx==2 is a case,
# nashdx==0 is a control, and we toss out any borderline
# nashdx where nashdx is neither 0 nor 2:
nashstat <- nashclin[rownames(nashgeno.onpath),"nashdx"]==2
nashstat[!nashclin[rownames(nashgeno.onpath),"nashdx"]%in%c(0,2)] <- NA

# toss out any samples with NA status befor computing
nashgeno.onpath <- nashgeno.onpath[!is.na(nashstat),]
nashstat <- nashstat[!is.na(nashstat)]


# pathBlock goes from 1 to 36
# each is a chunk of 40 paths
paths <- (1:40)+40*(pathBlock-1)
# but the last chunk should only go to 1433
paths <- paths[paths<=length(path2symb)]

# make the randomization reproducible and different across all runs
set.seed(pathBlock*i)


# i ranges from 0 [true pathways] to 
# 10 [10 runs of 10 resamplings = 1000 rand paths]
if (i==0) {
	# non-scrambled pathways
	outfile <- paste("PoDA-out/nashPoDA.pathblock_",pathBlock,".truePaths.Rdata",sep="")
	out <- podaDS(nashgeno.onpath, nashstat, path2symb[paths], symb2snp, NsamplePerms=100)
	save(out,file=outfile)	
} else {
	outfile <- paste("PoDA-out/nashPoDA.pathBlock_",pathBlock,".permRun_",i,".Rdata",sep="")
	out <- sapply(1:10, function(run){ # 10 random pathways
		runOut<-try(podaDS(nashgeno.onpath, nashstat, path2symb[paths], symb2snp, NsamplePerms=100, resamplePaths=TRUE))
		if(class(runOut)=="try-error"){runOut<-rep(NA,length(path2symb[paths]))}
		return(runOut)
	})
	save(out,file=outfile)
}

traceback()
