# load the generated data
nash.S <- NULL
nash.DS <- NULL
nash.randDS <- NULL
for (pathBlock in 1:36){  # we had 36 chunks
  # load the "true" data and push it into nash.S and nash.DS from file:
  trueName <- paste("PoDA-out/nashPoDA.pathblock_", pathBlock,".truePaths.Rdata", sep="")
  # put in a failsafe in case it's missing
  npaths <- 40-7*(pathBlock==36)
  out <- list(S=matrix(NA,nr=239,nc=npaths),DS=rep(NA,npaths))
  try(load(trueName))
  nash.S <- cbind(nash.S, out$S)
  nash.DS <- c(nash.DS, out$DS)
  # now load the "random pathway" data and push it into nash.randDS
  tmp <- NULL # (again, just easier)
  for (permRun in 1:100){
    permName <- paste("PoDA-out/nashPoDA.pathBlock_", pathBlock,".permRun_",permRun,".Rdata",sep="")
    # this is a failsafe if one is missing
    out <- matrix(NA,nr=npaths,nc=10) #in case it fails to load
    try(load(permName),silent=T)
    tmp <- cbind(tmp,out)
  }
  nash.randDS <- rbind(nash.randDS,tmp)
}

# clear out some memory
rm(npaths,out,pathBlock,permName,permRun,tmp,trueName);gc()

# get the p-values!
nash.pDS <- rowSums(nash.DS<nash.randDS,na.rm=T)/rowSums(!is.na(nash.randDS))

# make a nice table

load("nash-dataForPoDA.RData")
invlist <- function(x){
	split(rep(names(x),sapply(x,length)),unlist(x))
}
# next 3 lines get the # of genes in each pathway that are represented in 
# the NASH data
ng <- unique(unlist(invlist(symb2snp)[colnames(nashgeno.onpath)]))
nash.pl <- sapply(path2symb[names(nash.DS)],function(genes){sum(genes%in%ng,na.rm=T)})
rm(ng);gc()


nash.tab <- data.frame(
	pathName=names(nash.DS), # path name
	l=nash.pl, # path length 
	DS=nash.DS, # PoDA distinction score
	p.DS=nash.pDS, # PoDA p-value
	stringsAsFactors=FALSE
)

write.table(nash.tab, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, file="nash.PoDA_table.txt")
save(nash.tab,nash.S,file="nash.PoDA_table_and_S_vals.Rdata")


