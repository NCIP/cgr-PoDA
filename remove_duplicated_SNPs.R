rm(list=ls())

load("CGEMS-and-HCC-for-PoDA.RData");

gene2snps.liv2 <- list();
the.names <- names(gene2snps.liv);
for (i in 1:length(the.names)){
  n       <- the.names[i];
  n.v     <- (gene2snps.liv)[[i]];
  snp.num <- length(n.v);
  snp.u   <- unique(n.v);
  snp.uni <- length(unique(n.v));
  if (snp.num > 0){
    if (snp.num == snp.uni){
      gene2snps.liv2[[i]] <- n.v;
    } else {
      snp.o               <- n.v[!duplicated(n.v)];
      gene2snps.liv2[[i]] <- snp.o;
    }
  } else {
    gene2snps.liv2[[i]] <- n.v;
  }
}

names(gene2snps.liv2) <- names(gene2snps.liv);
out.list <- "CGEMS-and-HCC-for-PoDA2.RData";
save(gene2snps.cgems, gene2snps.liv, gene2snps.liv2, geno.cgems, 
     geno.liv, path2genes.cgems, path2genes.liv, stat.cgems, stat.liv, file=out.list);
