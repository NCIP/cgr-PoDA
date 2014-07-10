rm(list=ls())

load("CGEMS-and-HCC-for-PoDA.RData");

out.s <- c();
the.names <- names(gene2snps.liv);
for (i in 1:length(the.names)){
  n       <- the.names[i];
  n.v     <- (gene2snps.liv)[[i]];
  snp.num <- length(n.v);
  snp.uni <- length(unique(n.v));
  if (snp.num > 0){
    if (snp.num == snp.uni){
      out.s <- rbind(out.s, c(n, snp.num, snp.uni, "eq"));
    } else {
      out.s <- rbind(out.s, c(n, snp.num, snp.uni, "ne"));
    }
  } else {
    out.s <- rbind(out.s, c(n, 0, NA, NA));
  }
}

colnames(out.s) <- c("GeneID", "SNP.number", "SNP.unique.number", "two.number");
out.i  <- which(out.s[,4]=="ne");
out.s[out.i,]; 