#!/usr/bin/env Rscript

# plotSFS.R <SFS file> <output prefix> <fold>
# if fold = 0 the input SFS is unfolded, otherwise if fold = 1 the input SFS is assumed to be folded


# parse inputs
fold <- 1
args <- commandArgs(trailingOnly=TRUE)

sfs <- scan(args[1])
outprefix <- args[2]
fold <- as.numeric(args[3])
n <- (length(sfs)-1)/2
xlabel <- NULL
idx <- NULL
name.n <- NULL

if (fold == 1) {
   xlabel = "MAF"
   idx = n+1
   name.n = n
   if (sum(sfs[(n+2):(length(sfs))]) > 0) cat("WARNING: Assuming folded SFS input despite that it looks unfolded...\n")
} else {
   idx = 2*n
   xlabel = "Derived allele frequency"
   name.n = idx-1
}

# plot
pdf(file=paste0(outprefix,".pdf"),width=14,height=7)
barplot(sfs[2:idx], xlab=xlabel, ylab="Number SNPs", names=1:name.n, cex.names=0.8, cex.axis=1.2, cex.lab=1.2) # plot SFS without fixed categories
invisible(dev.off())
