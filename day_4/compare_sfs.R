#!/usr/bin/env Rscript

# compare_sfs.R <SFS file> <output prefix> <use folded (0|1)> <plot type: "bar"|"scatter">

# parse inputs
args <- commandArgs(trailingOnly=TRUE)

obs <- scan(args[1])
outprefix <- args[2]
fold <- as.numeric(args[3]) # is the SFS folded (1 = yes, 0 = no)
plottype <- args[4]

obs <- obs[-1] # remove fixed zero category
nchr = length(obs) # number chromosomes
n = nchr/2
xlabel = NULL
if (fold) {
	obs <- obs[1:n]
	xlabel = "MAF"
} else {
	obs <- obs[-length(obs)] # format observed SFS
	xlabel = "Derived allele frequency"
}
s =  sum(obs) # number segregating sites
nlab = ifelse(fold,n,nchr-1)

# calculate expected unfolded SFS counts
a = sum(1/1:(nchr-1))
sfs.expected <- sapply(1:(nchr-1),function(x,segsites,a){segsites/(x*a)}, segsites=s, a=a)

# fold SFS if necessary
if (fold) {
	sfs.fold = sapply(1:(n-1),function(x,sfs,nchr){sfs[x]+sfs[nchr-x]},nchr=nchr,sfs=sfs.expected)
	sfs.fold = c(sfs.fold,sfs.expected[n])
	sfs.expected = sfs.fold
}

if (plottype == "bar") {
# plot observed vs expected as a barplot
	pdf(file=paste0(outprefix,"_barplot.pdf"),width=14,height=7)
	sfs.mat <- t(as.matrix(data.frame(obs=obs, expected=sfs.expected)))
	barplot(sfs.mat, beside=TRUE, names=c(1:nlab), cex.names=0.7, col=rep(c("grey80","grey10"),nlab),xlab=xlabel,ylab="Number SNPs", cex.lab=1.2, cex.axis=1.2)
	legend('topright',c("observed","expected"),fill=c("grey80","grey10"),bty='n',cex=1.2)
	invisible(dev.off())
} else if (plottype == "scatter") {
	# plot expected vs observed scatterplot
	pdf(file=paste0(outprefix,"_scatter.pdf"))
	sfs.obs.p <- obs/sum(obs)
	sfs.e.p <- sfs.expected/sum(sfs.expected)
	px = -log10(sfs.e.p)
        py = -log10(sfs.obs.p)
        plot(y=py, x=-log10(sfs.e.p), ylab=paste0("-log10(observed ", xlabel, " probability)"), xlab=paste0("-log10(expected ", xlabel, " probability)"), 
        xlim=c(min(c(px,py)),max(c(px,py))),cex.lab=1.2, cex.axis=1.2)
	abline(0,1,col="red",lty=2)
	invisible(dev.off())
} else stop("Unrecognized plot type")
