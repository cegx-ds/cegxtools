#!/usr/bin/Rscript

thisdir = getwd()
tempdir <- function() thisdir
unlockBinding("tempdir", baseenv())
assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())

library(methylKit)
library(tools)
args<-commandArgs(TRUE)

# readBismarkCytosineReport function
devtools::source_gist("4839e615e2401d73fe51")

file.list = list(args[1],args[2],args[3],args[4])

# This reads straight from the output of bismark, works for both strands
myobj = readBismarkCytosineReport(
        file.list,
          sample.id=list(
  	  strsplit(basename(args[1]), '[_]')[[1]][1],
          strsplit(basename(args[2]), '[_]')[[1]][1],
  	  strsplit(basename(args[3]), '[_]')[[1]][1],
  	  strsplit(basename(args[4]), '[_]')[[1]][1]
  	       ),
  	assembly="hg38",
  	treatment=c(1,1,0,0),
  	min.cov=args[5])

pdf(file="Methylkit.01.MethylationStats.Sample01.pdf",width=12,height=12) ; getMethylationStats(myobj[[1]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.MethylationStats.Sample02.pdf",width=12,height=12) ; getMethylationStats(myobj[[2]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.MethylationStats.Sample03.pdf",width=12,height=12) ; getMethylationStats(myobj[[3]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.MethylationStats.Sample04.pdf",width=12,height=12) ; getMethylationStats(myobj[[4]],plot=T,both.strands=T); dev.off()

pdf(file="Methylkit.01.CoverageStats.Sample01.pdf",width=12,height=12) ; getCoverageStats(myobj[[1]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.CoverageStats.Sample02.pdf",width=12,height=12) ; getCoverageStats(myobj[[2]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.CoverageStats.Sample03.pdf",width=12,height=12) ; getCoverageStats(myobj[[3]],plot=T,both.strands=T); dev.off()
pdf(file="Methylkit.01.CoverageStats.Sample04.pdf",width=12,height=12) ; getCoverageStats(myobj[[4]],plot=T,both.strands=T); dev.off()

meth = unite(myobj)
pdf(file="Methylkit.02.Meth_CorrelationPlot.pdf",width=12,height=12); getCorrelation(meth,plot=T); dev.off()

# cluster all samples using correlation distance and return a tree object for plclust
hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

# cluster all samples using correlation distance and plot hiarachical clustering
pdf(file="Methylkit.03.PCASamples.ctreeplot.pdf",width=12,height=12); clusterSamples(meth, dist="correlation", method="ward", plot=TRUE); dev.off()

# screeplot of principal component analysis.
pdf(file="Methylkit.03.PCASamples.screeplot.pdf",width=12,height=12); PCASamples(meth, screeplot=TRUE); dev.off()

# principal component anlaysis of all samples.
pdf(file="Methylkit.03.PCASamples.pcaxyplot.pdf",width=12,height=12); PCASamples(meth); dev.off()

myDiff=calculateDiffMeth(meth)
write.table(myDiff, file="Methylkit.04.DiffMeth.tsv", sep='\t', quote=FALSE)
myDiff25p.hyper = get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
write.table(myDiff25p.hyper,"Methylkit.04.hyper_methylated.txt",sep='\t', quote=FALSE)
myDiff25p.hypo = get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
write.table(myDiff25p.hypo,"Methylkit.04.hypo_methylated.txt",sep='\t', quote=FALSE)
myDiff25p = get.methylDiff(myDiff,difference=25,qvalue=0.01)
write.table(myDiff25p,"Methylkit.04.differentialy_methylated.txt",sep='\t', quote=FALSE)
pdf("Methylkit.04.DiffMethPerChr.pdf",width=12,height=12); diffMethPerChr(myDiff,plot=TRUE,qvalue.cutoff=0.01,meth.cutoff=25); dev.off()

gene.obj=read.transcript.features("/bi/group/cegx/UCSC_database/refseq.hg38.bed.txt")

diffAnn = annotate.WithGenicParts(myDiff25p, gene.obj)
write.table(getAssociationWithTSS(diffAnn),"Methylkit.05.AssociationWithTSS.txt", sep='\t', quote=FALSE)

pdf("Methylkit.05.TargetAnnotation.piechart1.pdf",width=12,height=12); plotTargetAnnotation(diffAnn, precedence = TRUE, main ="differential methylation annotation"); dev.off()

write.table(getFeatsWithTargetsStats(diffAnn, percentage = TRUE),"Methylkit.05.FeatsWithStargetsStats.txt",sep='\t', quote=FALSE)

tmpfiles <- dir(path=getwd(), pattern="*.CpG_report.txt")
file.remove(tmpfiles)

1
