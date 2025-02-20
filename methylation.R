library(methylKit)
setwd("/path/to/9spine_wgbs/pca")

directory="./mbe_gz"

filenames=list.files(path=directory, full.names=TRUE)
filenames=as.list(filenames)
names=list.files(path=directory)
names=gsub(".cov.gz", "", names)
names=as.list(names)

# Metadata for all samples
info=read.csv("./sample.csv")
info$Sex=NULL
library(dplyr)
info=info %>%
  slice(match(names, Label))

my.methRaw=methRead(location = filenames,
                    sample.id = names,
                    assembly = "9spine",
                    pipeline = 'bismarkCoverage',
                    context = "CpG",
                    treatment = c(rep(1,12), 
                                  rep(0,14)), # marine is 1, freshwater is 0
                    mincov = 5)

filtered.my.methRaw = filterByCoverage(my.methRaw, 
                                       lo.count = 5,
                                       lo.perc = NULL,
                                       hi.count = NULL,
                                       hi.perc = 99.9)

# normalize read coverages between samples to avoid bias introduced by systematically more sequenced sameples
normalized.myobj=normalizeCoverage(my.methRaw, method="median")

# merging samples (step should be saved)------
meth.all=unite(normalized.myobj, destrand = F, mc.cores = 8)

save(meth.all, info, file = "./meth_9spine.RData")

# extract filtered CpG in DSS from methylkit
# dss_filtered_cpg=read.table("~/Downloads/9spine_wgbs/dss/filtered_CpG_dss.txt")
# load("./meth_9spine.RData")
##########################################
#Exclude SNPs for methylation analyses####
##########################################
bed_to_granges = function(file){
  df = read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df = df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header = c('chr','start','end','id','score','strand')
  names(df) = header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand = gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr = with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr = with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

# The C/T and A/G location is produced before purging high LD sites, using parental samples
# load CT SNPs and duplicate column and save as bedfile, result_pat_C_T_maf.bed
CT_maf= read.csv(file="./result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./result_pat_G_A_mafSTARTEND.bed")

# interesect bedfile and unite-file
#### create overlap --> these positions are corrected for CT SNPs
unite_norm_10x_GRanges <- as(meth.all, "GRanges")
Overlap_CT=unite_norm_10x_GRanges[countOverlaps(unite_norm_10x_GRanges, blacklist_CT_bed) == 0L]
meth.1=makeMethylDB(meth.all,"methylBaseDB")
unite_norm_10x_CT <- selectByOverlap(meth.1, granges(Overlap_CT))

# make methylkit database
objDB_CT=makeMethylDB(unite_norm_10x_CT,"methylBaseDB")

##### create overlap --> these positions are corrected for GA SNPs
unite_norm_10x_CT_GRanges <- as(unite_norm_10x_CT,"GRanges")
Overlap_GA=unite_norm_10x_CT_GRanges[countOverlaps(unite_norm_10x_CT_GRanges, blacklist_GA_bed) == 0L]
unite_norm_10x_CT_GA <- selectByOverlap(objDB_CT, Overlap_GA)

unite_norm_10x_CT_GA_DB=makeMethylDB(unite_norm_10x_CT_GA,"methylBaseDB")

# generate methylation percentage for RDA analysis of gene expression and alternative splicing
meth.rda=regionCounts(meth.all, regions = as(unite_norm_10x_CT_GA, "GRanges"))
perc.meth.rda=percMethylation(meth.rda)
perc.meth.rda.t=t(perc.meth.rda)
write.table(perc.meth.rda.t, "./meth.rda.txt", sep = "\t", col.names = F, quote = F)

#save.image("/Volumes/Juntao/9spine_wgbs/pca/pca_9spine.RData")

# PCA of methylation----
library(ggbiplot)

# check if the order of meth count table is the same as info
colnames(perc.meth.rda)==as.character(info$Label)

pca=prcomp(t(perc.meth.rda), center = T)
summary(pca)

write.table(pca$x, "pca_summary_methylation.txt", sep = "\t", quote = F)

# pca of methylation in brain
perc.meth.rda.brain=perc.meth.rda[,colnames(perc.meth.rda) %in% info[info$Tissue=="brain",]$Label]
colnames(perc.meth.rda.brain)

pca.brain=prcomp(t(perc.meth.rda.brain), center = T)
summary(pca.brain)

write.table(pca.brain$x, "pca_summary_methylation.brain.txt", sep = "\t", quote = F)

# pca of methylation in brain
perc.meth.rda.liver=perc.meth.rda[,colnames(perc.meth.rda) %in% info[info$Tissue=="liver",]$Label]
colnames(perc.meth.rda.liver)

pca.liver=prcomp(t(perc.meth.rda.liver), center = T)
summary(pca.liver)

write.table(pca.liver$x, "pca_summary_methylation.liver.txt", sep = "\t", quote = F)


# Define factors
ecotype=as.factor(info$Ecotype)
tissue=as.factor(info$Tissue)

# visulize PCA results
g1=ggbiplot(pca, 
            obs.scale = 1, 
            var.scale = 1,
            varname.size = 0,
            # labels = info$Label,
            var.axes = FALSE)
g2=g1+geom_point(aes(color=ecotype, fill = ecotype, shape=tissue), size=2)+
  scale_color_manual(values = c("light blue", "orange"))+
  scale_fill_manual(values = c("light blue", "orange"))+
  scale_shape_manual(values = c(21,24))+
  labs(x=paste0("PC1 (", "22.26", "%)"), y=paste0("PC2 (", "6.56", "%)"))
g3=g2 + theme_classic()+theme(text = element_text(size = 15))
print(g3)

## RDA of methylation-----
# Constructing tiling windows
# tiles=tileMethylCounts(unite_norm_10x_CT_GA, win.size = 1000, step.size = 1000) 
# nrow(tiles)
# 275598 windows (regions)

# Calculate Euclidean distance matrix using the tiling windows
# get methylation percentage for all CpG sites
# meth.regions=regionCounts(meth.all, regions = as(tiles, "GRanges"))
# perc.meth.regions=percMethylation(meth.regions)
# perc.meth.regions.t=t(perc.meth.regions)

# Combine tissue and ecotype info into the methylation data
# rownames(perc.meth.regions.t)==info$Label
df=cbind(info[,c(2:3)], perc.meth.rda.t)
str(df)

library(vegan)
library(cluster)
dist=dist(df[ ,3:ncol(df)], method = "euclidean")
library(ape)
meth.pcoa=pcoa(dist)
meth.pcoa$values # select eigenvalues > 2.75%

#Produce variable selection of the db-RDA (var1 = Tissue, var2 = Ecotype)
ordistep(rda(meth.pcoa$vectors[,1:15]~df$Tissue+df$Ecotype, scale=F))

# Results show both tissue and ecotype are significant factors (P < 0.005)

#Global RDA with only significant factors according to ordistep
rda_genet=rda(formula = meth.pcoa$vectors[,1:15]~df$Tissue+df$Ecotype, scale=F)
rda_genet
# RDA1   RDA2   
# 84361026 23185646 
# RDA1 explains 84361026/2.885e+08=0.2924126, RDA2 explains 23185646/2.885e+08=0.08036619

summary(rda_genet)

#test significance of global model
anova(rda_genet, step=1000) # P = 0.001

#Test significance for each variable
anova(rda_genet, step=1000, by='margin')
#             Df Variance       F   Pr(>F)    
#  df$Tissue   1  84239598 10.7083  0.001 ***
#  df$Ecotype  1  23307074  2.9627  0.009 ** 
#  Residual   23 180934840    

#Compute Rsquare     
RsquareAdj(rda_genet) 
# $adj.r.squared
# [1] 0.3182637

Tissue=as.factor(df$Tissue)
Ecotype=as.factor(df$Ecotype)
#partial db-RDA for tissue, and control for ecotype
rda_genet2=varpart(meth.pcoa$vectors[,1:15], Tissue, Ecotype)
rda_genet2 # X1 is tissue

anova.cca(rda(meth.pcoa$vectors[,1:15], Tissue, Ecotype), step=1000)
# $adj.r.squared for Tissue=0.26251, F=10.708, P<0.001

#partial db-RDA for ecotype, and control for tissue
rda_genet3=varpart(meth.pcoa$vectors[,1:15], Ecotype, Tissue)
rda_genet3 # X1 is Ecotype

anova.cca(rda(meth.pcoa$vectors[,1:15], Ecotype, Tissue), step=1000)
# $adj.r.squared for Ecotype=0.04249, F=2.9627, P<0.001

# Plotting RDA -------------------
library(ggplot2)
sommaire = summary(rda_genet)
sommaire
df1 = data.frame(sommaire$sites[,1:2])       # PC1 and PC2
df2 = data.frame(sommaire$species[,1:2])

# prepare df1 with group info
df1=merge(df1, df[,1:2], by=0, all=TRUE)
rownames(df1)=df1$Row.names
df1$Row.names = NULL
colnames(df1)
df1$Tissue=as.factor(df1$Tissue)
df1$Ecotype=as.factor(df1$Ecotype)
# loadings for PC1 and PC2
library(ggplot2)
rda.plot.rna = ggplot(df1, aes(x=RDA1, y=RDA2, colour=Tissue)) + 
  geom_point(size=4, aes(shape = Ecotype)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()
  # scale_colour_manual(values = c("#009E73", "#D55E00", "#999999"))
rda.plot.rna 

df2subset=df2[c("Axis.1","Axis.2"),]
title_lab=expression(bold(paste(" adj. ",R^{2}," = 0.32; ", bolditalic(P), " = 0.001"))) # fill with you own stats in sommaire

theme1=theme(panel.grid.major = element_blank(),
              panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              axis.text.x=element_text(colour="black",size=16),
              axis.text.y=element_text(colour="black",size=16),
              axis.title.x=element_text(colour="black",size=16,face="bold"),
              axis.title.y=element_text(colour="black",size=16,face="bold"),
              axis.ticks=element_line(colour="black"),
              plot.margin=unit(c(1,1,1,1),"line"),
              plot.title = element_text(size=12,hjust=0),
              aspect.ratio=1,
              legend.background = element_rect(fill ="NA"),
              panel.border=element_rect(colour = "black", fill=NA, size=1),
              legend.position = "none")
# legend.text = element_text(size=12),
# legend.title = element_text(size=12,face="bold"),
# legend.justification = c(1, 0),
# legend.key = element_rect(fill=NA)) 

# couleurs=c("gray0", "dodgerblue2","tomato2")
legend_title=c("Tissue ***", "Ecotype **")

rda.biplot.rna = rda.plot.rna +
  geom_segment(data=df2subset, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  labs(x="RDA1 (29.2%)",y="RDA2 (8.0%)") +
  ggtitle(title_lab) + theme1 
rda.biplot.rna

RDA_expr=rda.biplot.rna +
  theme1+ 
  geom_text(data = df2subset, 
            aes(x=1.2*df2subset[ ,1], 
                y=1.2*df2subset[ ,2], 
                label=legend_title), 
                size=5, vjust=0, hjust = -0.05,
                nudge_x = 0.05, colour="black", 
                fontface=c("bold"))
RDA_expr
pdf(file = "./pca/fig.1.pdf", width = 8, height = 8)
RDA_expr
dev.off()

# DMCs in brain-------
brain_subset=info[info$Tissue=="brain",]
brain_subset$Label=as.character(brain_subset$Label)
# Order is important, adjust order of sample id based on original sample id in meth
ind.brain=unite_norm_10x_CT_GA_DB@sample.ids[seq(from = 1, to = 2*nrow(brain_subset), by = 2)]
library(dplyr)
brain_subset=brain_subset %>%
  slice(match(ind.brain, Label))
meth.brain=reorganize(unite_norm_10x_CT_GA,
                    sample.ids = as.character(brain_subset$Label),
                    treatment = as.numeric(as.factor(brain_subset$Ecotype)))
meth.brain@treatment=c(rep(1, 6), rep(0, 7)) # correct treatment to marine 1, freshwater 0

# DMCs in brain
myDiff.brain=calculateDiffMeth(meth.brain, mc.cores = 4, overdispersion = "MN")
myDiff.sig.brain=getMethylDiff(myDiff.brain, difference = 0, qvalue = 0.05)
nrow(myDiff.sig.brain) # 1433

myDiff.sig.brain.hyper=getMethylDiff(myDiff.brain,  difference = 0, qvalue = 0.05, type = "hyper")
nrow(myDiff.sig.brain.hyper) # 722
myDiff.sig.brain.hypo=getMethylDiff(myDiff.brain,  difference = 0, qvalue = 0.05, type = "hypo")
nrow(myDiff.sig.brain.hypo) # 711

library(DescTools)
expected=c(0.5,0.5)
observed=c(722,711)
GTest(x=observed,
      p=expected,
      correct="none")
# G = 0.084439, X-squared df = 1, p-value = 0.7714

# heatmap of dmcs in brain
library(pheatmap)
meth.brain.dmcs.count=regionCounts(meth.brain, regions = as(myDiff.sig.brain, "GRanges"))
perc.meth.brain.dmcs=percMethylation(meth.brain.dmcs.count)
sampleinfo.brain=brain_subset[,3,drop=F]
colnames(sampleinfo.brain)="ecotype"
rownames(sampleinfo.brain)=colnames(perc.meth.brain.dmcs)
ann_colors.brain=list(ecotype=c(marine="orange", freshwater="light blue"))
pheatmap(perc.meth.brain.dmcs, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         show_colnames = FALSE,
         cluster_cols=TRUE,
         border_color = F,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sampleinfo.brain,
         annotation_colors = ann_colors.brain,
         annotation_names_col = F,
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4, width = 6, filename = "./dmc_brain.tiff")

# DMCs in liver-------
liver_subset=info[info$Tissue=="liver",]
liver_subset$Label=as.character(liver_subset$Label)
# Order is important, adjust order of sample id based on original sample id in meth
ind.liver=unite_norm_10x_CT_GA_DB@sample.ids[seq(from = 2, to = 2*nrow(liver_subset), by = 2)]
library(dplyr)
liver_subset=liver_subset %>%
  slice(match(ind.liver, Label))
meth.liver=reorganize(unite_norm_10x_CT_GA,
                      sample.ids = as.character(liver_subset$Label),
                      treatment = as.numeric(as.factor(liver_subset$Ecotype)))
meth.liver@treatment=c(rep(1, 6), rep(0, 7)) # correct treatment to marine 1, freshwater 0

# DMCs in liver
myDiff.liver=calculateDiffMeth(meth.liver, mc.cores = 4, overdispersion = "MN")
myDiff.sig.liver=getMethylDiff(myDiff.liver, difference = 0, qvalue = 0.05)
nrow(myDiff.sig.liver) # 1732

myDiff.sig.liver.hyper=getMethylDiff(myDiff.liver,  difference = 0, qvalue = 0.05, type = "hyper")
nrow(myDiff.sig.liver.hyper) # 609
myDiff.sig.liver.hypo=getMethylDiff(myDiff.liver,  difference = 0, qvalue = 0.05, type = "hypo")
nrow(myDiff.sig.liver.hypo) # 1123

library(DescTools)
expected=c(0.5,0.5)
observed=c(609,1123)
GTest(x=observed,
      p=expected,
      correct="none")
# G = 154.86, X-squared df = 1, p-value < 2.2e-16

# heatmap of dmcs in liver
library(pheatmap)
meth.liver.dmcs.count=regionCounts(meth.liver, regions = as(myDiff.sig.liver, "GRanges"))
perc.meth.liver.dmcs=percMethylation(meth.liver.dmcs.count)
sampleinfo.liver=liver_subset[,3,drop=F]
colnames(sampleinfo.liver)="ecotype"
rownames(sampleinfo.liver)=colnames(perc.meth.liver.dmcs)
ann_colors.liver=list(ecotype=c(marine="orange", freshwater="light blue"))
pheatmap(perc.meth.liver.dmcs, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         show_colnames = FALSE,
         cluster_cols=TRUE,
         border_color = F,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = sampleinfo.liver,
         annotation_colors = ann_colors.liver,
         annotation_names_col = F,
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4, width = 6, filename = "./dmc_liver.tiff")

## overlap of dmcs between brain and liver
# brain.dmcs.loc=getData(myDiff.sig.brain)[,1:2]
# brain.dmcs.loc=paste(brain.dmcs.loc$chr, brain.dmcs.loc$start, sep = "-")
# liver.dmcs.loc=getData(myDiff.sig.liver)[,1:2]
# liver.dmcs.loc=paste(liver.dmcs.loc$chr, liver.dmcs.loc$start, sep = "-")

# library(ggVennDiagram)
# brain.liver.overlap=list(brain=brain.dmcs.loc, liver=liver.dmcs.loc)
# brain.liver.overlap.venn=ggVennDiagram(brain.liver.overlap)
# brain.liver.overlap.venn

# annotated dmcs in brain and liver--------
# annotate genes associated with the above snps
library(GenomicFeatures)
library(ChIPpeakAnno)
stickle=makeTxDbFromGFF("/path/to/9spine_wgbs/9spine_genome/liftover/ncbi/v7.phased.gff")

stickle_gene=genes(stickle)
stickle_gene_df=as.data.frame(stickle_gene)
gr2_stickle = toGRanges(stickle, format="GTF", header=FALSE, feature="gene")
stickle_promoter=promoters(gr2_stickle, upstream = 1000, downstream = 1000)
stickle_promoter_df=as.data.frame(stickle_promoter)
stickle_promoter_df$gene_id=rownames(stickle_promoter_df)

brain.dmcs.loc.grange=getData(myDiff.sig.brain)[,1:3]
colnames(brain.dmcs.loc.grange)=c("chr", "start", "end")

promoter_brain=annotatePeakInBatch(as(brain.dmcs.loc.grange, "GRanges"), AnnotationData = as(stickle_promoter_df, "GRanges"), output = "inside")
promoter_brain=data.frame(promoter_brain)
promoter_brain=promoter_brain[!is.na(promoter_brain$fromOverlappingOrNearest),]
promoter_brain$loc=paste(promoter_brain$seqnames, promoter_brain$start, sep = "_")

gene_brain=annotatePeakInBatch(as(brain.dmcs.loc.grange, "GRanges"), AnnotationData = stickle_gene, output = "inside")
gene_brain=data.frame(gene_brain)
gene_brain=gene_brain[!is.na(gene_brain$fromOverlappingOrNearest),]
gene_brain$loc=paste(gene_brain$seqnames, gene_brain$start, sep = "_")
gene_brain=gene_brain[!gene_brain$loc %in% promoter_brain$loc,]

gene_brain.anno=rbind(promoter_brain, gene_brain)
length(unique(gene_brain.anno$feature)) # 714

liver.dmcs.loc.grange=getData(myDiff.sig.liver)[,1:3]
colnames(liver.dmcs.loc.grange)=c("chr", "start", "end")

promoter_liver=annotatePeakInBatch(as(liver.dmcs.loc.grange, "GRanges"), AnnotationData = as(stickle_promoter_df, "GRanges"), output = "inside")
promoter_liver=data.frame(promoter_liver)
promoter_liver=promoter_liver[!is.na(promoter_liver$fromOverlappingOrNearest),]
promoter_liver$loc=paste(promoter_liver$seqnames, promoter_liver$start, sep = "_")

gene_liver=annotatePeakInBatch(as(liver.dmcs.loc.grange, "GRanges"), AnnotationData = stickle_gene, output = "inside")
gene_liver=data.frame(gene_liver)
gene_liver=gene_liver[!is.na(gene_liver$fromOverlappingOrNearest),]
gene_liver$loc=paste(gene_liver$seqnames, gene_liver$start, sep = "_")
gene_liver=gene_liver[!gene_liver$loc %in% promoter_liver$loc,]

gene_liver.anno=rbind(promoter_liver, gene_liver)
length(unique(gene_liver.anno$feature)) # 815

# check functions of genes annotated with DMCs in brain and liver-------
library(topGO)
library(biomaRt)
mart=useDataset("gaculeatus_gene_ensembl", useMart("ensembl"))
go_gene_brain.anno=getBM(filters = "external_gene_name", 
                         attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description", "name_1006"), 
                         values = unique(gene_brain.anno$feature), 
                         mart = mart,
                         useCache=F)
go_gene_liver.anno=getBM(filters = "external_gene_name", 
                         attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "description", "name_1006"), 
                         values = unique(gene_liver.anno$feature), 
                         mart = mart,
                         useCache=F)

# overlap between DEG and DMC-associated genes
deg.brain=read.csv("brain_DEG.csv")
colnames(deg.brain)[1]="feature"
length(intersect(unique(gene_brain.anno$feature), unique(deg.brain$feature))) # 70
# 70/length(deg.brain$feature)=0.0335249 0f 2088 DEG
# 70/length(unique(gene_brain.anno$feature))=0.09803922 of DMC-associated genes
dsg.brain=read.csv("DSG_brain.csv", row.names = 1)
colnames(dsg.brain)[1]="feature"
dsg.brain$feature=gsub("gene-", "", dsg.brain$feature)
length(intersect(unique(gene_brain.anno$feature), unique(dsg.brain$feature))) # 61
# 61/length(dsg.brain$feature)=0.03385128 0f 1802 DSG
# 61/length(unique(gene_liver.anno$feature))=0.07484663 of DMC-associated genes

deg.liver=read.csv("liver_DEG.csv")
colnames(deg.liver)[1]="feature"
length(intersect(unique(gene_liver.anno$feature), unique(deg.liver$feature))) # 67
# 67/length(deg.liver$feature)=0.04019196 of 1667 DEG
dsg.liver=read.csv("DSG_liver.csv", row.names = 1)
colnames(dsg.liver)[1]="feature"
dsg.liver$feature=gsub("gene-", "", dsg.liver$feature)
length(intersect(unique(gene_liver.anno$feature), unique(dsg.liver$feature))) # 40
# 40/length(dsg.liver$feature)=0.03338898 0f 1198 DSG

# CpG island enrichment-------
# Use TaJoCGI https://github.com/lucasnell/TaJoCGI to extract CpG island
# Note that to change 'import cyFuns as cgi' to 'import pyFuns as cgi' in in TaJoCGI.py, see the website for more details
# ./TaJoCGI.py -c 8 -o 9spine_CGI GCA_902500615.3.fasta
cgi=read.table("9spine_CGI.bed")

colnames(cgi)=c("seqname", "start", "end")
cgi$name=c(rep("island", nrow(cgi)))
cgi$start=cgi$start+1 # convert 0 based start in bedfile into 1 based start

cgi=as(cgi, "GRanges")

###############################################################
#             Extract CpG island shores
###############################################################
# extract the shore defined by 2000 bp upstream of cpg islands
shore1=flank(cgi, 2000)
# extract the shore defined by 2000 bp downstream of cpg islands
shore2=flank(cgi,2000,FALSE)
# perform intersection and combine the shores where they overlap
shore1_2=reduce(c(shore1,shore2))
# extract the features (ranges) that are present in shores only and not in cgi (ie., shores not overlapping islands)
cpgi_shores=setdiff(shore1_2, cgi)
cpgi_shores$name="shore"
###############################################################
#             Extract CpG island shelves
###############################################################
# extract the shore defined by 4000 bp upstream of cpg islands
shelves1=flank(cgi, 4000)
# extract the shore defined by 2000 bp downstream of cpg islands
shelves2=flank(cgi,4000,FALSE)
# perform intersection and combine the shelves where they overlap
shelves1_2=reduce(c(shelves1,shelves2))
# create a set of ranges consisting CpG Islands, Shores
island_shores=c(cgi, cpgi_shores)
# extract the features (ranges) that are present in shelves only and not in cgi or shores(ie., shelves not overlapping islands or shores)
cpgi_shelves=setdiff(shelves1_2, island_shores)
cpgi_shelves$name="shelf"

library(ChIPpeakAnno)
# dmcs in brain
island_dmcs.brain=annotatePeakInBatch(as(brain.dmcs.loc.grange, "GRanges"), AnnotationData = cgi,
                                      output = "inside")
island_dmcs.brain_df=as.data.frame(island_dmcs.brain)
island_dmcs.brain_df=island_dmcs.brain_df[!is.na(island_dmcs.brain_df$insideFeature),]
nrow(island_dmcs.brain_df) #47
write.csv(island_dmcs.brain_df, "island_dmcs.brain.csv", row.names=F)

shore_dmcs.brain=annotatePeakInBatch(as(brain.dmcs.loc.grange, "GRanges"), AnnotationData = cpgi_shores,
                                     output = "inside")
shore_dmcs.brain_df=as.data.frame(shore_dmcs.brain)
shore_dmcs.brain_df=shore_dmcs.brain_df[!is.na(shore_dmcs.brain_df$insideFeature),]
nrow(shore_dmcs.brain_df) # 354
write.csv(shore_dmcs.brain, "shore_dmcs.brain.csv")

shelf_dmcs.brain=annotatePeakInBatch(as(brain.dmcs.loc.grange, "GRanges"), AnnotationData = cpgi_shelves,
                                     output = "inside")
shelf_dmcs.brain_df=as.data.frame(shelf_dmcs.brain)
shelf_dmcs.brain_df=shelf_dmcs.brain_df[!is.na(shelf_dmcs.brain_df$insideFeature),]
nrow(shelf_dmcs.brain_df) # 206
write.csv(shelf_dmcs.brain, "shelves_dmcs.brain.csv")

sum(nrow(island_dmcs.brain_df), nrow(shore_dmcs.brain_df), nrow(shelf_dmcs.brain_df))
# 607 dmcs in brain locate within island
# nrow(brain.dmcs.loc.grange)-607 = 826 dmcs locate outside island (open sea)

# dmcs in liver
island_dmcs.liver=annotatePeakInBatch(as(liver.dmcs.loc.grange, "GRanges"), AnnotationData = cgi,
                                      output = "inside")
island_dmcs.liver_df=as.data.frame(island_dmcs.liver)
island_dmcs.liver_df=island_dmcs.liver_df[!is.na(island_dmcs.liver_df$insideFeature),]
nrow(island_dmcs.liver_df) # 38
write.csv(island_dmcs.liver_df, "island_dmcs.liver.csv", row.names=F)

shore_dmcs.liver=annotatePeakInBatch(as(liver.dmcs.loc.grange, "GRanges"), AnnotationData = cpgi_shores,
                                     output = "inside")
shore_dmcs.liver_df=as.data.frame(shore_dmcs.liver)
shore_dmcs.liver_df=shore_dmcs.liver_df[!is.na(shore_dmcs.liver_df$insideFeature),]
nrow(shore_dmcs.liver_df) # 404
write.csv(shore_dmcs.liver, "shore_dmcs.liver.csv")

shelf_dmcs.liver=annotatePeakInBatch(as(liver.dmcs.loc.grange, "GRanges"), AnnotationData = cpgi_shelves,
                                     output = "inside")
shelf_dmcs.liver_df=as.data.frame(shelf_dmcs.liver)
shelf_dmcs.liver_df=shelf_dmcs.liver_df[!is.na(shelf_dmcs.liver_df$insideFeature),]
nrow(shelf_dmcs.liver_df) # 253
write.csv(shelf_dmcs.liver, "shelves_dmcs.liver.csv")

sum(nrow(island_dmcs.liver_df), nrow(shore_dmcs.liver_df), nrow(shelf_dmcs.liver_df))
# 695 dmcs in liver locate within island
# nrow(liver.dmcs.loc.grange)-695 = 1037 dmcs locate outside island (open sea)

# plot enrichment of cpgi in brain and liver-------
# plot hyper vs. hypo ecotype dmcs in brain and liver-------
library(ggplot2)
cpgi=data.frame(status=c("CpGi", "CpGi shores", "CpGi shelves", "Open sea",
                         "CpGi", "CpGi shores", "CpGi shelves", "Open sea"),
                tissue=c("brain", "brain", "brain", "brain", "liver", "liver", "liver", "liver"),
                value=c(47, 354, 206, 826,
                        38, 404, 253, 1037))
stack.cpgi=ggplot(cpgi, aes(x=tissue, y=value, fill=status))+
  geom_bar(position="fill", stat="identity", width = 0.5)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position = "top",
        axis.text.x=element_text(colour = "black", size = 12),
        # axis.text.y = element_text(colour = "black"),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(colour = "black"),
        panel.border = element_blank(),
        legend.title = element_blank())+
  scale_fill_manual(values = c("#CCCCCC", "#999999", "#666666", "#333333"))+
  scale_y_continuous(expand = c(0, 0), labels = scales::percent)+
  labs(y="% of DMCs", x="")

stack.cpgi

# plot genomic distribution of ecotype dmcs in brain and liver-------
# for brain
dmcs.brain.1=data.frame(myDiff.sig.brain)
dmcs.brain.1$value=ifelse(dmcs.brain.1$meth.diff>0, "hyper", "hypo")
dmcs.brain.1$value=as.factor(dmcs.brain.1$value)

dmcs.brain.1$chr=as.character(dmcs.brain.1$chr)
library(tidyr)
dmcs.brain.1=separate(data = dmcs.brain.1, col = chr, into = c("1", "2", "chr"), sep = "\\|")
dmcs.brain.1=dmcs.brain.1[,-c(1:2)]

dmcs.brain.1$chr=gsub("CABVRH030000001.1", "LG1", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000002.1", "LG2", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000003.1", "LG3", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000004.1", "LG4", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000005.1", "LG5", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000006.1", "LG6", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000007.1", "LG7", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000008.1", "LG8", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000009.1", "LG9", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000010.1", "LG10", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000011.1", "LG11", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000012.1", "LG12", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000013.1", "LG13", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000014.1", "LG14", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000015.1", "LG15", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000016.1", "LG16", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000017.1", "LG17", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000018.1", "LG18", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000019.1", "LG19", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000020.1", "LG20", dmcs.brain.1$chr)
dmcs.brain.1$chr=gsub("CABVRH030000021.1", "LG21", dmcs.brain.1$chr)

dmcs.brain.1.chr=dmcs.brain.1[grep("LG", dmcs.brain.1$chr),] # 1396 dmcs in brain locate on chromosome, 1433-1396=37 dmcs locate on non-chromosome

library(ggplot2)
scatter_dmcs.brain=ggplot(dmcs.brain.1.chr, aes(colour=chr, x=start, y=meth.diff))+
  geom_point(size=2)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position = "none",
        axis.text.x=element_text(angle = 90),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(colour = "black"),
        panel.border = element_blank())+
  # scale_y_continuous(breaks = seq(-75, 75, 25))+
  geom_hline(yintercept = c(0), linetype="dashed")+
  # geom_jitter()+
  labs(x="",y="change in % methylation")

scatter_dmcs.brain.1=scatter_dmcs.brain+
  facet_grid(.~factor(chr, levels = paste("LG", c(1:21), sep = "")), scales = 'free_x', space = "free_x", switch = 'x')+
  theme(axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(angle = 90),
        panel.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        panel.spacing = unit(0.2,"lines"))+
  scale_x_continuous(expand = c(0, 0))

scatter_dmcs.brain.1

# for liver
dmcs.liver.1=data.frame(myDiff.sig.liver)
dmcs.liver.1$value=ifelse(dmcs.liver.1$meth.diff>0, "hyper", "hypo")
dmcs.liver.1$value=as.factor(dmcs.liver.1$value)

dmcs.liver.1$chr=as.character(dmcs.liver.1$chr)
library(tidyr)
dmcs.liver.1=separate(data = dmcs.liver.1, col = chr, into = c("1", "2", "chr"), sep = "\\|")
dmcs.liver.1=dmcs.liver.1[,-c(1:2)]

dmcs.liver.1$chr=gsub("CABVRH030000001.1", "LG1", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000002.1", "LG2", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000003.1", "LG3", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000004.1", "LG4", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000005.1", "LG5", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000006.1", "LG6", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000007.1", "LG7", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000008.1", "LG8", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000009.1", "LG9", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000010.1", "LG10", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000011.1", "LG11", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000012.1", "LG12", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000013.1", "LG13", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000014.1", "LG14", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000015.1", "LG15", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000016.1", "LG16", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000017.1", "LG17", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000018.1", "LG18", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000019.1", "LG19", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000020.1", "LG20", dmcs.liver.1$chr)
dmcs.liver.1$chr=gsub("CABVRH030000021.1", "LG21", dmcs.liver.1$chr)

dmcs.liver.1.chr=dmcs.liver.1[grep("LG", dmcs.liver.1$chr),] # 1732 dmcs in liver locate on chromosome, 1732-1638=94 dmcs locate on non-chromosome

library(ggplot2)
scatter_dmcs.liver=ggplot(dmcs.liver.1.chr, aes(colour=chr, x=start, y=meth.diff))+
  geom_point(size=2)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position = "none",
        axis.text.x=element_text(angle = 90),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(colour = "black"),
        panel.border = element_blank())+
  # scale_y_continuous(breaks = seq(-75, 75, 25))+
  geom_hline(yintercept = c(0), linetype="dashed")+
  # geom_jitter()+
  labs(x="",y="change in % methylation")

scatter_dmcs.liver.1=scatter_dmcs.liver+
  facet_grid(.~factor(chr, levels = paste("LG", c(1:21), sep = "")), scales = 'free_x', space = "free_x", switch = 'x')+
  theme(axis.text.x = element_blank(),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(angle = 90),
        panel.background = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        panel.spacing = unit(0.2,"lines"))+
  scale_x_continuous(expand = c(0, 0))

scatter_dmcs.liver.1

### analysis enrichment of genomic features of dmcs in brain and liver---------
library(genomation)
# for v7 genome, because liftoff did not output phase info for CDS,
# I used 'gt gff3 -sort -tidy -retainids v7.gff > v7.phased.gff' to add phase info to this assembly
# thus use './gff3ToGenePred v7.phased.gff v7.phased.generpred' and './genePredToBed v7.phased.generpred v7.phased.bed' to convert gff3 to bed12 format
anno=readTranscriptFeatures("/Volumes/Juntao/9spine_wgbs/9spine_genome/liftover/ncbi/v7.phased.bed")
# null genomic feature distribution for all CpGs passed filtering
annotateWithGeneParts(as(unite_norm_10x_CT_GA, "GRanges"),anno)

# percentage of target features overlapping with annotation:
#  (with promoter > exon > intron precedence):
#  promoter   exon      intron   intergenic 
#  10.57      13.46      43.95      32.02 

# for dmcs in brain
annotateWithGeneParts(as(brain.dmcs.loc.grange, "GRanges"),anno)

# percentage of target features overlapping with annotation:
#  (with promoter > exon > intron precedence):
#  promoter    exon     intron    intergenic 
#  15.21       9.98      40.47      34.33 

# G-test
library(DescTools)
observed=c(34.33 , 100-34.33)
expected=c(32.02/100, 1-32.02/100)

GTest(x=observed,
      p=expected)

# promoter: G = 2.0359, X-squared df = 1, p-value = 0.1536
# exon: G = 1.1273, X-squared df = 1, p-value = 0.2884
# intron: G = 0.49485, X-squared df = 1, p-value = 0.4818
# intergenic: G = 0.24218, X-squared df = 1, p-value = 0.6226

# # for dmcs in liver
annotateWithGeneParts(as(liver.dmcs.loc.grange, "GRanges"),anno)

# percentage of target features overlapping with annotation:
#  (with promoter > exon > intron precedence):
#  promoter     exon     intron   intergenic 
#  12.41       8.78      39.78      39.03 

# G-test
library(DescTools)
observed=c(39.03, 100-39.03)
expected=c(32.02/100, 1-32.02/100)

GTest(x=observed,
      p=expected)

# promoter: G = 0.3413, X-squared df = 1, p-value = 0.5591
# exon: G = 2.1062, X-squared df = 1, p-value = 0.1467
# intron: G = 0.71161, X-squared df = 1, p-value = 0.3989
# intergenic: G = 2.1826, X-squared df = 1, p-value = 0.1396

####### codes below are not used at all#############
# Combine all sex, conversion rate, habitat info into the snp data
# snp=read.table("pca/final.012")
# snp=snp[,-1]

# ind=read.table("pca/final.012.indv")
# rownames(snp)=ind$V1

# pos=read.table("pca/final.012.pos")
# colnames(snp)=paste(pos$V1, pos$V2, sep = "_")

# snp_info=data.frame(ecotype=c(rep("marine", 5), rep("marine", 5), rep("freshwater", 4), rep("freshwater", 5)),
#                    population=c(rep("HEL", 5), rep("BOL", 5), rep("ABB", 4), rep("BYN", 5)))
# rownames(snp_info)=ind$V1

# rownames(snp)==rownames(snp_info)
# df_snp=cbind(snp_info, snp)
# str(df_snp)

# df_snp[df_snp==-1]=NA

# snp=snp[!rownames(snp)=="ID_515_ABB",]
#snp_info=snp_info[!rownames(snp_info)=="ID_515_ABB",]

# pca for snp----
# pca.snp=prcomp(snp, center = T)
# summary(pca.snp)

#Define factors
# ecotype.snp=as.factor(snp_info$ecotype)
# population.snp=as.factor(snp_info$population)

#visulize pca.snp results
# g1.snp=ggbiplot(pca.snp, 
#                 obs.scale = 1, 
#                 var.scale = 1,
#                 varname.size = 0,
                # labels = info$ID,
#                 var.axes = FALSE)
# g2.snp=g1.snp+geom_point(aes(colour=ecotype.snp, fill=ecotype.snp, shape=population.snp), size=2)+
#   scale_shape_manual(values = c(21,22,23,24))
# to channge the color +scale_color_manual(values = c("red", "blue"))
# g3.snp=g2.snp + theme_bw()
# print(g3.snp)

### association between snps and methylation in shared dmls between brain and liver----
# using brain as representative due to high correlation 
# in direction of methylation change between brain and liver in shared dmls
meth.count=regionCounts(meth.all, regions = as(unite_norm_5x_dss, "GRanges"))
perc.meth.all=percMethylation(meth.count)

# check if the order of meth count table is the same as info
# colnames(perc.meth.all)==as.character(info$Label)

colnames(perc.meth.all)==unite_norm_5x_dss@sample.ids

cpg=data.frame(perc.meth.all)
rownames(cpg)=paste(unite_norm_5x_dss$chr, unite_norm_5x_dss$start, sep = "_")
colnames(cpg)=unite_norm_5x_dss@sample.ids

same_direction_dmls=read.table("./pca/same_direction_dmls_brain_liver.txt")
same_direction_dmls$loc=paste(same_direction_dmls$V1, same_direction_dmls$V2, sep = "_")

ind_brain=paste(ind$V1, "-BB", sep = "")

cpg_same_dir=cpg[rownames(cpg) %in% same_direction_dmls$loc,colnames(cpg) %in% ind_brain]
colnames(cpg_same_dir)==ind_brain
cpg_same_dir=cpg_same_dir/100

cpg_same_dir$CpGid=rownames(cpg_same_dir)
cpg_same_dir=cpg_same_dir[,c(ncol(cpg_same_dir),1:(ncol(cpg_same_dir)-1))]

library(tidyr)
cpg_same_dir=separate(data = cpg_same_dir, col = CpGid, into = c("1", "2", "CpGid"), sep = "\\|")
cpg_same_dir=cpg_same_dir[,-c(1:2)]

cpg_same_dir$CpGid=gsub("CABVRH030000001.1", "LG1", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000002.1", "LG2", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000003.1", "LG3", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000004.1", "LG4", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000005.1", "LG5", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000006.1", "LG6", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000007.1", "LG7", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000008.1", "LG8", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000009.1", "LG9", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000010.1", "LG10", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000011.1", "LG11", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000013.1", "LG13", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000014.1", "LG14", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000015.1", "LG15", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000016.1", "LG16", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000017.1", "LG17", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000018.1", "LG18", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000019.1", "LG19", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000020.1", "LG20", cpg_same_dir$CpGid)
cpg_same_dir$CpGid=gsub("CABVRH030000021.1", "LG21", cpg_same_dir$CpGid)

write.table(cpg_same_dir, "./pca/cpg.txt", row.names = F, sep = "\t", quote = F)

# exclude snps on sex chromosome, LG12
snps=read.table("./pca/final.012", row.names = 1) # 1796 SNPs

snps_pos=read.table("./pca/final.012.pos")
library(tidyr)
snps_pos=separate(data = snps_pos, col = V1, into = c("1", "2", "chr"), sep = "\\|")
snps_pos=snps_pos[,-c(1:2)]

snps_pos$chr=gsub("CABVRH030000001.1", "LG1", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000002.1", "LG2", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000003.1", "LG3", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000004.1", "LG4", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000005.1", "LG5", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000006.1", "LG6", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000007.1", "LG7", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000008.1", "LG8", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000009.1", "LG9", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000010.1", "LG10", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000011.1", "LG11", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000012.1", "LG12", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000013.1", "LG13", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000014.1", "LG14", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000015.1", "LG15", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000016.1", "LG16", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000017.1", "LG17", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000018.1", "LG18", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000019.1", "LG19", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000020.1", "LG20", snps_pos$chr)
snps_pos$chr=gsub("CABVRH030000021.1", "LG21", snps_pos$chr)

sexchr=grep("LG12", as.character(snps_pos$chr))
snps=snps[,-sexchr]

snps_pos=snps_pos[-sexchr,]
snps_pos=paste(snps_pos$chr, snps_pos$V2, sep = "_")

snps_ind=read.table("./pca/final.012.indv")
ind_brain==paste(snps_ind$V1, "-BB", sep = "")
snps_ind$V1=as.character(snps_ind$V1)
rownames(snps)=snps_ind$V1
colnames(snps)=snps_pos

info_brain=info[info$Tissue=="brain",]

# make sure all orders are same in info, cpg, and snp tables
rownames(snps)==info_brain$ID
paste(rownames(snps), "-BB", sep = "") ==colnames(cpg_same_dir)[2:ncol(cpg_same_dir)]

snps=data.frame(t(snps))
snps$snpid=rownames(snps)
snps=snps[,c(20, 1:19)]
snps[snps==-1]=NA

write.table(snps, "./pca/snps.txt", row.names = F, sep = "\t", quote = F) # 1687 SNPs

library(MatrixEQTL)
base.dir = "./pca/"
useModel = modelLINEAR
SNP_file_name = paste(base.dir, "snps.txt", sep="")
meth_file_name = paste(base.dir, "cpg.txt", sep="")
# covariates_file_name = paste(base.dir, "covariate.txt", sep="") 
output_file_name = paste(base.dir, "output.txt", sep="") 
pvOutputThreshold = 1e-5
errorCovariance = as.numeric()

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;     # read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name)    # 1687 SNPs

meth.qtl = SlicedData$new();
meth.qtl$fileDelimiter = "\t";      # the TAB character
meth.qtl$fileOmitCharacters = "NA"; # denote missing values;
meth.qtl$fileSkipRows = 1;          # one row of column labels
meth.qtl$fileSkipColumns = 1;       # one column of row labels
meth.qtl$fileSliceSize = 1000000;   # read file in pieces of 1,000,000 rows
meth.qtl$LoadFile(meth_file_name)   # 4709 CpGs

me = Matrix_eQTL_engine(
  snps = snps,
  gene = meth.qtl,
  # cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  verbose = TRUE,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = T) 

plot(me)

output=read.table("./pca/output.txt", sep = "\t", skip = 1)
colnames(output)=c("SNP", "CpG", "beta", "t-stat", "p-value")
output$BH=p.adjust(output$`p-value`, method = "BH", n = length(output$`p-value`))
p.adjust=5e-8/4709
output.sig=output[output$BH<p.adjust,] 
nrow(output.sig) # 62 sig mQTL

write.table(output.sig, "./pca/sig.meqtl.txt", sep = "\t", quote = F)

library(tidyr)
output.sig=separate(data = output.sig, col = SNP, into = c("snp_chr", "snp_start"), sep = "_")
output.sig=separate(data = output.sig, col = CpG, into = c("cpg_chr", "cpg_start"), sep = "_")

output.sig$CpG=paste(output.sig$cpg_chr, output.sig$cpg_start, sep = "_")
length(unique(output.sig$CpG)) # 19 unique CpGs in sig meQTLs
output.sig$SNP=paste(output.sig$snp_chr, output.sig$snp_start, sep = "_")
length(unique(output.sig$SNP)) # 57 unique SNPs in sig meQTLs

output.sig$cis_trans=ifelse(output.sig$snp_chr==output.sig$cpg_chr & abs(as.numeric(output.sig$snp_start)-as.numeric(output.sig$cpg_start))<1000000, "cis", "trans")
nrow(output.sig[output.sig$cis_trans=="cis",]) # 2 cis
nrow(output.sig[output.sig$cis_trans=="trans",]) # 60 trans

library(DescTools)
observed=c(2, 60)
expected=c(0.5, 0.5)

GTest(x=observed,
      p=expected)

# G = 68.28, X-squared df = 1, p-value < 2.2e-16

output.sig.snp=output.sig[,1:2]
output.sig.snp$end=output.sig.snp$snp_start
colnames(output.sig.snp)[1:2]=c("chr", "start")
output.sig.snp=dplyr::distinct(output.sig.snp)
output.sig.snp.grange=as(output.sig.snp, "GRanges")

# import liftover annotation from 3spine
liftover.gtf=rtracklayer::import("~/Downloads/9spine_wgbs/9spine_genome/liftover/gaculeatus/v7_3spine.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)

library(tidyr)
liftover.gtf_df=separate(liftover.gtf_df, seqnames, into = c("1", "2", "seqnames"), sep = "\\|")
liftover.gtf_df=liftover.gtf_df[,-c(1:2)]

# modify seqnames name
liftover.gtf_df$seqnames=gsub("CABVRH030000001.1", "LG1", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000002.1", "LG2", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000003.1", "LG3", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000004.1", "LG4", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000005.1", "LG5", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000006.1", "LG6", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000007.1", "LG7", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000008.1", "LG8", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000009.1", "LG9", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000010.1", "LG10", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000011.1", "LG11", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000012.1", "LG12", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000013.1", "LG13", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000014.1", "LG14", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000015.1", "LG15", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000016.1", "LG16", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000017.1", "LG17", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000018.1", "LG18", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000019.1", "LG19", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000020.1", "LG20", liftover.gtf_df$seqnames)
liftover.gtf_df$seqnames=gsub("CABVRH030000021.1", "LG21", liftover.gtf_df$seqnames)

liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",] # 20131 genes
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

promoter_liftover.gtf=promoters(as(liftover.gtf_df_gene, "GRanges"), upstream = 1000, downstream = 1000)
promoter_liftover.gtf_df=as.data.frame(promoter_liftover.gtf)
rownames(promoter_liftover.gtf_df)=c(1:nrow(promoter_liftover.gtf_df))
promoter_liftover.gtf_df_grange=as(promoter_liftover.gtf_df[,1:3], "GRanges")

# find overlaps between all meqtl annotated with genes and liftover annotation
library(GenomicRanges)
pool.promoter=findOverlaps(output.sig.snp.grange, promoter_liftover.gtf_df_grange, type = c("within"), ignore.strand=T)
pool.promoter_df=data.frame(pool.promoter)
pool.promoter_df=promoter_liftover.gtf_df$gene_id[pool.promoter_df$subjectHits]

pool.gene=findOverlaps(output.sig.snp.grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
pool.gene_df=data.frame(pool.gene)
pool.gene_df=liftover.gtf_df_gene$gene_id[pool.gene_df$subjectHits]

meqtl_gene=c(pool.promoter_df, pool.gene_df)
meqtl_gene=unique(meqtl_gene)

# extract gene names
library(biomaRt)
# if encounter error 'Error in bmRequest(request = request, verbose = verbose) : 
# Internal Server Error (HTTP 500).', change host to 'useast.ensembl.org', 'uswest.ensembl.org' or 'asia.ensembl.org'
mart=useDataset("gaculeatus_gene_ensembl", useMart("ensembl"))
meqtl_description=getBM(filters = "ensembl_gene_id", 
              attributes = c("ensembl_gene_id","go_id", "name_1006", "description", "external_gene_name"), 
              values = meqtl_gene, 
              mart = mart,
              useCache = F)

# genes of genomic regions with 3spine ecotype divergence----
# library(tidyr)
# jones=read.csv("./pca/Jones.csv")
# jones=separate(data = jones, col = region, into = c("chr", "regions"), sep = ":")
# jones=separate(data = jones, col = regions, into = c("start", "end"), sep = "-")
# jones$chr=sub("chr", "group", jones$chr)

# hohenlohe=read.csv("./pca/hohenlohe.csv")
# terekhanova=read.csv("./pca/terekhanova.csv")

# divergence_region=rbind(jones, hohenlohe, terekhanova)
# divergence_region=unique(divergence_region) # excluding duplicated regions
# divergence_region_granges=as(divergence_region, "GRanges")

# annotate genes associated with the above snps
library(GenomicFeatures)
library(ChIPpeakAnno)
stickle=makeTxDbFromGFF("~/Downloads/9spine_wgbs/9spine_genome/liftover/gaculeatus/Gasterosteus_aculeatus.BROADS1.105.gtf")

stickle_gene=genes(stickle)
stickle_gene_df=as.data.frame(stickle_gene)
gr2_stickle = toGRanges(stickle, format="GTF", header=FALSE, feature="gene")
stickle_promoter=promoters(gr2_stickle, upstream = 1000, downstream = 1000)
stickle_promoter_df=as.data.frame(stickle_promoter)
stickle_promoter_df$gene_id=rownames(stickle_promoter_df)

# overlaps.promoter_3spine_divergence=annotatePeakInBatch(divergence_region_granges, AnnotationData =  as(stickle_promoter_df, "GRanges"), output = "overlapping")
# overlaps.promoter_3spine_divergence_df=as.data.frame(overlaps.promoter_3spine_divergence)
# overlaps.promoter_3spine_divergence_df=overlaps.promoter_3spine_divergence_df[!is.na(overlaps.promoter_3spine_divergence_df$insideFeature),]

# overlaps.gene_3spine_divergence=annotatePeakInBatch(divergence_region_granges, AnnotationData =  as(stickle_gene_df, "GRanges"), output = "overlapping")
# overlaps.gene_3spine_divergence_df=as.data.frame(overlaps.gene_3spine_divergence)
# overlaps.gene_3spine_divergence_df=overlaps.gene_3spine_divergence_df[!is.na(overlaps.gene_3spine_divergence_df$insideFeature),]

# give a precedence of promoter to gene, and combine info
# remove loci from overlaps.gene_3spine_divergence_df if they are inside promoters
# overlaps.gene_3spine_divergence_df=overlaps.gene_3spine_divergence_df[!overlaps.gene_3spine_divergence_df$loc %in% overlaps.promoter_3spine_divergence_df$loc,]
# divergence_region.anno=c(overlaps.promoter_3spine_divergence_df$feature, overlaps.gene_3spine_divergence_df$feature)
# divergence_region.anno=unique(divergence_region.anno)

# overlap_9spine_3spine=meqtl_gene[meqtl_gene%in%divergence_region.anno]
# length(overlap_9spine_3spine) # 0

# meqtls associated with methylation ecotype divergence in 3spine from Hu et al. 2021 Genetics----
threespine_meqtl=read.table("3spine_meqtl.txt") # 536 SNPs (sig meqtls)
library(tidyr)
threespine_meqtl=separate(data = threespine_meqtl, col = V1, into = c("chr", "start"), sep = ",")
threespine_meqtl$end=threespine_meqtl$start
threespine_meqtl_grange=as(threespine_meqtl, "GRanges")

overlaps.promoter_3spine_meqtl_divergence=annotatePeakInBatch(threespine_meqtl_grange, AnnotationData =  as(stickle_promoter_df, "GRanges"), output = "inside")
overlaps.promoter_3spine_meqtl_divergence_df=as.data.frame(overlaps.promoter_3spine_meqtl_divergence)
overlaps.promoter_3spine_meqtl_divergence_df=overlaps.promoter_3spine_meqtl_divergence_df[!is.na(overlaps.promoter_3spine_meqtl_divergence_df$insideFeature),]

overlaps.gene_3spine_meqtl_divergence=annotatePeakInBatch(divergence_region_granges, AnnotationData =  as(stickle_gene_df, "GRanges"), output = "inside")
overlaps.gene_3spine_meqtl_divergence_df=as.data.frame(overlaps.gene_3spine_meqtl_divergence)
overlaps.gene_3spine_meqtl_divergence_df=overlaps.gene_3spine_meqtl_divergence_df[!is.na(overlaps.gene_3spine_meqtl_divergence_df$insideFeature),]

divergence_3spine_meqtl_region.anno=c(overlaps.promoter_3spine_meqtl_divergence_df$feature, overlaps.gene_3spine_meqtl_divergence_df$feature)
divergence_3spine_meqtl_region.anno=unique(divergence_region.anno)

overlap_9spine_3spine_meqtl=meqtl_gene[meqtl_gene%in%divergence_3spine_meqtl_region.anno]
length(overlap_9spine_3spine_meqtl) # 0


save.image("~/Downloads/9spine_wgbs/pca/pca_9spine.RData")
