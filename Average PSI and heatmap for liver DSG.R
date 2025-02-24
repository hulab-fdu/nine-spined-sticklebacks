
setwd("/path/to/AS/dsg/liver")
library(reshape2)
#A3SS
A3SS=read.table("A3SS.MATS.JC.txt",header=F)
A3SS=A3SS[-1,]
colnames(A3SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC1","SJC1","IJC2","SJC2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
A3SS=transform(A3SS,IJC1=colsplit(IJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A3SS=transform(A3SS,IJC2=colsplit(IJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A3SS=transform(A3SS,SJC1=colsplit(SJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A3SS=transform(A3SS,SJC2=colsplit(SJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A3SS$IJC_sum=rowSums(A3SS$IJC1)+rowSums(A3SS$IJC2)
A3SS$SJC_sum=rowSums(A3SS$SJC1)+rowSums(A3SS$SJC2)
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
A3SS_fdr=A3SS_count_over_20[A3SS_count_over_20$FDR<0.05,]
write.csv(A3SS_count_over_20,file="A3SS_all_count_over_20.csv")

#A5SS
A5SS=read.table("A5SS.MATS.JC.txt",header=F)
A5SS=A5SS[-1,]
colnames(A5SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC1","SJC1","IJC2","SJC2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
A5SS=transform(A5SS,IJC1=colsplit(IJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A5SS=transform(A5SS,IJC2=colsplit(IJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A5SS=transform(A5SS,SJC1=colsplit(SJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A5SS=transform(A5SS,SJC2=colsplit(SJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
A5SS$IJC_sum=rowSums(A5SS$IJC1)+rowSums(A5SS$IJC2)
A5SS$SJC_sum=rowSums(A5SS$SJC1)+rowSums(A5SS$SJC2)
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
A5SS_fdr=A5SS_count_over_20[A5SS_count_over_20$FDR<0.05,]
write.csv(A5SS_count_over_20,file="A5SS_all_count_over_20.csv")

#SE
SE=read.table("SE.MATS.JC.txt",header=F)
SE=SE[-1,]
colnames(SE)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC1","SJC1","IJC2","SJC2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
SE=transform(SE,IJC1=colsplit(IJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
SE=transform(SE,IJC2=colsplit(IJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
SE=transform(SE,SJC1=colsplit(SJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
SE=transform(SE,SJC2=colsplit(SJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
SE=transform(SE,IncLevel1=colsplit(IncLevel1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
SE$IJC_sum=rowSums(SE$IJC1)+rowSums(SE$IJC2)
SE$SJC_sum=rowSums(SE$SJC1)+rowSums(SE$SJC2)
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
SE_fdr=SE_count_over_20[SE_count_over_20$FDR<0.05,]
write.csv(SE_count_over_20,file="SE_all_count_over_20.csv")

#MXE
MXE=read.table("MXE.MATS.JC.txt",header=F)
MXE=MXE[-1,]
colnames(MXE)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2stExonStart_0base","2stExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC1","SJC1","IJC2","SJC2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
MXE=transform(MXE,IJC1=colsplit(IJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
MXE=transform(MXE,IJC2=colsplit(IJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
MXE=transform(MXE,SJC1=colsplit(SJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
MXE=transform(MXE,SJC2=colsplit(SJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
MXE$IJC_sum=rowSums(MXE$IJC1)+rowSums(MXE$IJC2)
MXE$SJC_sum=rowSums(MXE$SJC1)+rowSums(MXE$SJC2)
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
MXE_fdr=MXE_count_over_20[MXE_count_over_20$FDR<0.05,]
write.csv(MXE_count_over_20,file="MXE_all_count_over_20.csv")

#RI
RI=read.table("RI.MATS.JC.txt",header=F)
RI=RI[-1,]
colnames(RI)=c("ID","GeneID","geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC1","SJC1","IJC2","SJC2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
RI=transform(RI,IJC1=colsplit(IJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
RI=transform(RI,IJC2=colsplit(IJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
RI=transform(RI,SJC1=colsplit(SJC1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
RI=transform(RI,SJC2=colsplit(SJC2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
RI=transform(RI,IncLevel1=colsplit(IncLevel1,",",names = c('ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA','ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,",",names = c('ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA','ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA')))
RI$IJC_sum=rowSums(RI$IJC1)+rowSums(RI$IJC2)
RI$SJC_sum=rowSums(RI$SJC1)+rowSums(RI$SJC2)
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
RI_fdr=RI_count_over_20[RI_count_over_20$FDR<0.05,]
write.csv(RI_count_over_20,file="RI_all_count_over_20.csv")

A3SS_fdr$Type=rep("A3SS",nrow(A3SS_fdr))
A5SS_fdr$Type=rep("A5SS",nrow(A5SS_fdr))
MXE_fdr$Type=rep("MXE",nrow(MXE_fdr))
SE_fdr$Type=rep("SE",nrow(SE_fdr))
RI_fdr$Type=rep("RI",nrow(RI_fdr))


#id,fdr,lncleveldiff,type

risebind=rbind(as.matrix(RI_fdr[,c(2,20,23,26)]),as.matrix(SE_fdr[,c(2,20,23,26)]))
risebind_df=as.data.frame(risebind)
a35bind=rbind(as.matrix(A3SS_fdr[,c(2,20,23,26)]),as.matrix(A5SS_fdr[,c(2,20,23,26)]))
DSG=rbind(risebind,a35bind)
DSG=as.data.frame(rbind(DSG,as.matrix(MXE_fdr[,c(2,22,25,28)])))
DSG$absIncleveldiff=abs(as.numeric(DSG$IncLevelDifference))

# The DS event with the largest absIncleveldifference of a DSG represent the differential splicing level of the gene
DSG_unique_gene=DSG[order(-DSG[,'absIncleveldiff']),]
DSG_unique_gene_du=DSG_unique_gene[!duplicated(DSG_unique_gene$GeneID),]
write.csv(DSG_unique_gene_du,file = "DSG_liver.csv")


#splicing events by genes (no fdr, only count over 20)
risebind2=rbind(as.matrix(RI_count_over_20[,2]),as.matrix(SE_count_over_20[,2]))
risebind2_df=as.data.frame(risebind2)
a35bind2=rbind(as.matrix(A3SS_count_over_20[,2]),as.matrix(A5SS_count_over_20[,2]))
splicing_events=rbind(risebind2,a35bind2)
splicing_events=as.data.frame(rbind(splicing_events,as.matrix(MXE_count_over_20[,2])))
write.csv(splicing_events,"liver_events.csv")
SG=splicing_events[!duplicated(splicing_events$V1),]




#barplot of DS events
library(ggplot2)
df_events=read.csv("plot_liver_events.csv")#summary table of each event: number of event( count over 20 and FDR < 0.05 )
colnames(df_events)=c("event","count_over20","fdr")
#add 'mapping=' if error in 'mapping must created by aes()'

library(ggforce)
ggplot()+
  geom_bar(df_events,mapping=aes(event,count_over20),stat="identity",fill="grey")+geom_bar(df_events,mapping=aes(event,fdr,fill=event),stat="identity",color=NA)+
  scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02"))+
  scale_x_discrete(label=c("A3SS","A5SS","MXE","RI","SE"))+
  labs(x="Splicing events",y="Number of events")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 15))+
  guides(fill=F)+
  facet_zoom(ylim=c(0,1200),
             zoom.size = 1)
ggsave("liver_events.tiff",width = 6,height = 4,dpi = 300)


#figs2
RI_boxplot=cbind(RI_fdr$IncLevel1,RI_fdr$IncLevel2)
RI_boxplot$ave1=rowMeans(RI_fdr$IncLevel1)
RI_boxplot$ave2=rowMeans(RI_fdr$IncLevel2)
RI_boxplot_ave=RI_boxplot[,c("ave1","ave2")]
library(reshape2)
RI_boxplot_ave=melt(RI_boxplot_ave)
RI_boxplot_ave$ecotype=rep(c("marine","freshwater"),c(1030,1030))
RI_boxplot_ave$variable=NULL

library(ggplot2)
ggplot(RI_boxplot_ave,aes(x=ecotype,y=value,fill=ecotype))+geom_boxplot()+
  labs(x=NULL,y="Average PSI")+guides(fill=F)+
  scale_fill_manual(values = c("light blue","yellow"))+
  theme(axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15))+theme_bw()
ggsave("D:/rmats/liver/RI_liver_boxplot.tiff",width = 4,height = 3,dpi = 300)

t.test(RI_boxplot$ave1,RI_boxplot$ave2,paired = T)#p-value=0.932
kruskal.test(value ~ ecotype, data = RI_boxplot_ave)#Kruskal-Wallis chi-squared = 4.6762, df = 1, p-value = 0.03058
summary(as.numeric(RI$IncLevelDifference))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.000000 -0.018000  0.001000  0.001724  0.027000  1.000000

# round 1 conment: luo revised to geom_violin()####
violin=ggplot(RI_boxplot_ave,aes(x=ecotype,y=value,fill=ecotype))+geom_violin()+
  #geom_jitter(shape=21,color="grey",aes(fill=ecotype),position = position_jitter(width = 0.1),alpha = 0.2)+
  geom_boxplot(outlier.shape=NA,width=0.1,position=position_dodge(0.2),alpha = 0.1,color="grey30")+
  labs(x=NULL,y="Average PSI")+guides(fill=F)+
  scale_fill_manual(values = c("light blue","yellow"))+theme_bw()+
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 15))


tiff("RI_liver_boxplot_violin.tiff", width = 6, height = 4, units = "in", res = 300)
print(violin)
dev.off()

#dsg heatmap
samples=read.csv("/path/to/dsg/samples.csv", stringsAsFactors = F)
count=read.table("/path/to/dsg/matrix_gene.txt",header = F,sep = "\t")
colnames(count)=c("sample","gene_id","count")
library(reshape2)
data=dcast(count,gene_id~sample,value.var = "count")#make count table
rownames(data)=data$gene_id
data=data[,-1]
noint = rownames(data) %in% c("__alignment_not_unique","__ambiguous","__no_feature", "__not_aligned","__too_low_aQual")
counts=data[!noint,]

#liver dsg read counts
liver_dsg=substring(DSG_unique_gene_du$GeneID,6)
liver_dsg_count=counts[liver_dsg,]
LA=grep("LA",colnames(liver_dsg_count))
liver_dsg_count=liver_dsg_count[,LA]
library(DESeq2)
dds=DESeqDataSetFromMatrix(countData = counts,
                           colData = samples,
                           design = ~ecotype+tissue)
dds_liver=grep("LA",colnames(dds))
dds_liver=dds[,dds_liver]
colData(dds_liver)
design(dds_liver)=formula(~ecotype)
dds_liver=DESeq(dds_liver)

select_liver <- rownames(liver_dsg_count)
nt_liver <- normTransform(dds_liver, f = log2, pc = 1) 
log2.norm.liver <- assay(nt_liver)
log2.norm.liver=log2.norm.liver[select_liver,]
#heatmap liver
col_liver=data.frame(sample=rep(c("marine","freshwater"),c(6,7)))
row.names(col_liver)=colnames(log2.norm.liver)
colnames(col_liver)=c("ecotype")
ann_colors=list(ecotype=c(marine="orange",freshwater="light blue"))
library(pheatmap)
library("RColorBrewer")
pheatmap(log2.norm.liver, 
         annotation_col=col_liver,
         angle_col = 45,
         cluster_rows=T, 
         cluster_cols=TRUE,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",#ͼ��Ҫд
         clustering_method="complete",
         scale = "row",
         show_colnames = F,
         show_rownames = F,
         annotation_colors = ann_colors,
         border_color = F,
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4,width = 6,filename = "D:/rmats/liver/liver_dsg_heatmap.tiff"
)
