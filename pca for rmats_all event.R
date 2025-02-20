setwd("path/to/rmats_all")
library(reshape2)
library(tidyr)
#A3SS
A3SS=read.table("A3SS.MATS.JC.txt",header=F)
colnames(A3SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A3SS=transform(A3SS,IJC=colsplit(IJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A3SS=transform(A3SS,SJC=colsplit(SJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A3SS=transform(A3SS,IncLevel=colsplit(IncLevel,",",
                                      names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                                'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                                'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                                'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                                'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                                'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                                'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                                'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))

A3SS$IJC_sum=rowSums(A3SS$IJC)
A3SS$SJC_sum=rowSums(A3SS$SJC)
A3SS$IncLevel_mean=rowMeans(A3SS$IncLevel,na.rm = TRUE)
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]

#A5SS
A5SS=read.table("A5SS.MATS.JC.txt",header=F)
colnames(A5SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A5SS=transform(A5SS,IJC=colsplit(IJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A5SS=transform(A5SS,SJC=colsplit(SJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
A5SS=transform(A5SS,IncLevel=colsplit(IncLevel,",",
                                      names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                                'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                                'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                                'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                                'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                                'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                                'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                                'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))

A5SS$IJC_sum=rowSums(A5SS$IJC)
A5SS$SJC_sum=rowSums(A5SS$SJC)
A5SS$IncLevel_mean=rowMeans(A5SS$IncLevel,na.rm = TRUE)
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]

#SE
SE=read.table("SE.MATS.JC.txt",header=F)
colnames(SE)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
SE=transform(SE,IJC=colsplit(IJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
SE=transform(SE,SJC=colsplit(SJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
SE=transform(SE,IncLevel=colsplit(IncLevel,",",
                                      names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                                'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                                'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                                'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                                'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                                'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                                'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                                'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))

SE$IJC_sum=rowSums(SE$IJC)
SE$SJC_sum=rowSums(SE$SJC)
SE$IncLevel_mean=rowMeans(SE$IncLevel,na.rm = TRUE)
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]

#MXE
MXE=read.table("MXE.MATS.JC.txt",header=F)
colnames(MXE)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2stExonStart_0base","2stExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
MXE=transform(MXE,IJC=colsplit(IJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
MXE=transform(MXE,SJC=colsplit(SJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
MXE=transform(MXE,IncLevel=colsplit(IncLevel,",",
                                      names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                                'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                                'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                                'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                                'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                                'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                                'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                                'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))

MXE$IJC_sum=rowSums(MXE$IJC)
MXE$SJC_sum=rowSums(MXE$SJC)
MXE$IncLevel_mean=rowMeans(MXE$IncLevel,na.rm = TRUE)
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]

#RI
RI=read.table("RI.MATS.JC.txt",header=F)
colnames(RI)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
RI=transform(RI,IJC=colsplit(IJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
RI=transform(RI,SJC=colsplit(SJC,",",
                                 names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                           'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                           'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                           'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                           'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                           'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                           'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                           'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))
RI=transform(RI,IncLevel=colsplit(IncLevel,",",
                                      names = c('ID_579_BYN-BA','ID_675_BYN-BA','ID_683_BYN-BA','ID_697_BYN-BA',
                                                'ID_451_ABB-BA','ID_493_ABB-BA','ID_494_ABB-BA',
                                                'ID_209_HEL-BA','ID_225_HEL-BA','ID_271_HEL-BA',
                                                'ID_305_BOL-BA','ID_307_BOL-BA','ID_315_BOL-BA',
                                                'ID_579_BYN-LA','ID_675_BYN-LA','ID_683_BYN-LA','ID_697_BYN-LA',
                                                'ID_451_ABB-LA','ID_493_ABB-LA','ID_494_ABB-LA',
                                                'ID_209_HEL-LA','ID_225_HEL-LA','ID_271_HEL-LA',
                                                'ID_305_BOL-LA','ID_307_BOL-LA','ID_315_BOL-LA')))

RI$IJC_sum=rowSums(RI$IJC)
RI$SJC_sum=rowSums(RI$SJC)
RI$IncLevel_mean=rowMeans(RI$IncLevel,na.rm = TRUE)
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]

genes_A3SS=unique(as.data.frame(A3SS_count_over_20[,2]))
colnames(genes_A3SS)="gene"
genes_A5SS=unique(as.data.frame(A5SS_count_over_20[,2]))
colnames(genes_A5SS)="gene"
genes_SE=unique(as.data.frame(SE_count_over_20[,2]))
colnames(genes_SE)="gene"
genes_MXE=unique(as.data.frame(MXE_count_over_20[,2]))
colnames(genes_MXE)="gene"
genes_RI=unique(as.data.frame(RI_count_over_20[,2]))
colnames(genes_RI)="gene"
genes_all_clean=rbind(genes_A3SS,genes_A5SS)
genes_all_clean=rbind(genes_all_clean,genes_MXE)
genes_all_clean=rbind(genes_all_clean,genes_SE)                 
genes_all_clean=rbind(genes_all_clean,genes_RI) 
genes_all_clean=unique(genes_all_clean)#10467 genes,57325 events

write.csv(A3SS_count_over_20,file = "A3SS_count_over_20.csv")
write.csv(A5SS_count_over_20,file = "A5SS_count_over_20.csv")
write.csv(SE_count_over_20,file = "SE_count_over_20.csv")
write.csv(MXE_count_over_20,file = "MXE_count_over_20.csv")
write.csv(RI_count_over_20,file = "RI_count_over_20.csv")

A3SS_count_over20=read.csv(file = "A3SS_count_over_20.csv")
A5SS_count_over20=read.csv(file = "A5SS_count_over_20.csv")
SE_count_over20=read.csv(file = "SE_count_over_20.csv")
MXE_count_over20=read.csv(file = "MXE_count_over_20.csv")
RI_count_over20=read.csv(file = "RI_count_over_20.csv")

A3SS_count_over20$GeneID <- paste("A3SS", substring(A3SS_count_over20$GeneID,6), sep="_")
A5SS_count_over20$GeneID <- paste("A5SS", substring(A5SS_count_over20$GeneID,6), sep="_")
MXE_count_over20$GeneID <- paste("MXE", substring(MXE_count_over20$GeneID,6), sep="_")
RI_count_over20$GeneID <- paste("RI", substring(RI_count_over20$GeneID,6), sep="_")
SE_count_over20$GeneID <- paste("SE", substring(SE_count_over20$GeneID,6), sep="_")
all_count_over20=rbind(A3SS_count_over20[70:95],A5SS_count_over20[70:95],MXE_count_over20[72:97],SE_count_over20[70:95],RI_count_over20[70:95])
colnames(all_count_over20)=substring(colnames(all_count_over20),10)
all_count_over20_nona=all_count_over20[complete.cases(all_count_over20),]

#pca
splicing_pca=prcomp(t(all_count_over20_nona))
summ=summary(splicing_pca)
splicing_pca_df=as.data.frame(splicing_pca$x)
write.csv(splicing_pca_df,"/path/to/splicing_pca.csv")
splicing_pca_df$ecotype=c(rep("freshwater",7),rep("marine",6),rep("freshwater",7),rep("marine",6))
splicing_pca_df$tissue=c(rep("brain",13),rep("liver",13))
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")
library(ggplot2)
#without label
rmats_pca=ggplot(splicing_pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue))+ geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))
ggsave("rmats_pca.tiff",width = 5,height = 4,dpi=300)

tiff("rmats_pca_.tiff", width = 6, height = 5, units = "in", res = 300)
print(rmats_pca)
dev.off()

write.csv(splicing_pca_df,file = "splicing_pca_df.csv")#luo for QDA analysis

#with label
ggplot(splicing_pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue,label=rownames(splicing_pca_df)))+ geom_point()+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))+
  geom_text_repel(box.padding = 0.6,
                  max.overlaps = 30,
                  size=2,
                  segment.size=0.3)
ggsave("rmats_pca_label.tiff",width = 5,height = 4,dpi=300)

###analysis for screen plot
#
explained_variance_ratio <- splicing_pca$sdev^2 / sum(splicing_pca$sdev^2)

explained=qplot(x = 1:length(explained_variance_ratio), 
      y = explained_variance_ratio, 
      geom = c("line", "point")) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained") +
  scale_x_continuous(breaks = seq(1, length(explained_variance_ratio), by = 2)) + 
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) + 
  #scale_y_continuous(breaks = seq(0, max(explained_variance_ratio), by = 0.1)) +
  theme_minimal()+
  theme_classic()+theme(text = element_text(size = 15))

tiff("Screenplot_rmat.tiff", width = 6, height = 5, units = "in", res = 300)
print(explained)
dev.off()

#