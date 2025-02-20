#setwd("E:/mirna")
setwd("/path/to/mirna")

samples=read.csv("samples.csv", stringsAsFactors = F)
data=read.table("miRNAs_expressed.count.txt",header = F,sep = "\t")
colnames(data)=c("miRNA","ID_305_BOL-BA","ID_305_BOL-LA",
                 "ID_307_BOL-BA","ID_307_BOL-LA",
                 "ID_209_HEL-BA","ID_209_HEL-LA",
                 "ID_315_BOL-BA","ID_315_BOL-LA",
                 "ID_225_HEL-BA","ID_225_HEL-LA",
                 "ID_451_ABB-BA","ID_451_ABB-LA",
                 "ID_271_HEL-BA","ID_271_HEL-LA",
                 "ID_675_BYN-BA","ID_675_BYN-LA",
                 "ID_579_BYN-BA","ID_579_BYN-LA",
                 "ID_683_BYN-BA","ID_683_BYN-LA",
                 "ID_493_ABB-BA","ID_493_ABB-LA",
                 "ID_494_ABB-BA","ID_494_ABB-LA",
                 "ID_697_BYN-BA","ID_697_BYN-LA")
row.names(data)=data$miRNA
data$miRNA=NULL
#order
order(colnames(data))
data=data[,order(colnames(data))]
#filter
data_filter=subset(data,apply(data,1, function(x) sum(x > 1))>= 26)

#118 miRNA for general miranda####
write.csv(t(data_filter),"./118miRNAs/miRNA_filter.csv",col.names = T,row.names = T)

data_filter=t(data_filter)
#dds all
library(DESeq2)
dds_all=DESeqDataSetFromMatrix(countData = t(data_filter), 
                               colData = samples,
                               design = ~ecotype+tissue)
dds_all_rlog=rlog(dds_all,blind = F)
dds_all_rlog_df=data.frame(dds_all_rlog@assays@data@listData)
write.csv(dds_all_rlog_df,"mirna_rlog.csv",col.names = T,row.names = T)
dds_all_rlog_df=t(dds_all_rlog_df)

pca=prcomp(dds_all_rlog_df,center = T)
pca_df=data.frame(pca$x,ecotype=samples$ecotype,tissue=samples$tissue)
summ=summary(pca)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")


#
##add analysis for screen plot
explained_variance_ratio <- pca$sdev^2 / sum(pca$sdev^2)

explained=qplot(x = 1:length(explained_variance_ratio), 
      y = explained_variance_ratio, 
      geom = c("line", "point")) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained") +
  scale_x_continuous(breaks = seq(1, length(explained_variance_ratio), by = 2)) + 
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) + 
  #scale_y_continuous(breaks = seq(0, max(explained_variance_ratio), by = 0.1)) +
  theme_minimal()+
  theme_classic()+theme(text = element_text(size = 15))
#

tiff("Screenplot_miRNA_luo.tiff", width = 6, height = 5, units = "in", res = 300)
print(explained)
dev.off()


library(ggplot2)
#without label
mirna_pca=ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue))+ geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))
ggsave("mirna_pca.tiff",width = 5,height = 4,dpi=300)


tiff("mirna_pca_luo.tiff", width = 6, height = 5, units = "in", res = 300)
print(mirna_pca)
dev.off()

write.csv(pca_df,file = "miRNA_pca_df.csv")# luo for QDA analysis

#with label
ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue,label=rownames(pca_df)))+ geom_point()+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))+
  geom_text_repel(box.padding = 0.5,
                  max.overlaps = 30,
                  size=2,
                  segment.size=0.3)
ggsave("mirna_pca_label.tiff",width = 5,height = 4,dpi=300)

dds_brain=grep("BA",colnames(dds_all))
dds_brain=dds_all[,dds_brain]
colData(dds_brain)
design(dds_brain)=formula(~ecotype)
dds_brain=DESeq(dds_brain)
res_brain=results(dds_brain,contrast = c("ecotype","marine","freshwater"))
resSig_brain=res_brain[which(res_brain$padj<0.05),]

select_brain <- rownames(resSig_brain)
nt_brain <- normTransform(dds_brain, f = log2, pc = 1) 
log2.norm.brain <- assay(nt_brain)
log2.norm.brain=log2.norm.brain[select_brain,]


col_brain=data.frame(sample=rep(c("marine","freshwater"),c(6,7)))
row.names(col_brain)=colnames(log2.norm.brain)
colnames(col_brain)=c("ecotype")
ann_colors=list(ecotype=c(marine="orange",freshwater="light blue"))

library(pheatmap)
library("RColorBrewer")
png("./brain_mirna_heatmap.png", res=300, height=4, width=8, units="in")
pheatmap(log2.norm.brain, 
         annotation_col=col_brain,
         angle_col = 45,
         cluster_rows=T, 
         cluster_cols=TRUE,
         scale = "row",
         show_colnames = F,
         show_rownames = F,
         annotation_colors = ann_colors,
         border_color = F,
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4,width = 6,filename = "E:/mirna/brain_mirna_heatmap.tiff")

#liver
dds_liver=grep("LA",colnames(dds_all))
dds_liver=dds_all[,dds_liver]
colData(dds_liver)
design(dds_liver)=formula(~ecotype)
dds_liver=DESeq(dds_liver)
res_liver=results(dds_liver,contrast = c("ecotype","marine","freshwater"))
resSig_liver=res_liver[which(res_liver$padj<0.05),]

select_liver <- rownames(resSig_liver)
nt_liver <- normTransform(dds_liver, f = log2, pc = 1) 
log2.norm.liver <- assay(nt_liver)
log2.norm.liver=log2.norm.liver[select_liver,]

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
         scale = "row",
         show_colnames = F,
         show_rownames = F,
         annotation_colors = ann_colors,
         border_color = F,
         clustering_method = "complete",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4,width = 6,filename = "E:/mirna/liver_mirna_heatmap.tiff")

