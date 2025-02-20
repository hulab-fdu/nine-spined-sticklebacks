setwd("/path/to/deg")
library(DESeq2)
library(reshape2)

samples=read.csv("samples.csv", stringsAsFactors = F)
count=read.table("matrix_gene.txt",header = F,sep = "\t")
colnames(count)=c("sample","gene_id","count")
data=dcast(count,gene_id~sample,value.var = "count")#make count table
rownames(data)=data$gene_id
data=data[,-1]
noint = rownames(data) %in% c("__alignment_not_unique","__ambiguous","__no_feature", "__not_aligned","__too_low_aQual")
counts=data[!noint,]

#counts_raw=counts
#counts_raw$genename=rownames(counts_raw)

#dds_raw count: used to analysis the pearson correlation between miRNA and mRNA in pearson.R ####
dds_raw=DESeqDataSetFromMatrix(countData = counts,
                               colData = samples,
                               design = ~ecotype+tissue)
dds_raw_rlog=rlog(dds_raw,blind = F)
dds_raw_rlog_df=data.frame(dds_raw_rlog@assays@data@listData)
write.csv(dds_raw_rlog_df,"raw_gene_rlog.csv",col.names = T,row.names = T)


#filter
counts_filter=subset(counts,apply(counts,1, function(x) sum(x > 1))>= 26)
write.csv(counts_filter,"count_filter.csv")

counts_filter_liver=counts_filter[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]#luo
library(edgeR)
test1=DGEList(count=counts_filter_liver)
test2=calcNormFactors(test1)
#dds_all
dds_all=DESeqDataSetFromMatrix(countData = counts_filter,
                               colData = samples,
                               design = ~ecotype+tissue)
dds_all_rlog=rlog(dds_all,blind = F)
dds_all_rlog_df=data.frame(dds_all_rlog@assays@data@listData)
write.csv(dds_all_rlog_df,"gene_rlog.csv",col.names = T,row.names = T)
dds_all_rlog_df=t(dds_all_rlog_df)

#liver_pca
dds_liver_rlog_df=dds_all_rlog_df[c(2,4,6,8,10,12,14,16,18,20,22,24,26),]
pca_l=prcomp(dds_liver_rlog_df,center = T)
pca_l_df=data.frame(pca_l$x,ecotype=rep(c("marine","freshwater"),c(6,7)))
summ=summary(pca_l)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

library(ggplot2)
library(ggrepel)
#pca_l without  label
liver_pca=ggplot(pca_l_df,aes(x=PC1,y=PC2,color=ecotype))+ 
  geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+
  scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))


tiff("liver_pca.tiff", width = 6, height = 5, units = "in", res = 300)
print(liver_pca)
dev.off()

#brain_pca
dds_brain_rlog_df=dds_all_rlog_df[c(1,3,5,7,9,11,13,15,17,19,21,23,25),]
pca_b=prcomp(dds_brain_rlog_df,center = T)
pca_b_df=data.frame(pca_b$x,ecotype=rep(c("marine","freshwater"),c(6,7)))
summ=summary(pca_b)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

#pca_b without label
brain_pca=ggplot(pca_b_df,aes(x=PC1,y=PC2,color=ecotype))+ 
  geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+
  scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))


tiff("brain_pca.tiff", width = 6, height = 5, units = "in", res = 300)
print(brain_pca)
dev.off()


#pca for all
pca=prcomp(dds_all_rlog_df,center = T)
pca_df=data.frame(pca$x,ecotype=samples$ecotype,tissue=samples$tissue)
summ=summary(pca)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

##round 1 comment:  luo add analysis for screen plot

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

tiff("Screenplot_gene_luo.tiff", width = 6, height = 5, units = "in", res = 300)
print(explained)
dev.off()


library(ggplot2)
library(ggrepel)
#pca for two tissues without label
gene_pca=
  ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue))+ geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+guides(shape=guide_legend(override.aes = list(shape=c(1:2))))

tiff("gene_pca_luo.tiff", width = 6, height = 5, units = "in", res = 300)
print(gene_pca)
dev.off()

#pca with label in supp
ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue,label=rownames(pca_df)))+ geom_point()+labs(x=xlab,y=ylab)+scale_color_manual(values = c("light blue","orange"))+
  theme_classic()+theme(text = element_text(size = 15))+guides(shape=guide_legend(override.aes = list(shape=c(1:2))))+
  geom_text_repel(box.padding = 0.4,
                  max.overlaps = 30,
                  size=2,
                  segment.size=0.3)
ggsave("gene_pca_label.tiff",width = 5,height = 4,dpi=300)

write.csv(pca_df,file="pca_df.csv")#  for QDA

#brain
dds_brain=grep("BA",colnames(dds_all))
dds_brain=dds_all[,dds_brain]
colData(dds_brain)
design(dds_brain)=formula(~ecotype)
dds_brain=DESeq(dds_brain)

#significant:brain
res_brain=results(dds_brain,contrast = c("ecotype","marine","freshwater"))
resSig_brain=res_brain[which(res_brain$padj<0.05),]
resSig_brain=data.frame(resSig_brain)
write.csv(resSig_brain,file = "./png/brain_DEG.csv")

write.csv(resSig_brain,file = "./brain_DEG.csv")#luo

select_brain <- rownames(resSig_brain)
nt_brain <- normTransform(dds_brain, f = log2, pc = 1) 
log2.norm.brain <- assay(nt_brain)
log2.norm.brain=log2.norm.brain[select_brain,]

#heatmap brain
col_brain=data.frame(sample=rep(c("marine","freshwater"),c(6,7)))
row.names(col_brain)=colnames(log2.norm.brain)
colnames(col_brain)=c("ecotype")
ann_colors=list(ecotype=c(marine="orange",freshwater="light blue"))
library(pheatmap)
library("RColorBrewer")
pheatmap(log2.norm.brain, 
         annotation_col=col_brain,
         angle_col = 45,
         cluster_rows=T, 
         cluster_cols=TRUE,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         scale = "row",
         show_colnames = T,
         show_rownames = F,
         annotation_colors = ann_colors,
         border_color = F,
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4,width = 6)#,filename = "brain_DEG_heatmap.tiff")

#liver
dds_liver=grep("LA",colnames(dds_all))
dds_liver=dds_all[,dds_liver]
colData(dds_liver)
design(dds_liver)=formula(~ecotype)
dds_liver=DESeq(dds_liver)

#significant:liver
res_liver=results(dds_liver,contrast = c("ecotype","marine","freshwater"))
resSig_liver=res_liver[which(res_liver$padj<0.05),]
resSig_liver=data.frame(resSig_liver)
write.csv(resSig_liver,file = "./png/liver_DEG.csv")

write.csv(resSig_liver,file = "./liver_DEG.csv")#luo

select_liver <- rownames(resSig_liver)
nt_liver <- normTransform(dds_liver, f = log2, pc = 1) 
log2.norm.liver <- assay(nt_liver)
log2.norm.liver=log2.norm.liver[select_liver,]

#heatmap
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
         clustering_distance_cols="euclidean",#Í¼ÀýÒªÐ´
         clustering_method="complete",
         scale = "row",
         show_colnames = T,
         show_rownames = F,
         annotation_colors = ann_colors,
         border_color = F,
         color = colorRampPalette(colors = c("#4169E1","#FEFBFB","#CD2626"))(8),
         height = 4,width = 6)#,filename = "liver_DEG_heatmap.tiff")



#fig3 venn
#deg,dsg,demir target genes(not DEG),dmc gene
brain_dmc=read.csv("./dmc_associated_genes_brain.csv",header = T)
brain_target=read.csv("/path/to/mirna/pearson/b_pearson.csv",header = T)
brain_venn_list=list(brain_DSG$GeneID,brain_DEG$X,brain_target$mR,brain_dmc$x)
library(VennDiagram)
library(venn)
library(ggplot2)
venn.diagram(brain_venn_list,
             category.names = c("","","",""),
             cat.cex=1.3,
             cex=1.5,
             filename = "brain_four.tiff")

liver_dmc=read.csv("./dmc_associated_genes_liver.csv",header = T)
liver_target=read.csv("/path/to/mirna/pearson/l_pearson.csv",header = T)
liver_venn_list=list(liver_DSG$GeneID,liver_DEG$X,liver_target$mR,liver_dmc$x)
venn.diagram(liver_venn_list,
     category.names = c("","","",""),
     cat.cex=1.3,
     cex=1.5,
     filename = "liver_four.tiff")

