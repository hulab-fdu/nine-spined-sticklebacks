
setwd("/path/to/dexseq-pca")

library(DEXSeq)
inDir=setwd("/path/to/dexseq/all")
countFiles=list.files(pattern = ".txt$",full.names = T)#inDir, input the count files
basename(countFiles)
flattenedFile=list.files(pattern = "gff$",full.names = T)#inDir,
#sampletable
samples=read.csv("samples.csv", stringsAsFactors = F)

dxd=DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = samples,
  design = ~sample+exon+ecotype:exon+tissue:exon,
  flattenedfile=flattenedFile
)
colData(dxd)
sampleAnnotation(dxd)

dxd_est = estimateSizeFactors( dxd )
dxd_dis = estimateDispersions( dxd_est )#norm

df=data.frame(dxd_dis@assays@data@listData$mu)
colnames(df)[1:26]=samples$samplename

norm_df=df[,1:26]
norm_df_t=t(norm_df)
data1 <- norm_df_t[ , which(apply(norm_df_t, 2, var) != 0)]#去掉方差=0的列

#pca
pca=prcomp(data1,center = T)
pca_df=data.frame(pca$x,ecotype=samples$ecotype,tissue=samples$tissue)
summ=summary(pca)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

library(ggplot2)
plot=ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue))+ geom_point(size=4,alpha = 0.6)+labs(x=xlab,y=ylab)+
  scale_color_manual(values = c("light blue","orange"))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))+
  theme_classic()+theme(text = element_text(size = 15))

ggsave("./dex_pca.tiff",width = 6,height = 4,dpi=300)

tiff("dex_pca_luo.tiff", width = 6, height = 4, units = "in", res = 300)
print(plot)
dev.off()

