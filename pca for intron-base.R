
setwd("/path/to/AS/leafcutter-pca")
all=read.table("./pca/test_all_perind.counts.gz",header = T,stringsAsFactors = F)
row.names(all)=all$chrom
all$chrom=NULL

all_ratio=sapply(X = all, FUN = function(v) {
  sapply(X = v,
         FUN = function(w) eval(parse(text=w)))
}
)
row.names(all_ratio)=rownames(all)
all_df=data.frame(all_ratio)

samples=read.csv("samples_leafcutter_dex.csv", stringsAsFactors = F)

all_df_nona=na.omit(all_df)
#order
order(colnames(all_df_nona))
all_df_nona=all_df_nona[,order(colnames(all_df_nona))]


library(pcaMethods)
pca=prcomp(all_df_nona,center = T)
pca_df=data.frame(pca$x,ecotype=samples$ecotype,tissue=samples$tissue)

library(ggplot2)
summ=summary(pca)
xlab=paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
ylab=paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")
plot=ggplot(pca_df,aes(x=PC1,y=PC2,color=ecotype,shape=tissue))+ geom_point(size=4,alpha = 0.6)+ labs(x=xlab,y=ylab)+
  scale_color_manual(values = c("light blue","orange"))+
  guides(shape=guide_legend(override.aes = list(shape=c(1:2))))+
  theme_classic()+theme(text = element_text(size = 15))

#ggsave("./all/all_intron_pca.tiff",width = 6,height = 4,dpi=300)

tiff("all_intron_pca.tiff", width = 6, height = 4, units = "in", res = 300)
print(plot)
dev.off()
