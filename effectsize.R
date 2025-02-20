
setwd("/path/to/effectsize")

library(reshape2)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(VennDiagram)
#library(effectsize)
library(rstatix)
library(car)
library(ggbiplot)
library(pheatmap)
library(dplyr)
library(tidyr)
library(patchwork)

library(ggplot2)  
library(reshape2)
library(ggpubr)


##1 DEG:repeat ANOVA for each DEgene in brain or liver####
group=read.table("samples_location.csv",header = T,sep = ",",row.names = 1)
group$sample=rownames(group)
group_brain=group[which(group$tissue=="brain"),]
group_brain$tissue=NULL
group_liver=group[which(group$tissue=="liver"),]
group_liver$tissue=NULL
str(group_brain)
as.factor(group_brain$ecotype)
as.factor(group_brain$location)

#1.1 brain DEG####
parallel_genes=read.table("brain_deg_counttable.csv",header = T,sep = ",",row.names = 1)
genelist_brain <- parallel_genes
genelist_brain <-t(genelist_brain)
#genelist_brain$sample=rownames(genelist_brain)
rownames(genelist_brain)==group_brain$sample
genelist_brain=merge(genelist_brain,group_brain,by.x=0,by.y="sample")
rownames(genelist_brain)=genelist_brain$Row.names
genelist_brain$Row.names=NULL


effectsizegenes <- data.frame(gene = colnames(genelist_brain)[1:(ncol(genelist_brain)-3)], 
                              eta_eco = rep(NA, ncol(genelist_brain)-3))

### calculate effect size of ecotype for each deg, put into table

for(i in 1:nrow(effectsizegenes)){
  effectsizegenes[i,"eta_eco"] <- eta_squared(lm(genelist_brain[,i] ~ ecotype+location+location*ecotype, genelist_brain))[1]  #[1,2]
  
} 

#1.1 liver DEG####

parallel_genes_liver=read.table("liver_deg_counttable.csv",header = T,sep = ",",row.names = 1)
genelist_liver <- parallel_genes_liver
genelist_liver <-t(genelist_liver)
#genelist_liver$sample=rownames(genelist_liver)
rownames(genelist_liver)==group_liver$sample
genelist_liver=merge(genelist_liver,group_liver,by.x=0,by.y="sample")
rownames(genelist_liver)=genelist_liver$Row.names
genelist_liver$Row.names=NULL


effectsizegenes_liver <- data.frame(gene = colnames(genelist_liver)[1:(ncol(genelist_liver)-3)], 
                              eta_eco = rep(NA, ncol(genelist_liver)-3))


### calculate effect size of ecotype for each deg, put into table

for(i in 1:nrow(effectsizegenes_liver)){
  effectsizegenes_liver[i,"eta_eco"] <- eta_squared(lm(genelist_liver[,i] ~ ecotype+location+location*ecotype, genelist_liver))[1]  #[1,2]
  } 



##2 DSgenes :repeat ANOVA for each DSgene in brain or liver####
setwd("/path/to/effectsize")
# genes here
DSgenes.dir <- "PSI"

# read in gene names and locations; these are only the CDSs
setwd(DSgenes.dir)
#2.1 brain DSGs####
DS_brain <- list.files(pattern = "brain_")
DSgene_brain <- c()
for (i in 1:length(DS_brain)) {
  temp_brain <- read.table(DS_brain[i],header = T,sep = ",",check.names = F)
  DSgene_brain=rbind(DSgene_brain,temp_brain)
}

DSgene_brain=separate(data=DSgene_brain, col="GeneID", into=c("gene","geneID"), sep = "-", remove = TRUE)
DSgene_brain$gene=NULL
DSgene_brain_gene=unique(DSgene_brain$geneID)

group_ds=read.table("samples.csv",header = T,sep = ",",row.names = 1)
group_ds_brain=group_ds[which(group_ds$tissue=="brain"),]
group_ds_liver=group_ds[which(group_ds$tissue=="liver"),]


brain_rlog=read.table("brainrlog_all.csv",header = T,sep = ",",check.names = F,row.names = 1)
DSgene_brain_count=brain_rlog[DSgene_brain_gene,]


genelist_brain_ds <- DSgene_brain_count
genelist_brain_ds<-na.omit(genelist_brain_ds)
genelist_brain_ds <-t(genelist_brain_ds)
#genelist_brain_ds$sample=rownames(genelist_brain_ds)
rownames(genelist_brain_ds)==rownames(group_ds_brain)
genelist_brain_ds=merge(genelist_brain_ds,group_ds_brain,by.x=0,by.y=0)
rownames(genelist_brain_ds)=genelist_brain_ds$Row.names
genelist_brain_ds$Row.names=NULL


effectsizegenes_DSbrain <- data.frame(gene = colnames(genelist_brain_ds)[1:(ncol(genelist_brain_ds)-3)], 
                              eta_eco = rep(NA, ncol(genelist_brain_ds)-3))


### calculate effect size of ecotype for each dsg, put into table

for(i in 1:nrow(effectsizegenes_DSbrain)){
  effectsizegenes_DSbrain[i,"eta_eco"] <- eta_squared(lm(genelist_brain_ds[,i] ~ ecotype+location+location*ecotype, genelist_brain_ds))[1]  #[1,2]
  } 


#2.2 liver DSGs####
DS_liver <- list.files(pattern = "liver_")
DSgene_liver <- c()
for (i in 1:length(DS_liver)) {
  temp_liver <- read.table(DS_liver[i],header = T,sep = ",",check.names = F)
  DSgene_liver=rbind(DSgene_liver,temp_liver)
}

DSgene_liver=separate(data=DSgene_liver, col="GeneID", into=c("gene","geneID"), sep = "-", remove = TRUE)
DSgene_liver$gene=NULL

DSgene_liver_gene=unique(DSgene_liver$geneID)

group_ds_liver=group_ds[which(group_ds$tissue=="liver"),]

liver_rlog=read.table("liverrlog_all.csv",header = T,sep = ",",check.names = F,row.names = 1)
DSgene_liver_count=liver_rlog[DSgene_liver_gene,]


genelist_liver_ds <- DSgene_liver_count
genelist_liver_ds<-na.omit(genelist_liver_ds)
genelist_liver_ds <-t(genelist_liver_ds)
#genelist_liver_ds$sample=rownames(genelist_liver_ds)
rownames(genelist_liver_ds)==rownames(group_ds_liver)
genelist_liver_ds=merge(genelist_liver_ds,group_ds_liver,by.x=0,by.y=0)
rownames(genelist_liver_ds)=genelist_liver_ds$Row.names
genelist_liver_ds$Row.names=NULL


effectsizegenes_DSliver <- data.frame(gene = colnames(genelist_liver_ds)[1:(ncol(genelist_liver_ds)-3)], 
                                      eta_eco = rep(NA, ncol(genelist_liver_ds)-3))



### calculate effect size of ecotype for each dsg, put into table

for(i in 1:nrow(effectsizegenes_DSliver)){
  effectsizegenes_DSliver[i,"eta_eco"] <- eta_squared(lm(genelist_liver_ds[,i] ~ ecotype+location+location*ecotype, genelist_liver_ds))[1]  #[1,2]
  } 




#3 DMC ####
# genes here
setwd("/path/to/effectsize")
DMC.dir <- "DMC"
# read in gene names and locations; these are only the CDSs
setwd(DMC.dir)

#3.1 brain DMC####
DMC_brain=read.table("eta.dmc.brain.csv",header = T,sep = ",",check.names = F,row.names = 1)


effectsizegenes_DMC_brain <- data.frame(gene = colnames(DMC_brain)[1:(ncol(DMC_brain)-3)], 
                                        eta_eco = rep(NA, ncol(DMC_brain)-3))


### calculate effect size of ecotype for each dmc, put into table

for(i in 1:nrow(effectsizegenes_DMC_brain)){
  effectsizegenes_DMC_brain[i,"eta_eco"] <- eta_squared(lm(DMC_brain[,i] ~ ecotype+location+location*ecotype, DMC_brain))[1]  #[1,2]
 } 




#3.2 liver DMC####
DMC_liver=read.table("eta.dmc.liver.csv",header = T,sep = ",",check.names = F,row.names = 1)


effectsizegenes_DMC_liver <- data.frame(gene = colnames(DMC_liver)[1:(ncol(DMC_liver)-3)], 
                                        eta_eco = rep(NA, ncol(DMC_liver)-3))


### calculate effect size of ecotype for each dmc, put into table

for(i in 1:nrow(effectsizegenes_DMC_liver)){
  effectsizegenes_DMC_liver[i,"eta_eco"] <- eta_squared(lm(DMC_liver[,i] ~ ecotype+location+location*ecotype, DMC_liver))[1]  #[1,2]
  } 




#4 DEmiRNA ####
setwd("/path/to/effectsize")
DEmiR.dir <- "miRNA"

setwd(DEmiR.dir)

#4.1 brain DEmiR####
DEmiR_brain_rlog=read.table("brain_mirna_rlog.csv",header = T,sep = ",",check.names = F,row.names = 1)
  
  DEmiR_brain_rlog <-t(DEmiR_brain_rlog)
  #genelist_brain_ds$sample=rownames(genelist_brain_ds)
  rownames(DEmiR_brain_rlog)==rownames(group_ds_brain)
  DEmiR_brain=merge(DEmiR_brain_rlog,group_ds_brain,by.x=0,by.y=0)
  rownames(DEmiR_brain)=DEmiR_brain$Row.names
  DEmiR_brain$Row.names=NULL
  
  effectsizegenes_DEmiR_brain <- data.frame(gene = colnames(DEmiR_brain)[1:(ncol(DEmiR_brain)-3)], 
                                            eta_eco = rep(NA, ncol(DEmiR_brain)-3))
  

  
  ### calculate effect size of ecotype for each DEmiR, put into table
  
  for(i in 1:nrow(effectsizegenes_DEmiR_brain)){
    effectsizegenes_DEmiR_brain[i,"eta_eco"] <- eta_squared(lm(DEmiR_brain[,i] ~ ecotype+location+location*ecotype, DEmiR_brain))[1]  #[1,2]
    } 
  


#4.1 liver DEmiR####
DEmiR_liver_rlog=read.table("liver_mirna_rlog.csv",header = T,sep = ",",check.names = F,row.names = 1)

DEmiR_liver_rlog <-t(DEmiR_liver_rlog)
#genelist_liver_ds$sample=rownames(genelist_liver_ds)
rownames(DEmiR_liver_rlog)==rownames(group_ds_liver)
DEmiR_liver=merge(DEmiR_liver_rlog,group_ds_liver,by.x=0,by.y=0)
rownames(DEmiR_liver)=DEmiR_liver$Row.names
DEmiR_liver$Row.names=NULL

effectsizegenes_DEmiR_liver <- data.frame(gene = colnames(DEmiR_liver)[1:(ncol(DEmiR_liver)-3)], 
                                          eta_eco = rep(NA, ncol(DEmiR_liver)-3))


### calculate effect size of ecotype for each DEmiR, put into table

for(i in 1:nrow(effectsizegenes_DEmiR_liver)){
  effectsizegenes_DEmiR_liver[i,"eta_eco"] <- eta_squared(lm(DEmiR_liver[,i] ~ ecotype+location+location*ecotype, DEmiR_liver))[1]  #[1,2]
 } 


# 5 final results£ºecotype effect size for DEG/DSG/DEmiR/DMC in brain and liver ####

#5.1 brain ####
violindata_brain=list("DEGs"=effectsizegenes$eta_eco,
                "DSGs"=effectsizegenes_DSbrain$eta_eco,
                "DEmiRs"=effectsizegenes_DEmiR_brain$eta_eco,
                "DMCs"=effectsizegenes_DMC_brain$eta_eco)

summary(effectsizegenes$eta_eco) 
summary(effectsizegenes_DSbrain$eta_eco)
summary(effectsizegenes_DEmiR_brain$eta_eco)
summary(effectsizegenes_DMC_brain$eta_eco)


violindata_brain=do.call(cbind, lapply(lapply(violindata_brain, unlist), `length<-`, max(lengths(violindata_brain))))

library(ggplot2)  
library(reshape2)
library(ggpubr)#

data_melt_brain <- melt(violindata_brain)

my_comparison=list(c("DEGs","DSGs"),
                   c("DEmiRs","DEGs"),
                   c("DSGs","DEmiRs"),
                   c("DSGs","DMCs"),
                   c("DEGs","DMCs"),
                   c("DEmiRs","DMCs"))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "NS"))


p_brain=
  ggplot(data_melt_brain, aes(x = Var2 , y = value )) +#, fill = Var2
  geom_violin()+
  geom_jitter(shape=21,color="grey",aes(fill=Var2),position = position_jitter(width = 0.1),alpha = 0.2)+
  geom_boxplot(outlier.shape=NA,width=0.4,position=position_dodge(0.2),alpha = 0.1,color="grey30")+
  #stat_summary(fun.y = mean,geom="point",shape=21,color="grey",fill="#8B040A",alpha = 1,size=5)+#
  labs(x = "", y= bquote("ecotype "~eta^2))+
  scale_fill_manual(values=c("#A8C8BD","#FFC28F", "#AAA9CB","#F2AED4"))+#"#C5F6FA"  "#66c2a5" "#ca0020  #'brown2','cornflowerblue', "#b06cc5","grey"
  theme_bw()+theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  ggtitle("brain")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons=my_comparison,
                     #label = "p.signif",#p.format #p.signif
                     tip.length = 0,
                     method = "wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                         symbols = c("****", "***", "**", "*", "NS")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=14,family = "sans"),#family = "serif"
        axis.text.y=element_text(vjust=1,size=14,family = "sans"),
        axis.title.y=element_text(vjust=1,size=14,family = "sans"),
        axis.title.x=element_text(vjust=1,size=14,family = "sans"))+ theme(text=element_text("sans"))

ggsave("DEG_DSG_DEmiR_DMC_brain.png", device="png",
       #path = "xxx",
       width = 6, height = 6, units="in",
       dpi=660)

data_mult_brain=data_melt_brain[,c("value","Var2")]
data_mult_brain$Var2 <- ordered(data_mult_brain$Var2,
                          levels = c("DEGs", "DSGs", "DEmiRs","DMCs"))
colnames(data_mult_brain)=c("value","group")


pairwise.wilcox.test(data_mult_brain$value, data_mult_brain$group,
                     p.adjust.method = "BH")


#5.2 liver####
violindata_liver=list("DEGs"=effectsizegenes_liver$eta_eco,
                      "DSGs"=effectsizegenes_DSliver$eta_eco,
                      "DEmiRs"=effectsizegenes_DEmiR_liver$eta_eco,
                      "DMCs"=effectsizegenes_DMC_liver$eta_eco)

summary(effectsizegenes_liver$eta_eco) 
summary(effectsizegenes_DSliver$eta_eco)
summary(effectsizegenes_DEmiR_liver$eta_eco)
summary(effectsizegenes_DMC_liver$eta_eco)#


violindata_liver=do.call(cbind, lapply(lapply(violindata_liver, unlist), `length<-`, max(lengths(violindata_liver))))

data_melt_liver <- melt(violindata_liver)

p_liver=
  ggplot(data_melt_liver, aes(x = Var2 , y = value )) +#, fill = Var2
  geom_violin()+
  geom_jitter(shape=21,color="grey",aes(fill=Var2),position = position_jitter(width = 0.1),alpha = 0.2)+
  geom_boxplot(outlier.shape=NA,width=0.4,position=position_dodge(0.2),alpha = 0.1,color="grey30")+
  #stat_summary(fun.y = mean,geom="point",shape=21,color="grey",fill="#8B040A",alpha = 1,size=5)+#
  labs(x = "", y= bquote("ecotype "~eta^2))+
  scale_fill_manual(values=c("#A8C8BD","#FFC28F", "#AAA9CB","#F2AED4"))+#"#C5F6FA"  "#66c2a5" "#ca0020  #'brown2','cornflowerblue', "#b06cc5","grey"
  theme_bw()+theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(comparisons=my_comparison,
                     label = "p.signif",#p.format #p.signif
                     tip.length = 0,
                     method = "wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                      symbols = c("****", "***", "**", "*", "NS")))+
  ggtitle("liver")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=14,family = "sans"),#family = "serif"
        axis.text.y=element_text(vjust=1,size=14,family = "sans"),
        axis.title.y=element_text(vjust=1,size=14,family = "sans"),
        axis.title.x=element_text(vjust=1,size=14,family = "sans"))+ theme(text=element_text("sans"))

ggsave("DEG_DSG_DEmiR_DMC_liver.png", device="png",
       #path = "xxx",
       width = 6, height = 6, units="in",
       dpi=660)


patch_eta_ecptype <- (p_brain|p_liver)
#patch_eta_ecptype + plot_annotation(tag_prefix='(',tag_levels = 'a',tag_suffix=')')

ggsave("DEG_DSG_DEmiR_DMC_brain_liver.png", device="png",
       #path = "xxx",
       width = 12, height = 6, units="in",
       dpi=880)


data_mult=data_melt_liver[,c("value","Var2")]
data_mult$Var2 <- ordered(data_mult$Var2,
                         levels = c("DEGs", "DSGs", "DEmiRs","DMCs"))
colnames(data_mult)=c("value","group")

pairwise.wilcox.test(data_mult$value, data_mult$group,
                     p.adjust.method = "BH")





