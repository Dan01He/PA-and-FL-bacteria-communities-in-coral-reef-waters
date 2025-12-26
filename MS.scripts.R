# Copyright (C) Dan He
#
# This file is part of the code accompanying:
# "Regional environmental heterogeneity under contrasting anthropogenic pressures
#  has differential effects on particle-attached than free-living bacteria
#  communities in coral reef waters."
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# For commercial or non-GPL licensing, please contact:
#   Dan He <hugh_20@163.com>


library(vegan)
library(ggplot2)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(picante)
library(ape)
library(parallel)
library(dplyr)
library(ggpubr)
library(reshape2)

###### data preparation
otu_table <- read.csv(file='otu_tb.csv',row.names=1,head=T)
metadata <- read.csv(file='metadata.csv',row.names=1,head=T)
env_data <- read.csv(file='env.csv',row.names=1,head=T)
otu_tax <- read.csv(file='otu_tax.csv',row.names=1,head=T)
tree <- read.tree('otu.tree')
otu_tb=sqrt(otu_table)
lonlat=read.csv(file='lonlat.csv',row.names=1,head=T)



######################################## Fig. 1A --sampling map
library(ggmap)
library(ggspatial)
ggmap::register_stadiamaps('86a2cc1b-f9f0-43f7-9fa5-1d8dce4189ac',write=T)

lonlat2=cbind(lonlat, metadata[rownames(lonlat), ])
colnames(lonlat2)[1:2]=c('lat','lon')
sbbox <- make_bbox(lon = lonlat2$lon, lat = lonlat2$lat, f = 0.1)
sites.data = data.frame(lonlat2)

bigbox=c(105,5,115,24);names(bigbox)=c('left','bottom','right','top') ## maychange
maptype = "stamen_terrain" 
####maptype = "stamen_toner"
####maptype = "stamen_watercolor"

map.base<- get_map(location = sbbox, maptype = maptype, source='stadia')
map.bigbase<- get_stadiamap(bbox = bigbox,zoom=8)

map.scale2 <- ggmap(map.base,extent = "panel",maprange = FALSE)  +
  geom_point(data = subset(lonlat2),size = 4, shape=2, mapping = aes(x = lon, y = lat,color=Region2)) +
  geom_text(data = subset(lonlat2), aes(x = lon+0.15, y = lat-0.15, label = Site,color=Region2), size = 3.5,alpha=70)

map.scale2+annotation_north_arrow(location = "tr", which_north = "true",  
pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"), style = north_arrow_fancy_orienteering)+
annotation_scale(location = "bl", width_hint = 0.5) + 
coord_sf(crs = 4326)
ggsave('Fig.1A.sitemap.pdf', height= 9, width=7)

######################################## Fig. 1B --richness and evenness
com=t(otu_tb)

alphaf=data.frame(Richness=specnumber(com),Evenness=diversity(com)/log(specnumber(com)),metadata[3:6])
alphaf.m=melt(alphaf, id=c(3:6))
alphaf.m$Region2=factor(alphaf.m$Region2, levels=c('Hainan','ZhongXisha','Nansha'))

p1=ggplot(subset(alphaf.m, variable=='Richness'), aes(x=Region2, y=value)) +geom_boxplot()+theme_bw()+facet_wrap(~Niche)
p2=ggplot(subset(alphaf.m, variable=='Evenness'), aes(x=Region2, y=value)) +geom_boxplot()+theme_bw()+facet_wrap(~Niche)

plt <- ggarrange(p1, p2, nrow = 2, ncol = 1)
ggsave(plt,file='Fig. 1B.pdf', width=8, height=6)  # basic plot

ggplot(subset(alphaf.m), aes(x=Region2, y=value)) +geom_boxplot(aes(fill=Niche))+theme_bw()+facet_wrap(~variable, scales='free')



######################################## Fig.1C NMDS plot 
com=t(otu_tb)
bray <- vegdist(com, method = "bray")

mds.bray=metaMDS(bray)
scr.bMDS=as.data.frame(cbind(scores(mds.bray),metadata))
ggplot(data = scr.bMDS, mapping = aes(x = NMDS1, y = NMDS2))+geom_point(aes(colour = Region2, shape= Niche),size=3)+theme_bw()+
annotate("text",x=max(scr.bMDS$NMDS1)*0.88, y=max(scr.bMDS$NMDS2)*0.99, size=2.5,alpha=0.9, label=paste("Stress:",round(mds.bray$stress,2)),fontface='italic')+
annotate("text",x=0.5*(max(scr.bMDS$NMDS1)+min(scr.bMDS$NMDS1)),y=1.05*max(scr.bMDS$NMDS2), label="Eukaryotes Bray-Curtis",size=4,color="darkblue")
ggsave("Fig.1C.pdf", width = 8, height = 6)  # generating the basic plot that may need some outlook modification 



######################################## Fig.1D phyla barplot 
library(doBy)
library(ggpubr)
com=t(otu_table)

# total community
taxa_table <- otu_tax[colnames(com),]
colnames(taxa_table)[7]="subphylum"
subphylum <- taxa_table[, "subphylum"]
phylum_abundance <- aggregate(. ~ subphylum, data = data.frame(subphylum, t(com)), sum)
phylum_abundance[,1]->rownames(phylum_abundance)
phylum_abundance[,1]=NULL
phylum_abundance=phylum_abundance[order(rowSums(phylum_abundance),decreasing=T),]
phylum_rla=t(phylum_abundance)
phylum_rla=vegan::decostand(phylum_rla,'tot',1)*100
phylum_rlaf=cbind(phylum_rla, metadata[rownames(phylum_rla),])
phylum_abundance_top15 <- rbind(phylum_abundance[1:15,],colSums(phylum_abundance[-(1:15),])) 
rownames(phylum_abundance_top15)[16]='Minor groups'

taxa_forpl=t(phylum_abundance_top15)
taxa_forpl=vegan::decostand(taxa_forpl,'tot',1)*100
taxa_forpl=cbind(taxa_forpl, metadata[rownames(taxa_forpl),])
taxa.mlt=reshape2::melt(taxa_forpl,id=c((17:(dim(taxa_forpl)[2]))))
taxa.mlt$Region2=factor(taxa.mlt$Region2, levels=c('Hainan','ZhongXisha','Nansha'))

# barplots
color_32=c('#C1E8E9','#1F91FF','#FF915F','#B29FFF','#FFCD00','#00E97E','#0DC7C8','#CC90CB','#739C2F','#EEDFA7','#5912F4','#C71585','#CCCCCC','#00293F','#E2F100','#C93800','#87CEFA','#0915EC','#A52A2A','#00FF00','#FFA500','#A020F0','#008B00','#FFC0CB','#8B5A00','#8B8B00','#8299EA','#C6FF3F','#FF3BA1','#E62617','#006A66','#2B2B2B')
ggplot(taxa.mlt, aes(x=Region2, y=value)) +geom_bar(stat="summary",fun='mean',aes(fill=variable),position="stack")+scale_fill_manual(values=color_32[5:20])+theme_bw()+coord_flip()+facet_wrap(~Niche)+theme(legend.text = element_text(size=8),legend.key.size = unit(0.3,'cm'))
ggsave(file="Fig. 1D.pdf", width = 9, height = 3)  # generating the basic plot that may need some outlook modification 



########################################  Fig.2  env barplot
library(reshape2)
library(ggplot2)
Envs=cbind(env_data, metadata)

ENV=Envs[Envs$Niche=='Freeliving', c(1:8,12:14)]
ENV.m=melt(ENV,id=c(9:11))
ENV$Region2=factor(ENV$Region2, levels=c('Hainan','ZhongXisha','Nansha'))

envmc <- bear::summarySE(ENV.m, measurevar="value", groupvars=c("Region2",'Layer',"variable"),na.rm=T, conf.interval=.95)
envmc$Region2=factor(envmc$Region2, levels=c('Hainan','ZhongXisha','Nansha'))
dodge <- position_dodge(width=0.9)
envs.bar=ggplot(envmc, aes(x=Layer, y=value,fill=Region2)) +
geom_bar(stat="identity",position=dodge) +
geom_errorbar(width=.5, aes(ymin=value-se, ymax=value+se),position=dodge)+facet_wrap(~variable,scales="free")+theme_bw()+
theme(line = element_blank())

ggsave(envs.bar,file='Fig. 2.pdf', width=8, height=6)



########################################  Fig. 3  differential ASVs heatmap and FC bars

# 安装依赖包
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

com=t(otu_table)
com.sel=com[metadata$Niche=='PA_attach',]  # uncomment this line for PA communities
# com.sel=com[metadata$Niche=='Freeliving',]  # uncomment this line for FL communities
com.sel=com.sel[,colSums(com.sel)>0]
metadata_sel=metadata[rownames(com.sel),]

otu=t(com.sel); group=metadata_sel
group <- group %>% mutate(Anthr=ifelse( Region2=='Hainan', 'HighAnthr', 'LowAnthr') ) 
group$Anthr=factor(group$Anthr)

# DESeq
dds <- DESeqDataSetFromMatrix(countData = otu, colData = data.frame(group), design = ~ Anthr)
keep <- (rowSums(counts(dds) >= 10) >= 3) & (rowSums(otu) >= 0.00001 * sum(otu))
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, contrast=c("Anthr", "HighAnthr", "LowAnthr"))
res_sorted <- res[order(res$pvalue), ]

# limma
library(limma)
library(edgeR)
dge <- DGEList(counts=otu)

keep <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)
design <- model.matrix(~group$Anthr)
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, n=Inf, sort.by="P")

#deseq
res_df <- as.data.frame(res_sorted)
res_df$significant <- ifelse(res_df$padj<0.001 & abs(res_df$log2FoldChange)>=5, "Y", "N")
# res_df -> FL.deseq.res  # uncomment this line for FL communities
res_df -> PA.deseq.res  # uncomment this line for PA communities

#limma
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$adj.P.Val<0.001 & abs(res_df$logFC)>=5, "Y", "N")
# res_df -> FL.limma.res  # uncomment this line for FL communities
res_df -> PA.limma.res  # uncomment this line for PA communities

## pheatmap the differentiated ASVs by the two methods
FL.diffASVs=intersect(rownames(subset(FL.deseq.res, significant=='Y')), rownames(subset(FL.limma.res, significant=='Y'))) 
PA.diffASVs=intersect(rownames(subset(PA.deseq.res, significant=='Y')), rownames(subset(PA.limma.res, significant=='Y'))) 

FL.diffFC=data.frame(ASVs=FL.diffASVs, logFC=0.5*(FL.deseq.res[FL.diffASVs, 'log2FoldChange']- FL.limma.res[FL.diffASVs, 'logFC']), taxid=paste0(FL.diffASVs, ' ', otu_tax[FL.diffASVs,'Genus'],'(',otu_tax[FL.diffASVs,'taxon'],')'))
PA.diffFC=data.frame(ASVs=PA.diffASVs, logFC=0.5*(PA.deseq.res[PA.diffASVs, 'log2FoldChange']- PA.limma.res[PA.diffASVs, 'logFC']), taxid=paste0(PA.diffASVs, ' ',otu_tax[PA.diffASVs,'Genus'],'(',otu_tax[PA.diffASVs,'taxon'],')'))
rownames(FL.diffFC)=FL.diffFC$ASVs
rownames(PA.diffFC)=PA.diffFC$ASVs

FL.diffFC2 <- FL.diffFC %>%
  mutate(absolute_value = abs(FL.diffFC$logFC)) %>%  # 创建绝对值列
  group_by(sign = sign(FL.diffFC$logFC)) %>%         # 按正负号分组
  slice_max(absolute_value, n = 10) %>%                          # 取每组前10个
  ungroup() %>%                                   # 取消分组
  select(-absolute_value, -sign)     

PA.diffFC2 <- PA.diffFC %>%
  mutate(absolute_value = abs(PA.diffFC$logFC)) %>%  # 创建绝对值列
  group_by(sign = sign(PA.diffFC$logFC)) %>%         # 按正负号分组
  slice_max(absolute_value, n = 10) %>%                          # 取每组前10个
  ungroup() %>%                                   # 取消分组
  select(-absolute_value, -sign)     


com=t(otu_table)
com.sel=com[metadata$Niche=='PA_attach',]  # uncomment this line for PA communities
# com.sel=com[metadata$Niche=='Freeliving',]  # uncomment this line for FL communities

com.sel=com.sel[,colSums(com.sel)>0]
metadata_sel=metadata[rownames(com.sel),]
com.sel=decostand(com.sel, 'tot', MARGIN=1)*100

com.sel=com.sel[,PA.diffFC2$ASVs]  # uncomment this line for PA communities
# com.sel=com.sel[,FL.diffFC2$ASVs]   # uncomment this line for FL communities

mid=t(com.sel)

rownames(mid)=PA.diffFC2$taxid  # uncomment this line for PA communities
# rownames(mid)=FL.diffFC2$taxid  # uncomment this line for FL communities

annotation_col <- metadata_sel[c(6)]
p1=pheatmap::pheatmap(mid, show_rownames = T, show_colnames=F, scale="row",
 cluster_cols = F, cluster_rows= T,
annotation_col = annotation_col ) 
clustered_row_names <- rownames(mid)[p1$tree_row$order]
clustered_rowID=sapply(strsplit(clustered_row_names, " "), function(x) x[1])

PA.diffFC2_=PA.diffFC2[match(clustered_rowID, PA.diffFC2$ASVs),];PA.diffFC2_$taxid <- factor(PA.diffFC2_$taxid, levels = rev(unique(PA.diffFC2_$taxid)))  # uncomment this line for PA communities
# FL.diffFC2_=FL.diffFC2[match(clustered_rowID, FL.diffFC2$ASVs),];FL.diffFC2_$taxid <- factor(FL.diffFC2_$taxid, levels = rev(unique(FL.diffFC2_$taxid)))  # uncomment this line for FL communities

pdf(file='PA.2diff.heatmap.pdf', width=8, height=8)  # uncomment this line for PA communities 
# pdf(file='FL.2diff.heatmap.pdf', width=8, height=8)  # uncomment this line for FL communities
p1
dev.off()

ggplot(data=data.frame(PA.diffFC2_, sign=ifelse(PA.diffFC2_$logFC>0, 'pos', 'neg')),aes(x=abs(logFC), y=taxid))+geom_bar(width = 1, stat = "identity", position="stack",aes(fill=sign))+theme_bw()
ggsave(file='PA.2diff.bar.pdf', width=8, height=8) 

ggplot(data=data.frame(FL.diffFC2_, sign=ifelse(FL.diffFC2_$logFC>0, 'pos', 'neg')),aes(x=abs(logFC), y=taxid))+geom_bar(width = 1, stat = "identity", position="stack",aes(fill=sign))+theme_bw()
ggsave(file='FL.2diff.bar.pdf', width=8, height=8)



########################################  Fig. 4  randomForest importance plots and VPA plots

######################### Fig. 4 (A,C,E,G) randomForest importance plot
library(vegan)
library(randomForest)
library(ggplot2)
library(dplyr)

mat2cols = function(m1,val1="value"){ # transfer distance matrix to columns
  m1=as.matrix(m1)
  lt = lower.tri(m1)  
  res = data.frame(row = row(m1,as.factor = T)[lt],  
                   col = col(m1,as.factor = T)[lt],  
                   val1 = m1[lt])
  names(res)[3] = c(val1) 
  return(res)
}

hageodist=function(L1, phi1, L2, phi2){ # calculating geographical distance based on the longitude and latitude
 a = 6378.14
 f = 1/298.257
 F = (phi1+phi2)/2
 G = (phi1 - phi2)/2
 ramda <- (L1 - L2)/2
 
 S = (sin(G*pi/180)^2)*(cos(ramda*pi/180)^2) + (cos(F*pi/180)^2)*(sin(ramda*pi/180)^2)
 C= (cos(G*pi/180)^2)*(cos(ramda*pi/180)^2) + (sin(F*pi/180)^2)*(sin(ramda*pi/180)^2)
 
 omega = atan(sqrt(S/C))
 R = sqrt(S*C)/omega
 D = 2*omega*a
 
 H1 = (3*R-1)/(2*C)
 H2 = (3*R+1)/(2*S)
 
 res = D*(1 + f*H1*(sin(F*pi/180)^2)*(cos(G*pi/180)^2) - f*H2*(cos(F*pi/180)^2)*(sin(G*pi/180)^2))
 return(round(res,3))
}


com=t(otu_tb)
com.sel=com[metadata$Niche=='Freeliving',]
# com.sel=com[metadata$Niche=='PA_attach',] # uncomment this line for PA communities
com.sel=com.sel[,colSums(com.sel)>0]

Envs.sel=Envs[rownames(com.sel),]
Envs.sel[1:8] <- Envs.sel[1:8] %>% mutate(across(-pH,scale))

dis=matrix(nrow=length(Envs.sel$Lon),ncol=length(Envs.sel$Lat))
rownames(dis)<-rownames(Envs.sel)
colnames(dis)<-rownames(Envs.sel)
attach(Envs.sel)
for (i in 1:length(Lon)){
for (j in 1:length(Lat)){
if(i>=j)
{dis[i,j]=hageodist(Lon[i],Lat[i],Lon[j],Lat[j])}

}}
dis=as.dist(dis)
dis[dis=='NaN']=0.005 ## for surface and bottom water samples from the same site
geo_dis=mat2cols(dis); pcnm=pcnm(dis)
detach(Envs.sel)
Envs.sel=cbind(Envs.sel, pcnm$vectors[,1:2])
envs.sel <- Envs.sel %>% select(-c(9:13))

metadata.sel=metadata[rownames(com.sel), ]
metadata.sel <- apply(metadata.sel[c(2,4,6)], 2, as.factor)

alpha_diversity <- diversity(com.sel, index = "shannon")
beta_diversity <- vegdist(com.sel, method = "bray")

###### RF alpha
rf_alpha_data <- data.frame(
  alpha = alpha_diversity,
  envs.sel,
  metadata.sel
)

set.seed(123)  
rf_alpha_model <- randomForest(
  alpha ~ ., 
  data = rf_alpha_data,
  importance = TRUE,  
  ntree = 1000,       
  mtry = sqrt(ncol(rf_alpha_data) - 1)  
)
print(rf_alpha_model)
importance_alpha <- as.data.frame(importance(rf_alpha_model)) 
rownames(importance_alpha) -> importance_alpha$variable; colnames(importance_alpha)[1]='importance'


ggplot(subset(importance_alpha,importance>0), aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  # labs(title = "FL alpha", x = "", y = "%IncMSE") + # uncomment this line for FL communities
  labs(title = "PA alpha", x = "", y = "%IncMSE") + # uncomment this line for PA communities
  theme_minimal()
# ggsave(file='FL alpha rft.pdf',width=4,height=6) # uncomment this line for FL communities
ggsave(file='PA alpha rf.pdf',width=4,height=6)  # uncomment this line for PA communities


###### RF beta
pcoa <- cmdscale(beta_diversity, k = 5, eig = TRUE)
beta_scores <- as.data.frame(pcoa$points)
colnames(beta_scores) <- paste0("PCo", 1:ncol(beta_scores))
rf_beta_data <- data.frame(
  beta_scores,
  envs.sel,
  metadata.sel
)

# randomForest analysis for each main coordinate 
beta_models <- list()
beta_importance <- list()
for (i in 1:ncol(beta_scores)) {
  formula <- as.formula(paste0("PCo", i, " ~ ."))
  beta_models[[i]] <- randomForest(
    formula, 
    data = rf_beta_data[,c(i,6:(dim(rf_beta_data)[2]))],  
    importance = TRUE,
    ntree = 1000
  )
  beta_importance[[i]] <- importance(beta_models[[i]])
}

combined_importance <- do.call(cbind, lapply(beta_importance, function(x) x[, 1]))
rownames(combined_importance) <- rownames(beta_importance[[1]])
colnames(combined_importance) <- paste0("PCo", 1:ncol(beta_scores))

# averaging the importances for different coordinates
avg_importance <- rowMeans(combined_importance)
avg_importance_df <- data.frame(
  variable = names(avg_importance),
  importance = avg_importance
) %>% arrange(desc(importance))

ggplot(subset(avg_importance_df, importance>0), aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "FL beta", x = "", y = "%IncMSE") +  # uncomment this line for FL communities
  # labs(title = "PA beta", x = "", y = "%IncMSE") +  # uncomment this line for PA communities
  theme_minimal()
 ggsave(file='FL beta rf.pdf',width=4,height=6)
 # ggsave(file='PA beta rf.pdf',width=4,height=6)  # uncomment this line for PA communities



#########################  Fig.4 (B,D,F,H)  VPA plot
Envs=cbind(env_data, lonlat, metadata[c(3,4,6)])


for (ni in c('Freeliving','PA_attach')) {
com=t(otu_tb)
com.sel=com[metadata$Niche==ni,]
com.sel=com.sel[,colSums(com.sel)>0]

Envs.sel=Envs[rownames(com.sel),]
Envs.sel[1:8] <- Envs.sel[1:8] %>% mutate(across(-pH,scale))

attach(Envs.sel)
#diagol format
dis=matrix(nrow=length(Envs.sel$Lon),ncol=length(Envs.sel$Lat))
rownames(dis)<-rownames(Envs.sel)
colnames(dis)<-rownames(Envs.sel)
for (i in 1:length(Lon)){
for (j in 1:length(Lat)){
if(i>=j)
{dis[i,j]=hageodist(Lon[i],Lat[i],Lon[j],Lat[j])}

}}
dis=as.dist(dis)
dis[dis=='NaN']=0.005 ## for surface and bottom water samples from the same site

geo_dis=mat2cols(dis); pcnm=pcnm(dis)
assign(paste0(ni,'_geodis'), geo_dis); assign(paste0(ni,'_pcnm'), pcnm)
detach(Envs.sel)

ENVs.sel=cbind(Envs.sel, pcnm$vectors[,1:3])
VP1=varpart(dist(diversity(com.sel)), ~Region2, ENVs.sel[c(1:8)], ~PCNM1+PCNM2+PCNM3, ~ Layer, data=ENVs.sel)
VP2=varpart(vegdist(com.sel), ~Region2, ENVs.sel[c(1:8)], ~PCNM1+PCNM2+PCNM3, ~ Layer, data=ENVs.sel); 

assign( paste0(ni,'_VP1'), VP1)
assign( paste0(ni,'_VP2'), VP2)
} 


pdf(file="Fig. VPA alpha.pdf", width = 8, height = 5)  # generating the basic plot that may need some outlook modification  
par(mfrow = c(1, 2))
plot( Xnames=c("Region","Envrion.","PCNM", "Layer"),Freeliving_VP1 )
plot( Xnames=c("Region","Envrion.","PCNM", "Layer"),PA_attach_VP1 )
dev.off()

pdf(file="Fig. VPA beta.pdf", width = 8, height = 5)  # generating the basic plot that may need some outlook modification  
par(mfrow = c(1, 2))
plot( Xnames=c("Region","Envrion.","PCNM", "Layer"),Freeliving_VP2 )
plot( Xnames=c("Region","Envrion.","PCNM", "Layer"),PA_attach_VP2 )
dev.off()





###################################  Fig.4 -- network plots  and  zipi plots

## 1, generate niche(lifestyle) specific otu tables for network analysis
wd=getwd()
niche=c('FL','PA')
for (ni in niche) 
{
wdni=paste0(wd,'/',ni)
dir.create(wdni)

com=t(otu_table)
if (ni== 'FL') { nicheCOM=com[s.factor$Niche=='Freeliving',]; s.factor_ni=s.factor[s.factor$Niche=='Freeliving', ] }
else { nicheCOM=com[s.factor$Niche=='PA_attach',]; s.factor_ni=s.factor[s.factor$Niche=='PA_attach', ] }

nicheCOM=nicheCOM[, colSums(nicheCOM) > 0]
nicheCOM -> netcom
netcom=netcom[,(order(colSums(netcom), decreasing=T))[1:2000]] 
print(paste('dim for', ni))
print(dim(netcom))
nettax=otu_tax[match(colnames(netcom),rownames(otu_tax)),]
df=as.data.frame(t(netcom))
# copy rownames to the first column starting with the "#OTU ID" 
df_with_otu <- cbind('#OTU ID' = rownames(df), df)
rownames(df_with_otu) <- NULL

# write table, starting with the "#OTU ID" 
write.table(
  df_with_otu,
  file = paste0(wdni, '/', ni, '_otbtp2k.tsv'),
  sep = "\t",
  row.names = FALSE,  
  col.names = TRUE,   
  quote = FALSE       
)
}



## 2, in linux , run the "fastspar.ni.sh", then copy the files and folders into the R work diretory


## 3, in R, load the correlation and pvalues tsv files generated by the fastspar.sh in step 2
library(pacman)
p_load(igraph, ggplot2,reshape2)
p_load(RMThreshold)
p_load(WGCNA)
p_load(multtest)
p_load(vegan)
p_load(info.centrality)
allowWGCNAThreads(nThreads =3 )

rm(g)

connectivity_entropy <- function(graph) {  
  components_info <- igraph::components(graph)  
  nodes_per_component <- components_info$csize  
  total_nodes <- vcount(graph)  
  proportions <- nodes_per_component / total_nodes  
  entropy <- -sum(proportions * log(proportions))  
  
  return(entropy)  
} 
 
connectance <- function(g) {  
  components_info <- igraph::components(g)  
  reachable_pairs_count <- 0  
  total_pairs_count <- 0  
  component_sizes <- components_info$csize  
  for (size in component_sizes) {  
    if (size > 1) {  
      total_pairs_count <- total_pairs_count + choose(size, 2)  
    }  
  }  
  total_node_pairs <- choose(vcount(g), 2)  
  W <- total_node_pairs - total_pairs_count  
  v=vcount(g)
  C=1-(2*W/(v*(v-1)))
  return(C)  
}  


regionsub_sparcc.node_traits=NULL
regionsub_sparcc.basic_traits=NULL
taxon16=factor(rownames(phylum_abundance_top15),
 levels=c('Gammaproteobacteria',
 'Alphaproteobacteria',
 'Bacteroidota',
 'Cyanobacteria',
 'Actinobacteriota',
 'Planctomycetota',
 'Verrucomicrobiota',
 'Bdellovibrionota',
 'Desulfobacterota',
 'Marinimicrobia_(SAR406_clade)',
 'Proteobacteria_unclassified',
 'Myxococcota',
 'Bacteria_unclassified',
 'Firmicutes',
 'SAR324_clade(Marine_group_B)',
 'Minor groups'
))

region=c('Hainan','ZhongXisha','Nansha')
niche=c('FL','PA')
wd=getwd()

for (ni in niche) 
{
for (re in region)
{
wdi=paste0(wd, '/', ni)
setwd(wdi)
nm=paste0(ni,'_Sparcc.net')

corr=read.delim(file='median_correlation.tsv', row.names=1 )
pval=read.delim(file='pvalues.tsv', row.names=1)

col.rawp=mat2cols(pval)
col.corr=mat2cols(corr)

cols.corwithp=data.frame(obj1=col.rawp$row, obj2=col.rawp$col, corr=col.corr$value, rawp=col.rawp$value)
thre=0.7  # may change
Pval=0.01 # may change
subset(cols.corwithp, (abs(corr)>=thre & rawp<=Pval), obj1:rawp) -> cols.corwithp.sel #0.7 choose from SparCC
#save(cols.corwithp.sel, file='corwithp.RData')

cols.corwithp.sel[c(1,2,3)]->col.SparCCcorr
igraph::graph_from_data_frame(col.SparCCcorr,directed=F)-> SparCC.net
cleaned.matrix <- as.matrix(igraph::as_adjacency_matrix(SparCC.net,attr='corr'))

diag_forND=vegan::decostand(abs(cleaned.matrix),'pa')
 diag(diag_forND)=0
 isSymmetric.matrix(diag_forND)
 write.csv(diag_forND, file="diag_forND.csv")

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\86158\\AppData\\Local\\Programs\\Python\\Python312")  # may change path
file.copy(paste0(wd,'/ND_a.py'), 'ND_a.py')
source_python('ND_a.py') ###  network deconvolution
file.remove('ND_a.py')


ND_out=read.csv('ND_out.csv',head=T,row.names=1)
 diag(ND_out)=1
cleaned.matrix2 <- rm.denoise.mat(as.matrix(ND_out), threshold = thre) 
cleaned.matrix2 <- rm.discard.zeros(cleaned.matrix2) 
mat2cols(cleaned.matrix2)->col.ND_out; subset(col.ND_out,value!=0)->col.ND_out

cols.corwithp[paste(cols.corwithp[,1],cols.corwithp[,2])%in%paste(col.ND_out[,1],col.ND_out[,2])|
paste(cols.corwithp[,2],cols.corwithp[,1])%in%paste(col.ND_out[,1],col.ND_out[,2]),]-> ND.edges

graph_from_data_frame(ND.edges[c(1,2,3)],directed=F)-> SparCC_ND.net
SparCC_ND.nettax=otu_tax[match(as_ids(V(SparCC_ND.net)),rownames(otu_tax)),]; SparCC_ND.nettax=na.omit(SparCC_ND.nettax)
SparCC_ND.nettax=data.frame(OTUid=rownames(SparCC_ND.nettax),SparCC_ND.nettax)
ND.edges$Relationship=NA
ND.edges$Relationship[ND.edges$corr>0]='Positive'
ND.edges$Relationship[ND.edges$corr<0]='Negative'
graph_from_data_frame(ND.edges,directed=F,vertices=SparCC_ND.nettax)-> SparCC_ND.net


rm(cols.corwithp, col.corr)


SparCC_ND.net=delete.vertices(SparCC_ND.net, which(degree(SparCC_ND.net)<1))
assign(nm, SparCC_ND.net)
com=t(otu_tb^2)
if (ni== 'FL') { netcom=com[s.factor$Niche=='Freeliving',]; netcom=netcom[,(order(colSums(netcom), decreasing=T))[1:2000]]; s.factor_ni=s.factor[s.factor$Niche=='Freeliving', ] }
else { netcom=com[s.factor$Niche=='PA_attach',]; netcom=netcom[,(order(colSums(netcom), decreasing=T))[1:2000]]; s.factor_ni=s.factor[s.factor$Niche=='PA_attach', ] }

region_netcom=netcom[s.factor_ni$Region2==re,]
region_netcom=region_netcom[, colSums(region_netcom) > 0]
region_nettax=otu_tax[match(colnames(region_netcom),rownames(otu_tax)),]

net=get(nm)
colnames(region_netcom)->vids
which(as_ids(V(net))%in%vids)->vidsindex
nm2=paste(ni,re,'net',sep='_')

g=induced_subgraph(net, vidsindex)
g=delete.vertices(g, which(degree(g)<1))


# plot region specific nets
linecol_pal <- c("tomato1","skyblue")
edge_cat=as.numeric(ifelse(edge.attributes(g)$Relationship=='Negative',1,2))

sqrt(colSums(region_netcom[,vertex.attributes(g)$name]))->vertex.attributes(g)$Size
vertex.attributes(g)$taxon=factor(taxa_table[vertex.attributes(g)$name,]$subphylum, levels=levels(taxon16))
 vertex.attributes(g)$taxon[is.na(vertex.attributes(g)$taxon)]='Minor groups'

node_cat=as.numeric(vertex.attributes(g)$taxon)
V(g)$color <-(color_32[5:20])[node_cat] #


V(g)$size=plotrix::rescale(c(vertex.attributes(g)$Size, 1, 550 ),c(2,30))[1:vcount(g)] # rescale size
V(g)$frame.color=NA; V(g)$label=NA
E(g)$color= c("tomato1","grey10")[edge_cat]; E(g)$width=0.002; E(g)$curved=T
E(g)$name=as_ids(E(g))
lo=layout_with_fr(g, niter=1000, start.temp=500)
# g <- set_g_attr(g, "layout", layout_nicely(g))
edge.attr=as.data.frame(edge.attributes(g))
edge.attr <- tidyr::separate(edge.attr, 'name', c("Source", "Target"),  "\\|", F, F) %>% dplyr::select(c(7,8),everything())  ## for Gephi

assign(nm2,g)
node_attr <- as.data.frame(vertex.attributes(g)) %>% dplyr::rename(Id=name, Ori_Size=Size)
write.csv(edge.attr, paste0(nm2,'.edges.csv'),row.names=F)  ## for gephi
write.table(node_attr , paste0(nm2,'.nodes.tsv'), row.names=FALSE, sep='\t') ## for gephi
write.graph(g, paste0(nm2,'.edgelist'),"edgelist")

##### 4 load the edges.csv and nodes.csv into gephi and generate the network plots that were combined to Fig. 4 (A-F)


########################  Fig. 4 (G-L) zi pi plots
write.csv(as_adj(g,sparse=F,type="both"),file="net.adja.csv")
p_load(reticulate)
Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\86158\\AppData\\Local\\Programs\\Python\\Python312")

file.copy(paste0(wd,'/zp_cal_a.py'), 'zp_cal_a.py')
source_python('zp_cal_a.py')
file.remove('zp_cal_a.py')

zp=read.csv("zipi.csv",head=F,stringsAsFactors=F)
colnames(zp)=c('name','Zi','Pi','module')

net.vertraits=as.data.frame(vertex.attributes(g))
zp=cbind(zp, net.vertraits[zp$name,])
zp[,1]=NULL
zp$taxpl=paste0(zp$Genus,' (', zp$taxon,')')
zp$taxpl2=paste0(zp$name,' (', zp$Genus,')')

assign(paste0(nm2,'.zp'),zp)
assign(paste0(nm2,'.keynodes.zp'),subset(zp,Pi>=0.82|Zi>=2.5))
write.csv(get(paste0(nm2,'.keynodes.zp')),file=paste0(nm2,'.keynodes.zp.csv'))

library(ggplot2)
zppl=ggplot(data=zp,aes(x=Pi,y=Zi))+geom_point(data=subset(zp,Pi>=0.82|Zi>=2.5),color=subset(zp,Pi>=0.82|Zi>=2.5)$color,size=1)+geom_point(data=subset(zp,Pi<0.82&Zi<2.5),color=subset(zp,Pi<0.82&Zi<2.5)$color, alpha=0.2,size=0.5)+scale_x_continuous(limits=c(0,0.8))+
scale_y_continuous(limits=c(-2,4.5))+geom_vline(aes(xintercept=0.82), linetype="dashed")+
geom_hline(aes(yintercept=2.5), linetype="dashed")+
geom_text(data=subset(zp,Pi>=0.82|Zi>=2.5),aes(label=taxpl2,x=Pi,y=Zi),color=subset(zp,Pi>=0.82|Zi>=2.5)$color, nudge_y=0.08,size=2,alpha=20)+
ggtitle(paste0(nm2,'.net',' ZiPi plot'))+theme_bw()+
theme(line=element_blank(),legend.position='none') 

# generate the zipi plots that were combined to Fig. 5 
ggsave(zppl, file=paste0(nm2,'.net.zp.pdf') ,dpi=600, width=4, height=4 )  
#zipi plots done

# calculate region-net  network traits
posneg2=function(net){
require(igraph)
Vids=as_ids(V(net))
n=length(Vids)
df=data.frame(as_edgelist(net),Relationship=edge.attributes(net)$Relationship) # check if it is relationship or interactionType
mat=data.frame(matrix(nrow=n,ncol=4))
for ( i in (1:n )) {
df2=df[grepl(Vids[i],df[,1])|grepl(Vids[i],df[,2]),]
pos=sum(grepl("Positive",df2[,3]))  # check if it is Positive/Negative or others 
neg=sum(grepl("Negative",df2[,3]))
mat[i,1]=Vids[i]; mat[i,2]=pos; mat[i,3]=neg; mat[i,4]=round(pos/(neg+pos)*100,2)
}
mat=as.data.frame(mat,stringsAsFactors=F)
colnames(mat)=c("OTUid","posdegree","negdegree","PPE")
return(mat)
}

# node level traits
betweenness(g, directed=F, weights=NA,normalized=T)->v1
closeness(g, mode="all",normalized=T)->v2
degree(g, mode="all",normalized=F)->v3
transitivity(g, type="local")->v4
knn(g)$knn->v5
posneg2(g)$PPE ->v6
info.centrality.vertex(g)->v7
eigen_centrality(g)$vector -> v8
rep(re, vcount(g)) -> v9   #### may change as the sample group info
rep(ni, vcount(g)) -> v10   #### may change as the sample group info
V(g)$taxon -> v11

# net level traits
mean_distance(g)->g1
mean(degree(g))->g2
centr_betw(g)$centralization->g3
centr_clo(g)$centralization->g4
edge_density(g)->g5
assortativity_degree(g)->g6
transitivity(g, type="global")->g7
max(membership(cluster_fast_greedy(g)))->g8
modularity(cluster_fast_greedy(g))->g9
vcount(g)->g10
ecount(g)->g11
diameter(g)->g12
g13<-round(sum(edge.attributes(g)$Relationship=='Positive')/
length(edge.attributes(g)$Relationship)*100,2)
g14<-info.centrality.network(g)
g15<-network.efficiency(g)
g16 <- connectance(g)
g17 <- connectivity_entropy(g)
g18 <- nm2

mid1=as.data.frame(matrix(ncol=18,nrow=1))# ## mid should be classifed before the loop
for (j in 1:18) mid1[,j]=get(paste0('g',j))

colnames(mid1)=c('mean_distance',
'mean_degree',
'centr_betw',
'centr_clo',
'edge_density',
'assortativity_degree',
'transitivity',
'cluster_number',
'modularity',
'vcount',
'ecount',
'diameter',
'PPE',
'infocentrality', 
'efficiency',
'connectance',
'connectivity_entropy',
'group')
rm(list=paste0('g',1:18))
assign(paste0(nm2,'.basic_traits'),mid1)

mid2=as.data.frame(matrix(ncol=11,nrow=vcount(g)))# mid should be classifed before the loop
for (k in 1:11) {
mid2[,k]=get(paste0('v',k))
}
colnames(mid2)=c('betweenness',
'closeness',
'degree',
'transitivity',
'knn',
'PPE',
'info.centrality',
'eigen_centrality',
'Region',
'Niche',
'subphylum')

rm(list=paste0('v',1:11))
assign(paste0(nm2,'.node.traits'),mid2)

if(ni==niche[1] & re==region[1] ) { regionsub_sparcc.node_traits=get(paste0(nm2,'.node.traits'))
regionsub_sparcc.basic_traits=get(paste0(nm2,'.basic_traits')) }
else {
regionsub_sparcc.node_traits=rbind(regionsub_sparcc.node_traits,get(paste0(nm2,'.node.traits')))
regionsub_sparcc.basic_traits=rbind(regionsub_sparcc.basic_traits,get(paste0(nm2,'.basic_traits'))) }

}
}
pdf(file='net.legend.pdf') # legend plot for zipi plots
plot.new()
legend("topright",legend=c("Negative","Positive"), col= c("tomato1","skyblue"), lty=1, lwd=2, seg.len=1, bty="n",
title="Relationship")
legend("bottomleft",legend=levels(taxon16), col=color_32[5:20], pch=19, pt.cex=1, bty="n",title="Subphylum")
dev.off()



####################################### Fig. 5  robustness plots

### 1.  random remove 
library(igraph)
library(dplyr)

datComb_Rb=NULL
for(ni in c('PA','FL'))
{
 for (re in c('Hainan','ZhongXisha','Nansha')) 
 {

if (ni =='PA') niche='PA_attach' 
else niche='Freeliving'

g=get(paste0(ni, '_', re,'_net'))
adj_matrix <- as_adjacency_matrix(g, attr = "corr", sparse = FALSE) ## should have attr 'corr' or 'weitht'

if (!isSymmetric(adj_matrix)) {
  warning("The adjacency matrix is asymmetric. Forcing symmetry now.")
  adj_matrix <- (adj_matrix + t(adj_matrix)) / 2
}

diag(adj_matrix) <- 0

# remove isolated nodes
nonzero_nodes <- which(rowSums(abs(adj_matrix)) > 0)
network_raw <- adj_matrix[nonzero_nodes, nonzero_nodes]

sp_netcom=COM[s.factor$Region2==re & s.factor$Niche==niche , rownames(network_raw)]
sp_ra=colMeans(sp_netcom)/20813
sum(rownames(network_raw)!=names(sp_ra)) # check names match, should be zero

prefix=paste0(ni,'_',re)

# use network_raw as netRaw in the function "rmsimu" 

Weighted.simu <- rmsimu(netRaw = network_raw, 
                        rm.p.list = seq(0.05, 1, by = 0.05), 
                        sp.ra = sp_ra, 
                        abundance.weighted = TRUE, 
                        nperm = 999, prefix=prefix)
Unweighted.simu <- rmsimu(netRaw = network_raw, 
                        rm.p.list = seq(0.05, 1, by = 0.05), 
                        sp.ra = sp_ra, 
                        abundance.weighted = FALSE, 
                        nperm = 999, prefix=prefix)
dat_Rb<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Region=rep(re,40),Niche=rep(ni,40))

assign(paste0('dat_Rb_',ni,re), dat_Rb) 
if(ni==  'PA' & re =='Hainan') datComb_Rb=dat_Rb  else datComb_Rb=rbind(datComb_Rb, dat_Rb)
}
}


##### 2 . target remove
library(igraph)
library(dplyr)

datComb_Rb_tgt=NULL
for(ni in c('PA','FL'))
{
 for (re in c('Hainan','ZhongXisha','Nansha')) 
 {

if (ni =='PA') niche='PA_attach' 
else niche='Freeliving'

g=get(paste0(ni, '_', re,'_net'))
adj_matrix <- as_adjacency_matrix(g, attr = "corr", sparse = FALSE) ## should have attr 'corr'

if (!isSymmetric(adj_matrix)) {
  warning("The adjacency matrix is asymmetric. Forcing symmetry now.")
  adj_matrix <- (adj_matrix + t(adj_matrix)) / 2
}

diag(adj_matrix) <- 0

# remove isolated nodes
nonzero_nodes <- which(rowSums(abs(adj_matrix)) > 0)
network_raw <- adj_matrix[nonzero_nodes, nonzero_nodes]

sp_netcom=COM[s.factor$Region2==re & s.factor$Niche==niche , rownames(network_raw)]
sp_ra=colMeans(sp_netcom)/20813
sum(rownames(network_raw)!=names(sp_ra)) # check names match, should be zero

zp <- get(paste0(paste0(ni, '_', re,'_net.zp')))  %>% arrange(desc(Zi),desc(Pi))
# keyzp <- subset(zp, Zi > 2.5 | Pi > 0.62) ###### some nets have very less keyzps like this
# module.hub<-rownames(keyzp)
# module.hub<-rownames(zp)[1:30] ###### uncommnet this to use the top 30 nodes for Zi and Pi values

## choose degree 
V(g)$degree=degree(g)
V(g)$betw=betweenness(g, directed=F, weights=NA,normalized=T)
node.attri<-as.data.frame(vertex.attributes(g)) %>% arrange(desc(degree))######## uncommnet this to target remove by nodes' degree
# node.attri<-as.data.frame(vertex.attributes(g)) %>% arrange(desc(betw))######## uncommnet this to target remove by nodes' betweenness
module.hub<-node.attri$name[1:30] 

prefix=paste0(ni,'_',re)
# use network_raw as netRaw in the function "rmsimu" 

Weighted.simu <- rmsimu2(netRaw = network_raw, 
                        rm.p.list=1:length(module.hub), 
                        keystonelist=module.hub,
                        sp.ra = sp_ra, 
                        abundance.weighted = TRUE, 
                        nperm = 999, prefix=prefix)
Unweighted.simu <- rmsimu2(netRaw = network_raw, 
                        rm.p.list=1:length(module.hub), 
                        keystonelist=module.hub,
                        sp.ra = sp_ra, 
                        abundance.weighted = FALSE, 
                        nperm = 999, prefix=prefix)

dat_Rb_tgt<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=length(module.hub)),
                 Region=rep(re,2*length(module.hub)),Niche=rep(ni,2*length(module.hub)))

assign(paste0('dat_Rb_tgt_',ni,re), dat_Rb_tgt) 
if(ni==  'PA' & re =='Hainan') datComb_Rb_tgt=dat_Rb_tgt  else datComb_Rb_tgt=rbind(datComb_Rb_tgt, dat_Rb_tgt)
}
}



##### 3. plot robustness boxplots at one specific step 
require(dplyr)
require(ggplot2)

for (ni in c('FL', 'PA')) {
for (re in c('Hainan', 'ZhongXisha', 'Nansha')) {
tgt_remainW=read.csv(paste0(getwd(),'/tgt_outputs/',ni,'_',re,'15_tgt.remains.W.csv')) ## choose remove 15 top degree nodes
tgt_remainUW=read.csv(paste0(getwd(),'/tgt_outputs/',ni,'_',re,'15_tgt.remains.UW.csv'))

rdm_remainW=read.csv(paste0(getwd(),'/rdm_outputs/',ni,'_',re,'0.5_rdm.remains.W.csv')) ## choose remove 50% nodes
rdm_remainUW=read.csv(paste0(getwd(),'/rdm_outputs/',ni,'_',re,'0.5_rdm.remains.UW.csv'))

tgt_remainW <- tgt_remainW %>% mutate(Niche=rep(ni, nrow(tgt_remainW)), Region=rep(re, nrow(tgt_remainW)), weight= rep('Weighted', nrow(tgt_remainW)) )
tgt_remainUW <- tgt_remainUW %>% mutate(Niche=rep(ni, nrow(tgt_remainUW)), Region=rep(re, nrow(tgt_remainUW)), weight= rep('Unweighted', nrow(tgt_remainUW)) )

rdm_remainW <- rdm_remainW %>% mutate(Niche=rep(ni, nrow(rdm_remainW)), Region=rep(re, nrow(rdm_remainW)), weight= rep('Weighted', nrow(rdm_remainW)) )
rdm_remainUW <- rdm_remainUW %>% mutate(Niche=rep(ni, nrow(rdm_remainUW)), Region=rep(re, nrow(rdm_remainUW)), weight= rep('Unweighted', nrow(rdm_remainUW)) )

if (ni == 'FL' & re == 'Hainan') { tgt_remains=rbind(tgt_remainW, tgt_remainUW)
rdm_remains=rbind(rdm_remainW, rdm_remainUW) } else  { tgt_remains= rbind(tgt_remains, tgt_remainW, tgt_remainUW)
rdm_remains= rbind(rdm_remains, rdm_remainW, rdm_remainUW) }
}
}

ggplot(subset(rdm_remains, weight='Weighted'), aes(x=Region, y=x)) +geom_boxplot(aes(fill=Niche))+theme_bw()+facet_wrap(~weight, scales='free')
ggsave(file='Fig. 5A.pdf',width=7,height=3)

ggplot(subset(tgt_remains, weight='Weighted'), aes(x=Region, y=x)) +geom_boxplot(aes(fill=Niche))+theme_bw()+facet_wrap(~weight, scales='free')
ggsave(file='Fig. 5B.pdf',width=7,height=3)









