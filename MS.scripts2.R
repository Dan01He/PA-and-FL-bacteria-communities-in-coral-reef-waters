
######################################## Fig. S1 --ggvenn
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("ggvenn")) install.packages("ggvenn")
library(ggplot2)
library(ggvenn)

com=t(otu_table)
Hainan_FLcomm=as.data.frame( com[metadata$Region2=='Hainan' & metadata$Niche=='Freeliving',] ); Hainan_FLcomm=Hainan_FLcomm[,colSums(Hainan_FLcomm)>0]
Hainan_PAcomm=as.data.frame( com[metadata$Region2=='Hainan' & metadata$Niche=='PA_attach',] ); Hainan_PAcomm=Hainan_PAcomm[,colSums(Hainan_PAcomm)>0]

ZhongXisha_FLcomm=as.data.frame( com[metadata$Region2=='ZhongXisha' & metadata$Niche=='Freeliving',] ); ZhongXisha_FLcomm=ZhongXisha_FLcomm[,colSums(ZhongXisha_FLcomm)>0]
ZhongXisha_PAcomm=as.data.frame( com[metadata$Region2=='ZhongXisha' & metadata$Niche=='PA_attach',] ); ZhongXisha_PAcomm=ZhongXisha_PAcomm[,colSums(ZhongXisha_PAcomm)>0]

Nansha_FLcomm=as.data.frame( com[metadata$Region2=='Nansha' & metadata$Niche=='Freeliving',] ); Nansha_FLcomm=Nansha_FLcomm[,colSums(Nansha_FLcomm)>0]
Nansha_PAcomm=as.data.frame( com[metadata$Region2=='Nansha' & metadata$Niche=='PA_attach',] ); Nansha_PAcomm=Nansha_PAcomm[,colSums(Nansha_PAcomm)>0]

# region and niche specific 
niche_data1 <- list(Hainan_FL = colnames(Hainan_FLcomm),Hainan_PA = colnames(Hainan_PAcomm))
niche_data2 <- list(ZhongXisha_FL = colnames(ZhongXisha_FLcomm),ZhongXisha_PA = colnames(ZhongXisha_PAcomm))
niche_data3 <- list(Nansha_FL = colnames(Nansha_FLcomm),Nansha_PA = colnames(Nansha_PAcomm))

p1=ggvenn(niche_data1, fill_color = c("lightblue", "lightgreen"), stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE) 
p2=ggvenn(niche_data2, fill_color = c("lightblue", "lightgreen"), stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE) 
p3=ggvenn(niche_data3, fill_color = c("lightblue", "lightgreen"), stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE) 
plt <- ggarrange( p1, p2, p3, nrow = 3,ncol=2)
ggsave(plt, file="Fig. S1.pdf", width = 8, height = 7)  # generating the basic plot that may need some outlook modification  



######################################## Fig. S2 -- Phylum Env correlation plot
PA_attach_taxa_forpl=subset(taxa_forpl,Niche=='PA_attach')
Freeliving_taxa_forpl=subset(taxa_forpl,Niche=='Freeliving')
for (ni in c('PA_attach', 'Freeliving')) {
taxa_forpl=get(paste0(ni,'_taxa_forpl'))
PhylumRL.sel=taxa_forpl[,1:16] 
env_data_sel=env_data[rownames(PhylumRL.sel),]

cor_results <- data.frame() 
# correlation analyses for the env_data and PhylumRL.sel
for (env_col in colnames(env_data_sel)) { 
 for (abun_col in colnames(PhylumRL.sel)) { 
 cor_test_result <- cor.test(env_data_sel[, env_col], PhylumRL.sel[, abun_col],method='spearman') 
 cor_results <- rbind(cor_results, data.frame( 
 env_col = env_col, 
 abun_col = abun_col, 
 correlation = cor_test_result$estimate, 
 p.value = cor_test_result$p.value 
 )) 
 } 
} 

cor_results <- cor_results %>% filter(env_col != abun_col) 

p1=ggplot(cor_results, aes(x = abun_col, y = env_col, fill = correlation)) + 
 geom_tile() + 
 scale_fill_gradient2( 
 name = "Correlation Coefficient", 
 limits = c(-1, 1), 
 low = "blue", 
 mid = "white", 
 high = "red", 
 midpoint = 0 
 ) + 
 geom_text( 
 data = cor_results %>% filter(p.value < 0.05), 
 aes(label = round(correlation, 2)), 
 color = "black" 
 ) + 
 theme_minimal() + 
 theme( 
 axis.text.x = element_text(angle = 45, hjust = 1), 
 axis.title = element_blank(), 
 legend.position = "bottom" 
 ) + 
 labs(fill = "Correlation") 
assign( paste0(ni,'.PhylumEnv.corplt'),p1)
}
plt <- ggarrange(PA_attach.PhylumEnv.corplt, Freeliving.PhylumEnv.corplt, nrow = 2,ncol=1)
ggsave(plt, file="Fig. S2.pdf", width = 5, height = 7)  # generating the basic plot that may need some outlook modification 




######################################## Fig. S3  --assembly processes
library(vegan)
library(picante)
library(iCAMP)

wd <- getwd()  
region=c('Hainan', 'ZhongXisha', 'Nansha')
style=c('PA_attach','Freeliving')
imp_icamp=NULL

for (re in region) {
for (st in style ) {
com=t(otu_table)
com.sel=com[metadata$Region2== re & metadata$Niche==st,]
com.sel=com.sel[,colSums(com.sel) >= 2] # 过滤低丰度OTU


save.wd <- paste0(wd, "/icamp_out/", re, st )
if(! dir.exists(save.wd)) dir.create(save.wd, recursive=T)
tre=tree
matched <- match.phylo.comm2(tre, com.sel)
comm <- sqrt(matched$comm)
phylo <- matched$phy

clas <- otu_tax[colnames(comm),]
treat <- metadata[rownames(comm),]
env <- env_data[rownames(comm),]


# iCAMP big function

pd.big <- pdist.big(tree = phylo, wd = save.wd, nworker = 20)
icres <- icamp.big(comm = comm, pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd, 
                   rand = 1000, tree = phylo, prefix = paste0(re,'_',st), ds = 0.2, pd.cut = NA, sp.check = TRUE, 
                   phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all", 
                   phylo.metric = "bMPD", sig.index = "Confidence", bin.size.limit = 5, 
                   nworker = 20, rtree.save = FALSE, detail.save = TRUE, 
                   qp.save = TRUE, detail.null = FALSE, ignore.zero = TRUE, 
                   output.wd = save.wd, correct.special = TRUE, unit.sum = rowSums(comm), 
                   special.method = "depend", ses.cut = 1.96, rc.cut = 0.95, 
                   conf.cut = 0.975, omit.option = "no", meta.ab = NULL) 


# iCAMP bin 
icbin <- icamp.bins(icamp.detail = icres$detail, treat = treat[4], clas = clas, silent = FALSE, boot = TRUE, 
                    rand.time = 1000, between.group = FALSE) 

imp_icbin <- icbin$Pt
imp_icbin$Group <- rownames(imp_icbin)
imp_icbin$group=paste0(re,st)

if(re=='Hainan' & st =='PA_attach' ) imp_icamp = imp_icbin
   else imp_icamp=rbind(imp_icamp, imp_icbin)
}
}


library(pacman)
p_load(ggplot2, reshape2, dplyr)

imp_icamp[, 2:3]=NULL
imp_icamp.m=melt(imp_icamp, id= c(1,7:9))
imp_icamp.m$value=as.numeric(imp_icamp.m$value)

ggplot(subset(imp_icamp.m,Method=='CbMPDiCbraya'), aes(x = 2, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity", position="stack") +
  coord_polar("y", start = 0) +geom_text(aes(y = value, label = round(100*value,2)), color = "white")+
  xlim(0.5, 2.5)+theme_void()+
  facet_grid(re~st)
ggsave(file="Fig. S3.pdf", width = 7, height = 4.5)  # generating the basic plot that may need some outlook modification 




######################################## Fig. S4 .  top bins contributing to the HoS process and correlating with water properties
sort(subset(HainanPA_attach_icbin$BRPtk, Process=='HoS')[-(1:4)], decreasing=T)
BRPtk_index=order(subset(HainanPA_attach_icbin$BRPtk, Process=='HoS')[-(1:4)], decreasing=T)
names(sort(subset(HainanPA_attach_icbin$BRPtk, Process=='HoS')[-(1:4)], decreasing=T)) -> binnames
gsub('bin', 'Bin', binnames) -> DecreasBinames; rm(binnames)

HainanPA_attach_icbin$Bin.TopClass->HainanPA_BinTopClass
HainanPA_BinTopClass[BRPtk_index[1:5],]->HainanPA_BinTop5Hos
HainanPA_BinTopClass[BRPtk_index[1:10],]->HainanPA_BinTop10Hos

com=t(otu_table)
com.sel=com[metadata$Niche=='PA_attach'  ,]
# com.sel=com[metadata$Niche=='PA_attach' &  metadata$Region2=='Hainan',]
com.sel=com.sel[,colSums(com.sel)>0]

taxa_table <- BAC.tax[colnames(com.sel),]
Genusbin <- taxa_table[ taxa_table$Genus%in%HainanPA_BinTop5Hos$TopTaxon.Genus[c(1:3,5)] | (taxa_table$Genus=='uncultured' &taxa_table$Family=='Cryomorphaceae'),]
com.sel=com.sel[,rownames(Genusbin)]
com.sel=com.sel[,colSums(com.sel)>0]

Genusbin_abundance <- aggregate(. ~ Genus, data = data.frame(Genus=Genusbin$Genus, t(com.sel)), sum)
Genusbin_abundance[,1]->rownames(Genusbin_abundance)
Genusbin_abundance[,1]=NULL
Genusbin_abundance=Genusbin_abundance[order(rowSums(Genusbin_abundance),decreasing=T),]
Genusbin_abundance=t(Genusbin_abundance)
Genusbin_rla=Genusbin_abundance/20813*100
env_data_PA=env_data[rownames(Genusbin_abundance),]

library(corrplot)
corBinPA=cor(cbind(Genusbin_rla,env_data_PA),method="spearman")

cor.mtest.spearman <- function(mat, method="spearman"){
 mat <- as.matrix(mat)
 n <- ncol(mat)
 p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
 diag(p.mat) <- 0
 diag(lowCI.mat) <- diag(uppCI.mat) <- 1
 for(i in 1:(n-1)){
 for(j in (i+1):n){
 tmp <- cor.test(mat[,i], mat[,j], method=method)
 p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
 lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$p.value
 uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$p.value
 }
 }
 return(list(p.mat, lowCI.mat, uppCI.mat))
}
resPA<- cor.mtest.spearman(cbind(Genusbin_rla,env_data_PA))

col3 <- colorRampPalette(c("green", "white", "blue"))
rownames(corBinPA)
corM=corBinPA[5:12, 1:4]

pmat=resPA[[1]][5:12, 1:4]; dimnames(corM)->dimnames(pmat)
pdf(file='Fig. S4.pdf',width=6,height=8)
corrplot(corM, p.mat = pmat, col=col3(100), method="square", sig.level = 0.05, addCoef.col="grey")
dev.off()



######################################## Fig. S5 node-level network traits
setwd(wd)
nodetr=regionsub_sparcc.node_traits
nodetr2 <- nodetr
nodetr2$transitivity[is.na(nodetr2$transitivity)]=0
round(Ca <- cor(as.matrix(scale(nodetr2[c(1:8)]))), 2) # the 4th column transitivity has NA values
heatmap(Ca, symm = TRUE, margins = c(6,6)) # with reorder()
c(1:8)[(heatmap(Ca, symm = TRUE, margins = c(6,6)))$rowInd]->nodetr2_Ind

nodelvl.m=reshape2::melt(nodetr2[c(nodetr2_Ind,9:11)],id=c(9:11)) # 9:11, and 8:10 may be changed
head(nodelvl.m)

nodelvl.mse=data.frame( aggregate(data=nodelvl.m,value~variable+Region+Niche,mean),
se=(aggregate(data=nodelvl.m,value~variable+Region+Niche,function(x){se=sd(x)/sqrt(length(x))}))$value )

nodelvl.mse3=data.frame( aggregate(data=nodelvl.m,value~variable+Region+Niche+subphylum,mean),
se=(aggregate(data=nodelvl.m,value~variable+Region+Niche+subphylum,function(x){se=sd(x)/sqrt(length(x))}))$value )

nodelvl.mse$Region=factor(nodelvl.mse$Region,levels=c('Hainan','ZhongXisha','Nansha'))
pl.nodelvl.m=ggplot(subset(nodelvl.mse,variable!='info.centrality' & variable!='transitivity'), aes(x=Region, y=value,fill=Niche)) +
	geom_bar(stat="identity",position=position_dodge(width=0.9)) +
	geom_errorbar( aes(ymin=value-se, ymax=value+se),position=position_dodge(width=0.9),width=0.5)+
	facet_wrap(~variable,scales="free")+theme_bw()+
	theme(line = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(pl.nodelvl.m,file='Fig. S5.pdf.pdf',width=10,height=7,dpi=600)




######################################## Fig. S6 -- link env with network traits
## each sample
library(pacman)
p_load(igraph)
p_load(vegan)
p_load(info.centrality)
p_load(linkET)

comFL=read.delim('FL_otbtp2k.tsv'); rownames(comFL)=comFL[,1]; comFL=t(comFL[,-1])
comPA=read.delim('FL_otbtp2k.tsv'); rownames(comPA)=comPA[,1]; comPA=t(comPA[,-1])

sparcc_samplenet.basic_traits=NULL
for (ni in c('FL','PA')) {
comni=get(paste0('com',ni))
netni=get(paste0(ni,'_Sparcc.net'))
for (i in rownames(comni)) {
colnames(comni)[comni[i,]>0]->vids
which(as_ids(V(netni))%in%vids)->vidsindex
net=paste(i,'net',sep='_')
assign(net,subgraph(netni,vidsindex))

g=get(net)

# net-level traits
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
g18 <- i


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
assign(paste0(i,'.basic_traits'),mid1)


if(ni==niche[1] & i==rownames(comni)[1]) {
sparcc_samplenet.basic_traits=get(paste0(i,'.basic_traits')) }
else {
sparcc_samplenet.basic_traits=rbind(sparcc_samplenet.basic_traits,get(paste0(i,'.basic_traits'))) }

}
}
rownames(sparcc_samplenet.basic_traits)=sparcc_samplenet.basic_traits$group
sparcc_samplenet.basic_traits$group=NULL
sparcc_samplenet.basic_traits=sparcc_samplenet.basic_traits[rownames(metadata),]
sparcc_samplenet_FLtraits=sparcc_samplenet.basic_traits[metadata$Niche=='Freeliving', -4]
sparcc_samplenet_PAtraits=sparcc_samplenet.basic_traits[metadata$Niche=='PA_attach', -4]
env_data_FL=env_data[rownames(sparcc_samplenet_FLtraits),]
env_data_PA=env_data[rownames(sparcc_samplenet_PAtraits),]

## keystone species composition and network traits correlate with env
com=t(otu_table)
unique(c(rownames(FL_Nansha_net.keynodes.zp), rownames(FL_Hainan_net.keynodes.zp), 
rownames(FL_ZhongXisha_net.keynodes.zp))) -> FL.keynodesid
FLcom_keynode=com[metadata$Niche=='Freeliving',FL.keynodesid]

unique(c(rownames(PA_Nansha_net.keynodes.zp), rownames(PA_Hainan_net.keynodes.zp), 
rownames(PA_ZhongXisha_net.keynodes.zp))) -> PA.keynodesid
PAcom_keynode=com[metadata$Niche=='PA_attach',PA.keynodesid]

# dim(FLcom_keynode)
# dim(sparcc_samplenet_FLtraits)
# anyNA(vegan::decostand(FLcom_keynode,'chi.square'))
# anyNA(scale(sparcc_samplenet_FLtraits))

 for (ni in c('FL','PA')) {
 NIcom_keynode=get(paste0(ni,'com_keynode'))
 env_dataNI=get(paste0('env_data_',ni))
 nettraitsNI=get(paste0("sparcc_samplenet_",ni,"traits"))

mantel_com <- mantel_test(sqrt(NIcom_keynode), env_dataNI
 ,spec_select = list(keystone_species = 1:7), 
 spec_dist = dist_func(.FUN = "vegdist", method = "bray")) %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
 labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")), 
 pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), 
 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  ) 

mantel_net <- mantel_test(scale(nettraitsNI), env_dataNI
 ,spec_select = list(nettraits = 1:4), 
 spec_dist = dist_func(.FUN = "vegdist", method = "euclidean")) %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
 labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")), 
 pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), 
 labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  ) 
 
mantel_comb=rbind(mantel_com, mantel_net)

p1=linkET::qcorrplot(correlate(scale(env_dataNI)), type = "lower", diag = FALSE) +
  geom_square() +linkET::geom_mark(sep = '\n',size = 1.8, sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05,color="white") +
  linkET::geom_couple(aes(colour = pd, size = rd), 
              data = mantel_comb, 
			  # data = mantel_gower, 
              curvature = linkET::nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = linkET::color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(color = "black"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
assign(paste0('corplt',ni), p1)
}
plt <- ggarrange(corpltFL, corpltPA, nrow = 2,ncol=1, labels='AUTO')
ggsave(plt, file="Fig. S6.pdf", width = 6, height = 8)  # generating the basic plot that may need some outlook modification 



write.csv(df_vbedge, file='df_vbedge.csv')
