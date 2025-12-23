library(RVAideMemoire)
library(lmPerm)
library(dplyr)
library(agricolae)
library(multcompView)


## 1, ENV data comparisons


aovp.env=vector('list',8)
for (i in 1:8)
{
ancova_model=aovp(ENV[,i] ~ Region2*Layer, ENV)
pairwise_results_adj1 <- pairwise.perm.t.test(subset(ENV, Layer=='Bottom')[,i],subset(ENV, Layer=='Bottom')$Region2)
pairwise_results_adj2 <- pairwise.perm.t.test(subset(ENV, Layer=='Surface')[,i],subset(ENV, Layer=='Surface')$Region2)
aovp.env[[i]]=list(summary(ancova_model),bottom=pairwise_results_adj1, surface=pairwise_results_adj2)
names(aovp.env)[i]=colnames(ENV)[i]
}
print(aovp.env)


aovp.env2=vector('list',length(colnames(ENV))-5)
for (i in 1:8)
{
ancova_model=aovp(ENV[,i] ~ Region2*Layer, ENV)
pairwise_results_adj1 <- pairwise.perm.t.test(subset(ENV, Region2=='Hainan')[,i],subset(ENV, Region2=='Hainan')$Layer)
pairwise_results_adj2 <- pairwise.perm.t.test(subset(ENV, Region2=='ZhongXisha')[,i],subset(ENV, Region2=='ZhongXisha')$Layer)
pairwise_results_adj3 <- pairwise.perm.t.test(subset(ENV, Region2=='Nansha')[,i],subset(ENV, Region2=='Nansha')$Layer)

aovp.env2[[i]]=list(summary(ancova_model),Hainan=pairwise_results_adj1, ZhongXisha=pairwise_results_adj2, Nansha=pairwise_results_adj3)
names(aovp.env2)[i]=colnames(ENV)[i]
}
print(aovp.env2)



## 2, alpha  diversity comparisons
aovp.alpha=vector('list',2)
for (i in 1:2)
{
ancova_model=aovp(alphaf[,i] ~ Region2+Niche+Layer, alphaf)
subset1.1=subset(alphaf, Layer=='Bottom' & Niche=='Freeliving')
subset1.2=subset(alphaf, Layer=='Surface' & Niche=='Freeliving')
pairwise_results_adj1.1 <- pairwise.perm.t.test(subset1.1[,i], paste0(subset1.1$Region2,'_', subset1.1$Layer))
pairwise_results_adj1.2 <- pairwise.perm.t.test(subset1.2[,i], paste0(subset1.2$Region2,'_', subset1.2$Layer))

subset2.1=subset(alphaf, Layer=='Bottom' & Niche=='PA_attach')
subset2.2=subset(alphaf, Layer=='Surface' & Niche=='PA_attach')
pairwise_results_adj2.1 <- pairwise.perm.t.test(subset2.1[,i], paste0(subset2.1$Region2,'_', subset2.1$Layer))
pairwise_results_adj2.2 <- pairwise.perm.t.test(subset2.2[,i], paste0(subset2.2$Region2,'_', subset2.2$Layer))

subset3.1=subset(alphaf, Layer=='Bottom' & Region2=='Hainan')
subset3.2=subset(alphaf, Layer=='Surface' & Region2=='Hainan')

subset3.3=subset(alphaf, Layer=='Bottom' & Region2=='ZhongXisha')
subset3.4=subset(alphaf, Layer=='Surface' & Region2=='ZhongXisha')

subset3.5=subset(alphaf, Layer=='Bottom' & Region2=='Nansha')
subset3.6=subset(alphaf, Layer=='Surface' & Region2=='Nansha')


pairwise_results_adj3.1 <- pairwise.perm.t.test(subset3.1[,i], paste0(subset3.1$Niche))
pairwise_results_adj3.2 <- pairwise.perm.t.test(subset3.2[,i], paste0(subset3.2$Niche))
pairwise_results_adj3.3 <- pairwise.perm.t.test(subset3.3[,i], paste0(subset3.3$Niche))
pairwise_results_adj3.4 <- pairwise.perm.t.test(subset3.4[,i], paste0(subset3.4$Niche))
pairwise_results_adj3.5 <- pairwise.perm.t.test(subset3.5[,i], paste0(subset3.5$Niche))
pairwise_results_adj3.6 <- pairwise.perm.t.test(subset3.6[,i], paste0(subset3.6$Niche))



aovp.alpha[[i]]=list(summary(ancova_model), FL.b=pairwise_results_adj1.1, FL.s=pairwise_results_adj1.2, PA.b=pairwise_results_adj2.1, PA.s=pairwise_results_adj2.2, Hainan.b=pairwise_results_adj3.1,  Hainan.s=pairwise_results_adj3.2, ZhongXisha.b=pairwise_results_adj3.3, ZhongXisha.s=pairwise_results_adj3.4, Nansha.b=pairwise_results_adj3.5, Nansha.s=pairwise_results_adj3.6)
names(aovp.alpha)[i]=colnames(alphaf)[i]
}
print(aovp.alpha)


aovp.alpha2=vector('list',2)
for (i in 1:2)
{
ancova_model=aovp(alphaf[,i] ~ Region2*Layer, alphaf)
pairwise_results_adj1.1 <- pairwise.perm.t.test(subset(alphaf, Region2=='Hainan'& Niche=="Freeliving")[,i],subset(alphaf, Region2=='Hainan'& Niche=="Freeliving")$Layer)
pairwise_results_adj1.2 <- pairwise.perm.t.test(subset(alphaf, Region2=='Hainan'& Niche=="PA_attach")[,i],subset(alphaf, Region2=='Hainan'& Niche=="PA_attach")$Layer)
pairwise_results_adj2.1 <- pairwise.perm.t.test(subset(alphaf, Region2=='ZhongXisha'& Niche=="Freeliving")[,i],subset(alphaf, Region2=='ZhongXisha'& Niche=="Freeliving")$Layer)
pairwise_results_adj2.2 <- pairwise.perm.t.test(subset(alphaf, Region2=='ZhongXisha'& Niche=="PA_attach")[,i],subset(alphaf, Region2=='ZhongXisha'& Niche=="PA_attach")$Layer)
pairwise_results_adj3.1 <- pairwise.perm.t.test(subset(alphaf, Region2=='Nansha'& Niche=="Freeliving")[,i],subset(alphaf, Region2=='Nansha'& Niche=="Freeliving")$Layer)
pairwise_results_adj3.2 <- pairwise.perm.t.test(subset(alphaf, Region2=='Nansha'& Niche=="PA_attach")[,i],subset(alphaf, Region2=='Nansha'& Niche=="PA_attach")$Layer)

aovp.alpha2[[i]]=list(summary(ancova_model),Hainan.fl=pairwise_results_adj1.1, Hainan.pa=pairwise_results_adj1.2, ZhongXisha.fl=pairwise_results_adj2.1, ZhongXisha.pa=pairwise_results_adj2.2, Nansha.fl=pairwise_results_adj3.1, Nansha.pa=pairwise_results_adj3.2)
names(aovp.alpha2)[i]=colnames(alphaf)[i]
}
print(aovp.alpha2)




## 3, beta  diversity comparisons
library(vegan)
com=t(otu_tb)
adonis2((com)~Region2+Niche+Layer, data=metadata,by='m')
adonis2((com[metadata$Niche=='Freeliving', ])~ Region2+Layer, data=metadata[metadata$Niche=='Freeliving', ],by='m')
adonis2((com[metadata$Niche=='PA_attach', ])~ Region2+Layer, data=metadata[metadata$Niche=='PA_attach', ],by='m')

adonis2((com[metadata$Niche=='Freeliving', ])~ Niche+Layer, data=metadata[metadata$Niche=='Freeliving', ],by='m')
adonis2((com[metadata$Niche=='PA_attach', ])~ Region2+Layer, data=metadata[metadata$Niche=='PA_attach', ],by='m')


library(veganEx)
adonis.pairwise(
  (com[metadata$Niche=='Freeliving', ]),
  metadata[metadata$Niche=='Freeliving', 'Region2'],
  sim.function = "vegdist",
  sim.method = "bray",
  p.adjust.m = "fdr",
  perm = 999
)

adonis.pairwise(
  (com[metadata$Niche=='PA_attach', ]),
  metadata[metadata$Niche=='PA_attach', 'Region2'],
  sim.function = "vegdist",
  sim.method = "bray",
  p.adjust.m = "fdr",
  perm = 999
)


## 4, phyla abundance comparisons
aovp.phyla=vector('list', (length(taxa_forpl)-6) )
for (i in 1:(length(taxa_forpl)-6) )
{
ancova_model=aovp(taxa_forpl[,i] ~ Region2+Niche+Layer, taxa_forpl)
subset1=subset(taxa_forpl, Niche=='Freeliving')
pairwise_results_adj1 <- pairwise.perm.t.test(subset1[,i], paste(subset1$Region2))

subset2=subset(taxa_forpl, Niche=='PA_attach')
pairwise_results_adj2 <- pairwise.perm.t.test(subset2[,i], paste(subset2$Region2))

subset3=subset(taxa_forpl, Region2=='Hainan')
pairwise_results_adj3 <- pairwise.perm.t.test(subset3[,i], paste(subset3$Niche))

subset4=subset(taxa_forpl, Region2=='ZhongXisha')
pairwise_results_adj4 <- pairwise.perm.t.test(subset4[,i], paste(subset4$Niche))

subset5=subset(taxa_forpl, Region2=='Nansha')
pairwise_results_adj5 <- pairwise.perm.t.test(subset5[,i], paste(subset5$Niche))

aovp.phyla[[i]]=list(summary(ancova_model), FL=pairwise_results_adj1, PA=pairwise_results_adj2,
Hainan=pairwise_results_adj3, ZhongXisha=pairwise_results_adj4, Nansha=pairwise_results_adj5)
names(aovp.phyla)[i]=colnames(taxa_forpl)[i]
}
print(aovp.phyla)


## 5, network node-level traits comparisons
aovp.netraits=vector('list')
for (i in 1:8 )
{
ancova_model=aovp(nodetr2[,i] ~ Region*Niche, nodetr2)
subset1=subset(nodetr2, Niche=='FL')
pairwise_results_adj1 <- pairwise.perm.t.test(subset1[,i], paste(subset1$Region))
subset2=subset(nodetr2, Niche=='PA')
pairwise_results_adj2 <- pairwise.perm.t.test(subset2[,i], paste(subset2$Region))
aovp.netraits[[i]]=list(summary(ancova_model), FL=pairwise_results_adj1, PA=pairwise_results_adj2)
names(aovp.netraits)[i]=colnames(nodetr2)[i]
}


aovp.netraits2=vector('list')
for (i in 1:8 )
{
ancova_model=aovp(nodetr2[,i] ~ Region*Niche, nodetr2)
subset1=subset(nodetr2, Region=='Hainan')
pairwise_results_adj1 <- pairwise.perm.t.test(subset1[,i], paste(subset1$Niche))
subset2=subset(nodetr2, Region=='ZhongXisha')
pairwise_results_adj2 <- pairwise.perm.t.test(subset2[,i], paste(subset2$Niche))
subset3=subset(nodetr2, Region=='Nansha')
pairwise_results_adj3 <- pairwise.perm.t.test(subset3[,i], paste(subset3$Niche))
aovp.netraits2[[i]]=list(summary(ancova_model), Hainan=pairwise_results_adj1, ZhongXisha=pairwise_results_adj2, Nansha=pairwise_results_adj3)
names(aovp.netraits2)[i]=colnames(nodetr2)[i]
}




## 6, network robustness comparisons

# random remove
subset(rdm_remains, weight=='Weighted') -> rdm_remains_W

ancova_model=aovp(x ~ Region*Niche, rdm_remains_W)
subset1=subset(rdm_remains_W, Niche=='FL')
pairwise_results_adj1 <- pairwise.perm.t.test(subset1[,'x'], paste(subset1$Region))
subset2=subset(rdm_remains_W, Niche=='PA')
pairwise_results_adj2 <- pairwise.perm.t.test(subset2[,'x'], paste(subset2$Region))

subset3=subset(rdm_remains_W, Region=='Hainan')
pairwise_results_adj3 <- pairwise.perm.t.test(subset3[,'x'], paste(subset3$Niche))
subset4=subset(rdm_remains_W, Region=='ZhongXisha')
pairwise_results_adj4 <- pairwise.perm.t.test(subset4[,'x'], paste(subset4$Niche))
subset5=subset(rdm_remains_W, Region=='Nansha')
pairwise_results_adj5 <- pairwise.perm.t.test(subset5[,'x'], paste(subset5$Niche))

rdm_robust_stat=list(summary(ancova_model), FL=pairwise_results_adj1, PA=pairwise_results_adj2,
Hainan=pairwise_results_adj3, ZhongXisha=pairwise_results_adj4, Nansha=pairwise_results_adj5)
names(rdm_robust_stat)=c('anova', 'FL', 'PA', 'Hainan','ZhongXisha','Nansha')

# target remove
subset(tgt_remains, weight=='Weighted') -> tgt_remains_W

ancova_model=aovp(x ~ Region*Niche, tgt_remains_W)
subset1=subset(tgt_remains_W, Niche=='FL')
pairwise_results_adj1 <- pairwise.perm.t.test(subset1[,'x'], paste(subset1$Region))
subset2=subset(tgt_remains_W, Niche=='PA')
pairwise_results_adj2 <- pairwise.perm.t.test(subset2[,'x'], paste(subset2$Region))

subset3=subset(tgt_remains_W, Region=='Hainan')
pairwise_results_adj3 <- pairwise.perm.t.test(subset3[,'x'], paste(subset3$Niche))
subset4=subset(tgt_remains_W, Region=='ZhongXisha')
pairwise_results_adj4 <- pairwise.perm.t.test(subset4[,'x'], paste(subset4$Niche))
subset5=subset(tgt_remains_W, Region=='Nansha')
pairwise_results_adj5 <- pairwise.perm.t.test(subset5[,'x'], paste(subset5$Niche))

tgt_robust_stat=list(summary(ancova_model), FL=pairwise_results_adj1, PA=pairwise_results_adj2,
Hainan=pairwise_results_adj3, ZhongXisha=pairwise_results_adj4, Nansha=pairwise_results_adj5)
names(tgt_robust_stat)=c('anova', 'FL', 'PA', 'Hainan','ZhongXisha','Nansha')


