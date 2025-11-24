library(Seurat)
library(ggplot2)
library(reshape2)

load("DEG_sets.Rdata")

alpha_up[alpha_up %in% feature_names]
alpha_up_present <- alpha_up[alpha_up %in% feature_names]
allmerge <- AddModuleScore(object=allmerge,features=list(alpha_up_present),name="BulkAlphaUp")

beta_up_present <- beta_up[beta_up %in% feature_names]
allmerge <- AddModuleScore(object=allmerge,features=list(beta_up_present),name="BulkBetaUp")

gamma_up_present <- gamma_up[gamma_up %in% feature_names]
allmerge <- AddModuleScore(object=allmerge,features=list(gamma_up_present),name="BulkGammaUp")

FeaturePlot(allmerge,features=c("BulkAlphaUp1","BulkBetaUp1","BulkGammaUp1"),cols=jonny_pal,order=TRUE,pt.size=1,split.by="Genotype")
ggsave("rod_genes_featplots_withbulk_genosplit.png",device="png")
ggsave("rod_genes_featplots_withbulk_genosplit.pdf",device="pdf")


library(ggplot2)
library(reshape2)

abg_genotype_table <- allmerge@meta.data[, c("Genotype", "BulkAlphaUp1", "BulkBetaUp1", "BulkGammaUp1")]
abg_genotype_table$Cell <- row.names(abg_genotype_table)

abg_melt <- melt(abg_genotype_table,id.vars = c("Cell","Genotype"),variable.name = "ABG",value.name = "ModuleScore")
abg_melt$ABG <- factor(abg_melt$ABG,levels=levels(abg_melt$ABG)[c(1,3,2)])
#shorten it
levels(abg_melt$ABG) <- c("Alpha","Beta","Gamma")

npscols <- c("#6E8C6B","#AE062D","#6DA0E7","#4E3BA3")

box_no_out <- ggplot(abg_melt) +
  geom_boxplot(aes(Genotype,ModuleScore,fill=Genotype),outlier.shape=NA) +
  facet_grid(~ABG) +
  guides(fill=FALSE) +
  scale_fill_manual(values=npscols) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  coord_cartesian(ylim = c(-0.25,0.75)) +
  labs(y="Gene Module Score", x=NULL)


box_all <- ggplot(abg_melt) +
  geom_boxplot(aes(Genotype,ModuleScore,fill=Genotype),outlier.color = "gray25") +
  facet_grid(~ABG) +
  guides(fill=FALSE) +
  scale_fill_manual(values=npscols) +
  theme_bw() +
  theme(panel.grid = element_blank())+
  labs(y="Gene Module Score", x=NULL)

library(patchwork)
box_all/box_no_out +  plot_layout(heights = c(2,1))

#675x650
ggsave("geno_split_abg_boxplots.pdf",device = "pdf", width = 6.75, height=6.50)

#stats
atable <-abg_genotype_table[c(1,2)]
btable <-abg_genotype_table[c(1,4)]
gtable <-abg_genotype_table[c(1,3)]

amodel <- lm(BulkAlphaUp1 ~ Genotype, data=atable)
bmodel <- lm(BulkBetaUp1 ~ Genotype, data=btable)
gmodel <- lm(BulkGammaUp1 ~ Genotype, data=gtable)

anova(amodel)
anova(bmodel)
anova(gmodel)
#all significant

library(lsmeans)
library(multcompView)
library(multcomp)

a_ls <- lsmeans(amodel,
                pairwise ~ Genotype,
                adjust = "tukey")

b_ls <- lsmeans(bmodel,
                pairwise ~ Genotype,
                adjust = "tukey")

g_ls <- lsmeans(gmodel,
                pairwise ~ Genotype,
                adjust = "tukey")


