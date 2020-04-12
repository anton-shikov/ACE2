library(ggplot2)
library(reshape2)
library(wesanderson)
library("ggsci")
names(wes_palettes)
library(dplyr)

af_data = read.table('/home/anton/ACE2_gnomad_all_extracted.tsv', header=T, sep='\t', stringsAsFactors = F)
mis_dat <- af_data[af_data$Type=='missense_variant',]
cont_dat_1 <- af_data[af_data$Type!='intron_variant',]
cont_dat_2 <- cont_dat_1[cont_dat_1$Type!='downstream_gene_variant',]
cont_dat_3 <- cont_dat_2[cont_dat_2$Type!='3_prime_UTR_variant',]
cont_dat_4 <- cont_dat_3[cont_dat_3$Type!='5_prime_UTR_variant',]
cont_dat_5 <- cont_dat_4[cont_dat_4$Type!='upstream_gene_variant',]
cont_dat_6 <- cont_dat_5[cont_dat_5$Type!='intron_variant&non_coding_transcript_variant',]
cont_dat_7 <- cont_dat_6[cont_dat_6$Type=='missense_variant',]


table(af_data$Effect)
table(cont_dat_6$Effect)
table(cont_dat_6$Type)
table(af_data$Type)

AF_sums_dat=data.frame(AC_sum=colSums(af_data[,grepl('AN', colnames(af_data))]),
                       AC_mean=colMeans(af_data[,grepl('AN', colnames(af_data))]),
                       AN_count=colSums(af_data[,grepl('AN', colnames(af_data))]>0),
                       AC_count=colSums(af_data[,grepl('AC', colnames(af_data))]>0),
                       AC_mis=colSums(mis_dat[,grepl('AC', colnames(mis_dat))]),
                       AC_cont_sum=colSums(cont_dat_6[,grepl('AN', colnames(cont_dat_6))]),
                       AC_cont_mean=colMeans(cont_dat_6[,grepl('AN', colnames(cont_dat_6))]),
                       AN_cont_count=colSums(cont_dat_6[,grepl('AN', colnames(af_data))]>0))

AF_sums_dat$mis_proportion <- AF_sums_dat$AC_mis/AF_sums_dat$AC_cont_sum
AF_sums_dat <- AF_sums_dat[AF_sums_dat$AN_count>619,]


ggplot(cont_dat_6, aes(x=Effect, fill=Effect))+ 
  geom_bar(stat="count",col='black') +scale_fill_lancet()+
  theme_bw() +guides(fill=guide_legend(title="Variant effect"))+ ylab('Number of variants')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = c(0.57, 0.86),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))


ggplot(AF_sums_dat, aes(x=AC_cont_mean, y=rownames(AF_sums_dat),fill=rownames(AF_sums_dat))) + 
    geom_bar(stat="identity",col='black') + scale_fill_lancet(labels=c('AFR','EAS','FIN','NFE','NFE_EST','NFE_NWE','NFE_ONF','NFE_SEU'))+
   theme_bw() +guides(fill=guide_legend(title="Population"))+ xlab('Mean allele count')+
  theme( axis.text.y = element_blank(),
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = c(0.77, 0.75),
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))


ggplot(AF_sums_dat, aes(x=mis_proportion, y=rownames(AF_sums_dat),fill=rownames(AF_sums_dat))) + 
  geom_bar(stat="identity",col='black') + scale_fill_lancet()+
  theme_bw() + xlab('Missence proportion')+
  theme( axis.text.y = element_blank(),
         axis.title.y=element_blank(),
         axis.ticks.y=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = 'None',
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                    size=14))


af_data_filt = af_data[,c(2,9,12,15,18,30,33,36,39)] 
af_data_mis_filt = cont_dat_6[,c(2,9,12,15,18,30,33,36,39)]
af_data_mis_filt = cont_dat_7[,c(2,9,12,15,18,30,33,36,39)]
#af_data_mis_filt = cont_dat_6[,c(2,9,12,15,18,30,36,39)]
#af_data_mis_filt = af_data_mis_filt[rowSums(af_data_mis_filt[, c(2:8)] > 0) > 3, ]

af_data_mis_filt_seu = af_data %>% filter(AN_nfe_seu>1000) %>% filter(AF_eas >0) %>% filter(AF_afr >0) #60
af_data_mis_filt_seu = af_data_mis_filt_seu[,c(2,9,12,15,18,30,33,36,39)]

pcaout = prcomp(as.matrix(t(af_data_mis_filt[, c(2:9)])))
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-PC1, y=PC2, fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() + xlab('PC1 - 78% of variance') +guides(fill=guide_legend(title="Population"))+ ylab('PC2 - 16% of variance')+scale_fill_lancet(labels=c('AFR','EAS','FIN','NFE','NFE_EST','NFE_NWE','NFE_ONF','NFE_SEU'))+
  theme( axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14))


# Same on variants observed in all populations
rs_filt <- af_data[af_data$rsID %in% af_data_common$rsID,]
af_rs_filt <- af_data[af_data$AF_nfe %in% af_data_common$AF_nfe,]
af_data_common = af_data_filt[rowSums(af_data_filt[, c(2:9)] > 0) > 7, ]

af_data_common_for_graph = af_data[rowSums(af_data_filt[, c(2:9)] > 0) > 7, ]
#af_data_common_filt=af_data_filt[which(af_rs_filt$AN_nfe_seu>73),]

#af_data_seu_filt=af_data_filt[which(af_data$AN_nfe_seu>1000),]
#af_data_seu_filt = af_data_seu_filt[rowSums(af_data_seu_filt[, c(2,3,6:9)] > 0) > 3, ]
pcaout = prcomp(as.matrix(t(af_data_common[, c(2:9)])))
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-PC1, y=PC2,fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() + xlab('PC1 - 77% of variance') + ylab('PC2 - 17% of variance') +scale_fill_lancet()+
  theme( axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = 'None',
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14))

fisher.test(rbind(c(33,1082),c(86,1964)),
            alternative="two.sided")$p.value

fisher.test(rbind(c(33,1082),c(774,19045)),
            alternative="two.sided")$p.value


rus_data = read.table('/home/anton/ACE2_gnomad_rus.tsv', header=T, sep='\t', stringsAsFactors = F)
rus_data_filt = rus_data[,c(2,9,12,15,18,30,33,36,39,42)] 
write.table(rus_data_filt,'russian_al_freq.tsv',sep='\t')

pcaout = prcomp(as.matrix(t(rus_data_filt[, c(2:10)])))
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-PC1, y=PC2,fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() + xlab('PC1 - 92% of variance') + ylab('PC2 - 6% of variance') + scale_fill_lancet(labels=c('AFR','EAS','FIN','NFE','NFE_EST','NFE_NWE','NFE_ONF','NFE_SEU', 'RUS'))+guides(fill=guide_legend(title="Population"))+
  theme( axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14))

rus_data_fish = rus_data[,c(2,7,8,10,11,13,14,16,17,28,29,31,32,34,35,37,38,40,41)] 
colnames(rus_data_filt)
ncol(rus_data_fish)
pop_vec=c('nfe','fin','eas','afr','nfe_nwe','nfe_est','nfe_seu','nfe_onf','rus')
pvals <- data.frame(comp=c(1:40),pval=c(1:40),var=c(1:40))

ind=1
for (i in 1:nrow(rus_data_fish)){
  g=1
  for (j in seq(from=2, to=16, by=2)){
    n1=rus_data_fish[i,18]
    c1=rus_data_fish[i,19]
    c2=rus_data_fish[i,j+1]
    n2=rus_data_fish[i,j]
    pvals[ind,1]=pop_vec[g]
    pvals[ind,2]=fisher.test(rbind(c(c1,n1),c(c2,n2)),alternative="two.sided")$p.value
    pvals[ind,3]=rus_data_fish[i,1]
    ind=ind+1
    g=g+1
  }
}

pvals$p_adj=p.adjust(pvals$pval, method = 'fdr')

p_filt <- pvals[pvals$p_adj<0.05,]

write.table(p_filt,'russian_pvals.tsv',sep='\t')

