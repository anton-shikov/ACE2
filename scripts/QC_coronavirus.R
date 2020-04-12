
library(ggplot2)
mean_cov <- read.csv('/home/anton/samples_metrics_filtered.tsv', header = F, sep = "\t")

DP <- read.csv('/home/anton/DP.tsv', header = F, sep = "\t")

VAF <- read.csv('/home/anton/VAF.tsv', header = F, sep = "\t")


dodge <- position_dodge(width = 0.4)
median(mean_cov$V1)
mean(mean_cov$V1)
mean(mean_cov$V1[mean_cov$V2=='wes'])

ggplot(mean_cov , aes(y=V1,x=V2, fill=V2)) + 
  theme_bw()+
  geom_violin(scale = "width", trim = F, adjust=0.7, alpha=0.75)+
  geom_hline(yintercept=44, linetype="dashed", color = "green")+
  scale_fill_lancet(labels = c("CES", "WES"))+
  guides(fill=guide_legend(title="Platform"))+
  ylab('Mean coverage')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = c(0.89, 0.85),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14), 
         legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12))+ylim(c(0,200))

mean(DP$V1)

ggplot(DP[DP$V1>7,], aes(y=V1,x=V2, fill=V2)) + 
  theme_bw()+
  geom_violin(scale = "width", trim = F, adjust=0.5, alpha=0.9)+
  geom_hline(yintercept=22, linetype="dashed", color = "green")+
  scale_fill_lancet()+
  guides(fill=guide_legend(title="Variant"))+
  ylab('DP')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = c(0.79, 0.77),
         axis.text.y = element_text(color='black', 
                                                                   size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14), 
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))+ylim(c(0,250))

median(VAF$V1)
ggplot(VAF[VAF$V1>0,], aes(y=V1,x=V2, fill=V2)) + 
  theme_bw()+
  geom_violin(scale = "width", trim = F, adjust=0.7, alpha=0.9)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "green")+
  scale_fill_lancet()+
  ylab('VAF')+guides(fill=guide_legend(title="Variant"))+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.y = element_text(color='black', 
                                                             size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14))+ylim(c(0,1))

eqtls <- read.csv('/home/anton/ACE2_eqtls_gnomad.tsv', header = T, sep = "\t")

eqtl_sums_dat=data.frame(AC_sum=colSums(eqtls[,grepl('AN', colnames(eqtls))]),
                       AC_mean=colMeans(eqtls[,grepl('AN', colnames(eqtls))]),
                       AN_count=colSums(eqtls[,grepl('AN', colnames(eqtls))]>0),
                       AC_count=colSums(eqtls[,grepl('AC', colnames(eqtls))]>0),
                       AF_mean=colMeans(eqtls[,grepl('AF', colnames(eqtls))]))



eqtl_sums_dat= eqtl_sums_dat[eqtl_sums_dat$AN_count>0,]

ggplot(eqtl_sums_dat, aes(x=AF_mean, y=rownames(eqtl_sums_dat),fill=rownames(eqtl_sums_dat))) + 
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
         legend.text=element_text(size=12))+
  geom_vline(xintercept=0.6498567, linetype="dashed", color = "green")



max(eqtls$CROM)-min(eqtls$CROM)
eqtls$pval=-log10(eqtls$pval)

ggplot(eqtls, aes(POS, pval,fill = tissue)) +
  geom_point(alpha = 0.6,shape=21,size=4) +
  theme_bw()+scale_fill_lancet()+
  labs(y = "-log10P") + 
  labs(x = "Genomic coordinate")+
  geom_vline(xintercept=15580031, linetype="dashed", color = "green")+
  geom_vline(xintercept=15580136, linetype="dashed", color = "green")+
  geom_vline(xintercept=15582147, linetype="dashed", color = "green")+
  geom_vline(xintercept=15582341, linetype="dashed", color = "green")+
  geom_vline(xintercept=15584376, linetype="dashed", color = "green")+
  geom_vline(xintercept=15584492, linetype="dashed", color = "green")+
  geom_vline(xintercept=15585849, linetype="dashed", color = "green")+
  geom_vline(xintercept=15585949, linetype="dashed", color = "green")+
  geom_vline(xintercept=15588418, linetype="dashed", color = "green")+
  geom_vline(xintercept=15588476, linetype="dashed", color = "green")+
  geom_vline(xintercept=15589747, linetype="dashed", color = "green")+
  geom_vline(xintercept=15589919, linetype="dashed", color = "green")+
  geom_vline(xintercept=15590324, linetype="dashed", color = "green")+
  geom_vline(xintercept=15590446, linetype="dashed", color = "green")+
  geom_vline(xintercept=15591490, linetype="dashed", color = "green")+
  geom_vline(xintercept=15591588, linetype="dashed", color = "green")+
  geom_vline(xintercept=15593789, linetype="dashed", color = "green")+
  geom_vline(xintercept=15593933, linetype="dashed", color = "green")+
  geom_vline(xintercept=15596212, linetype="dashed", color = "green")+
  geom_vline(xintercept=15596438, linetype="dashed", color = "green")+
  geom_vline(xintercept=15599344, linetype="dashed", color = "green")+
  geom_vline(xintercept=15599513, linetype="dashed", color = "green")+
  geom_vline(xintercept=15603598, linetype="dashed", color = "green")+
  geom_vline(xintercept=15603695, linetype="dashed", color = "green")+
  geom_vline(xintercept=15605876, linetype="dashed", color = "green")+
  geom_vline(xintercept=15605981, linetype="dashed", color = "green")+
  geom_vline(xintercept=15607467, linetype="dashed", color = "green")+
  geom_vline(xintercept=15607579, linetype="dashed", color = "green")+
  geom_vline(xintercept=15609836, linetype="dashed", color = "green")+
  geom_vline(xintercept=15609979, linetype="dashed", color = "green")+
  geom_vline(xintercept=15610352, linetype="dashed", color = "green")+
  geom_vline(xintercept=15610445, linetype="dashed", color = "green")+
  geom_vline(xintercept=15612968, linetype="dashed", color = "green")+
  geom_vline(xintercept=15613126, linetype="dashed", color = "green")+
  geom_vline(xintercept=15618849, linetype="dashed", color = "green")+
  geom_vline(xintercept=15619034, linetype="dashed", color = "green")+
  geom_vline(xintercept=15579156, linetype="dashed", color = "red")+
  geom_vline(xintercept=15620271, linetype="dashed", color = "red")+
  theme( axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x=element_text(face="bold", color="black", 
                                  size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12))


intervals$z=intervals$V1-min(intervals$V1)+1

intervals$bin=1

s=1
for(i in seq(from=501, to=52217, by=500))
  {
  for (j in 1:length(intervals$z)){
  if (intervals$z[j] <= i && intervals$z[j] >= i-500){
    intervals$bin[j]=i
  }
}
}

intervals$bin=factor(intervals$bin)
15599940-15599613 #327
15597330-15596022 #1308
15583151-15582092 #1059
15590829-15590263 #566

327+1308+1059+566




