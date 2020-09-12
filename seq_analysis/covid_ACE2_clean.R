setwd("./")

library(ggplot2)
library(reshape2)
library(qqman)
library(cowplot)

mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

covid_data = read.table('genotypes.txt', sep='\t', header=T, stringsAsFactors = F)
covid_data$Controls = as.numeric(covid_data$Controls)
covid_data$Light = as.numeric(covid_data$Light)
covid_data$Severe = as.numeric(covid_data$Severe)
head(covid_data)

hist(covid_data$chisq_pval)
qq(covid_data$chisq_pval)

candidates = covid_data[covid_data$chisq_pval < 0.05, ]
candidates_2 = covid_data[covid_data$Light == 0 & covid_data$Severe > 1, ]
candidates_3 = covid_data[covid_data$Light == 0 & covid_data$Severe > covid_data$Light * 2, ]

table(covid_data$Controls + covid_data$Light + covid_data$Severe == 1)

covid_data[covid_data$Controls + covid_data$Light + covid_data$Severe == 1, 1:6]

counts = matrix(c(3, 21, 1, 19, 9, 18), nrow=2)

chisq.test(counts)

# Plot ACs
ace2_count = read.table('./genotypes.txt', sep='\t', header=T)
colnames(ace2_count)
plot_snp = melt(ace2_count, id.vars='SNP_ID', 
                measure.vars=c('Controls', 'Light', 'Severe'))

ggplot(plot_snp, aes(x=reorder(SNP_ID, value), y=value, fill=variable)) + 
  geom_bar(stat='identity', col='black', position='dodge') + theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1)) + coord_flip() +
  facet_wrap(~variable) + xlab('') + ylab('Alternative allele count') + 
  guides(fill=F) + scale_fill_manual(values=c(mycol1, mycol2, mycol3))

# Explore rare sites in gnomAD

gnomad_rare = read.table('ACE2_gnomad_rare_vars.tsv', header=T, sep="\t")

sum(ifelse(gnomad_rare$AF_afr <= 0.01, gnomad_rare$AF_afr, 0))
sum(ifelse(gnomad_rare$AF_fin <= 0.01, gnomad_rare$AF_fin, 0))
sum(ifelse(gnomad_rare$AF_nfe_seu <= 0.01, gnomad_rare$AF_nfe_seu, 0))
sum(ifelse(gnomad_rare$AF_nfe_est <= 0.01, gnomad_rare$AF_nfe_est, 0))
sum(ifelse(gnomad_rare$AF_fin <= 0.01, gnomad_rare$AF_fin, 0))


######### eQTL models ##############################################

eqtl_data = read.table('eqtl.tsv', sep='\t', header=T, dec=',')
af_cols = grepl('AF_', colnames(eqtl_data))
selector = af_cols
selector[2] = TRUE
eqtl_to_reg = melt(eqtl_data[, selector], id.vars='rsID')
eqtl_to_reg$value = as.numeric(eqtl_to_reg$value)
wqtl_to_reg = eqtl_to_reg[eqtl_to_reg$value > 0, ]
my_fit = lm(value~variable+rsID, eqtl_to_reg)
summary(my_fit)

# Accounting for sample size

ac_cols = grepl('AC_', colnames(eqtl_data))
selector = ac_cols
selector[2] = TRUE
eqtl_to_reg_AC = melt(eqtl_data[, selector], id.vars='rsID')
eqtl_to_reg_AC$value = as.numeric(eqtl_to_reg_AC$value)

an_cols = grepl('AN_', colnames(eqtl_data))
selector = an_cols
selector[2] = TRUE
eqtl_to_reg_AN = melt(eqtl_data[, selector], id.vars='rsID')
eqtl_to_reg_AN$value = as.numeric(eqtl_to_reg_AN$value)

eqtl_to_reg_AC$AN = eqtl_to_reg_AN$value
eqtl_to_reg_AC = eqtl_to_reg_AC[eqtl_to_reg_AC$value > 0, ]

my_fit = lm(value~variable*AN + rsID, eqtl_to_reg_AC)
summary(my_fit)


###################### Biochemical features ##################################

as_data = read.table('biochem.tsv', sep='\t', header=T)
snippet <- as_data[20:37, ]

wilcox.test(as_data$СРБ..количественный.~as_data$carrier_rare)

# CRP
wilcox.test(snippet$СРБ..количественный.~snippet$carrier_rare)
aggregate(СРБ..количественный.~carrier_rare, snippet, mean)

as_data$group = factor(ifelse(as_data$form == 'тяжелая',
                       ifelse(as_data$carrier_rare, 
                              'severe\ncarrier', 
                              'severe\nnon-carrier'), 'mild'),
                       levels=c('mild', 'severe\nnon-carrier', 'severe\ncarrier'))

ggplot(as_data, aes(x=group, y=СРБ..количественный., fill=group)) +
  geom_boxplot() + geom_point() + theme_bw() +
  ylab('C-reactive protein, mg/l') + xlab('COVID-19 patient group') +
  scale_fill_manual(values=c(mycol2, mycol3, mycol6)) + guides(fill=F)

#https://www.medrxiv.org/content/10.1101/2020.04.03.20047977v3

# LDG
wilcox.test(snippet$ЛДГ~snippet$carrier_rare)
aggregate(ЛДГ~carrier_rare, snippet, mean)

# Fibrinogen
wilcox.test(snippet$Фибриноген~snippet$carrier_rare)
aggregate(Фибриноген~carrier_rare, snippet, mean)

# D-dimer
wilcox.test(snippet$Д.димер~snippet$carrier_rare)
aggregate(Д.димер~carrier_rare, snippet, mean)

# Ferritin
wilcox.test(snippet$Ферритин~snippet$carrier_rare)
aggregate(Ферритин~carrier_rare, snippet, mean)

# Lymphocytes
wilcox.test(snippet$Лимфоциты..LYMPH.~snippet$carrier_rare)
wilcox.test(snippet$Лимфоциты....LYMPH..~snippet$carrier_rare)
  
