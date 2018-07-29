library(data.table)
library(ggplot2)

devtools::document('~/suvrelR')
devtools::load_all('~/suvrelR')

sample_info <- 
  fread('/home/jonatas/suvrelR/extdata/geuvadis_sample_info.tsv')

gene_wide <- 
  fread('/home/jonatas/suvrelR/extdata/geuvadis_quantification.txt') 

gene_long <- 
  melt(gene_wide, measure.vars = 5:dim(gene_wide)[2], 
       variable.name = 'subject', 
       value.name = 'fpkm')

datum <- 
  gene_long[sample_info[,.(name, pop)], 
            on = c(subject = 'name')
            ][pop != 'CEU' && !is.na(Chr)]

datum[, Chr := as.integer(Chr)]
datum <- datum[!is.na(Chr)]

setnames(datum, 'TargetID', 'gene_id')

erg <- suvrel_diag(datum, 
                   value_var = 'fpkm', 
                   group = 'pop', 
                   feature = 'gene_id')
erg[,summary(g)]

new_datum <- datum[erg, on = .(gene_id)]
new_datum[, fpkm_g := fpkm * sqrt(g)]

X_old <- dcast(new_datum, subject + pop ~ gene_id, value.var = 'fpkm')
X_new <- dcast(new_datum, subject + pop ~ gene_id, value.var = 'fpkm_g')
n <- dim(X_old)[2]

pc_old <- prcomp(X_old[,3:n], center = TRUE)
po <- as.data.table(pc_old$rotation)[, 1:4]

pc_new <- prcomp(X_new[,3:n], center = TRUE)
pn <- as.data.table(pc_new$rotation)[, 1:4]
X_new[,1:5]

ggplot(pn) + 
  geom_point(aes(x = PC1, y = PC2))
