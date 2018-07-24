library(data.table)
library(ggplot2)

devtools::document('~/suvrelR')
devtools::load_all('~/suvrelR')

sample_info <- 
  fread('/home/jonatas/suvrelR/extdata/geuvadis_sample_info.tsv')

gene_long <- 
  fread('/home/jonatas/suvrelR/extdata/gene_quantifications.tsv') 

datum <- 
  gene_long[sample_info[,.(name, pop)], 
            on = c(subject = 'name')]

erg <- suvrel_diag(datum, 
                   value_var = 'tpm', 
                   group = 'pop', 
                   feature = 'gene_id')

erg[,summary(g)]
new_datum <- datum[erg, on = .(gene_id)]
new_datum[, tpm_g := tpm * sqrt(g)]

X_old <- dcast(new_datum, subject + pop ~ gene_id, value.var = 'tpm')
X_new <- dcast(new_datum, subject + pop ~ gene_id, value.var = 'tpm_g')
n <- dim(X_old)[2]

pc_old <- prcomp(X_old[,3:n], center = TRUE)
po <- as.data.table(pc_old$rotation)[, 1:4]

pc_new <- prcomp(X_new[,3:n], center = TRUE)
pn <- as.data.table(pc_new$rotation)[, 1:4]

ggplot(pn) + 
  geom_point(aes(x = PC1, y = PC2))


