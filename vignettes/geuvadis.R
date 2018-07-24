library(data.table)
library(ggplot2)

devtools::document('~/suvrelR')
devtools::load_all('~/suvrelR')

sample_info <- 
  fread('/home/jonatas/suvrelR/extdata/geuvadis_sample_info.tsv')

gene_long <- 
  fread('/home/jonatas/suvrelR/extdata/gene_quantifications.tsv') 

sample_info

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

X_old <- dcast(new_datum, gene_id ~ subject, value.var = 'tpm')
X_new <- dcast(new_datum, gene_id ~ subject, value.var = 'tpm_g')

