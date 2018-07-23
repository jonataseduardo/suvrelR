library(data.table)
library(ggplot2)

devtools::document('~/suvrelR')
devtools::load_all('~/suvrelR')

#fname <- paste0(system.file("extdata", package = "diverdt"), 
#                '/', pop_name, '_Illumina')
sample_info <- 
  fread('/raid/genevol/users/vitor/hlaexpression/geuvadis_reanalysis/data/sample_info/geuvadis_sample_info.tsv')

gene_long <- 
  fread('/raid/genevol/users/vitor/hlaexpression/geuvadis_reanalysis/expression/3-map_to_transcriptome/hla_personalized/gene_quantifications.tsv') 

sample_info
datum <- 
  gene_long[sample_info[,.(name, pop)], 
            on = c(subject = 'name')]

devtools::load_all('~/suvrelR')
erg <- suvrel_diag(datum, value_var = 'tpm', group = 'pop', feature = 'gene_id')

erg[,summary(g)]
new_datum <- datum[erg, on = .(gene_id)]
new_datum

