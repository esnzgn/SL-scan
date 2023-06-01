rm(list = ls())

# Load required packages
# library(DataEditR)

# CMDB class A (litrature based) cancer metabolic gene list
# class_a = read.csv('../MetabolicSLInput/data/class_a_entrezid.csv')

#### expression matrix preparation ####
# read ccle data for 22 celllines then then match them with hgnc:ids in recon.202
# read ccle celllines expression data
ccle_GE = read.csv('../MetabolicSLInput/data/CCLE_expression.csv',header = T,row.names = 1)

for (i in 1:ncol(ccle_GE)){
  colnames(ccle_GE)[i] = unlist(strsplit(colnames(ccle_GE)[i],'[.][.]'))[1]
}

all_genes = unique(colnames(ccle_GE))
# write.csv(all_genes,'../MetabolicSLInput/data/all_genes.csv')

sample_info = read.csv('../MetabolicSLInput/data/sample_info.csv',header = T)

# read hgnc:ids list of recon.204
recon_entrz = read.csv('../MetabolicSLInput/data/Recon2v4_genes.csv',header = F)
# hgnc_ls = read.csv('./data/hgnc_complete_set.txt', header = T, sep = "")

# find each entrz_id gymbol through appling get_sym
source('./code/func__tans_dom__get_sym.R')
# met_genes = apply(recon_entrz,MARGIN = 1,FUN = get_sym)
# met_genes0 = cbind(met_genes, recon_entrz$V1)
met_genes0 = NULL
for (i in 1:length(recon_entrz$V1)){
  tmp = get_sym(recon_entrz$V1[i])
  if (!is.na(tmp)){
    met_genes0 = rbind(met_genes0, cbind(tmp , recon_entrz$V1[i]))
  }
}
met_genes0 = data.frame(met_genes0)
colnames(met_genes0) = c('gymbol', 'entrz')
 
# slice ccle expression matrix for metabolic genes
# grep('FPGT', colnames(ccle_GE))
ccle_col = colnames(ccle_GE) %in% met_genes0$gymbol
met_ccle_GE = ccle_GE[,ccle_col]

# add primary disease name beside cellline depmap ids
p_disease_ls = sample_info$primary_disease[match(rownames(met_ccle_GE),sample_info$DepMap_ID)]

# change all slash '/' in colon and endometrial cancer to under score '_' in disease list
p_disease_ls = gsub(p_disease_ls, pattern = '/', replacement = '_')

# cbind p_disease_ls with met_ccle_GE
met_ccle_GE = cbind(p_disease_ls,met_ccle_GE)
id_rpmi = grep('RPMI', sample_info$culture_medium, ignore.case = T)
rid = rownames(met_ccle_GE) %in% sample_info$DepMap_ID[id_rpmi]
met_ccle_GE = met_ccle_GE[rid,]

# number of cell lines per cancer barplot for manuscript ####
can_freq = table(met_ccle_GE$p_disease_ls)
id_can = which(can_freq >= 10)
can_cell_freq = can_freq[id_can][order(-can_freq[id_can])]
can_cell_freq = data.frame(can_cell_freq)

color.function <- colorRampPalette(c("blue", "white", "red"))
stepSize = (1/nrow(can_cell_freq))
col.seq <- round(seq(min(can_cell_freq$Freq)-stepSize
                     , max(can_cell_freq$Freq)+stepSize, stepSize), 1)
colors <- colorRampPalette(c("blue", "white", "red"))(length(col.seq)) 

x <- can_cell_freq$Freq
x.map <- round(x, 1)
x.index <- match(x.map,col.seq)
barCol = colors[x.index]
barCol = rep("blue",length(can_cell_freq$Freq))
# pdf(paste0("../MetabolicSLOutput/cancer_RPMI_celllines_freq.pdf"), width = 12, height = 8)
png(paste0("../MetabolicSLOutput/cancer_RPMI_celllines_freq.png"), width = 600, height = 400)
# create the plot with ggplot2
ggplot(can_cell_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = barCol) +
  labs(
    title = "Bar plot of frequency of each cancer cell line cultured in RPMI medium",
    y = "Frequency of cell lines"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 2*rel(1)),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 0.5, 1, 1), "cm")
  )

dev.off()

# save environment ####
# save|load env
# save.image('./env/raw_data_expression_cancer_selection.rdata')
# load('./env/raw_data_expression_cancer_selection.rdata')

# select celllines with more than 20 frequency
frq_can = which(table(met_ccle_GE$p_disease_ls) >= 10)
# length(frq_can)
expr0 = subset(met_ccle_GE, (met_ccle_GE$p_disease_ls %in% names(frq_can)))
expr0 = expr0[order(expr0$p_disease_ls),]
# transpose
expr = t(expr0)
# length(is.na(met_genes0$entrz))
expr1 = cbind(met_genes0$entrz[match(rownames(expr), met_genes0$gymbol)], expr)
expr1[,1] = paste0(expr1[,1],'.1')
expr1[1,1] = NA
rownames(expr1)[1] = NA
# colnames(expr1) = gsub(x = colnames(expr1), pattern = '-', replacement = '_')
write.csv(expr1, file = '../MetabolicSLInput/data/rnaData_gymbol_entrz_cel.csv', row.names = T)
can = unique(expr1[1,])[!is.na(unique(expr1[1,]))]
cantyp = data.frame(matrix(NA, ncol = 3, nrow = length(can)))
colnames(cantyp) = c('cancer', 'start', 'end')
  
colnames(expr) = unname(expr[1,])
expr = expr[-1,]
# expr = expr[,order(colnames(expr))]
# prepare for gene expression geymbol primary disease matrix gymbol
# grep("*[.]1",rownames(expr), value=TRUE)
# expr00 = expr[-c(grep("*[.]1",rownames(expr))),]
# grep("*[.]1",rownames(expr00), value=TRUE)
# write.csv(expr00, file = '../MetabolicSLInput/data/rnaData_gymbol.csv', row.names = T)

# preparing expr for iMAT matlab entrz 
expr000 = expr00
rownames(expr000) = met_genes0[match(rownames(expr00),met_genes0[,1]),2]
rownames(expr000) = paste0(rownames(expr000),'.1')

write.csv(expr000, file = '../MetabolicSLInput/data/rnaData_entrz.csv', row.names = T)
canty = unique(colnames(expr000))
# cell lines with only rpmi culture medium
write.csv(expr000, file = '../MetabolicSLInput/data/rnaData_entrz.csv', row.names = T)

# making a coordination file for rnaData_entrz.csv to be able to map each cancer to its related expression data
x = data.frame(matrix(data = NA, nrow = length(canty),ncol = 3))
for (i in 1:length(canty)){
  x[i,1] = canty[i]
  x[i,2] = min(which(canty[i] == colnames(expr000)))
  x[i,3] = max(which(canty[i] == colnames(expr000)))
}

write.csv(x, file = '../MetabolicSLInput/data/cantyp.csv', row.names = F)


# FastCore core set of genes for different cancers 
end_cntr = length(unique(colnames(expr000)))
cor_can = data.frame(matrix(NA, nrow = nrow(class_a), ncol = end_cntr))
cor_can_var = NULL

for (i in 1:end_cntr){
  can_tmp = unique(colnames(expr000))[i]
  dis_expr = expr000[,which(can_tmp == colnames(expr000))]
  # print(paste0(colnames(dis_expr)[1],' _ ncol: ',ncol(dis_expr),' _ start: ',strd,' _ end: ',endd,'\n'))
  cs = core_set(class_a$ENTREZID,dis_expr,0.6,0.4,'.1')
  cor_can[1:length(cs),i] = cs
  colnm = paste(unlist(strsplit(can_tmp,split = ' ')),collapse = '_')
  # colnm = paste(unlist(strsplit(can_tmp,split = '/')),collapse = '_')
  colnames(cor_can)[i] = colnm
}

# View(sort(table(unlist(cor_can))))

cor_can_var = cor_can + 0.1
# writing 22 cancers core set genes to find their corresponding rxns in matlab for fastcore algo
write.csv(cor_can_var, file = '../MetabolicSLInput/data/cancer_core_genes_entrz.csv', row.names = T)



