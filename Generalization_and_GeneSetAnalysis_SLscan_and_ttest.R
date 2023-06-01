rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(plyr)
library(dplyr)
library(stringr)
library("readxl")
library(igraph) # Load the igraph package
library("rTRM")
library(png)
library(stats)
library(msigdbr)

#### reading fdr corrected #### 
fname = list.files(pattern ='_ac_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
SL_achiles = NULL
SL_scan  = NULL
for (cn in 1:length(fname)){
  print(cn/length(fname))
  cantype = unlist(strsplit(fname[cn], '_ac_mn_w&t_tst_cons_fdr.csv')[1])
  tmp = read.csv(paste0('../MetabolicSLOutput/',fname[cn]))
  mn = which(as.numeric(tmp$ttst_metabolic_fdr) <= 0.05)
  ach = which(as.numeric(tmp$ttst_achiles_fdr) <= 0.05)
  SL_scan = rbind(SL_scan , tmp[mn,])
  SL_achiles = rbind(SL_achiles , tmp[ach,])
}; rm(tmp)
SL_achiles$ttst_achiles_stat = formatC(SL_achiles$ttst_achiles_stat, format = "e", digits = 4)
SL_scan$ttst_metabolic_stat = formatC(SL_scan$ttst_metabolic_stat, format = "e", digits = 4)
# View(SL_achiles[order(SL_achiles$ttst_achiles_fdr),])
write.csv(SL_achiles,'../MetabolicSLOutput/SL_Achiles_sig_fdr.csv')
write.csv(SL_scan,'../MetabolicSLOutput/SL_scan_sig_fdr.csv')

#### frequency of each mut gene and sl pairs in both methods ####
# find frequency of each single mutated and SL pairs of genes in different cancers for SLscan and Achiles ttest
SL_scan = SL_scan[which(SL_scan$ttst_metabolic_fdr < 0.05),]
freq_SLscan = data.frame(table(SL_scan$Mutated_gene))
top_fifteen_driver = freq_SLscan[order(freq_SLscan$Freq, decreasing = T),][1:15,]
library(ggplot2)
top_fifteen_driver$Var1 <- reorder(top_fifteen_driver$Var1, top_fifteen_driver$Freq)
png(paste0("../MetabolicSLOutput/freq_driver_genes_slscan.png"), width = 800, height = 1000)
ggplot(data = top_fifteen_driver, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", fill = "blue3") +
  coord_flip()  +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 25),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 45)
        ,plot.margin=unit(c(1,1,1,1),"cm"))
dev.off()


freq_SLscan = table(SL_scan$Mutated_gene, SL_scan$KO_gene)
freq_SLscan = as.data.frame(freq_SLscan)
freq_SLscan = subset(freq_SLscan, Freq > 0)
freq_SLscan = freq_SLscan[order(-freq_SLscan$Freq),]
colnames(freq_SLscan) = c("Mutated_gene","KO_gene", "count")
       
write.csv(SL_scan, "../MetabolicSLOutput/SLscan_sig.csv")

freq_SLscan =  cbind(freq_SLscan, Cancers_SLscan = NA, Cancers_CRISPR = NA)
for (i in 1:nrow(freq_SLscan)){
  freq_SLscan$Cancers_SLscan[i] = paste(SL_scan$cancer[intersect(which(freq_SLscan$Mutated_gene[i] == SL_scan$Mutated_gene) 
                   ,which(freq_SLscan$KO_gene[i] == SL_scan$KO_gene))], collapse = "; ")
  freq_SLscan$Cancers_CRISPR[i] = paste(SL_achiles$cancer[intersect(which(freq_SLscan$Mutated_gene[i] == SL_achiles$Mutated_gene) 
                                                          ,which(freq_SLscan$KO_gene[i] == SL_achiles$KO_gene))], collapse = "; ")
}


write.csv(freq_SLscan, "../MetabolicSLOutput/enrichment_results/SL_scan/freq_SLscan.csv")

freq_Single_genes_slscan =   table(SL_scan$Mutated_gene) %>%  data.frame()
freq_Single_genes_slscan = freq_Single_genes_slscan[order(-freq_Single_genes_slscan$Freq),]
colnames(freq_Single_genes_slscan) = c("Mutated_gene", "count")
  
# freq_SL_achiles
freq_SL_achiles = table(SL_achiles$Mutated_gene, SL_achiles$KO_gene) %>% data.frame()
freq_SL_achiles = subset(freq_SL_achiles, freq_SL_achiles$Freq > 0)
freq_SL_achiles = freq_SL_achiles[order(-freq_SL_achiles$Freq),]
colnames(freq_SL_achiles) = c("Mutated_gene","KO_gene", "count")

freq_SL_achiles =  cbind(freq_SL_achiles, cancers = NA)
for (i in 1:nrow(freq_SL_achiles)){
  freq_SL_achiles$cancers[i] = paste(SL_achiles$cancer[intersect(which(freq_SL_achiles$Mutated_gene[i] == SL_achiles$Mutated_gene) 
                                                          ,which(freq_SL_achiles$KO_gene[i] == SL_achiles$KO_gene))], collapse = "; ")
}
write.csv(freq_SL_achiles, "../MetabolicSLOutput/freq_cancer_fdrsig_achiles.csv")

# frequency of driver genes in achilles
freq_Single_genes_achiles = table(SL_achiles$Mutated_gene) %>% data.frame()
freq_Single_genes_achiles = freq_Single_genes_achiles[order(-freq_Single_genes_achiles$Freq),]
colnames(freq_Single_genes_achiles) = c("Mutated_gene", "count")

# library("dplyr")
# freq_achilles_driver = SL_achiles %>% 
#   group_by( SL_achiles$Mutated_gene ) %>% 
#   summarise(count = n()) %>% 
#   arrange( desc(count) )


# put all KO genes of each single frequent mutated genes among different SL-pairs and cancers 
# and send them for GO enrichment and another time find each one frequency and set their frequency as
# their score for GSEA

#### SLscan GSEA ####
aggregator_ko_genes_slscan = function(mut_gene, df){
  tmp = df$ttst_metabolic_stat[which(mut_gene == df$Mutated_gene)]
  if (length(grep(mut_gene, df$KO_gene)) != 0){
    ko_esti = mean(as.numeric(df$ttst_metabolic_stat[which(mut_gene == df$KO_gene)]))
    tmp1 = c(ko_esti,as.numeric(tmp))
    names(tmp1) = c(mut_gene, df$KO_gene[which(mut_gene == df$Mutated_gene)])
  } else {
    tmp1 = as.numeric(tmp);
    names(tmp1) = c(df$KO_gene[which(mut_gene == df$Mutated_gene)])
    }
  abs(tmp1)[order(-abs(tmp1))]
}

gene_set_SLscan = NULL
for (freq_i in 1:length(which(freq_Single_genes_slscan$count >= 2))){
  mut_gene = as.character(freq_Single_genes_slscan$Mutated_gene[freq_i])
  tmp = aggregator_ko_genes_slscan(mut_gene, SL_scan)
  gene_set_SLscan[[mut_gene]] = tmp
  rm(mut_gene, tmp)
}

dir.create("../MetabolicSLOutput/enrichment_results/SL_scan/", showWarnings = F, recursive = T)

for (sl_i in 1:length(gene_set_SLscan)){
  print(sl_i/length(gene_set_SLscan)*100)
  gene_list = gene_set_SLscan[[sl_i]]
  
  tryCatch({
    tryCatch({
      res_gene_set_SLscan = gseGO(gene_list
                                  , OrgDb = 'org.Hs.eg.db'
                                  ,ont = "BP" ,keyType = "SYMBOL"
                                  # ,eps = 1e-100
                                  ,minGSSize    = 10
                                  ,maxGSSize    = 1000
                                  ,scoreType = "pos"
                                  ,pvalueCutoff = 0.05, verbose = FALSE
                                  , use_internal_data = TRUE
                                  ,gsDatabase = "msigdb", geneSetType = "c2.cp.kegg"
                                  )
    }, error=function(e) {
      message(paste("Error occurred in iteration", sl_i, "with message:", e$message))
    })
    
    if (exists("res_gene_set_SLscan") && nrow(res_gene_set_SLscan@result) > 0) {
      mutgene = names(gene_set_SLscan[sl_i])
      write.csv(res_gene_set_SLscan@result, paste0("../MetabolicSLOutput/enrichment_results/SL_scan/GSEA_table_mutG_",mutgene,"_GSEA.csv"))
      pdf(paste0("../MetabolicSLOutput/enrichment_results/SL_scan/GSEA_RunnerPlot_mutG_",mutgene,"_GSEA.pdf"), width = 10, height = 6)
      for (i in 1:nrow(res_gene_set_SLscan@result)) {
        print(gseaplot(res_gene_set_SLscan, geneSetID = i, title = res_gene_set_SLscan$Description[i]))
      }
      dev.off()
    }
  },
  error=function(e) {
    message(paste("Error occurred in iteration", sl_i, "with message:", e$message))
  }
  )
  rm(res_gene_set_SLscan)
}

#### Achiles GSEA ####
aggregator_ko_genes_achilles = function(mut_gene, df){
  tmp = df$ttst_achiles_stat[which(mut_gene == df$Mutated_gene)]
  ko_esti = mean(df$ttst_achiles_stat[which(mut_gene == df$KO_gene)])
  tmp1 = c(ko_esti,tmp)
  names(tmp1) = c(mut_gene, df$KO_gene[which(mut_gene == df$Mutated_gene)])
  abs(tmp1)[order(-abs(tmp1))]
}

gene_set_achilles = NULL
for (freq_i in 1:length(which(freq_Single_genes_achiles$count >= 10))){
  mut_gene = as.character(freq_Single_genes_achiles$Mutated_gene[freq_i])
  tmp = aggregator_ko_genes_achilles(mut_gene, SL_achiles)
  gene_set_achilles[[mut_gene]] = tmp
  rm(mut_gene, tmp)
}

dir.create("../MetabolicSLOutput/enrichment_results/SL_achilles/", showWarnings = F, recursive = T)

for (sl_i in 1:length(gene_set_achilles)){
  print(sl_i/length(gene_set_achilles)*100)
  gene_list = gene_set_achilles[[sl_i]]
  
  tryCatch({
    tryCatch({
      res_gene_set_achilles = gseGO(gene_list
                                  , OrgDb = 'org.Hs.eg.db'
                                  ,ont = "BP" ,keyType = "SYMBOL"
                                  ,eps = 1e-100
                                  ,minGSSize    = 10
                                  ,maxGSSize    = 1000
                                  ,scoreType = "pos"
      )
    }, error=function(e) {
      message(paste("Error occurred in iteration", sl_i, "with message:", e$message))
    })
    
    if (exists("res_gene_set_achilles") && nrow(res_gene_set_achilles@result) > 0) {
      mutgene = names(gene_set_achilles[sl_i])
      write.csv(res_gene_set_achilles@result, paste0("../MetabolicSLOutput/enrichment_results/SL_achilles/GSEA_table_mutG_",mutgene,"_GSEA.csv"))
      pdf(paste0("../MetabolicSLOutput/enrichment_results/SL_achilles/GSEA_RunnerPlot_mutG_",mutgene,"_GSEA.pdf"), width = 10, height = 6)
      for (i in 1:nrow(res_gene_set_achilles@result)) {
        print(gseaplot(res_gene_set_achilles, geneSetID = i, title = res_gene_set_achilles$Description[i]))
      }
      dev.off()
    }
  },
  error=function(e) {
    message(paste("Error occurred in iteration", sl_i, "with message:", e$message))
  }
  )
  rm(res_gene_set_achilles)
}

#### BP freq SLscan ####
fname = list.files(pattern ='GSEA_table_mutG_', '../MetabolicSLOutput/enrichment_results/SL_scan/')

BP_slscan = NULL
for (mutg_i in 1:length(fname)){
  print(mutg_i/length(fname)*100)
  BP_tmp = read.csv(file = paste0('../MetabolicSLOutput/enrichment_results/SL_scan/', fname[mutg_i]))
  mut_gene = str_extract(fname[mutg_i], "(?<=GSEA_table_mutG_).*?(?=_GSEA.csv)")
  BP_slscan = rbind(BP_slscan, cbind(mut_gene, BP_tmp)); rm(BP_tmp, mut_gene)
}
freq_BP_slscan <- BP_slscan %>%
  group_by(Description) %>%
  summarize(count = n()) %>%
  arrange(desc(count))



write.csv(freq_BP_slscan
          ,"../MetabolicSLOutput/enrichment_results/SL_scan/SLscan_freq_BP.csv")

freq_BP_slscan = read.csv("../MetabolicSLOutput/enrichment_results/SL_scan/SLscan_freq_BP.csv")

# fname = list.files(pattern ='GSEA_table_mutG_', '../MetabolicSLOutput/enrichment_results/SL_achilles/')
# 
# BP_achiles = NULL
# for (mutg_i in 1:length(fname)){
#   print(mutg_i/length(fname)*100)
#   BP_tmp = read.csv(file = paste0('../MetabolicSLOutput/enrichment_results/SL_achilles/', fname[mutg_i]))
#   mut_gene = str_extract(fname[mutg_i], "(?<=GSEA_table_mutG_).*?(?=_GSEA.csv)")
#   BP_achiles = rbind(BP_achiles, cbind(mut_gene, BP_tmp)); rm(BP_tmp, mut_gene)
# }
# 
# #### intersect between slscan and achiles enriched BP and cortest of Pearson ####
# freq_BP_achiles = BP_achiles %>%
#   group_by(Description) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count))
# write.csv(freq_BP_achiles
#           ,"../MetabolicSLOutput/enrichment_results/SL_achilles/Achilles_freq_BP.csv")
# 
# BP_inter = intersect(freq_BP_achiles$Description, freq_BP_slscan$Description)
# x = freq_BP_achiles$count[match(BP_inter, freq_BP_achiles$Description)]
# y = freq_BP_slscan$count[match(BP_inter, freq_BP_slscan$Description)]
# corr = cor.test(x, y, alternative = "two.sided", method = "pearson")

# BP_slscan_achilles = data.frame(cbind(Biological_process = freq_BP_slscan$Description[match(BP_inter, freq_BP_slscan$Description)]
# ,freq_BP_achiles = as.numeric(freq_BP_achiles$count[match(BP_inter, freq_BP_achiles$Description)])
# ,freq_BP_slscan = as.numeric(freq_BP_slscan$count[match(BP_inter, freq_BP_slscan$Description)])))
# 
# write.csv(BP_slscan_achilles, "../MetabolicSLOutput/BP_slscan_achilles.csv")

# # waterfall plot of x and y ####
# zscor_norm = function(x){ 
#   sdx = sd(x); meanx = mean(x)
#   (x - meanx)/sdx
#   }
# x = zscor_norm(x); y = zscor_norm(y)
# df_waterfall = data.frame(Achilles = x, SL_scan = y)
# corr <- cor.test(df_waterfall$Achilles, df_waterfall$SL_scan)
library(ggplot2)
library(reshape2)
freq_BP_slscan = freq_BP_slscan[1:15,]
df_stacked = data.frame(x_axis_var = freq_BP_slscan$Description  
                        ,SL_scan = freq_BP_slscan$count)

# create stacked bar plot
library(ggplot2)
png(paste0("../MetabolicSLOutput/BP_Slscan_SL-scan_enriched_top_20.png"), width = 2000, height = 1000)
# par(mar = c(30, 4, 4, 10))

df_stacked$x_axis_var <- reorder(df_stacked$x_axis_var, df_stacked$SL_scan, decreasing = T)

ggplot(df_stacked, aes(x = x_axis_var, y = SL_scan)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "SL-scan enriched top 20 frequent biological process", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size = 45)
       ,plot.margin=unit(c(1,1,1,20),"cm"))
dev.off()
# substr(corr$p.value,1,6)
# color.function <- colorRampPalette(c("blue", "white", "red"))
# stepSize = (1/nrow(df))
# col.seq <- round(seq(min(df_waterfall$Achilles)-stepSize
#                      , max(df_waterfall$Achilles)+stepSize, stepSize), 1)
# colors <- colorRampPalette(c("blue", "white", "red"))(length(col.seq)) 
# colors = "blue"
# x <- df_waterfall$Achilles
# x.map <- round(x, 1)
# x.index <- match(x.map,col.seq)
# barCol = colors[x.index]
# 
# pdf(paste0("../MetabolicSLOutput/BP_achillesVsSl-scan_waterfall_corrTest.pdf"), width = 10, height = 6)
# barplot(df_waterfall$SL_scan,col = barCol,names.arg=""
# , main=paste0("WFP of Achilles and SL-scan Biological process frequency with Pearson correlation 
#               coefficient of "
#             , corr$estimate," and P-value of " ,corr$p.value)
# ,ylab=paste0(colnames(df_waterfall)[2]," BP frequency ")
# ,xlab=paste0(colnames(df_waterfall)[1]," BP frequency as bar colors "))
# 
# dev.off()

#### SLscan four mode network of BP, KO_genes, driver_genes, and cancers ####
#####  prepare BP_slscan for four mode network #### 
BP_gene_divider = function(df){-
  df_final = NULL
  for (dfRow_i in 1:nrow(df)){
    tmp = unlist(strsplit(df$core_enrichment[dfRow_i], split = "/"))
    df_final = rbind(df_final, cbind(df$mut_gene[dfRow_i], df$ID[dfRow_i]
                             , df$Description[dfRow_i], df$NES[dfRow_i], df$p.adjust[dfRow_i], tmp))
    rm(dfRow_i)
  }
  colnames(df_final) = c("mut_gene","ID","Description", "NES","p.adjust","KO_gene")
  data.frame(df_final)
}

BP_slscan_divided = BP_gene_divider(BP_slscan)
# suffix for cytoscape ####
# concatenating some suffixes for different nodes for cytoscape 4 mode network 
# BP_slscan_divided$mut_gene = paste0(BP_slscan_divided$mut_gene, "__DG")
# BP_slscan_divided$KO_gene = paste0(BP_slscan_divided$KO_gene, "__KG")
# BP_slscan_divided$Description = paste0(BP_slscan_divided$Description, "__BP")
# write.csv(BP_slscan_divided
#           ,"../MetabolicSLOutput/enrichment_results/SL_scan/SLscan_KoGenes_mutGene_vs_BP_.csv")

# freq description and kogenes ####
freq_BP_slscan_divided = BP_slscan_divided %>% 
  group_by(Description, KO_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

# select top 95% quantile count of ko genes in bps 
freq_BP_slscan_divided = subset(freq_BP_slscan_divided, count > 3)
BP_slscan_divided_samll = (BP_slscan_divided[((BP_slscan_divided$KO_gene %in% freq_BP_slscan_divided$KO_gene) 
  & (BP_slscan_divided$Description %in% freq_BP_slscan_divided$Description)),])

# save|load env
# save.image('./env/WS_MN_slscan_achilesTtest_GSEA_generalization.rdata')
# load('./env/WS_MN_slscan_achilesTtest_GSEA_generalization.rdata')

write.csv(BP_slscan_divided_samll
          ,"../MetabolicSLOutput/enrichment_results/SL_scan/SLscan_KoGenes_mutGene_vs_BP_small.csv")

freq_BP_slscan_divided_mutgene$count = BP_slscan_divided %>% 
  group_by(Description, mut_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))
inter_bp = intersect(freq_BP_slscan_divided_mutgene$Description, freq_BP_slscan_divided$Description)

tmp_freq_ko_bp = freq_BP_slscan_divided[freq_BP_slscan_divided$Description %in% inter_bp,]
tmp_freq_ko_bp$Description = paste(tmp_freq_ko_bp$Description, "_BP")
# colnames(tmp_freq_ko_bp) = c("fromNode", "toNode", "weight")
write.table(tmp_freq_ko_bp
            ,"../MetabolicSLOutput/enrichment_results/SL_scan/2_SLscan_KO_genes_vs_BPcommon_cytoscape.txt"
            , sep = "\t", row.names = FALSE)

tmp_freq_mut_bp = freq_BP_slscan_divided_mutgene[freq_BP_slscan_divided_mutgene$Description %in% inter_bp,]
tmp_freq_mut_bp$Description = paste(tmp_freq_mut_bp$Description, "_BP")
# colnames(tmp_freq_mut_bp) = c("fromNode", "toNode", "weight")
write.table(tmp_freq_mut_bp
          ,"../MetabolicSLOutput/enrichment_results/SL_scan/3_SLscan_mut_genes_vs_BPcommon_cytoscape.txt"
          , sep = "\t", row.names = FALSE)

##### find cancers SL-scan ####

id_mut = (SL_scan$Mutated_gene %in% tmp_freq_mut_bp$mut_gene)
id_ko = (SL_scan$KO_gene %in% tmp_freq_ko_bp$KO_gene)
inter_SLscan = (SL_scan[(id_ko & id_mut),c("cancer","Mutated_gene","KO_gene")])
colnames(inter_SLscan) = c("cancer" ,"mut_gene","KO_gene")
inter_SLscan$cancer = paste(inter_SLscan$cancer, "*")
write.table(inter_SLscan
            ,"../MetabolicSLOutput/enrichment_results/SL_scan/4_inter_SLscan_mut_KO_cancers_genes_vs_BPcommon_cytoscape.txt"
            , sep = "\t", row.names = FALSE)



#### SL-achilles four mode network of BP, KO_genes, driver_genes, and cancers ####
#####  prepare BP_achilles for four mode network #### 
BP_gene_divider = function(df){
  df_final = NULL
  for (dfRow_i in 1:nrow(df)){
    tmp = unlist(strsplit(df$core_enrichment[dfRow_i], split = "/"))
    df_final = rbind(df_final, cbind(df$mut_gene[dfRow_i], df$ID[dfRow_i]
                                     , df$Description[dfRow_i], df$NES[dfRow_i], df$p.adjust[dfRow_i], tmp))
    rm(dfRow_i)
  }
  colnames(df_final) = c("mut_gene","ID","Description", "NES","p.adjust","KO_gene")
  data.frame(df_final)
}

BP_slscan_divided = BP_gene_divider(BP_slscan)
write.csv(BP_slscan_divided
          ,"../MetabolicSLOutput/enrichment_results/SL_scan/SLscan_KoGenes_mutGene_vs_BP_.csv")


write.csv(freq_SL_achiles, "../MetabolicSLOutput/freq_cancer_fdrsig_achiles.csv")


freq_BP_slscan_divided = BP_slscan_divided %>% 
  group_by(Description, KO_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

# select top 95% quantile count of ko genes in bps 
freq_BP_slscan_divided = subset(freq_BP_slscan_divided, count 
                                > quantile(freq_BP_slscan_divided$count, 0.95))

freq_BP_slscan_divided_mutgene$count = BP_slscan_divided %>% 
  group_by(Description, mut_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))
inter_bp = intersect(freq_BP_slscan_divided_mutgene$Description, freq_BP_slscan_divided$Description)

tmp_freq_ko_bp = freq_BP_slscan_divided[freq_BP_slscan_divided$Description %in% inter_bp,]
tmp_freq_ko_bp$Description = paste(tmp_freq_ko_bp$Description, "_BP")
# colnames(tmp_freq_ko_bp) = c("fromNode", "toNode", "weight")
write.table(tmp_freq_ko_bp
            ,"../MetabolicSLOutput/enrichment_results/SL_scan/2_SLscan_KO_genes_vs_BPcommon_cytoscape.txt"
            , sep = "\t", row.names = FALSE)

tmp_freq_mut_bp = freq_BP_slscan_divided_mutgene[freq_BP_slscan_divided_mutgene$Description %in% inter_bp,]
tmp_freq_mut_bp$Description = paste(tmp_freq_mut_bp$Description, "_BP")
# colnames(tmp_freq_mut_bp) = c("fromNode", "toNode", "weight")
write.table(tmp_freq_mut_bp
            ,"../MetabolicSLOutput/enrichment_results/SL_scan/3_SLscan_mut_genes_vs_BPcommon_cytoscape.txt"
            , sep = "\t", row.names = FALSE)

##### find cancers SL-scan ####

id_mut = (SL_scan$Mutated_gene %in% tmp_freq_mut_bp$mut_gene)
id_ko = (SL_scan$KO_gene %in% tmp_freq_ko_bp$KO_gene)
inter_SLscan = (SL_scan[(id_ko & id_mut),c("cancer","Mutated_gene","KO_gene")])
colnames(inter_SLscan) = c("cancer" ,"mut_gene","KO_gene")
inter_SLscan$cancer = paste(inter_SLscan$cancer, "*")
write.table(inter_SLscan
            ,"../MetabolicSLOutput/enrichment_results/SL_scan/4_inter_SLscan_mut_KO_cancers_genes_vs_BPcommon_cytoscape.txt"
            , sep = "\t", row.names = FALSE)



# #### igraph ####
# can_ig_nodelist <- graph.data.frame(freq_BP_slscan_divided[,c("KO_gene", "Description")])
# 
# V(can_ig_nodelist)$type <- V(can_ig_nodelist)$name %in% freq_BP_slscan_divided$KO_gene
# 
# shapee <- shapes()[c(9,1)]
# 
# shape = shapee[as.numeric(V(can_ig_nodelist)$type+1)]
# #  color and shape mappings definiton.
# col <- c("lightblue", "blue")
# col = col[as.numeric(V(can_ig_nodelist)$type+1)]
# lnc_typ = V(can_ig_nodelist)$name[V(can_ig_nodelist)$type]
# 
# for (col_cntr in 1:length(lnc_typ)){
#   # V(can_ig_nodelist)$name[3]
#   # V(can_ig_nodelist)$type[col_cntr]
#   indx = which(V(can_ig_nodelist)$name[col_cntr]==freq_BP_slscan_divided$KO_gene)
#   if (sum(freq_BP_slscan_divided$count[indx])>0){
#     col[col_cntr] = "red"                 
#   }        
# }
# 
# edge_color = rep("blue3",nrow(freq_BP_slscan_divided))
# 
# x = freq_BP_slscan_divided$count
# edge_width = (x-min(x))/(max(x)-min(x)) #Normalized dg
# 
# for (j in 1:nrow(freq_BP_slscan_divided)){
#   tmp = freq_BP_slscan_divided[j, c("Description","KO_gene")]
#   ids = which((BP_slscan_divided$Description == tmp$Description) 
#         & (BP_slscan_divided$KO_gene == tmp$KO_gene))
#   if (mean(as.numeric(BP_slscan_divided$NES[ids])) > 0 ){
#     edge_color[j] <- "red3"
#   }
#   edge_width[j] =  (edge_width[j]+.5)^4
# }
# 
# pdf(paste0("../MetabolicSLOutput/enrichment_results/SL_scan/","1_SLscan_KO_genes_vs_BP.pdf"),width=30, height=30)
# 
# plot(can_ig_nodelist,
#      # layout=layout.circle,
#      # layout=layout.sphere,
#      # main="sphere",
#      layout = layout.fruchterman.reingold(can_ig_nodelist, niter=1000,area=5*vcount(can_ig_nodelist)^2),
#      # layout = l,
#      # label.color = "black",
#      # vertex.size = 10,
#      # vertex.label.angle = .5,
#      vertex.label.family="Helvetica",
#      # vertex.label.color="black",
#      vertex.label.size= 5 ,
#      edge.color = edge_color,
#      edge.width = edge_width,
#      edge.arrow.size = 1,
#      # arrow.width = .1 ,
#      # edge.arrow.color = "black" ,
#      vertex.label.color = "black",
#      # vertex.label.dist = 1,
#      vertex.color = col,
#      vertex.shape = shape,
#      vertex.size=15
#      # ,vertex.size2=15
# )
# mtext(paste0("SLscan_KO_genes_vs_BP"), side=1)
# 
# dev.off()


# save|load env
# save.image('./env/WS_MN_slscan_achilesTtest_GSEA_generalization.rdata')
# load('./env/WS_MN_slscan_achilesTtest_GSEA_generalization.rdata')

