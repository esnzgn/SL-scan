rm(list = ls())

library(dplyr)

# reviewer 1 experimental comarison ####
## comment 1) It is better to clearly list the number of metabolites and # of interactions on ####
# the metabolic networks used in different studies.

## comment 2) 2) It is important to list the number of available evaluated SLs in the existing datasets.  ####
# How many of them are validated; and how many are SLs (what is the definition of SL in these experiments, better than individual targets?
# do it using shRNA data of ataris ####
rm(list = ls())

library("dplyr")
library(ggplot2)
library(tidyr)

source('./code/func__tans_dom__get_sym.R')

# loading data
gen_ls = read.csv('../MetabolicSLInput/data/tb_gene_ls.csv')
KO_FBA = read.csv('../MetabolicSLOutput/KO_res_all.csv'); rownames(KO_FBA) = gen_ls$Var1;
KO_FBA = data.frame(t(KO_FBA))
rownames(KO_FBA) = gsub('_', '-',rownames(KO_FBA))
cantyp = read.csv('../MetabolicSLInput/data/cantyp.csv')
maf_met = read.csv('../MetabolicSLInput/data/maf_met_ge.csv')
maf_met$DepMap_ID = gsub('_', '-', maf_met$DepMap_ID)
sampleInfo <- read.csv("../MetabolicSLInput/data/sample_info.csv", stringsAsFactors = FALSE, header = TRUE)
# achiles = readRDS('../MetabolicSLInput/data/Achilles_gene_dependency.rds')
achiles = t(readRDS('../MetabolicSLInput/data/DRIVE_ATARiS_data.rds'))
# D2_DRIVE_gene_dep_scores.csv

rownames(achiles) = sampleInfo$DepMap_ID[match(toupper(rownames(achiles)), sampleInfo$CCLE_Name)]
achiles = achiles[!is.na(rownames(achiles)), ]

# subset achiles and KO_FBA for intersect met genes
# genes
intr_gen = intersect(colnames(achiles), colnames(KO_FBA))
achiles = achiles[,colnames(achiles) %in% intr_gen] 
achiles = achiles[,order(colnames(achiles))]
KO_FBA = KO_FBA[,colnames(KO_FBA) %in% intr_gen]; rm(intr_gen)
KO_FBA = KO_FBA[,order(colnames(KO_FBA))]
# cell lines
intr_cel = intersect(rownames(achiles), rownames(KO_FBA))
achiles = achiles[rownames(achiles) %in% intr_cel,] 
achiles = achiles[order(rownames(achiles)),]
KO_FBA = KO_FBA[rownames(KO_FBA) %in% intr_cel,]; rm(intr_cel)
KO_FBA = KO_FBA[order(rownames(KO_FBA)),]
# make map of cell lines names and diseases by a two col df
pd_cel = data.frame(cbind(primary_disease = sampleInfo$primary_disease[match(rownames(achiles)
                                                                             , sampleInfo$DepMap_ID)], DepMap_ID = rownames(achiles)))

# table(pri_fre_ac$pd)
ac_cantyp = ach_can_ord(pd_cel, cantyp$cancer);
# sampleInfo$primary_disease[match(ac_cantyp$`Bile Duct Cancer`, sampleInfo$DepMap_ID)]

altn = 'less'
thre = 2
# parallelization parameters
# install.packages('foreach')
# BiocManager::install('foreach')
library('foreach')
nrCores =  parallel::detectCores() - 1
# nrCores = nrCores[1]-1 #not to overload your computer
# nrCores = 4

if(Sys.info()['sysname'] != "Windows") {
  library(doMC)
  registerDoMC(nrCores)
  getDoParWorkers()
} else {
  # install.packages('doParallel')
  library(doParallel)
  cl <- makeCluster(nrCores)
  registerDoParallel(cl)
}
fin_can = NULL
fin_can = foreach(can_i= 1:(length(ac_cantyp)), .combine=rbind) %dopar%{
  # for (can_i in c(6,7,9,12,14)){
  # print(rep('###', 10));print(cantyp$X1[can_i])
  can = names(ac_cantyp[can_i]);
  can_cel = ac_cantyp[[can]]
  if (length(can_cel) > 1){
    # tmp_mut = mutbck_ls[, srt:end]
    tmp_mut = mutbackmkr(colnames(achiles), maf_met, can_cel)
    # head(rowSums(tmp_mut)[order(rowSums(tmp_mut), decreasing = T)])
    # slice achiles for the selected cell lines of the cancer
    id_cel = match(can_cel, pd_cel$DepMap_ID)
    KO_ac = achiles[(id_cel[!is.na(id_cel)]), ]
    KO_mn = KO_FBA[(id_cel[!is.na(id_cel)]), ]
    for (mn_i in 1:ncol(KO_mn)){ 
      KO_mn[is.na(KO_mn[,mn_i]),mn_i] = 1
    }
    rm(id_cel)
    can_tbl = NULL
    # find each cancer A and B genes
    id_GnA = which(rowSums(tmp_mut) >= thre)
    
    if (length(id_GnA) > 0){
      GnA = names(id_GnA)
      GnB = colnames(KO_ac)
      
      for(i in 1:length(GnA)){
        # print(i/length(GnA)*100)
        res_cel = celselector(GnA[i], tmp_mut)
        for (j in 1:length(GnB)){
          # ttest
          t_res = NULL
          t_res = tw_tst(res_cel, GnB[j], KO_ac, KO_mn, altn, thre)
          if (!is.na(t_res[1])){          
            if (!is.na(t_res$ttst_obs$p.value)){
              can_tbl = rbind(can_tbl ,  cbind('cancer' = can, 'Mutated_gene' = GnA[i], 'KO_gene' = GnB[j]
                                               , 'No_mut' = length(res_cel$mut), 'No_not_mut' = length(res_cel$not_mut)
                                               , 'ttst_obs_stat' =  t_res$ttst_obs$statistic, 'ttst_obs_pvl' = t_res$ttst_obs$p.value
                                               # , 'wtst_obs_stat' =  t_res$wtst_obs$statistic, 'wtst_obs_pvl' = t_res$wtst_obs$p.value
                                               , 'ttst_pre_stat' =  t_res$ttst_pre$statistic, 'ttst_pre_pvl' = t_res$ttst_pre$p.value
                                               # , 'wtst_pre_stat' =  t_res$wtst_pre$statistic, 'wtst_pre_pvl' = t_res$wtst_pre$p.value
              ))
            }
          }
        }
      }
    }
    write.csv(can_tbl, file = paste0('../MetabolicSLOutput/', can,'_at_mn_w&t_tst_cons.csv'))
    can_tbl 
  }
  # print(paste0(round(can_i / nrow(cantyp) * 100, 1), ' % current cancer'))
}

write.csv(fin_can, file = '../MetabolicSLOutput/fin_can_at_mn_w&ttst_cons_ataris.csv')
fin_can = data.frame(fin_can)

# read the file of each cancer adjust them and agian write them in the same place
fname = list.files(pattern ='_at_mn_w&t_tst_cons.csv', '../MetabolicSLOutput/')

#### fdr ####
for(f_i in 1:length(fname)){
  print(f_i / length(fname))
  tmp = data.frame(NULL)
  tryCatch({
    tmp = read.csv(paste0('../MetabolicSLOutput/',fname[f_i]))
  }, error = function(e) {
    if (grepl("first five rows are empty: giving up", e)) {
      tmp = data.frame()
      message("The file is empty.")
    } else {
      stop(e)
    }
  })
  if (nrow(tmp) != 0){
    tmp = tmp[,-1]
    colnames(tmp) = c("cancer","Mutated_gene","KO_gene","No_mut","No_not_mut","ttst_ataris_stat", "ttst_ataris_pvl"
                      # ,"wtst_ataris_stat", "wtst_ataris_pvl"
                      , "ttst_metabolic_stat", "ttst_metabolic_pvl"
                      # , "wtst_metabolic_stat", "wtst_pre_pvl"
    )
    tmp = cbind(tmp
                , ttst_ataris_fdr = p.adjust(as.numeric(tmp$ttst_ataris_pvl), method = 'fdr')
                , ttst_metabolic_fdr = p.adjust(as.numeric(tmp$ttst_metabolic_pvl), method = 'fdr'))
    # tmp = cbind(tmp
    #             , ttst_achiles_fdr = (as.numeric(tmp$ttst_achiles_pvl))
    #             , ttst_metabolic_fdr = (as.numeric(tmp$ttst_metabolic_pvl)))
    name = unlist(strsplit(fname[f_i], split = '_at_mn_w&t_tst_cons.csv'))[1]
    write.csv(tmp, paste0('../MetabolicSLOutput/',name, '_at_mn_w&t_tst_cons_fdr.csv'))
  }
}; #rm(tmp)
# adj_final_all_ac = subset(final_all_ac, final_all_ac$fdr_tt <= 0.05)

# hypergeometric test shRNA ####
#### SL-scan shRNA hypergeometric test ####
fname = list.files(pattern ='_at_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_SLscan_shRNA = data.frame(matrix(NA, nrow = length(fname), ncol = 6))
colnames(cantyp_hyp_SLscan_shRNA) = c('cancer', 'hypergeom', 'no_SLscan', 'no_SL_atr', 'overlp', 'total')
SL_scan = NULL
SL_scan_valid = NULL
for (cn in 1:length(fname)){
  print(cn/length(fname)*100)
  cantyp_hyp_SLscan_shRNA$cancer[cn] = unlist(strsplit(fname[cn], '_at_mn_w&t_tst_cons_fdr.csv')[1])
  tmp = read.csv(paste0('../MetabolicSLOutput/',fname[cn]))
  mn = which(as.numeric(tmp$ttst_metabolic_pvl) <= 0.05)
  atr = which(as.numeric(tmp$ttst_ataris_pvl) <= 0.05)
  SL_scan = rbind(SL_scan , tmp[(mn),])
  SL_scan_valid = rbind(SL_scan_valid , tmp[intersect(mn, atr),])
  cantyp_hyp_SLscan_shRNA$overlp[cn] = length(intersect(mn, atr))
  cantyp_hyp_SLscan_shRNA$no_SL_atr[cn] = length(atr)
  cantyp_hyp_SLscan_shRNA$no_SLscan[cn] = length(mn)
  cantyp_hyp_SLscan_shRNA$total[cn] = nrow(tmp)
  
  cantyp_hyp_SLscan_shRNA$hypergeom[cn] = phyper(cantyp_hyp_SLscan_shRNA$overlp[cn]-1, cantyp_hyp_SLscan_shRNA$no_SLscan[cn]
                                                 , cantyp_hyp_SLscan_shRNA$total[cn]-cantyp_hyp_SLscan_shRNA$no_SLscan[cn]
                                                 , cantyp_hyp_SLscan_shRNA$no_SL_atr[cn], lower.tail = FALSE)
}

cantyp_hyp_SLscan_shRNA = rbind(cantyp_hyp_SLscan_shRNA, NA)

cantyp_hyp_SLscan_shRNA$cancer[cn+1] = "Pan_cancer"
id_exclude = !is.na(cantyp_hyp_SLscan_shRNA$hypergeom)
cantyp_hyp_SLscan_shRNA$no_SL_atr[cn+1] = sum(cantyp_hyp_SLscan_shRNA$no_SL_atr[id_exclude])
cantyp_hyp_SLscan_shRNA$no_SLscan[cn+1] = sum(cantyp_hyp_SLscan_shRNA$no_SLscan[id_exclude])
cantyp_hyp_SLscan_shRNA$overlp[cn+1] = sum(cantyp_hyp_SLscan_shRNA$overlp[id_exclude])
cantyp_hyp_SLscan_shRNA$total[cn+1] = sum(cantyp_hyp_SLscan_shRNA$total[id_exclude])

# hypergeom of pvalue sig slscan and synlethdb records
cantyp_hyp_SLscan_shRNA$hypergeom[cn+1] = phyper(cantyp_hyp_SLscan_shRNA$overlp[cn+1]-1, cantyp_hyp_SLscan_shRNA$no_SLscan[cn+1]
                                                 , cantyp_hyp_SLscan_shRNA$total[cn+1]-cantyp_hyp_SLscan_shRNA$no_SLscan[cn+1]
                                                 , cantyp_hyp_SLscan_shRNA$no_SL_atr[cn+1], lower.tail = FALSE)

# write.csv(cantyp_hyp, '../MetabolicSLOutput/hyper_table_fdrcorrect_cons_ataris.csv')
# cortest_at_slscan = cor.test(cantyp_hyp$no_SL_MN, cantyp_hyp$no_SL_atr, method = "pearson")
# TP <- sum(cantyp_hyp$overlp)
# FP <- abs(sum(cantyp_hyp$no_SL_MN)-sum(cantyp_hyp$overlp))
# FN <- abs(sum(cantyp_hyp$no_SL_atr)-sum(cantyp_hyp$overlp))
# TN <- abs((sum(cantyp_hyp$total)-sum(cantyp_hyp$no_SL_atr))-sum(cantyp_hyp$overlp))
# 
# accuracy <- (TP + TN) / sum(TP,FP,TN,FN)
# precision <- TP / (TP + FP)
# recall <- TP / (TP + FN)
# f1_score <- 2 * precision * recall / (precision + recall)
# 
# 
# # Create a matrix for the confusion matrix
# confusion_matrix <- matrix(c(TP, FP, FN, TN), ncol = 2, byrow = TRUE)
# conf_matrix <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE,
#                       dimnames = list(c("True", "False"), c("Positive", "Negative")))

#### FastSL hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_FastSL' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_FastSL = read.csv("../MetabolicSLOutput/R24_iMAT_cons_FastSL.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_FastSL),2)){ can_ls[i] = unlist(strsplit(colnames(R24_iMAT_FastSL)[i], "[.|_]"))[1]}
fname = list.files(pattern ='_at_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_FastSL_shRNA = data.frame(matrix(NA, nrow = length(fname), ncol = 8))
colnames(cantyp_hyp_FastSL_shRNA) = c('cancer', 'hypergeom', 'no_SL_shRNA', 'no_SL_FastSL', 'overlp_FastSL'
                                      , 'total', "no_SLscan", 'no_overlap_SLscan')
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_at_mn_w&t_tst_cons_fdr.csv')[1])
  # rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  
  cantyp_hyp_FastSL_shRNA$cancer[i] = can
  if (length(srt) == 0){
    next
  }
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_FastSL[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_FastSL[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  
  tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[i]))
  # tmp_gene = intersect(unique(c(tmp_ac$Mutated_gene, tmp_ac$KO_gene)), tmp_gls)
  # total = expand.grid(Mutated_gene = tmp_gene , KO_gene = tmp_gene)
  tmp_mut = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
  tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
  tmp_ko = tmp_ko[!(tmp_ko %in% tmp_mut)]
  total = expand.grid(Mutated_gene = (tmp_mut) , KO_gene = tmp_ko)
  
  dup_sl = cbind(sl$col2, sl$col1)
  colnames(dup_sl) = colnames(sl)
  dup_sl = rbind(dup_sl,sl)
  dup_sl = cbind(dup_sl, Algo = "FastSL")
  colnames(dup_sl)[c(1,2)] = c("Mutated_gene", "KO_gene")
  
  cantyp_hyp_FastSL_shRNA$total[i] = nrow(total)
  if (nrow(total) < 1){
    next
  }
  # hyper_geom test
  cantyp_hyp_FastSL_shRNA$overlp_FastSL[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), dup_sl))
  cantyp_hyp_FastSL_shRNA$no_SL_shRNA[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), total))
  cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i] = nrow(dup_sl[!duplicated(dup_sl),])/2
  # cantyp_hyp_FastSL_shRNA$total[i] = nrow(total)
  # cantyp_hyp_FastSL_shRNA$no_overlap_SLscan[i] = length(intersect(which(chek_concordance$SLscan == "YES")
  #                                                         , which(chek_concordance$CRISPR == "YES")))
  cantyp_hyp_FastSL_shRNA$no_SLscan[i] = nrow(inner_join(subset(tmp_ac, ttst_metabolic_fdr <= 0.05), total))
  cantyp_hyp_FastSL_shRNA$no_overlap_SLscan[i] = nrow(inner_join(subset(tmp_ac,(ttst_metabolic_fdr <= 0.05 
                                                                                & ttst_ataris_fdr <= 0.05)), total))
  
  cantyp_hyp_FastSL_shRNA$hypergeom[i] = phyper(cantyp_hyp_FastSL_shRNA$overlp[i]-1, cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i]
                                                , cantyp_hyp_FastSL_shRNA$total[i]-cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i]
                                                , cantyp_hyp_FastSL_shRNA$no_SL_shRNA[i], lower.tail = FALSE)
  
  rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
}

cantyp_hyp_FastSL_shRNA = rbind(cantyp_hyp_FastSL_shRNA, NA)

cantyp_hyp_FastSL_shRNA$cancer[i+1] = "Pan_cancer"
id_exclude = !is.na(cantyp_hyp_FastSL_shRNA$hypergeom)
cantyp_hyp_FastSL_shRNA$no_SL_shRNA[i+1] = sum(cantyp_hyp_FastSL_shRNA$no_SL_shRNA[id_exclude])
cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i+1] = sum(cantyp_hyp_FastSL_shRNA$no_SL_FastSL[id_exclude])
cantyp_hyp_FastSL_shRNA$overlp_FastSL[i+1] = sum(cantyp_hyp_FastSL_shRNA$overlp_FastSL[id_exclude])
cantyp_hyp_FastSL_shRNA$total[i+1] = sum(cantyp_hyp_FastSL_shRNA$total[id_exclude])
cantyp_hyp_FastSL_shRNA$no_SLscan[i+1] = sum(cantyp_hyp_FastSL_shRNA$no_SLscan[id_exclude])
cantyp_hyp_FastSL_shRNA$no_overlap_SLscan[i+1] = sum(cantyp_hyp_FastSL_shRNA$no_overlap_SLscan[id_exclude])

# hypergeom of pvalue sig slscan and synlethdb records
cantyp_hyp_FastSL_shRNA$hypergeom[i+1] = phyper(cantyp_hyp_FastSL_shRNA$overlp_FastSL[i+1]-1, cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i+1]
                                                , cantyp_hyp_FastSL_shRNA$total[i+1]-cantyp_hyp_FastSL_shRNA$no_SL_FastSL[i+1]
                                                , cantyp_hyp_FastSL_shRNA$no_SL_shRNA[i+1], lower.tail = FALSE)

# cantyp_hyp_FastSL_shRNA = subset(cantyp_hyp_FastSL_shRNA, !is.na(hypergeom))


#### MCS hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_MCS' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_MCS = read.csv("../MetabolicSLOutput/R24_iMAT_cons_MCS.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_MCS),2)){ can_ls[i] = unlist(strsplit(colnames(R24_iMAT_MCS)[i], "[.|_]"))[1]}
fname = list.files(pattern ='_at_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_MCS_shRNA = data.frame(matrix(NA, nrow = length(fname), ncol = 8))
colnames(cantyp_hyp_MCS_shRNA) = c('cancer', 'hypergeom', 'no_SL_shRNA', 'no_SL_MCS', 'overlp_MCS'
                                   , 'total', "no_SLscan", 'no_overlap_SLscan')
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_at_mn_w&t_tst_cons_fdr.csv')[1])
  # rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  
  cantyp_hyp_MCS_shRNA$cancer[i] = can
  if (length(srt) == 0){
    next
  }
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_MCS[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_MCS[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  
  tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[i]))
  # tmp_gene = intersect(unique(c(tmp_ac$Mutated_gene, tmp_ac$KO_gene)), tmp_gls)
  # total = expand.grid(Mutated_gene = tmp_gene , KO_gene = tmp_gene)
  tmp_mut = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
  tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
  tmp_ko = tmp_ko[!(tmp_ko %in% tmp_mut)]
  total = expand.grid(Mutated_gene = (tmp_mut) , KO_gene = tmp_ko)
  
  dup_sl = cbind(sl$col2, sl$col1)
  colnames(dup_sl) = colnames(sl)
  dup_sl = rbind(dup_sl,sl)
  dup_sl = cbind(dup_sl, Algo = "MCS")
  colnames(dup_sl)[c(1,2)] = c("Mutated_gene", "KO_gene")
  
  cantyp_hyp_MCS_shRNA$total[i] = nrow(total)
  if (nrow(total) < 1){
    next
  }
  # hyper_geom test
  cantyp_hyp_MCS_shRNA$overlp_MCS[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), dup_sl))
  cantyp_hyp_MCS_shRNA$no_SL_shRNA[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), total))
  cantyp_hyp_MCS_shRNA$no_SL_MCS[i] = nrow(dup_sl[!duplicated(dup_sl),])/2
  # cantyp_hyp_MCS_shRNA$total[i] = nrow(total)
  # cantyp_hyp_MCS_shRNA$no_overlap_SLscan[i] = length(intersect(which(chek_concordance$SLscan == "YES")
  #                                                         , which(chek_concordance$CRISPR == "YES")))
  cantyp_hyp_MCS_shRNA$no_SLscan[i] = nrow(inner_join(subset(tmp_ac, ttst_metabolic_fdr <= 0.05), total))
  cantyp_hyp_MCS_shRNA$no_overlap_SLscan[i] = nrow(inner_join(subset(tmp_ac,(ttst_metabolic_fdr <= 0.05 
                                                                             & ttst_ataris_fdr <= 0.05)), total))
  
  cantyp_hyp_MCS_shRNA$hypergeom[i] = phyper(cantyp_hyp_MCS_shRNA$overlp[i]-1, cantyp_hyp_MCS_shRNA$no_SL_MCS[i]
                                             , cantyp_hyp_MCS_shRNA$total[i]-cantyp_hyp_MCS_shRNA$no_SL_MCS[i]
                                             , cantyp_hyp_MCS_shRNA$no_SL_shRNA[i], lower.tail = FALSE)
  
  rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
}

cantyp_hyp_MCS_shRNA = rbind(cantyp_hyp_MCS_shRNA, NA)

cantyp_hyp_MCS_shRNA$cancer[i+1] = "Pan_cancer"
id_exclude = !is.na(cantyp_hyp_MCS_shRNA$hypergeom)
cantyp_hyp_MCS_shRNA$no_SL_shRNA[i+1] = sum(cantyp_hyp_MCS_shRNA$no_SL_shRNA[id_exclude])
cantyp_hyp_MCS_shRNA$no_SL_MCS[i+1] = sum(cantyp_hyp_MCS_shRNA$no_SL_MCS[id_exclude])
cantyp_hyp_MCS_shRNA$overlp_MCS[i+1] = sum(cantyp_hyp_MCS_shRNA$overlp_MCS[id_exclude])
cantyp_hyp_MCS_shRNA$total[i+1] = sum(cantyp_hyp_MCS_shRNA$total[id_exclude])
cantyp_hyp_MCS_shRNA$no_SLscan[i+1] = sum(cantyp_hyp_MCS_shRNA$no_SLscan[id_exclude])
cantyp_hyp_MCS_shRNA$no_overlap_SLscan[i+1] = sum(cantyp_hyp_MCS_shRNA$no_overlap_SLscan[id_exclude])

# hypergeom of pvalue sig slscan and synlethdb records
cantyp_hyp_MCS_shRNA$hypergeom[i+1] = phyper(cantyp_hyp_MCS_shRNA$overlp_MCS[i+1]-1, cantyp_hyp_MCS_shRNA$no_SL_MCS[i+1]
                                             , cantyp_hyp_MCS_shRNA$total[i+1]-cantyp_hyp_MCS_shRNA$no_SL_MCS[i+1]
                                             , cantyp_hyp_MCS_shRNA$no_SL_shRNA[i+1], lower.tail = FALSE)

# cantyp_hyp_MCS_shRNA = subset(cantyp_hyp_MCS_shRNA, !is.na(hypergeom))

#### gMCS hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_gmcs' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_gmcs = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gmcs.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_gmcs),2)){ can_ls[i] = unlist(strsplit(colnames(R24_iMAT_gmcs)[i], "[.|_]"))[1]}
fname = list.files(pattern ='_at_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_gMCS_shRNA = data.frame(matrix(NA, nrow = length(fname), ncol = 8))
colnames(cantyp_hyp_gMCS_shRNA) = c('cancer', 'hypergeom', 'no_SL_shRNA', 'no_SL_gmcs', 'overlp_gmcs'
                                    , 'total', "no_SLscan", 'no_overlap_SLscan')
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_at_mn_w&t_tst_cons_fdr.csv')[1])
  # rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  
  cantyp_hyp_gMCS_shRNA$cancer[i] = can
  if (length(srt) == 0){
    next
  }
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_gmcs[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_gmcs[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  
  tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[i]))
  # tmp_gene = intersect(unique(c(tmp_ac$Mutated_gene, tmp_ac$KO_gene)), tmp_gls)
  # total = expand.grid(Mutated_gene = tmp_gene , KO_gene = tmp_gene)
  tmp_mut = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
  tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
  tmp_ko = tmp_ko[!(tmp_ko %in% tmp_mut)]
  total = expand.grid(Mutated_gene = (tmp_mut) , KO_gene = tmp_ko)
  
  dup_sl = cbind(sl$col2, sl$col1)
  colnames(dup_sl) = colnames(sl)
  dup_sl = rbind(dup_sl,sl)
  dup_sl = cbind(dup_sl, Algo = "gMCS")
  colnames(dup_sl)[c(1,2)] = c("Mutated_gene", "KO_gene")
  
  cantyp_hyp_gMCS_shRNA$total[i] = nrow(total)
  if (nrow(total) < 1){
    next
  }
  # hyper_geom test
  cantyp_hyp_gMCS_shRNA$overlp_gmcs[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), dup_sl))
  cantyp_hyp_gMCS_shRNA$no_SL_shRNA[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), total))
  cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i] = nrow(dup_sl[!duplicated(dup_sl),])/2
  # cantyp_hyp_gMCS_shRNA$total[i] = nrow(total)
  # cantyp_hyp_gMCS_shRNA$no_overlap_SLscan[i] = length(intersect(which(chek_concordance$SLscan == "YES")
  #                                                         , which(chek_concordance$CRISPR == "YES")))
  cantyp_hyp_gMCS_shRNA$no_SLscan[i] = nrow(inner_join(subset(tmp_ac, ttst_metabolic_fdr <= 0.05), total))
  cantyp_hyp_gMCS_shRNA$no_overlap_SLscan[i] = nrow(inner_join(subset(tmp_ac,(ttst_metabolic_fdr <= 0.05 
                                                                              & ttst_ataris_fdr <= 0.05)), total))
  
  cantyp_hyp_gMCS_shRNA$hypergeom[i] = phyper(cantyp_hyp_gMCS_shRNA$overlp[i]-1, cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i]
                                              , cantyp_hyp_gMCS_shRNA$total[i]-cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i]
                                              , cantyp_hyp_gMCS_shRNA$no_SL_shRNA[i], lower.tail = FALSE)
  
  rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
}

cantyp_hyp_gMCS_shRNA = rbind(cantyp_hyp_gMCS_shRNA, NA)

cantyp_hyp_gMCS_shRNA$cancer[i+1] = "Pan_cancer"
id_exclude = !is.na(cantyp_hyp_gMCS_shRNA$hypergeom)
cantyp_hyp_gMCS_shRNA$no_SL_shRNA[i+1] = sum(cantyp_hyp_gMCS_shRNA$no_SL_shRNA[id_exclude])
cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i+1] = sum(cantyp_hyp_gMCS_shRNA$no_SL_gmcs[id_exclude])
cantyp_hyp_gMCS_shRNA$overlp_gmcs[i+1] = sum(cantyp_hyp_gMCS_shRNA$overlp_gmcs[id_exclude])
cantyp_hyp_gMCS_shRNA$total[i+1] = sum(cantyp_hyp_gMCS_shRNA$total[id_exclude])
cantyp_hyp_gMCS_shRNA$no_SLscan[i+1] = sum(cantyp_hyp_gMCS_shRNA$no_SLscan[id_exclude])
cantyp_hyp_gMCS_shRNA$no_overlap_SLscan[i+1] = sum(cantyp_hyp_gMCS_shRNA$no_overlap_SLscan[id_exclude])

# hypergeom of pvalue sig slscan and synlethdb records
cantyp_hyp_gMCS_shRNA$hypergeom[i+1] = phyper(cantyp_hyp_gMCS_shRNA$overlp_gmcs[i+1]-1, cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i+1]
                                              , cantyp_hyp_gMCS_shRNA$total[i+1]-cantyp_hyp_gMCS_shRNA$no_SL_gmcs[i+1]
                                              , cantyp_hyp_gMCS_shRNA$no_SL_shRNA[i+1], lower.tail = FALSE)

# cantyp_hyp_gMCS_shRNA = subset(cantyp_hyp_gMCS_shRNA, !is.na(hypergeom))

#### ngMCS hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_ngMCS' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_ngMCS = read.csv("../MetabolicSLOutput/R24_iMAT_cons_ngMCS.csv")
can_ls = NULL
colnames(R24_iMAT_ngMCS) = colnames(R24_iMAT_gmcs)
for (i in seq(1,ncol(R24_iMAT_ngMCS),2)){ can_ls[i] = unlist(strsplit(colnames(R24_iMAT_ngMCS)[i], "[.|_]"))[1]}
fname = list.files(pattern ='_at_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_ngMCS_shRNA = data.frame(matrix(NA, nrow = length(fname), ncol = 8))
colnames(cantyp_hyp_ngMCS_shRNA) = c('cancer', 'hypergeom', 'no_SL_shRNA', 'no_SL_ngMCS', 'overlp_ngMCS'
                                     , 'total', "no_SLscan", 'no_overlap_SLscan')
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_at_mn_w&t_tst_cons_fdr.csv')[1])
  # rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  
  cantyp_hyp_ngMCS_shRNA$cancer[i] = can
  if (length(srt) == 0){
    next
  }
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_ngMCS[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_ngMCS[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  
  tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[i]))
  # tmp_gene = intersect(unique(c(tmp_ac$Mutated_gene, tmp_ac$KO_gene)), tmp_gls)
  # total = expand.grid(Mutated_gene = tmp_gene , KO_gene = tmp_gene)
  tmp_mut = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
  tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
  tmp_ko = tmp_ko[!(tmp_ko %in% tmp_mut)]
  total = expand.grid(Mutated_gene = (tmp_mut) , KO_gene = tmp_ko)
  
  dup_sl = cbind(sl$col2, sl$col1)
  colnames(dup_sl) = colnames(sl)
  dup_sl = rbind(dup_sl,sl)
  dup_sl = cbind(dup_sl, Algo = "ngMCS")
  colnames(dup_sl)[c(1,2)] = c("Mutated_gene", "KO_gene")
  
  cantyp_hyp_ngMCS_shRNA$total[i] = nrow(total)
  if (nrow(total) < 1){
    next
  }
  # hyper_geom test
  cantyp_hyp_ngMCS_shRNA$overlp_ngMCS[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), dup_sl))
  cantyp_hyp_ngMCS_shRNA$no_SL_shRNA[i] = nrow(inner_join(subset(tmp_ac, ttst_ataris_fdr <= 0.05), total))
  cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i] = nrow(dup_sl[!duplicated(dup_sl),])/2
  # cantyp_hyp_ngMCS_shRNA$total[i] = nrow(total)
  # cantyp_hyp_ngMCS_shRNA$no_overlap_SLscan[i] = length(intersect(which(chek_concordance$SLscan == "YES")
  #                                                         , which(chek_concordance$CRISPR == "YES")))
  cantyp_hyp_ngMCS_shRNA$no_SLscan[i] = nrow(inner_join(subset(tmp_ac, ttst_metabolic_fdr <= 0.05), total))
  cantyp_hyp_ngMCS_shRNA$no_overlap_SLscan[i] = nrow(inner_join(subset(tmp_ac,(ttst_metabolic_fdr <= 0.05 
                                                                               & ttst_ataris_fdr <= 0.05)), total))
  
  cantyp_hyp_ngMCS_shRNA$hypergeom[i] = phyper(cantyp_hyp_ngMCS_shRNA$overlp[i]-1, cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i]
                                               , cantyp_hyp_ngMCS_shRNA$total[i]-cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i]
                                               , cantyp_hyp_ngMCS_shRNA$no_SL_shRNA[i], lower.tail = FALSE)
  
  rm(thr_fdr, overlap, tmp_ac, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
}

cantyp_hyp_ngMCS_shRNA = rbind(cantyp_hyp_ngMCS_shRNA, NA)

cantyp_hyp_ngMCS_shRNA$cancer[i+1] = "Pan_cancer"
id_exclude = !is.na(cantyp_hyp_ngMCS_shRNA$hypergeom)
cantyp_hyp_ngMCS_shRNA$no_SL_shRNA[i+1] = sum(cantyp_hyp_ngMCS_shRNA$no_SL_shRNA[id_exclude])
cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i+1] = sum(cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[id_exclude])
cantyp_hyp_ngMCS_shRNA$overlp_ngMCS[i+1] = sum(cantyp_hyp_ngMCS_shRNA$overlp_ngMCS[id_exclude])
cantyp_hyp_ngMCS_shRNA$total[i+1] = sum(cantyp_hyp_ngMCS_shRNA$total[id_exclude])
cantyp_hyp_ngMCS_shRNA$no_SLscan[i+1] = sum(cantyp_hyp_ngMCS_shRNA$no_SLscan[id_exclude])
cantyp_hyp_ngMCS_shRNA$no_overlap_SLscan[i+1] = sum(cantyp_hyp_ngMCS_shRNA$no_overlap_SLscan[id_exclude])

# hypergeom of pvalue sig slscan and synlethdb records
cantyp_hyp_ngMCS_shRNA$hypergeom[i+1] = phyper(cantyp_hyp_ngMCS_shRNA$overlp_ngMCS[i+1]-1, cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i+1]
                                               , cantyp_hyp_ngMCS_shRNA$total[i+1]-cantyp_hyp_ngMCS_shRNA$no_SL_ngMCS[i+1]
                                               , cantyp_hyp_ngMCS_shRNA$no_SL_shRNA[i+1], lower.tail = FALSE)

# cantyp_hyp_ngMCS_shRNA = subset(cantyp_hyp_ngMCS_shRNA, !is.na(hypergeom))

# ## overlap ngMCS and gMCS ####
# gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
# list.files(pattern = 'R24_iMAT_cons_ngmcs_Wtanslated_EX.csv' ,'../MetabolicSLOutput/', ignore.case = T)
# R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
# R24_iMAT_gMCS = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gMCS.csv")
# R24_iMAT_ngMCS = read.csv("../MetabolicSLOutput/R24_iMAT_cons_ngmcs_Wtanslated_EX.csv")
# colnames(R24_iMAT_ngMCS) = colnames(R24_iMAT_gMCS)
# cantyp_hyp_ngMCS_gMCS = data.frame(matrix(NA, nrow = ncol(R24_iMAT_gMCS)/2, ncol = 6))
# can_ls = unique(sapply(colnames(R24_iMAT_ngMCS), function(x) unlist(strsplit(x, '[.]'))[1]))
# colnames(cantyp_hyp_ngMCS_gMCS) = c('cancer', 'hypergeom', 'no_SL_ngMCS', 'no_SL_gMCS', 'overlp', 'total')
# for(i in seq(can_ls)){
#   print(i/length(fname)*100)
#   can = can_ls[i]
#   rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
#   id_gls = grep(can, colnames(R24_iMAT_ngMCS))[1]
#   if (is.na(id_gls)){
#     message(paste0(">>>>>>>>>>>> no metaboilc gene list for: ", can))
#     cantyp_hyp_ngMCS_gMCS[i,] = NA
#     cantyp_hyp_ngMCS_gMCS$cancer[i] = can
#     next
#   }
#   # SL pairs of ngMCS
#   srt = grep(((can)), colnames(R24_iMAT_ngMCS))[1]; endd = srt + 1;
#   col1 = entrz2hgnc(paste0(R24_iMAT_ngMCS[,srt],'.1'), gen_ls)
#   col2 = entrz2hgnc(paste0(R24_iMAT_ngMCS[,endd],'.1'), gen_ls)
#   sl = data.frame(cbind(col1,col2))
#   sl = sl[!is.na(sl$col1),]
#   # sl_ngmcs = cbind(sl$col2, sl$col1)
#   # colnames(sl_ngmcs) = colnames(sl)
#   # sl_ngmcs = rbind(sl, sl_ngmcs)
#   sl_ngmcs = sl
#   
#   gls_WOdot = as.character(unique(lapply(as.character(R24_iMAT_gls[,i]), function(x) unlist(strsplit(x, '[.]'))[1])))
#   total_can = NULL
#   total_can = (length(!is.na(gls_WOdot)) * length(!is.na(gls_WOdot))-1)/2
#   
#   # SL pairs of ngMCS
#   srt_gmcs = grep(((can)), colnames(R24_iMAT_gMCS))[1]; endd_gmcs = srt_gmcs + 1;
#   col1_gmcs = entrz2hgnc(paste0(R24_iMAT_gMCS[,srt_gmcs],'.1'), gen_ls)
#   col2_gmcs = entrz2hgnc(paste0(R24_iMAT_gMCS[,endd_gmcs],'.1'), gen_ls)
#   sl_gmcs = data.frame(cbind(col1_gmcs,col2_gmcs))
#   sl_gmcs = sl_gmcs[!is.na(sl_gmcs$col1_gmcs),]
#   colnames(sl_gmcs) = colnames(sl_ngmcs)
#   
#   join_res = inner_join(sl_ngmcs, sl_gmcs)
#   
#   # hyper_geom test
#   cantyp_hyp_ngMCS_gMCS$cancer[i] = can
#   cantyp_hyp_ngMCS_gMCS$overlp[i] = nrow(join_res)
#   cantyp_hyp_ngMCS_gMCS$no_SL_ngMCS[i] = nrow(sl_ngmcs) # nrow(sl_ngmcs)/2
#   cantyp_hyp_ngMCS_gMCS$no_SL_gMCS[i] = nrow(sl_gmcs)
#   cantyp_hyp_ngMCS_gMCS$total[i] = total_can
#   
#   cantyp_hyp_ngMCS_gMCS$hypergeom[i] = phyper(cantyp_hyp_ngMCS_gMCS$overlp[i]-1, cantyp_hyp_ngMCS_gMCS$no_SL_ngMCS[i]
#                                                , cantyp_hyp_ngMCS_gMCS$total[i]-cantyp_hyp_ngMCS_gMCS$no_SL_ngMCS[i]
#                                                , cantyp_hyp_ngMCS_gMCS$no_SL_gMCS[i], lower.tail = FALSE)
#   
# }


rm(ls_shRNA_res)
for (i in 1:length(ls(pattern = c('_shRNA')))){
  tmpn = ls(pattern = '_shRNA')[i]  
  tmp = get(tmpn)
  write.csv(tmp, paste0('../MetabolicSLOutput/reviewer comments/cons_ataris_', tmpn,'.csv'))
}

df_SL_predict_algoshRNA = cbind("cancer" = cantyp_hyp_FastSL_shRNA$cancer
                                ,"#Total" = cantyp_hyp_FastSL_shRNA$total
                                , "#shRNA" = cantyp_hyp_FastSL_shRNA$no_SL_shRNA
                                , "#SLscan" = cantyp_hyp_FastSL_shRNA$no_SLscan
                                , "#overlap_SLscan" = cantyp_hyp_FastSL_shRNA$no_overlap_SLscan)
ls_shRNA_res = ls(pattern = '_shRNA')
ls_shRNA_res = ls_shRNA_res[!grepl('cantyp_hyp_SLscan_shRNA', ls_shRNA_res)]
for (i in seq(length(ls_shRNA_res))){
  tmp = get(ls_shRNA_res[i])
  df_SL_predict_algoshRNA = cbind(df_SL_predict_algoshRNA, tmp[,c(4, 5)])
  ncol_ids = ncol(df_SL_predict_algoshRNA)
  ncol_ids = c((ncol_ids -1), ncol_ids)
  colnames(df_SL_predict_algoshRNA)[ncol_ids] = colnames(tmp)[c(4, 5)]
}

df_SL_predict_algoshRNA = df_SL_predict_algoshRNA[!is.na(tmp$overlp_ngMCS),]

write.csv(df_SL_predict_algoshRNA, '../MetabolicSLOutput/reviewer comments/shRNA_analysis_comapre_algorithms.csv')

# save|load env ####
# save.image('./env/WS_MN_exh_ttest_cor_cons_ataris.rdata')
# load('./env/WS_MN_exh_ttest_cor_cons_ataris.rdata')