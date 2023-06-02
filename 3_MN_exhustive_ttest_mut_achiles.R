# adding drive attaris and prism for this analysis
# ko gene of achiles ataris and prism have to be limited to metaboic genes but for partnering genes in an additional version of data
# driven and even for our algorithm we can search for all genes with mutation in more than 2 cell lines in each cancer types and not
# only the metabolic ones, in addition for better results the mcmc sampling approach must be done

rm(list = ls())

library("dplyr")

# loading data
gen_ls = read.csv('../MetabolicSLInput/data/tb_gene_ls.csv')
# KO_FBA_n = list.files('../MetabolicSLOutput/can_KO/')
# KO_FBA = NULL
# for ( filenames in KO_FBA_n){
#   KO_FBA = cbind(KO_FBA, read.csv(paste0('../MetabolicSLOutput/can_KO/',filenames)))
# }; rm(KO_FBA_n)
# KO_FBA = read.csv('../MetabolicSLOutput/can_KO/ko_geninue_biomass__Bile Duct Cancer.csv'); rownames(KO_FBA) = gen_ls$Var1;
KO_FBA = read.csv('../MetabolicSLOutput/KO_res_all.csv'); rownames(KO_FBA) = gen_ls$Var1;
KO_FBA = data.frame(t(KO_FBA))
rownames(KO_FBA) = gsub('_', '-',rownames(KO_FBA))
cantyp = read.csv('../MetabolicSLInput/data/cantyp.csv')
maf_met = read.csv('../MetabolicSLInput/data/maf_met_ge.csv')
maf_met$DepMap_ID = gsub('_', '-', maf_met$DepMap_ID)
sampleInfo <- read.csv("../MetabolicSLInput/data/sample_info.csv", stringsAsFactors = FALSE, header = TRUE)
load('../MetabolicSLInput/data/Achilles_gene_effect.RData')
# save(achiles, file = '../MetabolicSLInput/data/Achilles_gene_effect.RData')
rownames(achiles) = achiles$DepMap_ID; achiles = achiles[,-match('DepMap_ID', colnames(achiles))]
# correct achiles colnames, rm .. and after that from colnames
for (i in 1:ncol(achiles)){
  colnames(achiles)[i] = unlist(strsplit(colnames(achiles)[i],'[.][.]'))[1]
}; rm(i)

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
maf_met = maf_met[maf_met$DepMap_ID %in% rownames(KO_FBA),]

#### ttest gold standard SL prediction #### 
source('./code/func__tans_dom__get_sym.R')
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
    write.csv(can_tbl, file = paste0('../MetabolicSLOutput/', can,'_ac_mn_w&t_tst_cons.csv'))
    can_tbl 
  }
  # print(paste0(round(can_i / nrow(cantyp) * 100, 1), ' % current cancer'))
}

write.csv(fin_can, file = '../MetabolicSLOutput/fin_can_ac_mn_w&ttst_cons.csv')
# fin_can = read.csv('../MetabolicSLOutput/fin_can_ac_mn_w&ttst_cons.csv')
fin_can = data.frame(fin_can)

# read the file of each cancer adjust them and agian write them in the same place
fname = list.files(pattern ='ac_mn_w&t_tst_cons.csv', '../MetabolicSLOutput/')
#### fdr ####
sl_scan_fdr = NULL
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
    colnames(tmp) = c("cancer","Mutated_gene","KO_gene","No_mut","No_not_mut","ttst_achiles_stat", "ttst_achiles_pvl"
                     # ,"wtst_achiles_stat", "wtst_achiles_pvl"
                     , "ttst_metabolic_stat", "ttst_metabolic_pvl"
                     # , "wtst_metabolic_stat", "wtst_pre_pvl"
                     )
    tmp = cbind(tmp
          , ttst_achiles_fdr = p.adjust(as.numeric(tmp$ttst_achiles_pvl), method = 'fdr')
          , ttst_metabolic_fdr = p.adjust(as.numeric(tmp$ttst_metabolic_pvl), method = 'fdr'))
    # tmp = cbind(tmp
    #             , ttst_achiles_fdr = (as.numeric(tmp$ttst_achiles_pvl))
    #             , ttst_metabolic_fdr = (as.numeric(tmp$ttst_metabolic_pvl)))
    name = gsub(".csv","",fname[f_i])
    write.csv(tmp, paste0("../MetabolicSLOutput/",name,"_fdr.csv"))
    sl_scan_fdr = rbind(sl_scan_fdr, tmp)
  }
}; #rm(tmp)
# adj_final_all_ac = subset(final_all_ac, final_all_ac$fdr_tt <= 0.05)
# View(sl_scan_fdr[intersect(match(SL_drug$Mutated_gene,sl_scan_fdr$Mutated_gene), match(SL_drug$KO_gene, sl_scan_fdr$KO_gene)),])

#### hypergeometric test SL-scan #### 
fname = list.files(pattern ='_ac_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp = data.frame(matrix(NA, nrow = length(fname), ncol = 6))
colnames(cantyp_hyp) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_ach', 'overlp', 'total')
SL_scan = NULL
SL_scan_valid = NULL
for (cn in 1:length(fname)){
  print(cn/length(fname))
  cantyp_hyp$cancer[cn] = unlist(strsplit(fname[cn], '_ac_mn_w&t_tst_cons_fdr.csv')[1])
  tmp = read.csv(paste0('../MetabolicSLOutput/',fname[cn]))
  mn = which(as.numeric(tmp$ttst_metabolic_fdr) <= 0.05)
  ach = which(as.numeric(tmp$ttst_achiles_fdr) <= 0.05)
  SL_scan = rbind(SL_scan , tmp[(mn),])
  SL_scan_valid = rbind(SL_scan_valid , tmp[intersect(mn, ach),])
  cantyp_hyp$overlp[cn] = length(intersect(mn, ach))
  cantyp_hyp$no_SL_ach[cn] = length(ach)
  cantyp_hyp$no_SL_MN[cn] = length(mn)
  # cantyp_hyp$total[cn] = ncol(KO_FBA)*(ncol(KO_FBA)-1)/2
  cantyp_hyp$total[cn] = nrow(tmp)
  
  cantyp_hyp$hypergeom[cn] = phyper(cantyp_hyp$overlp[cn]-1, cantyp_hyp$no_SL_MN[cn]
                                    , cantyp_hyp$total[cn]-cantyp_hyp$no_SL_MN[cn]
                                    , cantyp_hyp$no_SL_ach[cn], lower.tail = FALSE)
}
freq_slscan_driver = SL_scan %>% 
  group_by(SL_scan$Mutated_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

freq_slscan_KO = SL_scan %>% 
  group_by(SL_scan$KO_gene) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

write.csv(SL_scan_valid, "../MetabolicSLOutput/Pre_SL_SLScan_valid.csv")
write.csv(SL_scan, "../MetabolicSLOutput/Pre_SL_SLScan_all.csv")
table_top_SL_can_valid = can_SL(SL_scan_valid, 2, 4)
write.csv(table_top_SL_can_valid, "../MetabolicSLOutput/table_top_SL_can_valid_cons.csv")
# cortest_ac_slscan = cor.test(cantyp_hyp$no_SL_MN, cantyp_hyp$no_SL_ach, method = "pearson")

#### hypergeometric test fastSL results #### 
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_fastSL' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_fastSL = read.csv("../MetabolicSLOutput/R24_iMAT_cons_fastSL.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_fastSL),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_fastSL[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_ac_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_fsl = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_fsl) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_ach', 'overlp', 'total', 'no_SLscan')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_ac_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  col1 = entrz2hgnc(paste0(R24_iMAT_fastSL[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_fastSL[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  id = grep(unlist(strsplit(fname[i],'Cancer'))[1], fname)
  
  if (length(id) != 0){
    tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[id]))
    tmp_driver = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
    tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    
    cantyp_hyp_fsl$cancer[i] = can
    cantyp_hyp_fsl$total[i] = 0
    if (nrow(total)>0){
      chek_concordance = cbind(total, achilles = NA, MN = NA, slscan = NA)
      for (chec_i in 1:nrow(chek_concordance)){
        chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                              ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
        if (length(chck_prsm) != 0 ){
          if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
          if (any(tmp_ac$ttst_achiles_fdr[chck_prsm] <= 0.05)){chek_concordance$achilles[chec_i] = "YES"} else {chek_concordance$achilles[chec_i] = "NO"}
        }
        
        chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
        chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
        if (length(chck1) != 0 | length(chck2) != 0){
          chek_concordance$MN[chec_i] = "YES"
          # print("YES ............................")
        } 
      }
      
      chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
    
        # hyper_geom test
      cantyp_hyp_fsl$overlp[i] = length(intersect(which(chek_concordance$MN == "YES")
                                                  , which(chek_concordance$achilles == "YES")))
      cantyp_hyp_fsl$no_SL_ach[i] = length(which(chek_concordance$achilles == "YES"))
      cantyp_hyp_fsl$no_SL_MN[i] = length(which(chek_concordance$MN == "YES"))
      cantyp_hyp_fsl$total[i] = nrow(chek_concordance)
      cantyp_hyp_fsl$no_SLscan[i] = length(intersect(which(chek_concordance$slscan == "YES")
                                                     , which(chek_concordance$achilles == "YES")))
      
      cantyp_hyp_fsl$hypergeom[i] = phyper(cantyp_hyp_fsl$overlp[i]-1, cantyp_hyp_fsl$no_SL_MN[i]
                                           , cantyp_hyp_fsl$total[i]-cantyp_hyp_fsl$no_SL_MN[i]
                                           , cantyp_hyp_fsl$no_SL_ach[i], lower.tail = FALSE)
    }  
  # rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  }
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_achilles_FastSL_concordance_check.csv'))
# cortest_ac_fsl = cor.test(cantyp_hyp_fsl$no_SL_MN, cantyp_hyp_fsl$no_SL_ach, method = "pearson")

#### hypergeometric test MCS results ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_mcs' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_mcs = read.csv("../MetabolicSLOutput/R24_iMAT_cons_mcs.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_mcs),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_mcs[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_ac_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_mcs = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_mcs) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_ach', 'overlp', 'total', 'no_SLscan')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_ac_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_mcs[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_mcs[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  id = grep(unlist(strsplit(fname[i],'Cancer'))[1], fname)
  
  if (length(id) != 0){
    tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[id]))
    tmp_driver = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
    tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    
    cantyp_hyp_mcs$cancer[i] = can
    cantyp_hyp_mcs$total[i] = 0
    if (nrow(total)>0){
      chek_concordance = cbind(total, achilles = NA, MN = NA, slscan = NA)
      for (chec_i in 1:nrow(chek_concordance)){
        chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                              ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
        if (length(chck_prsm) != 0 ){
          if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
          if (any(tmp_ac$ttst_achiles_fdr[chck_prsm] <= 0.05)){chek_concordance$achilles[chec_i] = "YES"} else {chek_concordance$achilles[chec_i] = "NO"}
        }
        
        chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
        chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
        if (length(chck1) != 0 | length(chck2) != 0){
          chek_concordance$MN[chec_i] = "YES"
          # print("YES ............................")
        } 
      }
      
      chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
      
      # hyper_geom test
      cantyp_hyp_mcs$overlp[i] = length(intersect(which(chek_concordance$MN == "YES")
                                                  , which(chek_concordance$achilles == "YES")))
      cantyp_hyp_mcs$no_SL_ach[i] = length(which(chek_concordance$achilles == "YES"))
      cantyp_hyp_mcs$no_SL_MN[i] = length(which(chek_concordance$MN == "YES"))
      cantyp_hyp_mcs$total[i] = nrow(chek_concordance)
      cantyp_hyp_mcs$no_SLscan[i] = length(intersect(which(chek_concordance$slscan == "YES")
                                                     , which(chek_concordance$achilles == "YES")))
      
      cantyp_hyp_mcs$hypergeom[i] = phyper(cantyp_hyp_mcs$overlp[i]-1, cantyp_hyp_mcs$no_SL_MN[i]
                                           , cantyp_hyp_mcs$total[i]-cantyp_hyp_mcs$no_SL_MN[i]
                                           , cantyp_hyp_mcs$no_SL_ach[i], lower.tail = FALSE)
    }  
    # rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  }
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_achilles_MCS_concordance_check.csv'))

#### hypergeometric test gMCS results ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_gmcs' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_gmcs = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gmcs.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_gmcs),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_gmcs[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_ac_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_gmcs = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_gmcs) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_ach', 'overlp', 'total', 'no_SLscan')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i / length(fname) * 100)
  can = unlist(strsplit(fname[i], '_ac_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,(endd/2)]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(colnames(achiles), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  col1 = entrz2hgnc(paste0(R24_iMAT_gmcs[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_gmcs[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  id = grep(unlist(strsplit(fname[i],'Cancer'))[1], fname)
  
  if (length(id) != 0){
    tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[id]))
    tmp_driver = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
    tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    
    cantyp_hyp_gmcs$cancer[i] = can
    cantyp_hyp_gmcs$total[i] = 0
    if (nrow(total)>0){
      chek_concordance = cbind(total, achilles = NA, MN = NA, slscan = NA)
      for (chec_i in 1:nrow(chek_concordance)){
        chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                              ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
        if (length(chck_prsm) != 0 ){
          if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
          if (any(tmp_ac$ttst_achiles_fdr[chck_prsm] <= 0.05)){chek_concordance$achilles[chec_i] = "YES"} else {chek_concordance$achilles[chec_i] = "NO"}
        }
        
        chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
        chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
        if (length(chck1) != 0 | length(chck2) != 0){
          chek_concordance$MN[chec_i] = "YES"
          # print("YES ............................")
        } 
      }
      
      chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
      
      # hyper_geom test
      cantyp_hyp_gmcs$overlp[i] = length(intersect(which(chek_concordance$MN == "YES")
                                                   , which(chek_concordance$achilles == "YES")))
      cantyp_hyp_gmcs$no_SL_ach[i] = length(which(chek_concordance$achilles == "YES"))
      cantyp_hyp_gmcs$no_SL_MN[i] = length(which(chek_concordance$MN == "YES"))
      cantyp_hyp_gmcs$total[i] = nrow(chek_concordance)
      cantyp_hyp_gmcs$no_SLscan[i] = length((which(chek_concordance$slscan == "YES")))
      
      cantyp_hyp_gmcs$hypergeom[i] = phyper(cantyp_hyp_gmcs$overlp[i]-1, cantyp_hyp_gmcs$no_SL_MN[i]
                                            , cantyp_hyp_gmcs$total[i]-cantyp_hyp_gmcs$no_SL_MN[i]
                                            , cantyp_hyp_gmcs$no_SL_ach[i], lower.tail = FALSE)
    }  
    # rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  }
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_achilles_gMCS_concordance_check.csv'))

# cortest_ac_gmcs = cor.test(cantyp_hyp_gmcs$no_SL_MN, cantyp_hyp_gmcs$no_SL_ach, method = "pearson")
#### cancers model dimensions ####
# digit = 0
# for(i in 1:(ncol(R24_iMAT_mcs)/2)){
#   print(i)
#   endd = i * 2; srt = endd - 1
#   col1 = entrz2hgnc(paste0(R24_iMAT_mcs[,srt],'.1'), gen_ls)
#   col2 = entrz2hgnc(paste0(R24_iMAT_mcs[,endd],'.1'), gen_ls)
#   sl = data.frame(cbind(col1,col2))
#   print(colnames(R24_iMAT_mcs)[srt])
#   print(max(which(!is.na(col1))))
#   digit = digit +  max(which(!is.na(col1)))
# }
#### preparing and writing all hypergeometric test results ####


for (i in 1:length(ls(pattern = '_hyp'))){
  tmpn = ls(pattern = '_hyp')[i]  
  tmp = get(tmpn)
  tmp$hypergeom = formatC(tmp$hypergeom, format = "e", digits = 2)
  write.csv(tmp, paste0('../MetabolicSLOutput/cons_achiles_', tmpn,'.csv'))
}

# stackbarplot of algorithms of SL prediction ####
rm(stack_plotHyp)
namess = ls(pattern = "_hyp")
stack_plotHyp = NULL
for (i in 1:length(namess)){
  tmp_name = namess[i]
  tmp = get(tmp_name)
  if (i == 1 ){ stack_plotHyp = data.frame(cbind(cancer = tmp$cancer, SL_scan = tmp$overlp))
  } else {
    ids = match(stack_plotHyp$cancer, tmp$cancer)
    tmp = tmp$overlp[ids]; tmp[is.na(tmp)] = 0
    stack_plotHyp = cbind(stack_plotHyp, tmp)
    colnames(stack_plotHyp)[i+1] = namess[i]
  }
}
colnames(stack_plotHyp) = c("Cancers", "SL_scan", "FastSL", "gMCS", "MCS")
stack_plotHyp = stack_plotHyp[order(-as.numeric(stack_plotHyp$SL_scan)),]
colnames(stack_plotHyp)[colnames(stack_plotHyp) == "SL_scan"] = "SL-scan"
# stack_plotHyp$Cancers <- reorder(stack_plotHyp$Cancers, stack_plotHyp$SL_scan)

library(ggplot2)
longformat = reshape2::melt(stack_plotHyp, id.vars = "Cancers", variable.name = "Approaches", value.name = "Log_valid_SL_pairs")
longformat$Log_valid_SL_pairs[is.na(longformat$Log_valid_SL_pairs)] = 0
longformat$Log_valid_SL_pairs = log(as.numeric(longformat$Log_valid_SL_pairs))

cancer_order = unique(longformat$Cancers)
longformat$Cancers <- factor(longformat$Cancers, levels = cancer_order)

png(paste0("../MetabolicSLOutput/dodge_barplot_valid_SL_no.png"),width=1000, height=700)
ggplot(longformat, aes(fill = Approaches, y = exp(Log_valid_SL_pairs), x = Cancers)) +
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle("Grouped bar plot of SL predicted per cancer group using different approaches") +
  scale_y_log10() +
  labs(x = NULL, y = "Number of valid predicted SL pairs") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        plot.margin = unit(c(1, 1, 1, 2), "cm")) +
  scale_fill_manual(values = c("red2", "blue2", "green2", "yellow")) + # custom colors
  geom_text(aes(label = exp(Log_valid_SL_pairs)), position = position_dodge(width = 1), vjust = -0.5)
dev.off()


# save|load env
# save.image('./env/WS_MN_exh_ttest_cor_cons_achiles.rdata')
# load('./env/WS_MN_exh_ttest_cor_cons_achiles.rdata')

