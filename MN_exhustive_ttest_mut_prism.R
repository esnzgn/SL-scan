# adding drive attaris and prism for this analysis
# ko gene of achiles ataris and prism have to be limited to metaboic genes but for partnering genes in an additional version of data
# driven and even for our algorithm we can search for all genes with mutation in more than 2 cell lines in each cancer types and not
# only the metabolic ones, in addition for better results the mcmc sampling approach must be done

# importnat note is that I  have to consider those cell lines which were not constructed using iMAT, and all elemetns of their KO_FBA
# file vector elements are nan so they will translate to 1 when they should be eleminated

rm(list = ls())

library("dplyr")

# loading data
gen_ls = read.csv('../MetabolicSLInput/data/tb_gene_ls.csv')
KO_FBA = read.csv('../MetabolicSLOutput/KO_res_446_ttestV.csv'); rownames(KO_FBA) = gen_ls$Var1;
KO_FBA = data.frame(t(KO_FBA))
rownames(KO_FBA) = gsub('_', '-',rownames(KO_FBA))
cantyp = read.csv('../MetabolicSLInput/data/cantyp.csv')
maf_met = read.csv('../MetabolicSLInput/data/maf_met_ge.csv')
maf_met$DepMap_ID = gsub('_', '-', maf_met$DepMap_ID)
sampleInfo <- read.csv("../MetabolicSLInput/data/sample_info.csv", stringsAsFactors = FALSE, header = TRUE)
# achiles = readRDS('../MetabolicSLInput/data/Achilles_gene_dependency.rds')
# achiles = t(readRDS('../MetabolicSLInput/data/DRIVE_ATARiS_data.rds'))

# 
prism_info = read.delim('../MetabolicSLInput/data/primary-screen-replicate-collapsed-treatment-info.csv', stringsAsFactors = FALSE, header = TRUE,sep = ',')
# prism <- read.delim("../MetabolicSLInput/data/Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4.csv", stringsAsFactors = FALSE, header = TRUE,sep = ',')
prism <- read.delim("../MetabolicSLInput/data/primary-screen-replicate-collapsed-logfold-change.csv", stringsAsFactors = FALSE, header = TRUE,sep = ',')

prism_col = NULL
for (i in colnames(prism)){ prism_col = c(prism_col,unlist(strsplit(i,'..BRD.BRD.',fixed = T))[1])}

# prism has more than 1 target gene for each drug so I need to map drugs to each gene for selection of metabolic genes 
colnames(prism) = prism_col
prism_mapping = NULL
for (ps in 2:length(prism_col)){
  # print((ps)/length(prism_col)*100)
  (prism_col[ps])
  pattern = paste(unlist(strsplit(prism_col[ps], split = '\\.'))[1:5], collapse = "-")
  id_prs = which(toupper(pattern) == toupper(prism_info$broad_id))
  for (id_len in 1:length(id_prs)){
    maped_gen = unlist(strsplit(prism_info$target[id_prs[id_len]], ','))
    if (!is.null(maped_gen[1]) && !is.na(maped_gen[1])){
      if (length(maped_gen) == 1){
        add = cbind(ps, prism_info$name[id_prs[id_len]], prism_info$target[id_prs[id_len]], prism_info$dose[id_prs[id_len]])
        prism_mapping = rbind(prism_mapping, add)  
      } else if (length(maped_gen) > 1){
          len = length(maped_gen)
        for (ele_i in 1:len){
          ele = maped_gen[ele_i]
          add = cbind(ps, prism_info$name[id_prs[id_len]], ele, prism_info$dose[id_prs[id_len]])
          prism_mapping = rbind(prism_mapping, add) 
        }
      }
    }
  }
}
prism_mapping = data.frame(prism_mapping)
colnames(prism_mapping) = c('drug_prism_id', 'drug_info', 'target_info', 'dose')
# find prism metabolic gens for subseting prism based on it then subset prism data set of primary treatmet and extend it for target genes
prism_mapping_met = prism_mapping[prism_mapping$target_info %in% gen_ls$Var1, ]
# make all drug-gene combinations 
prism_mapping_met = prism_mapping_met[order(prism_mapping_met$target_info), ]
prism_met =  prism[, as.numeric(prism_mapping_met$drug_prism_id)]
colnames(prism_met) = paste(prism_mapping_met$drug_info, prism_mapping_met$target_info, sep = '==')
prism_mapping_met = cbind(prism_mapping_met, Prism_colnames = colnames(prism_met))
rownames(prism_met) = prism$X

# cell lines
intr_cel = intersect(rownames(prism_met), rownames(KO_FBA))
prism_met = prism_met[rownames(prism_met) %in% intr_cel,] 
prism_met = prism_met[order(rownames(prism_met)),]
KO_FBA = KO_FBA[rownames(KO_FBA) %in% intr_cel,]; rm(intr_cel)
KO_FBA = KO_FBA[order(rownames(KO_FBA)),]
# make map of cell lines names and diseases by a two col df
pd_cel = data.frame(cbind(primary_disease = sampleInfo$primary_disease[match(rownames(prism_met)
                                     , sampleInfo$DepMap_ID)], DepMap_ID = rownames(prism_met)))
write.csv(prism_mapping_met, '../MetabolicSLOutput/prism_mapping_met.csv')
maf_met = maf_met[maf_met$DepMap_ID %in% rownames(KO_FBA),]
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
nrCores =  parallel::detectCores() - 4
# nrCores = nrCores[1]-1 #not to overload your computer
# nrCores = 4

if(Sys.info()['sysname'] != "Windows") {
  # install.packages("doMC")
  library(doMC)
  registerDoMC(nrCores)
  getDoParWorkers()
} else {
  # install.packages('doParallel')
  library(doParallel)
  cl <- makeCluster(nrCores)
  registerDoParallel(cl)
}
# parallel PRISM ####
fin_can = NULL
fin_can = foreach(can_i= 1:(length(ac_cantyp)), .combine=rbind) %dopar%{
  # for (can_i in c(6,7,9,12,14)){
  # print(rep('###', 10));print(cantyp$X1[can_i])
  can = names(ac_cantyp[can_i]);
  can_cel = ac_cantyp[[can]]
  if (length(can_cel) >= 4){
    # tmp_mut = mutbck_ls[, srt:end]
    tmp_mut = mutbackmkr(unique(gen_ls$Var1), maf_met, can_cel)
    # head(rowSums(tmp_mut)[order(rowSums(tmp_mut), decreasing = T)])
    # slice achiles for the selected cell lines of the cancer
    id_cel = (pd_cel$DepMap_ID %in% can_cel)
    KO_PR = prism_met[id_cel, ]
    KO_mn = KO_FBA[id_cel, ]
    # find nan elements in fba ko matrix and put 1 instead of them
    for (mn_i in 1:ncol(KO_mn)){
      na_id = is.na(KO_mn[,mn_i])
      KO_mn[na_id ,mn_i] = 1
    }
    rm(id_cel)
    can_tbl = NULL;
    # find each cancer A and B genes
    id_GnA = which(rowSums(tmp_mut) >= thre)
    
    if (length(id_GnA) > 0){
      GnA = names(id_GnA)
      GnB = unique(colnames(KO_PR))
      
      for(i in 1:length(GnA)){
        # print(i/length(GnA)*100)
        # celselector function will find each gene mut and not_mut cell liens for statistical test in next step
        res_cel = celselector(GnA[i], tmp_mut)
        for (j in 1:length(GnB)){
          # ttest
          t_res = NULL
          t_res = tw_tst_prism(res_cel, GnB[j], KO_PR, KO_mn, altn, thre)
          if (!is.na(t_res[1])){
            if (!is.na(t_res$ttst_obs$p.value)){
              can_tbl = rbind(can_tbl ,  cbind('cancer' = can, 'Mutated_gene' = GnA[i], 'KO_gene' = unlist(strsplit(GnB[j], '=='))[2]
                                               , 'KO_gene_drug' =  GnB[j], 'No_mut' = length(res_cel$mut), 'No_not_mut' = length(res_cel$not_mut)
                                               , 'ttst_obs_stat' =  t_res$ttst_obs$statistic, 'ttst_obs_pvl' = t_res$ttst_obs$p.value
                                               # , 'wtst_obs_stat' =  t_res$wtst_obs$statistic, 'wtst_obs_pvl' = t_res$wtst_obs$p.value
                                               , 'ttst_pre_stat' =  t_res$ttst_pre$statistic, 'ttst_pre_pvl' = t_res$ttst_pre$p.value
                                               # , 'wtst_pre_stat' =  t_res$wtst_pre$statistic, 'wtst_pre_pvl' = t_res$wtst_pre$p.value
              )
            )
          } else {t_res = NULL}
        }
      }
    }
    write.csv(can_tbl, file = paste0('../MetabolicSLOutput/', can,'_PRISM_mn_w&t_tst_cons.csv'))
    can_tbl
  }
  # print(paste0(round(can_i / nrow(cantyp) * 100, 1), ' % current cancer'))
  }
}
# doMC::stopImplicitCluster()
stopCluster(cl)
registerDoSEQ()


write.csv(fin_can, file = '../MetabolicSLOutput/fin_can_at_mn_w&ttst_cons_prism.csv')
fin_can = data.frame(fin_can)
# read the file of each cancer adjust them and agian write them in the same place
fname = list.files(pattern ='_PRISM_mn_w&t_tst_cons.csv', '../MetabolicSLOutput/')
slscan_prism_fdr = NULL
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
    colnames(tmp) = c("cancer","Mutated_gene","KO_gene", 'KO_gene_drug',"No_mut","No_not_mut","ttst_PRISM_stat", "ttst_PRISM_pvl"
                      # ,"wtst_PRISM_stat", "wtst_PRISM_pvl", 
                      ,"ttst_metabolic_stat", "ttst_metabolic_pvl"
                      # , "wtst_metabolic_stat", "wtst_pre_pvl"
                      )
    tmp = cbind(tmp
                , ttst_PRISM_fdr = p.adjust(as.numeric(tmp$ttst_PRISM_pvl), method = 'fdr')
                , ttst_metabolic_fdr = p.adjust(as.numeric(tmp$ttst_metabolic_pvl), method = 'fdr'))
    # tmp = cbind(tmp
    #             , ttst_achiles_fdr = (as.numeric(tmp$ttst_achiles_pvl))
    #             , ttst_metabolic_fdr = (as.numeric(tmp$ttst_metabolic_pvl)))
    name = unlist(strsplit(fname[f_i], split = '_PRISM_mn_w&t_tst_cons.csv'))[1]
    write.csv(tmp, paste0('../MetabolicSLOutput/',name, '_PRISM_mn_w&t_tst_cons_fdr.csv'))
    slscan_prism_fdr = rbind(slscan_prism_fdr, tmp)
  }
}; #rm(tmp)
# adj_final_all_ac = subset(final_all_ac, final_all_ac$fdr_tt <= 0.05)

#### SL-scan hypergeom ####
gene_ls_PRISM = list()
fname = list.files(pattern ='_PRISM_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp = data.frame(matrix(NA, nrow = length(fname), ncol = 6))
colnames(cantyp_hyp) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_prism', 'overlp', 'total')
SL_drug = NULL
for (cn in 1:length(fname)){
  print(cn/length(fname)*100)
  cantyp_hyp$cancer[cn] = unlist(strsplit(fname[cn], '_PRISM_mn_w&t_tst_cons_fdr.csv')[1])
  tmp = read.csv(paste0('../MetabolicSLOutput/',fname[cn]))
  driver = unique(tmp$Mutated_gene); KO_gene = unique(tmp$KO_gene)
  gene_ls_PRISM[[cantyp_hyp$cancer[cn]]] = list(driver = driver, KO_gene = KO_gene)
  mn = which(as.numeric(tmp$ttst_metabolic_fdr) <= 0.05)
  prsm = which(as.numeric(tmp$ttst_PRISM_pvl) <= 0.05)
  SL_drug = rbind(SL_drug, tmp[intersect(mn, prsm),])
  cantyp_hyp$overlp[cn] = length(intersect(mn, prsm))
  # View(tmp[intersect(mn, prsm),])
  cantyp_hyp$no_SL_prism[cn] = length(prsm)
  cantyp_hyp$no_SL_MN[cn] = length(mn)
  cantyp_hyp$total[cn] = nrow(tmp)
  
  cantyp_hyp$hypergeom[cn] = phyper(cantyp_hyp$overlp[cn]-1, cantyp_hyp$no_SL_MN[cn]
                                    , cantyp_hyp$total[cn]-cantyp_hyp$no_SL_MN[cn]
                                    , cantyp_hyp$no_SL_prism[cn], lower.tail = FALSE)
}
unique_drg = unlist(strsplit(SL_drug$KO_gene_drug, "=="))
drug_ls_repur = NULL
for (i in seq(1,length(unique_drg),2)){
  drug_ls_repur = c(drug_ls_repur, unique_drg[i] )
}; drug_ls_repur = unique(drug_ls_repur)

SL_drug = (cbind(SL_drug[,-ncol(SL_drug)], ttst_metabolic_fdr = as.numeric(SL_drug$ttst_metabolic_fdr)))
write.csv(cantyp_hyp, '../MetabolicSLOutput/hyper_table_fdrcorrect_cons_PRISM.csv')
SL_drug$ttst_metabolic_fdr = formatC(SL_drug$ttst_metabolic_fdr, format = "e", 3); SL_drug$ttst_PRISM_pvl =  formatC(SL_drug$ttst_PRISM_pvl, format = "e", 3)
write.csv(SL_drug, '../MetabolicSLOutput/SLscan_drug_fdrcorrect_cons_PRISM.csv')

#### fastSL hypergeom ####
# library(tibble)
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_fastSL' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_fastSL = read.csv("../MetabolicSLOutput/R24_iMAT_cons_fastSL.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_fastSL),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_fastSL[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_PRISM_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_fsl = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_fsl) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_psrm', 'no_SLscan','overlp', 'total')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i)
  can = unlist(strsplit(fname[i], '_PRISM_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  tmp_g = gen_ls$Var1[gen_ls$Var2 %in% R24_iMAT_gls[,i]]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
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
    tmp_drg = unique(prism_mapping_met$drug_info[prism_mapping_met$target_info %in% tmp_ko])
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    # tmp_ac = length(unique(tmp_ac$KO_gene_drug)) *length(unique(tmp_ac$Mutated_gene))
    cantyp_hyp_fsl$cancer[i] = can
    cantyp_hyp_fsl$total[i] = 0
    if (nrow(total)>0){
          chek_concordance = cbind(total, prism = NA, MN = NA, slscan = NA)
    for (chec_i in 1:nrow(chek_concordance)){
      chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                            ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
      if (length(chck_prsm) != 0 ){
        if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
        if (any(tmp_ac$ttst_PRISM_pvl[chck_prsm] <= 0.05)){chek_concordance$prism[chec_i] = "YES"} else {chek_concordance$prism[chec_i] = "NO"}
      }
      
      
            chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
            chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
            if (length(chck1) != 0 | length(chck2) != 0){
              chek_concordance$MN[chec_i] = "YES"
              print("YES ............................")
          } 
      }

    chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
    # prism_sig = intersect((tmp_ac$Mutated_gene %in% tmp_driver),  (tmp_ac$KO_gene_drug, "==") %in% tmp_drg) 
    # hyper_geom test
    
    cantyp_hyp_fsl$overlp[i] = length(which(chek_concordance$MN == "YES"))
    cantyp_hyp_fsl$no_SL_psrm[i] = length(which(chek_concordance$prism == "YES"))
    cantyp_hyp_fsl$no_SLscan[i] = length(which(chek_concordance$slscan == "YES"))
    cantyp_hyp_fsl$no_SL_MN[i] = length(which(!is.na(chek_concordance$MN)))
    cantyp_hyp_fsl$total[i] = nrow(chek_concordance)
    
    cantyp_hyp_fsl$hypergeom[i] = phyper(cantyp_hyp_fsl$overlp[i]-1, cantyp_hyp_fsl$no_SL_MN[i]
                                      , cantyp_hyp_fsl$total[i]-cantyp_hyp_fsl$no_SL_MN[i]
                                      , cantyp_hyp_fsl$no_SL_psrm[i], lower.tail = FALSE)
    }

  }
  
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_prism_FastSL_concordance_check.csv'))

#### MCS hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_mcs.csv' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_mcs = read.csv("../MetabolicSLOutput/R24_iMAT_cons_mcs.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_mcs),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_mcs[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_PRISM_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_mcs = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_mcs) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_psrm', 'no_SLscan','overlp', 'total')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i)
  can = unlist(strsplit(fname[i], '_PRISM_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  tmp_g = gen_ls$Var1[match(R24_iMAT_gls[,i], gen_ls$Var2)]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(unique(prism_mapping_met$target_info), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  col1 = entrz2hgnc(paste0(R24_iMAT_mcs[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_mcs[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  id = grep(unlist(strsplit(fname[i],'Cancer'))[1], fname)
  
  if (length(id) != 0){
    tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[id]))
    tmp_driver = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
    tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
    tmp_drg = unique(prism_mapping_met$drug_info[prism_mapping_met$target_info %in% tmp_ko])
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    # tmp_ac = length(unique(tmp_ac$KO_gene_drug)) *length(unique(tmp_ac$Mutated_gene))
    # tmp_ac = length(unique(tmp_ac$KO_gene_drug)) *length(unique(tmp_ac$Mutated_gene))
    cantyp_hyp_mcs$cancer[i] = can
    cantyp_hyp_mcs$total[i] = 0
    if (nrow(total)>0){
      chek_concordance = cbind(total, prism = NA, MN = NA, slscan = NA)
      for (chec_i in 1:nrow(chek_concordance)){
        chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                              ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
        if (length(chck_prsm) != 0 ){
          if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
          if (any(tmp_ac$ttst_PRISM_pvl[chck_prsm] <= 0.05)){chek_concordance$prism[chec_i] = "YES"} else {chek_concordance$prism[chec_i] = "NO"}
        }
        
        
        chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
        chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
        if (length(chck1) != 0 | length(chck2) != 0){
          chek_concordance$MN[chec_i] = "YES"
          print("YES ............................")
        } 
      }
      
      chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
    # hyper_geom test
      cantyp_hyp_mcs$overlp[i] = length(which(chek_concordance$MN == "YES"))
      cantyp_hyp_mcs$no_SL_psrm[i] = length(which(chek_concordance$prism == "YES"))
      cantyp_hyp_mcs$no_SLscan[i] = length(which(chek_concordance$slscan == "YES"))
      cantyp_hyp_mcs$no_SL_MN[i] = length(which(!is.na(chek_concordance$MN)))
      cantyp_hyp_mcs$total[i] = nrow(chek_concordance)
      
      cantyp_hyp_mcs$hypergeom[i] = phyper(cantyp_hyp_fsl$overlp[i]-1, cantyp_hyp_fsl$no_SL_MN[i]
                                           , cantyp_hyp_fsl$total[i]-cantyp_hyp_fsl$no_SL_MN[i]
                                           , cantyp_hyp_fsl$no_SL_psrm[i], lower.tail = FALSE)
    
  }
  }
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_prism_MCS_concordance_check.csv'))

#### gmcs hypergeom ####
gen_ls = read.csv('../MetabolicSLOutput/tb_gene_ls.csv')
list.files(pattern = 'R24_iMAT_cons_gmcs.csv' ,'../MetabolicSLOutput/', ignore.case = T)
R24_iMAT_gls = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gls.csv")
R24_iMAT_gmcs = read.csv("../MetabolicSLOutput/R24_iMAT_cons_gmcs.csv")
can_ls = NULL
for (i in seq(1,ncol(R24_iMAT_gmcs),2)){ can_ls[i] = (unlist(strsplit(colnames(R24_iMAT_gmcs[i]), "[.|_]"))[1])}
fname = list.files(pattern ='_PRISM_mn_w&t_tst_cons_fdr.csv', '../MetabolicSLOutput/')
cantyp_hyp_gmcs = data.frame(matrix(NA, nrow = length(fname), ncol = 7))
colnames(cantyp_hyp_gmcs) = c('cancer', 'hypergeom', 'no_SL_MN', 'no_SL_psrm', 'no_SLscan','overlp', 'total')
chek_concordance_ttl = NULL
for(i in 1:length(fname)){
  print(i)
  can = unlist(strsplit(fname[i], '_PRISM_mn_w&t_tst_cons_fdr.csv')[1])
  rm(thr_fdr, overlap, tmp_ac, id, sl, col1, col2, srt , endd, total, tmp_g, tmp_gls)
  tmp_g = gen_ls$Var1[match(R24_iMAT_gls[,i], gen_ls$Var2)]
  tmp_gls = unique(tmp_g[!is.na(tmp_g)])
  # total = length(intersect(unique(prism_mapping_met$target_info), tmp_gls))
  # total = ((total * (total - 1))/2)
  # sl comparison of fastsl algo with achiles results
  srt = grep(unlist(strsplit(can, " "))[1], can_ls); endd = srt + 1
  col1 = entrz2hgnc(paste0(R24_iMAT_gmcs[,srt],'.1'), gen_ls)
  col2 = entrz2hgnc(paste0(R24_iMAT_gmcs[,endd],'.1'), gen_ls)
  sl = data.frame(cbind(col1,col2))
  sl = sl[!is.na(sl$col1), ]; sl = sl[!is.na(sl$col2), ]
  id = grep(unlist(strsplit(fname[i],'Cancer'))[1], fname)
  
  if (length(id) != 0){
    tmp_ac = read.csv(paste0('../MetabolicSLOutput/',fname[id]))
    tmp_driver = intersect(unique(tmp_ac$Mutated_gene), tmp_gls)
    tmp_ko = intersect(unique(tmp_ac$KO_gene), tmp_gls)
    tmp_drg = unique(prism_mapping_met$drug_info[prism_mapping_met$target_info %in% tmp_ko])
    total = expand.grid(driver = (tmp_driver) , ko_gene = tmp_ko)
    # tmp_ac = length(unique(tmp_ac$KO_gene_drug)) *length(unique(tmp_ac$Mutated_gene))
    # tmp_ac = length(unique(tmp_ac$KO_gene_drug)) *length(unique(tmp_ac$Mutated_gene))
    cantyp_hyp_gmcs$cancer[i] = can
    cantyp_hyp_gmcs$total[i] = 0
    if (nrow(total)>0){
      chek_concordance = cbind(total, prism = NA, MN = NA, slscan = NA)
      for (chec_i in 1:nrow(chek_concordance)){
        chck_prsm = intersect(which(chek_concordance$driver[chec_i] == tmp_ac$Mutated_gene) 
                              ,which(chek_concordance$ko_gene[chec_i] == tmp_ac$KO_gene))
        if (length(chck_prsm) != 0 ){
          if (any(tmp_ac$ttst_metabolic_fdr[chck_prsm] <= 0.05)){chek_concordance$slscan[chec_i] = "YES"} else {chek_concordance$slscan[chec_i] = "NO"}
          if (any(tmp_ac$ttst_PRISM_pvl[chck_prsm] <= 0.05)){chek_concordance$prism[chec_i] = "YES"} else {chek_concordance$prism[chec_i] = "NO"}
        }
        
        
        chck1 = intersect(which(chek_concordance$driver[chec_i] == sl$col1), which(chek_concordance$ko_gene[chec_i] == sl$col2))
        chck2 = intersect(which(chek_concordance$driver[chec_i] == sl$col2), which(chek_concordance$ko_gene[chec_i] == sl$col1))
        if (length(chck1) != 0 | length(chck2) != 0){
          chek_concordance$MN[chec_i] = "YES"
          print("YES ............................")
        } 
      }
      
      chek_concordance_ttl = rbind(chek_concordance_ttl, cbind(can, chek_concordance))
      # hyper_geom test
      cantyp_hyp_gmcs$overlp[i] = length(which(chek_concordance$MN == "YES"))
      cantyp_hyp_gmcs$no_SL_psrm[i] = length(which(chek_concordance$prism == "YES"))
      cantyp_hyp_gmcs$no_SLscan[i] = length(which(chek_concordance$slscan == "YES"))
      cantyp_hyp_gmcs$no_SL_MN[i] = length(which(!is.na(chek_concordance$MN)))
      cantyp_hyp_gmcs$total[i] = nrow(chek_concordance)
      
      cantyp_hyp_gmcs$hypergeom[i] = phyper(cantyp_hyp_fsl$overlp[i]-1, cantyp_hyp_fsl$no_SL_MN[i]
                                            , cantyp_hyp_fsl$total[i]-cantyp_hyp_fsl$no_SL_MN[i]
                                            , cantyp_hyp_fsl$no_SL_psrm[i], lower.tail = FALSE)
      
      cantyp_hyp_gmcs$hypergeom[i] = phyper(cantyp_hyp_gmcs$overlp[i]-1, cantyp_hyp_gmcs$no_SL_MN[i]
                                            , cantyp_hyp_gmcs$total[i]-cantyp_hyp_gmcs$no_SL_MN[i]
                                            , cantyp_hyp_gmcs$no_SL_psrm[i], lower.tail = FALSE)
    }
  }
}
write.csv(chek_concordance_ttl, paste0('../MetabolicSLOutput/cons_prism_gmcs_concordance_check.csv'))

hyp_test = ls(pattern = 'cantyp_hyp')
for (hy_i in 1:length(hyp_test)){
  tmp = get(hyp_test[hy_i])
  print(tmp$overlp)
  rm(tmp)
}


for (i in 1:length(ls(pattern = '_hyp'))){
  tmpn = ls(pattern = '_hyp')[i]  
  tmp = get(tmpn)
  write.csv(tmp, paste0('../MetabolicSLOutput/cons_prism_', tmpn,'.csv'))
}

write.csv(SL_drug, paste0('../MetabolicSLOutput/cons_prism_SL_scan_drugs_overlap.csv'))

# panel boxplots cancer drug ####
maf_met = cbind(maf_met, DepMap_ID =  sampleInfo$DepMap_ID[match(maf_met$ccle_name, sampleInfo$CCLE_Name)])
library(reshape2)
library(ggplot2)

for ( i in 1:nrow(SL_drug)){
  cancer = SL_drug$cancer[i]
  dep_id_cellline = (ac_cantyp[match(cancer, names(ac_cantyp))]) 
  prism_row = (rownames(prism_met) %in% unlist(unname(dep_id_cellline)))
  drug_target = SL_drug$KO_gene_drug[i]
  prism_col = which(drug_target==colnames(prism_met))[1]
  
  prism_tmp = prism_met[prism_row,prism_col]
  names(prism_tmp) = rownames(prism_met)[prism_row]
  
  mut_gene = SL_drug$Mutated_gene[i]
  id = intersect(which(mut_gene == maf_met$Hugo_Symbol)
                 ,which("damaging" == maf_met$Variant_annotation))
  id = intersect(id, which(cancer == maf_met$primary_disease))
  mut_cellline = maf_met$DepMap_ID[id]
  
  x = prism_tmp[names(prism_tmp) %in% mut_cellline]
  y = prism_tmp[!(names(prism_tmp) %in% mut_cellline)]
  altn = 'less'
  ttest = t.test(x, y, altn, var.equal = FALSE)
  
  row_id_l = (rownames(KO_FBA) %in% names(prism_tmp)) 
  KO_FBA_tmp = KO_FBA[row_id_l, (colnames(KO_FBA) %in% SL_drug$KO_gene[i])]
  names(KO_FBA_tmp) = rownames(KO_FBA)[row_id_l]
  KO_FBA_tmp[is.na(KO_FBA_tmp)] = 1
  a = KO_FBA_tmp[names(KO_FBA_tmp) %in% mut_cellline]
  b = KO_FBA_tmp[!names(KO_FBA_tmp) %in% mut_cellline]
  a = a + runif(length(a), min = 1e-12, max = 1e-11)
  b = b + runif(length(b), min = 1e-12, max = 1e-11)
  tttt = t.test(a, b, altn, var.equal = FALSE)
  
  if (tttt$p.value < 0.01){
    print(paste(i, cancer, ttest$statistic, mut_gene, drug_target, ttest$p.value, sep = " ** "))
    print(paste(tttt$statistic, tttt$p.value, sep = " ======= "))
    ttest$statistic = formatC(ttest$statistic, format = "e", 2)
    ttest$p.value = formatC(ttest$p.value, format = "e", 2)
    
    df = data.frame(variable = c(rep("Mutated_PRISM",length(x)), rep("Not_Mutated_PRISM",length(y))
                    ,c(rep("Mutated_SL-scan",length(a)), rep("Not_Mutated_SL-scan",length(b))))
                    , value =  c(x,y,a,b))
    pdf(paste0("../MetabolicSLOutput/SL_drug_BoxPlot_",cancer, "_MUt_",mut_gene
               , "_drugTarget_", drug_target , ".pdf"))
    plt = ggplot(df, aes(x = variable, y = value, fill = variable)) +
      geom_boxplot() +
      labs(x = NULL, y = "Data value")+
      ggtitle(paste0("t statistic: ",ttest$statistic, " P-value: "
         ,ttest$p.value, ", drug: ", unlist(strsplit(drug_target, "=="))[1])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(plt)
    dev.off()
  }
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

png(paste0("../MetabolicSLOutput/dodge_prism_barplot_valid_SL_no.png"),width=1000, height=700)
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
# save.image('./env/WS_MN_exh_ttest_cor_cons_prism.rdata')
# load('./env/WS_MN_exh_ttest_cor_cons_prism.rdata')




