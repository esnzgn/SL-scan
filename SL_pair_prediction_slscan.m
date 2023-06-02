
% #########  chunk 3 - SL-scan implementation using Recon2.v4 and iMAT for SL pair prediction #########  
clear

s1 = pwd;
s2 = 'D:\G\thesis\thesis\contextualization\MetabolicSL\code';
if ~strcmp(s1,s2)
    cd(s2);
end
clear s1 s2

global CBTDIR
load([CBTDIR filesep 'test' filesep 'models' filesep 'mat' filesep 'Recon2.v04.mat']);
modelR204.description = 'modelR204'; clear CBTDIR;

%changeCobraSolver ('ibm_cplex', 'all');
%changeCobraSolver ('gurobi', 'all');
epsilon = getCobraSolverParams('LP', 'feasTol');

% rxn 2173 must have the following grRules and rules but not mentioned in the
% model properly
idx_wrul = find(contains(modelR204.rules,'()','IgnoreCase',true));
modelR204.rules(idx_wrul) = {'(x(1039)) | (x(1040)) | (x(1041)) | (x(1042)) | (x(1043)) | (x(1044)) | (x(1045)) | (x(1046)) | (x(1047)) | (x(1048)) | (x(1049)) | (x(1050))'};
modelR204.grRules{idx_wrul, 1} = {'(761.1) or (771.1) or (377677.1) or (762.1) or (766.2) or (23632.1) or (766.1) or (771.2) or (759.1) or (768.1) or (760.1) or (765.1)'};

% model correction: genes '6241.1', '6242.1', '50484.1' all must be in grRules of 
% 3208, 3209, 3210, 3211 rxns as '(6240.1) and (6241.1 or 50484.1)'
find(ismember(modelR204.genes,'6240.1'));
find(ismember(modelR204.genes,'6241.1'));
find(ismember(modelR204.genes,'50484.1'));
genes = {'6241.1', '6242.1', '50484.1'};
[Reaclist] = findRxnsActiveWithGenes(modelR204, genes);
rxnns = find(ismember(modelR204.rxns,Reaclist)); clear Reaclist;
rxnns = rxnns(1:end-1);
modelR204.grRules(rxnns) = cellstr('(6240.1) and (6241.1 or 50484.1)');
modelR204.rules(rxnns) = cellstr('(x(1764)) & (x(1765) | x(1766))'); clear rxnns;

%changeCobraSolver ('ibm_cplex', 'all');
%changeCobraSolver ('glpk', 'all');
%changeCobraSolver ('gurobi', 'all');


% find modelR204 ATP and biomass rxns to add to contextualized models
atp_idx = find(contains(modelR204.rxnNames,'ATP demand','IgnoreCase',true));
atp_idx1 = find(contains(modelR204.rxns,'ATPM','IgnoreCase',true));
bio_idx = find(contains(modelR204.rxns,'biomass','IgnoreCase',true));
addi_idx = vertcat(bio_idx, atp_idx,atp_idx1);
clear atp_idx bio_idx atp_idx1 options idx_wrul ans genes

% reading expression table
% reading cancer core gene list for 22 cancer of fastcore 
T1 = readtable('..\..\MetabolicSLInput\data\rnaData_gymbol_entrz_cel.csv');
gene_ls = T1(2:end,1:2); 
tb_gene_ls = gene_ls;
writetable(tb_gene_ls,'../../MetabolicSLOutput/tb_gene_ls.csv');
%mapping of cancer and celllines names 
cancer_mapping = vertcat(T1.Properties.VariableNames, T1(1,:));
cancer_mapping(:,1:2) = [];

T1(:,1:2) = []; T1(1,:) = [];
expression.gene  = table2cell(gene_ls(:,2 ));
cantyp = table2cell(readtable('..\..\MetabolicSLInput\data\cantyp.csv'));
cel_all = T1.Properties.VariableNames';

mut = readtable('..\..\MetabolicSLInput\data\maf_met_ge.csv');
mut_dep = mut.DepMap_ID;

options.tol = epsilon;
options.core = modelR204.rxns(addi_idx);
options.solver = 'iMAT';

modelR204_cons = RPMI_Medium(modelR204);
gen_cons_FBA = optimizeCbModel(modelR204_cons); gen_cons_FBA.f

%lung cell lines with most sig ko gene CYBRD1 
%lng_cel = {'ACH_000335', 'ACH_000610', 'ACH_000844', 'ACH_000869', ...
   %                 'ACH_000888', 'ACH_000900', 'ACH_000924', 'ACH_000929'};
%lng_celid = find(ismember(T1.Properties.VariableNames, lng_cel(8)));
%T1.Properties.VariableNames(lng_celid )
% preparing KO score for each cell line
% ttl_expressionRxns_all_R24 = [];
%T0 = T1(:,lng_celid);
all_SL = {};
gen_cell_KO = [];
ko_res = cell(size(T1,1),size(T1,2));
%geneee = {'50484.1', '586.1', '54675.1', '686.1', '2628.1', '4706.1'}';
%ko_res_sg = cell(50,2);
%for i=1:size(T0,2)

for i=1:size(T1,2)
%for i=1:28
    disp('*********************************************') 
    disp('loop')
    disp(i)
    %}
    % reaction expression mapping
    %expression.value =  table2array(T0(:,i));
    expression.value =  table2array(T1(:,i));
    expression.value =  str2double(expression.value);
    [levels parsedGPRSP] = mapExpressionToReactions(modelR204, expression); % defualt: minMax
    
    options.expressionRxns = levels;
    res = {};
    for ub_i = 1:10
        for lb_i =5:19
            disp(((20 - ub_i) * 5))
            disp(((20 - lb_i) * 5))
            options.threshold_ub = prctile(levels(find(levels > -1)),(20 - ub_i) * 5);
            options.threshold_lb = prctile(levels(find(levels > -1)),(20 - lb_i) * 5);
            tmpl = createTissueSpecificModel(modelR204_cons, options);
            % ckeck for biomass and atp demand rxns activity
            id0 = find(contains(tmpl.rxns, options.core{1},'IgnoreCase',true));
            id1 = find(contains(tmpl.rxns, options.core{2},'IgnoreCase',true));
            tmp_fba= optimizeCbModel(tmpl);
            
            rxn_id = find(contains(tmpl.rxns, 'Biomass_reaction','IgnoreCase',true));
            %precursors = findReactionPrecursors(tmpl, rxn_id);
            %[prerxn, prerxn_id] = findReactionrxnPrecursors(model, reaction_id);
            res{ub_i,lb_i} = tmp_fba.f;
            if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.x([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
                break % break lb_i loop
            end
        end
        if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.x([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
            all_SL{1, i} = char(strcat( cel_all{i}));
            all_SL{10, i} = ((20 - ub_i) * 5);
            all_SL{11, i} = ((20 - lb_i) * 5);
            all_SL{3, i} = tmp_fba.f;
            all_SL{2, i} = tmpl;
            break % break ub_i loop
        else
            all_SL{1, i} = char(strcat( cel_all{i}));
            all_SL{10, i} = 'NA';
            all_SL{11, i} = 'NA';
            all_SL{3, i} = 0;
            all_SL{2, i} = 'NA';
        end
    end

    if ~eq(all_SL{3, i}, 0)
        tmp_fba = optimizeCbModel(all_SL{2, i}, 'max', 'one');
        rxn_id = (find(tmp_fba.x));
        genesList = {};
        for i_gr = 1:length(rxn_id)
            tmp = rxn_id(i_gr);
            if (~isempty(all_SL{2, i}.rules{tmp ,1}))
                genesList = [genesList, all_SL{2, i}.genes{find(all_SL{2, i}.rxnGeneMat(tmp,:))}];
            end
        end
        del_gen = genesList';
        del_gen = unique(del_gen);
        [genels_id, del_geid] = ismember(gene_ls.Var2, del_gen);
        del_gen = gene_ls.Var2(genels_id);
        clear rxn_id tmp genesList genels_id del_geid
        %findGenesFromRxns(model,rxnsList(1))        
        %disp(find(contains(all_SL{4, i}.rxns,'biomass','IgnoreCase',true)))
        %iddddd = find(tmp_fba.x);
        %iddddd = all_SL{2, i}.rxns(iddddd);
        % replace gene names because some of them are renamed and
        % following functions not work properly with modified gene
        % names
        % find FBA sol of each active genes to construct similar ataris
        % achiles matrix score for ttest and cortest
        %sample1 = optimizeCbModel(all_SL{4, i}); sample1 = sample1.x; 
        %tests = 'tTest';
        %{
        sampleFile = 'samplingf'; option.nStepsPerPoint = 10; option.nPointsReturned = 50; samplerName = 'ACHR';
        [modelSampling, samples] = sampleCbModel(all_SL{2, i}, sampleFile, samplerName, option);
        id10 = find(contains(modelSampling.rxns,'biomass','IgnoreCase',true));
        itr = 1;
        %}
        for d =1:length(del_gen)
            [gnd0, ~] = ismember(all_SL{2, i}.genes, del_gen(d));
            %[gnd0, ~] = ismember(geneee , del_gen(d));
            if ~isempty(find(gnd0))
                [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(all_SL{2, i}, 'FBA', del_gen(d));
              % if grRatio < .90
                    %disp(d)
                    %{
                    [modelDel,hasEffect,constrRxnNames] = deleteModelGenes(all_SL{2, i},del_gen(d));
                    sampleFileKO = 'samplingfKO'; option.nStepsPerPoint = 10; option.nPointsReturned = 50; samplerName = 'ACHR';                    
                    [modelSamplingKO, samplesKO] = sampleCbModel(modelDel, sampleFileKO, samplerName, option);
                    idKO = find(contains(modelSamplingKO.rxns,'biomass','IgnoreCase',true));
                    if isempty(idKO)
                        [ h, p, ci, stats] = ttest2(zeros(size(samplesKO,2),1)', samples(id10,:));
                    else
                        [ h, p, ci, stats] = ttest2(samplesKO(idKO,:), samples(id10,:));
                    end
                    %}
                    ko_gid = find(ismember(gene_ls.Var2, del_gen(d)));
                    %ko_res{ko_gid, i} = stats.tstat;  
                    %ko_res{itr,1} = stats.tstat
                    ko_res{ko_gid,i} = grRatio;
                    %itr = itr + 1;
                    %gen_ko_sm = horzcat(geneee, ko_res_sg);
                    %tb_gen_ko_sm = array2table(gen_ko_sm);
                    %writetable(tb_gen_ko_sm, '../../MetabolicSLOutput/tb_gen_ko_sm.csv');
                    %{
                    iiiiid= cell(150,2);
                    iiiiidi = 0;  
                    for iiii = 1:size(ko_res,1)
                        if ~isempty(ko_res{iiii, 1})
                            if ~isnan(ko_res{iiii, 1})
                                iiiiidi = iiiiidi + 1;
                                iiiiid{iiiiidi, 1} = gene_ls.Var1(iiii);
                                iiiiid{iiiiidi, 2} = ko_res{iiii, 1};
                            end        
                        end
                    end
                    tb_iiiiid = array2table(iiiiid);
                    writetable(tb_iiiiid, '../../MetabolicSLOutput/KO_res_rnd_smple.csv');
                    %}
              % end
            end
        end      
        % r = xcorr( ko_res_sg(1:9,1), ko_res_sg(1:9,2))
        clear del_gen d gnd0 grRatio grRateKO grRateWT hasEffect delRxns fluxSolution grRatio ko_gid 
    end
    % find exact mutational backgroung of each sample in comparison to the
    % GEM used to construct it
    all_SL{2, i} = [] ;
     if (size(all_SL,2) == 100 )
        save('../../MetabolicSLOutput/all_SL_save_100_ttestV.mat ', 'all_SL', '-v7.3');
        tb_ko_res = array2table(ko_res); % tb_gen_cell_KO
        tb_ko_res.Properties.VariableNames = cel_all(:, 1);
        writetable(tb_ko_res,'../../MetabolicSLOutput/KO_res_100_ttestV.csv');
        clear tb_ko_res;
     elseif (size(all_SL,2) == 200 )
        save('../../MetabolicSLOutput/all_SL_save_200_ttestV.mat ', 'all_SL', '-v7.3');
        tb_ko_res = array2table(ko_res); % tb_gen_cell_KO
        tb_ko_res.Properties.VariableNames = cel_all(:, 1);
        writetable(tb_ko_res,'../../MetabolicSLOutput/KO_res_200_ttestV.csv');
        clear tb_ko_res;
    elseif (size(all_SL,2)  == size(ko_res,2) )
        save('../../MetabolicSLOutput/all_SL_save_446_ttestV.mat ', 'all_SL', '-v7.3');
        tb_ko_res = array2table(ko_res);
        tb_ko_res.Properties.VariableNames = cel_all(:, 1);
        writetable(tb_ko_res,'../../MetabolicSLOutput/KO_res_446_ttestV.csv');
        clear tb_ko_res;
    end
end

% saved as ko_res.mat
% all_sl.mat
% merg two table of ko_res of 319 and 446 to make the total one
ko_mrg = table2cell(readtable('..\..\MetabolicSLoutput\KO_res_319_ttestV.csv'));
ko_mrg =  horzcat(ko_mrg(:,1:319), ko_res(:,320:446));

cel_all = T1.Properties.VariableNames';
mut = readtable('..\..\MetabolicSLInput\data\maf_met_ge.csv');

