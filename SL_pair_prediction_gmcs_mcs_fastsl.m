%%
% ######### chunk 1 - cancer modeling using Recon2.04 and iMAT and SL pair prediction using FastSL, MCS, and gMCS #########

clear

s1 = pwd;
s2 = 'D:\G\thesis\thesis\contextualization\MetabolicSL\code';
if ~strcmp(s1,s2)
    cd(s2);
end
clear s1 s2

global CBTDIR
load([CBTDIR filesep 'test' filesep 'models' filesep 'mat' filesep 'Recon2.v04.mat']);
modelR204.description = 'modelR204';

%changeCobraSolver ('ibm_cplex', 'all');
%changeCobraSolver ('gurobi', 'all');

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
rxnns = find(ismember(modelR204.rxns,Reaclist));
rxnns = rxnns(1:end-1);
modelR204.grRules(rxnns) = cellstr('(6240.1) and (6241.1 or 50484.1)');
modelR204.rules(rxnns) = cellstr('(x(1764)) & (x(1765) | x(1766))');

%changeCobraSolver ('ibm_cplex', 'all');
%changeCobraSolver ('glpk', 'all');
%changeCobraSolver ('gurobi', 'all');

s1 = pwd;
s2 = 'D:\G\thesis\thesis\contextualization\MetabolicSL\code';
if ~strcmp(s1,s2)
    cd(s2);
end
clear s1 s2 ans CBTDIR genes idx_wrul rxnns Reaclist 


% find modelR204 ATP and biomass rxns to add to contextualized models
atp_idx = find(contains(modelR204.rxnNames,'ATP demand','IgnoreCase',true));
atp_idx1 = find(contains(modelR204.rxns,'ATPM','IgnoreCase',true));
bio_idx = find(contains(modelR204.rxns,'biomass','IgnoreCase',true));
addi_idx = vertcat(atp_idx,bio_idx,atp_idx1);
clear atp_idx bio_idx atp_idx1 options

% reading expression table
% reading cancer core gene list for 22 cancer of fastcore 
T = readtable('..\..\MetabolicSLInput\data\rnaData_gymbol_entrz_cel.csv');
gene_ls = T(:,1:2);
gene_ls(1,:) = [];
T(:,1:2) = [];
expression.gene  = table2cell(gene_ls(:,2));
cantyp = table2cell(readtable('..\..\MetabolicSLInput\data\cantyp.csv'));
ttl_expressionRxns_R24 = [];
ttl_expressionRxns_samples_R24 = [];
ttl_expressionRxns_all_R24 = [];
for exp_i = 1:size(cantyp,1)    
    strt = cantyp{exp_i,3}-1;
    endd = cantyp{exp_i,4}-1;
    expression.value =  table2array(T(:,strt:endd));
    expression.value =  str2double(expression.value);
    [expressionRxns parsedGPR] = mapExpressionToReactions(modelR204, expression); % defualt: minMax
    %P = prctile(sort(expressionRxns),70)
    %res = expressionRxns > (mean(expressionRxns)+(2 * (std(expressionRxns))));
    ttl_expressionRxns_R24 = horzcat(ttl_expressionRxns_R24, expressionRxns);
    
    %{
	% for each sample of cancer type
    for exp_j =1:size(expression.value,2)
        expressionSp.value = expression.value(:,exp_j);
        [expressionRxnsSP parsedGPRSP] = mapExpressionToReactions(modelR204, expressionSp); % defualt: minMax
        ttl_expressionRxns_samples_R24 = horzcat(ttl_expressionRxns_samples_R24, expressionRxnsSP);
        disp(' sample reaction mapping remaining:') 
        disp(size(expression.value,2) - exp_j)
        disp(' out of:') 
        disp(size(expression.value,2))
    end
    ttl_expressionRxns_all_R24 = horzcat(ttl_expressionRxns_all_R24 , ttl_expressionRxns_samples_R24);
     %}
    
    disp('*** reaction mapping ***') 
    disp('remaining:')
    disp(size(cantyp,1) - exp_i)
    
end 
clear exp_i strt endd expression expressionRxns parsedGPR gene_ls T  ans
%length(find(res))
% here i am
core_g = struct([]);
epsilon = getCobraSolverParams('LP', 'feasTol');
options.tol = epsilon;

% using fastcc the flux consistent rxns will be determined for next
% sessions
%[A,orientation,V] = fastcc(modelR204,epsilon);

%ttl_expressionRxns_all; cantyp{:,1};
modelR204_cons = RPMI_Medium(modelR204);
gen_cons_FBA = optimizeCbModel(modelR204_cons); gen_cons_FBA.f
options.core = modelR204_cons.rxns(addi_idx);
options.solver = 'iMAT';

% applying medium culture


for i = 1:(size(cantyp,1))
    disp('*********') 
    disp('loop')
    disp(i)
    %%{
    %disp('sub_loop')
    %disp(1)
    
    endd = i; 
    levels = ttl_expressionRxns_R24(:,endd);
    core_g{5, i} = cantyp(endd,2);

    %disp(3)
    %funcModel = 0;
    
    options.expressionRxns = levels;
    res = {};
    for ub_i = 1:10
        for lb_i =10:19
            disp(((20 - ub_i) * 5))
            disp(((20 - lb_i) * 5))
            options.threshold_ub = prctile(levels(find(levels > -1)),((20 - ub_i) * 5));
            options.threshold_lb = prctile(levels(find(levels > -1)),((20 - lb_i) * 5));
            tmpl = createTissueSpecificModel(modelR204_cons, options);
            % ckeck for biomass andn atp demand rxns activity
            id0 = find(contains(tmpl.rxns, options.core{1},'IgnoreCase',true));
            id1 = find(contains(tmpl.rxns, options.core{2},'IgnoreCase',true));
            tmp_fba= optimizeCbModel(tmpl);
            res{ub_i,lb_i} = tmp_fba.f;
            if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.x([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
                break % break lb_i loop
            end
        end
        if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.x([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
            core_g{1, i} = res;
            core_g{2, i} = ((20 - ub_i) * 5);
            core_g{3, i} = ((20 - lb_i) * 5);
            core_g{4, i} = tmpl;
            core_g{4, i}.description = char(strcat('cons_R24_iMAT',core_g{5, i}));
            core_g{5, i} = char(strcat('cons_R24_iMAT_',core_g{5, i}));
            core_g{6, i} = tmp_fba.f;
            break % break ub_i loop
        end
    end
end
clear levels endd tmp_fba options i ttl_expressionRxns addi_idx cantyp levels endd lb_i ub_i gen_cons_FBA ...
options id0 id1 tmp_fba res tmpl ans ttl_expressionRxns_R24 ttl_expressionRxns_all_R24 ttl_expressionRxns_samples_R24

% mdl verfication
for modl_i = 1:size(core_g,2)
   tmp = core_g{4, modl_i};
   LPproblem = buildLPproblemFromModel(tmp);
    %initial check that all inputs are valid:
    if ~(verifyCobraProblem(LPproblem, [], [], false) == 1)
        warning(strcat('invalid problem with model:  ',string(modl_i)));
        return;
    end
    %a more strict check weather all models are valid or not:
    %core_g
    results = verifyModel(tmp,'fluxConsistency', true,'massBalance', true);
    %strcat('no of inconsistent rxns:  ',string(length(find(results.fluxConsistency.consistentReactionBool < 1))))
    %strcat('inconsistent rxns:  ')
    %results.fluxConsistency.consistentReactions{find(results.fluxConsistency.consistentReactionBool < 1),1}
    %4) do a simple check with a specified set of fields.
    results.simpleCheck = verifyModel(tmp,'simpleCheck', true);
    results.requiredFields = verifyModel(tmp,'requiredFields',{'S', 'lb', 'ub', 'c', 'rxns', 'mets', 'genes', 'rules'});
    core_g{7, modl_i} = results;
    clear results tmp
end
clear tmp modl_i LPproblem results ans

%[A,orientation,V] = fastcc(core_g{6, i},epsilon);
%optimizeCbModel(core_g{6, i})

%{%} test fba
for fbi = 1:size(core_g,2)
    atp_id = find(contains(core_g{4, fbi}.rxnNames,'ATP demand','IgnoreCase',true));
    disp(find(contains(core_g{4, fbi}.rxns,'biomass','IgnoreCase',true)))
    disp(core_g{4, fbi }.rxns(atp_id))
    disp(core_g{4, fbi }.rxns(find(core_g{4, fbi}.c)))
    %fastcc(core_g{6, fbi},epsilon,printLevel);
    aaa = optimizeCbModel(core_g{4,fbi});
    disp(aaa.full(find(core_g{4, fbi }.c)))
    disp(aaa.full(atp_id))
end
clear fbi atp_id aaa atp_id

% to test models functionalities
[totTest, numTests] = testFunctionality_old(modelR204_cons)
gnr_mdl_FBA = optimizeCbModel(modelR204_cons);
gnr_mdl_FBA.f

for j = 1:size(core_g,2)
    if ~isempty(core_g{4, j})
        [totTest, numTests] = testFunctionality_old(core_g{4, j});
        % percent of functionalities
        core_g{8, j} = totTest;
    end
end
clear gnr_mdl_FBA totTest numTests j spc_mdl_FBA

% thr functionality checking 100 iteration 17 min single core
% I wrote the following lines after a while so load the core strct again,
% if you are doing all sections sequential its not need to load it again
% so do not uncomment
% load('../../MetabolicSLOutput/R24_iMAT_core_g.mat'); core_g = iMAT_core_g;
ncor = size(core_g,2);
%M = 5;                     % M specifies maximum number of workers
%parfor  (j = 3:ncor,M)
tmp = {strcat('modelR204_cons','_iMAT')};
tl_msg = [];
seq_fun_fin = [];
for  j = 1:ncor
    disp(j)
    n_iter = 100;
    if ~isempty(core_g{4, j})
        [CS_totTest, seq_tst, top_1p, top_50p, Iqr, Max, seq_fun ] = functionality_thr(modelR204_cons, core_g{4, j}, n_iter);
        seq_fun_fin = vertcat(seq_fun_fin, seq_fun); 
        tl_msg = vertcat(tl_msg , horzcat(CS_totTest, top_1p, top_50p, Iqr, seq_tst));
    end
end
% cell to table conversion
tbl_msg = array2table((tl_msg));
tbl_msg .Properties.VariableNames([1:4]) = {'NoFuncs' 'top1p' 'top50p' 'IQR' };
writetable(tbl_msg , strcat('../../MetabolicSLOutput/func_test_ modelR204_cons_iMAT .csv'));
tl_seq_fun_fin = array2table((seq_fun_fin));
writetable(tl_seq_fun_fin,strcat('../../MetabolicSLOutput/func_test_seq_modelR204_cons_iMAT.csv'));

%{
% do FastSL on 22 can_specific models of core_g in core_g{7, i}
addpath('.\fastSL_Rxns');
%Exhange reactions Eliminate in parfastSL func for lethlaity analysis
parfor j = 1:size(core_g,2), parfastSL(core_g, j); end
delete(gcp('nocreate'))
clear j ans
% using func FastSLIO which find output of fastSL from HDD
core_g = FastSLIO(core_g,'R24_iMAT');
%save('../../MetabolicSLOutput/RD3iMAT_core_g.mat', 'core_g', '-v7.3')
%}

%%{
%Eliminate Exhange reactions for lethlaity analysis
for j = 1:size(core_g,2)
    if ~isempty(core_g{4, j})
         eliList  = {};
        for i_gr = 1:length(core_g{4, j}.rules)
            if (isempty(core_g{4, j}.rules{i_gr,1}))
                eliList = vertcat(eliList, core_g{4, j}.rxns{i_gr});
            end
        end
            %eliList = core_g{4, j}.rxns(find(findExcRxns(core_g{4, j}))); 
            fastSL(core_g{4, j},0.01,2,eliList); 
    end
    %core_g{7, j} = load('fastcoreM_ can name ... _Rxn_lethals.mat');
    %disp(core_g{7, j})
end
%}
% using func FastSLIO which find output of fastSL from HDD
core_g = FastSLIO(core_g);
save('../../MetabolicSLOutput/R24_cons_iMAT_core_g.mat', 'core_g', '-v7.3')

% here #########
% searching for SLs using gMCS 
%parpool(4)
%parfor j = 1:size(core_g,2)
for j = 1:size(core_g,2)
    if ~isempty(core_g{4, j})
        disp('gmcs:')
        disp(j)
        % searching for SLs using gMCS and MCS algorithms
        % searching for SLs using gMCS 
        model_name = core_g{4, j}.description;
         core_g{4, j} = buildRxnGeneMat(core_g{4, j});
        %model_name = 'Recon2v2';
        n_gMCS = 1000;
        max_len_gMCS = 2;
        %optional_inputs.KO = '6240';
        %optional_inputs.separate_transcript = '.';
        optional_inputs.timelimit = 4*60;
        optional_inputs.numWorkers = 4;

        % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
        % recommended if the G matrix has not been calculated yet.
        % parpool;
        [gmcs, ~] = calculateGeneMCS(model_name, core_g{4, j}, n_gMCS, max_len_gMCS, optional_inputs);
        core_g{9, j} = gmcs;
        filePattern    = fullfile(strcat('./G_',core_g{4, j}.description,'.mat'));
        delete(filePattern);
    end
end
clear gmcs gMCS_time optional_inputs max_len_gMCS n_gMCS model_name j filePattern

% searching for SLs using MCS 
%parpool(4)
for j = 1:size(core_g,2)
    if ~isempty(core_g{4, j})
            disp('mcs:')
    disp(j)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using MCS 
    %model_name = 'Recon2v2';
    n_MCS = 1000;
    max_len_MCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 4;
    rxnsList = {};
    for i_gr = 1:length(core_g{4, j}.rules)
        if (~isempty(core_g{4, j}.rules{i_gr,1}))
            rxnsList  = [rxnsList  ,core_g{4, j}.rxns{i_gr}];
        end
    end
    optional_inputs.rxn_set = rxnsList';
    [mcs, ~] = calculateMCS(core_g{4, j}, n_MCS, max_len_MCS, optional_inputs);
    core_g{10, j} = mcs;
    end
end
clear mcs mcs_time optional_inputs rxnsList i_gr max_len_MCS n_MCS j

% searching for SLs using doubleGeneDeletion 
%{
rxnsList = {};
j=1
for i_gr = 1:length(iMAT_core_g{6, j}.rules)
    if (~isempty(iMAT_core_g{6, j}.rules{i_gr,1}))
        rxnsList = [rxnsList ,iMAT_core_g{6, j}.rxns{i_gr}];
    end
end
rxnsList = rxnsList';
    
% genes = cellstr(model.genes(1:10));
len = length(rxnsList);
%ismember('ACACT1rm',rxnsList)
%ismember('ACACT1rm',iMAT_core_g{6,1}.rxns)
%gprs = findGPRFromRxns(iMAT_core_g{6,1}, cellstr(rxnsList));
%Updates the rules, genes, and rxnGeneMat fields based on the grRules field for reactions in the model.
%modelNew = removeUnusedGenes(iMAT_core_g{6,1});
%Update the model genes field (if it does not exist, generate it from the grRules field, if it exists, remove Unused genes and remove duplicate genes.
%[modelNew] = updateGenes(model);
   

for j = size(core_g,2)
    % searching for SLs using doubleGeneDeletion 
  
    optional_inputs.timelimit = 2*60;
    optional_inputs.numWorkers = 6;

    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    %core_g{6, j}.grRules = core_g{6, j}.rules;
    % parpool;
    [grRatioDble, grRateKO, grRateWT] = doubleGeneDeletion(core_g{6, j});
    core_g{10, j} = grRatioDble;
    core_g{11, j} = grRateKO;
end


%parpool(4)
for j = 1:size(core_g,2)
    % searching for SS (synthetic sick) and SL cases of doubleGeneDeletion
    %core_g{12, j} = core_g{10, j}<= 0.01;
    mat = core_g{10, j};
    %core_g{12, j} = zeros(size(core_g{10, j},1),2);
    res = [0,0];
    for i = 1:size(mat,1)
        for io = 1:size(mat,2)
            if i ~= io
                if (mat(i,io) / min(mat(i,i), mat(io,io))) < (epsilon*10)
                    idx_i = find(ismember(res(:,2),i));
                    if find(ismember(res(idx_i,1),io))
                    else
                        res = vertcat(res,[i, io]);   
                        %core_g{12, j} = vertcat(core_g{12, j},[i,io]);
                    end
                end
            end
        end
    end
    res(1,:) = [];
    core_g{12, j} = res;
        
end
%}

% gMCS SL res purifier
trans_sep = '.';
max_gmcs = 100;
iMAT_gmcs = [];
for modl_i = 1:size(core_g,2)
    if ~isempty(core_g{9, modl_i})
        gmc = gmcs_purifier(core_g{9, modl_i},trans_sep);
        %gmc.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
        output0 = vertcat(table2array(gmc) ,zeros(max_gmcs  - size(gmc,1),2));
        iMAT_gmcs  = [iMAT_gmcs  output0];
    end
end
clear trans_sep max_gmcs modl_i gmc output0

% fastSL jdl res purifier
trans_sep = '.';
max_jDL = 1000;
iMAT_jdl = [];
for modl_i = 1:size(core_g,2)
    if ~isempty(core_g{7, modl_i})
        output0 = Jsl_purifier(core_g{4, modl_i}, core_g{7, modl_i}.Jsl,trans_sep);
        output1 = Jdl_purifier(core_g{4, modl_i}, core_g{7, modl_i}.Jdl,trans_sep);
        output = vertcat(output0,output1);
       if isempty(output)
           output = table(([str2double('NA'), str2double('NA')]));
       end
        [~,idx] = unique(output,'rows','first');
        output = output(sort(idx),:);
        output2 = vertcat(table2array(output) ,zeros(max_jDL - size(output,1),2));
        iMAT_jdl  = [iMAT_jdl  output2];
    end
end
clear max_jDL modl_i JDL output Jdl model output0 trans_sep output output0 output1 output2 idx ans filePattern

% MCS order 2 res purifier
trans_sep = '.';
max_mcs = 1000;
iMAT_mcs = [];
for modl_i = 1:size(core_g,2)
    mcs = mcs_purifier(core_g{4, modl_i}, core_g{10, modl_i},trans_sep);
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    if isempty(mcs)
        output0 = zeros(max_mcs  - size(mcs,1),2);    
    else
        output0 = vertcat(table2array(mcs) ,zeros(max_mcs  - size(mcs,1),2));    
    end
    iMAT_mcs = [iMAT_mcs output0];
end
clear output0 mcs modl_i max_mcs trans_sep

%load('D:\G\thesis\thesis\contextualization\final_res\iMAT_Recon_204_can_specific_gMCS_MCS_fastSL.mat');

% find essential genes for comparison to CRISPR KO experiments
HasEff = [];
for esnGn_i =1:size(core_g,2)
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(core_g{6, esnGn_i});
    HasEff = horzcat(cellstr(core_g{6, esnGn_i}.genes(find(grRatio<.99))),string(grRatio(find(grRatio<.99))));
    core_g{12, esnGn_i} = HasEff;
end
clear HasEff grRatio grRateKO grRateWT hasEffect delRxns fluxSolution esnGn_i 

% find gene list of each model for dd method
max_gls = 3000;
% find each model gene list
iMAT_gls = [];
for modl_i = 1:size(core_g,2)
    if ~isempty(core_g{4, modl_i})
        disp(modl_i)
        iMAT_gls_tmp1 = core_g{4, modl_i}.grRules;
        iMAT_gls_tmp = [];
        for i_rxn = 1:length(iMAT_gls_tmp1)
            delimiter = {'or' 'and'};
            match = {'(' ')' ' '};
            bbb = strsplit(string(iMAT_gls_tmp1{i_rxn}),delimiter);
            bb = erase(bbb,match);
            if length(bb) == 1 && isempty(bb{1,1})
            else
                 iMAT_gls_tmp = [iMAT_gls_tmp bb];
            end
        end
        iMAT_gls_tmp = unique(iMAT_gls_tmp)';

        output0 = vertcat(iMAT_gls_tmp ,num2cell(zeros(max_gls  - length(iMAT_gls_tmp ),1)));    
        iMAT_gls = [iMAT_gls output0];
    end
end
clear iMAT_gls_tmp output0 bb bbb delimiter modl_i max_gls i_rxn match iMAT_gls_tmp1

% find each model essential gene list
iMAT_gls_esn = [];
max_gls_ens = 130;
for modl_i = 1:size(core_g,2)
    FASTCORE_gls_tmp = core_g{12, modl_i};
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    
	output0 = vertcat(string([core_g{5, modl_i} core_g{5, modl_i}]),FASTCORE_gls_tmp ,num2cell(zeros(max_gls_ens  - length(FASTCORE_gls_tmp),2)));    
    
    iMAT_gls_esn = [iMAT_gls_esn output0];
end
clear max_gls_ens output0 FASTCORE_gls_tmp modl_i 
%here##########
iMAT_core_g = core_g;
clear core_g
% fastcore model gene lists write as csv
iMAT_glsT = cell2table(cellstr(iMAT_gls));
writetable(iMAT_glsT,'../../MetabolicSLOutput/R24_iMAT_cons_gls.csv'); clear iMAT_glsT;
% fastcore model gmcs lists write as csv
iMAT_gmcsT = array2table(iMAT_gmcs);
writetable(iMAT_gmcsT,'../../MetabolicSLOutput/R24_iMAT_cons_gmcs.csv'); clear iMAT_gmcsT;
% fastcore model mcs lists write as csv
iMAT_mcsT = array2table(iMAT_mcs);
writetable(iMAT_mcsT,'../../MetabolicSLOutput/R24_iMAT_cons_mcs.csv'); clear iMAT_mcsT;
% fastcore model jdl lists write as csv
iMAT_fastSLT = array2table(iMAT_jdl);
writetable(iMAT_fastSLT,'../../MetabolicSLOutput/R24_iMAT_cons_fastSL.csv'); clear iMAT_fastSLT;
% fastcore model essential gene lists write as csv
%iMAT_gls_esn = array2table(iMAT_gls_esn);
%writetable(iMAT_gls_esn,'../../MetabolicSLOutput/R24_iMAT_cons_essen_genes.csv'); clear iMAT_gls_esn;
% fastcore model percent of functionalities and portion of contextualized model fba sim to generic model res lists write as csv
iMAT_funcSim_ls_esn = array2table(iMAT_core_g(8, 1:16));
writetable(iMAT_funcSim_ls_esn,'../../MetabolicSLOutput/R24_iMAT_cons_func_sim.csv'); clear iMAT_funcSim_ls_esn;
ab = iMAT_core_g(5, :);
writetable(array2table(ab),'../../MetabolicSLOutput/R24_iMAT_cons_canls.csv'); clear ab

save('../../MetabolicSLOutput/R24_iMAT_core_g_cons.mat ', 'iMAT_core_g', '-v7.3')
% load('R24_iMAT_core_g_cons.mat')

%cantype
imat_coreg_cantyp = array2table(iMAT_core_g(5, 1:16));
writetable(imat_coreg_cantyp ,'../../MetabolicSLOutput/R24_imat_cons_cantyp.csv'); clear ab
% percent of functionalities core_g{10, j}
% fraction of wild generic model fba sol core_g{11, j}
% essential genes core_g{12, j}
% '../../MetabolicSLOutput/fastcore_essen_genes.csv'

%%
% #########  chunk 2 - functionality checks and mets conversion to metNames based on Recon2.0v4 #########  
func_nm = [];
for j = 1:numel(bmMets)
            met = [bmMets{j}];
           id = find(ismember(modelR204.mets,met))
           func_nm = horzcat(func_nm, modelR204.metNames(id});
end

solgn_r24 = optimizeCbModel(modelR204);
solgn_r24.f
solgn_r3d = optimizeCbModel(Recon3D);
solgn_r3d.f
names = who('-regexp', 'core_g');
fba_ls = [];
label = {};
for i =1: numel(names)
    eval(['core_g = ',names{i},';']);
    solcs = [];
    for j = 1:size(core_g, 2)
        soltp = optimizeCbModel(core_g{6, j});
        solcs = horzcat(solcs, soltp.f);
        label = horzcat(label, core_g{5, j}{1, 1}  );
    end
    fba_ls = vertcat(fba_ls, solcs);
     %fba_ls{i,j+1} = names{i}; 
end
names
x = [1:22]; y = fba_ls(4,:);
bar(x, y)
%n = {'m1';'m2';'m3';'m4';'m5'};
set(gca,'xticklabel', label(1, :))

tb_fba_ls = array2table(fba_ls);
tb_fba_ls.Properties.VariableNames(1:22) = label(1, 1:22);
writetable(tb_fba_ls ,'../../MetabolicSLOutput/tb_fba_ls.csv'); 