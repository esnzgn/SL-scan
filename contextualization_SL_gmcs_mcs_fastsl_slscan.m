%%
% ######### chunk 1 Recon3D fastcore #########
clear
% model loading and correction based on lit and my observations
%
% go to lewis folder of MEM functions
% s0 = D:\Program Files\MATLAB\R2016b\bin
changeCobraSolver ('ibm_cplex', 'all');
%changeCobraSolver ('glpk', 'all');
%changeCobraSolver ('gurobi', 'all');


s1 = pwd;
s2 = 'D:\G\thesis\thesis\contextualization\MetabolicSL\code';
if ~strcmp(s1,s2)
    cd(s2);
end
clear s1 s2

Recon3D = readCbModel('..\..\MetabolicSLInput\data\Recon3D.xml');
Recon3D = creategrRulesField(Recon3D);

epsilon = getCobraSolverParams('LP', 'feasTol');

% find modelR204 ATP and biomass rxns to add to contextualized models
atp_idx = find(contains(Recon3D.rxnNames,'ATP demand','IgnoreCase',true));
atp_idx1 = find(contains(Recon3D.rxns,'ATPM','IgnoreCase',true));
bio_idx = find(contains(Recon3D.rxns,'biomass','IgnoreCase',true));
addi_idx = vertcat(atp_idx,bio_idx,atp_idx1);
clear atp_idx bio_idx atp_idx1

% using fastcc the fulx consistent rxns will be determined for next
% sessions
[A,orientation,V] = fastcc(Recon3D,epsilon);
clear orientation V

% fastcore models reconstruction

% reading cancer core gene list for 22 cancer of fastcore 
T = readtable('..\..\MetabolicSLInput\data\cancer_core_genes_entrzR3D.csv');
T(:,1) = [];

%T = readtable('cancer_core_genes_entrz.csv');
%T(:,1) = [];
%T = T1;

% making cell struct which contains all models of fastcore and next step
% will save at lower elements row of the struct
for i = 1:size(T,2)
    disp('loop')
    disp(i)
    disp('sub_loop')
    disp(1)
    %i = 1
    %T(:,i).Properties.VariableNames
    T0 = T{:,i};
    disp(2)
    TF = {T(find(~ismissing(T0,{'NA'})),i)};
    disp(3)
    core_g(1,i) = TF;
    disp(4)
    core_g(2,i) = {ones(size(TF{1,1},1),1)};
    disp(5)
    expression.gene  = string(table2cell(core_g{1, i}));
    disp(6)
    expression.value = core_g{2, i};
    disp(7)
    [expressionRxns parsedGPR] = mapExpressionToReactions(Recon3D, expression); % defualt: minMax
    disp(8)
    core_g{3, i} = find(expressionRxns>0);
    disp(9)
    core_g{4, i} = intersect((core_g{3, i}),A);
    % adding biomass and atp demand
    core_g{4, i} = vertcat(core_g{4, i},addi_idx);
    disp(9.1)
    core_g{5, i} = T(:,i).Properties.VariableNames;
    disp(10)
	%options.expressionCol = expressionRxns;
    %expressionCol = expressionRxns;
    disp(11)
    tol = 1e-4;
    %scaling = 1;
    options.epsilon = tol;  
    options.core = core_g{4,i};  
    %options.core = (10:30)';  
    options.printLevel = 2;
    options.solver = 'fastCore';
    core_g{6, i} = createTissueSpecificModel(Recon3D, options);
    core_g{6, i}.description = char(strcat('R3D_fastcoreM_',core_g{5, i}));
    %{
    core_g{6, i} = creategrRulesField(core_g{6, i});
    core_g{6, i}.rev = zeros(length(core_g{6, i}.rxns),1);
    core_g{6, i}.rev(find(core_g{6, i}.lb<0)) = 1;
    core_g{6, i} = buildRxnGeneMat(core_g{6, i});
    core_g{6, i}.rxnGeneMat = double(core_g{6, i}.rxnGeneMat);
    %}
    disp(12)
    %core_g{7, i} = fastcc(core_g{6, i},epsilon);
    %disp(13)
end
clear i  options tol expressionRxns TF T0 parsedGPR A addi_idx T ans expression epsilon


%{%}
cntr = 0;
for fbi = 1:size(core_g,2)
    %fastcc(core_g{6, fbi},epsilon,printLevel);
    disp('loop:')
    disp(fbi)
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxnNames,'ATP demand','IgnoreCase',true))))
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxns,'ATPM','IgnoreCase',true))))
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxns,'biomass','IgnoreCase',true))))
    aaa = optimizeCbModel(core_g{6,fbi});
    disp(core_g{6, fbi }.description)
    disp(core_g{6, fbi }.rxns(find(core_g{6, fbi }.c)))
    disp(aaa.f)
    %disp(aaa.v(find(contains(core_g{6, fbi}.rxnNames,'biomass','IgnoreCase',true))))
    idx = find(contains(core_g{6, fbi}.rxnNames,'ATP demand','IgnoreCase',true));
    disp(core_g{6, fbi }.rxns(idx))
    disp(aaa.v(idx ))
    if aaa.f > 0
        cntr = cntr + 1;
    end
    cntr
end
clear cntr aaa idx ans fbi

% do FastSL on 22 can_specific models of core_g in core_g{7, i}
addpath('.\fastSL_Rxns');
%Eliminate Exhange reactions for lethlaity analysis
for j = 1:size(core_g,2)
	rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList = [rxnsList ,core_g{6, j}.rxns{i_gr}];
        end
    end
    rxnsList = rxnsList';
    %findGenesFromRxns(model,rxnsList(1))
    eliList = rxnsList;

    fastSL(core_g{6, j},0.01,2,eliList); 
end
clear j rxnsList i_gr eliList

core_g = FastSLIO(core_g,'R3D_fastcoreM_');

mem = 'fastcoreM';
% searching for SLs using gMCS 
for j = 1:size(core_g,2)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using gMCS 
    model_name = core_g{6, j}.description;
    %model_name = 'Recon2v2';
    n_gMCS = 100;
    max_len_gMCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;

    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    % parpool;
    [gmcs, gMCS_time] = calculateGeneMCS(model_name, core_g{6, j}, n_gMCS, max_len_gMCS, optional_inputs);
    core_g{8, j} = gmcs;
    filePattern    = fullfile('./',strcat('G_',core_g{6, j}.description,'.mat'));
    delete(filePattern);
end
clear mem gmcs gMCS_time filePattern optional_inputs max_len_gMCS n_gMCS model_name

% searching for SLs using MCS 
for j = 1:size(core_g,2)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using MCS 
    %model_name = 'Recon2v2';
    n_MCS = 100;
    max_len_MCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;
    
    rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (~isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList  = [rxnsList  ,core_g{6, j}.rxns{i_gr}];
        end
    end
    optional_inputs.rxn_set = rxnsList';
    
    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    %core_g{6, j}.grRules = core_g{6, j}.rules;
    %core_g{6, j} = creategrRulesField(core_g{6, j});
    %core_g{6, j}.rev = zeros(length(core_g{6, j}.rxns),1);
    %core_g{6, j}.rev(find(core_g{6, j}.lb<0)) = 1;
    %core_g{6, j} = buildRxnGeneMat(core_g{6, j});
    %core_g{6, j}.rxnGeneMat = double(core_g{6, j}.rxnGeneMat);
    % parpool;
    [mcs, mcs_time] = calculateMCS(core_g{6, j}, n_MCS, max_len_MCS, optional_inputs);
    core_g{9, j} = mcs;
end
clear mcs mcs_time optional_inputs i_gr rxnsList j     n_MCS  max_len_MCS    

% to test each model functionalities
[totTest, numTests] = testFunctionality(Recon3D)
gnr_mdl_FBA = optimizeCbModel(Recon3D);

for j = 1:size(core_g,2)
    [totTest, numTests] = testFunctionality(core_g{6, j});
    % percent of functionalities
    core_g{10, j} = totTest/numTests;
    spc_mdl_FBA = optimizeCbModel(core_g{6, j}); 
    % portion of contextualized model fba sim to generic model
    core_g{11, j} = spc_mdl_FBA.f/gnr_mdl_FBA.f;
end
clear j spc_mdl_FBA gnr_mdl_FBA totTest numTests

% thr functionality checking 100 iteration 17 min single core
% I wrote the following lines after a while so load the core strct again,
% if you are doing all sections sequential its not need to load it again
% so do not uncomment
% load('../../MetabolicSLOutput/R3d_fastcore_core_g.mat'); core_g = fastcore_g;     
Gmodel = Recon3D;
%outp = {};
ncor = size(core_g,2);
%M = 5;                     % M specifies maximum number of workers
%parfor  (j = 3:ncor,M)
tmp = {strcat('Recon3D','fastcore')};
tl_msg = [];
for  j = 1:ncor
    % 22 coreg_models
    n_iter = 100;
    CS_model = core_g{6, j};
    [CS_totTest, seq_tst, top_1p, top_50p, Iqr, Max ] = functionality_thr(Gmodel, CS_model, n_iter);
     tl_msg = vertcat(tl_msg , horzcat(CS_totTest, top_1p, top_50p, Iqr, seq_tst));
    %outp = horzcat(outp ,msg);
end

% cell to table conversion
tl_msg = array2table(tl_msg);
tl_msg.Properties.VariableNames([1:4]) = {'NoFuncs' 'top1p' 'top50p' 'IQR' }
writetable(tl_msg,strcat('../../MetabolicSLOutput/func_test_', tmp{1, 1}  ,'.csv'));



% find essential genes for comparison to CRISPR KO experiments
HasEff = [];
for esnGn_i =1:size(core_g,2)
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(core_g{6, esnGn_i});
    HasEff = horzcat(cellstr(core_g{6, esnGn_i}.genes(find(grRatio<.99))),string(grRatio(find(grRatio<.99))));
    core_g{12, esnGn_i} = HasEff;
end
clear HasEff grRatio grRateKO grRateWT hasEffect delRxns fluxSolution esnGn_i ans


% fastSL jdl res purifier
max_jDL = 5000;
FASTCORE_jdl = [];
%generic model has suffix like recon3D
trans_sep  = '_AT';
for modl_i = 1:size(core_g,2)
    output0 = Jsl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jsl,trans_sep);
    output1 = Jdl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jdl,trans_sep);
    output = vertcat(output0,output1);
    [~,idx] = unique(output,'rows','first');
    output = output(sort(idx),:);

    output2 = vertcat(table2array(output) ,zeros(max_jDL - size(output,1),2));
    FASTCORE_jdl = horzcat(FASTCORE_jdl, output2);
end
clear max_jDL  output0 output JDL modl_i max_jDL modl_i JDL output Jdl model output0 trans_sep output output0 output1 output2 ans idx

% gMCS SL res purifier
max_gmcs = 500;
FASTCORE_gmcs = [];
trans_sep  = '_AT';
for modl_i = 1:size(core_g,2)
    gmc = gmcs_purifier(core_g{8, modl_i},trans_sep);
    %gmc.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    output0 = vertcat(table2array(gmc) ,zeros(max_gmcs  - size(gmc,1),2));
    FASTCORE_gmcs = [FASTCORE_gmcs output0];
end
clear max_gmcs output0 gmc modl_i gmc trans_sep

% MCS order 2 res purifier
trans_sep  = '_AT';
max_mcs = 5000;
FASTCORE_mcs = [];
for modl_i = 1:size(core_g,2)
    mcs = mcs_purifier(core_g{6, modl_i}, core_g{9, modl_i},trans_sep);
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    if isempty(mcs)
        output0 = zeros(max_mcs  - size(mcs,1),2);    
    else
        output0 = vertcat(table2array(mcs) ,zeros(max_mcs  - size(mcs,1),2));    
    end
    FASTCORE_mcs = [FASTCORE_mcs output0];
end
clear max_mcs modl_i mcs output0 trans_sep

% find gene list of each model for dd method
% find which model has most genes
max_gls = length(core_g{6, 1}.genes);
for modl_i = 1:size(core_g,2)
    gls_tmp = (length(core_g{6, modl_i}.genes));
    if max_gls < gls_tmp
        max_gls = gls_tmp;
        disp(gls_tmp)
        disp(modl_i)
    end
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
end
clear modl_i gls_tmp 

% find each model gene list
trans_sep = '_AT';
FASTCORE_gls = [];
for modl_i = 1:size(core_g,2)
    FASTCORE_gls_tmp = core_g{6, modl_i}.genes;
    %{
    FASTCORE_gls_tmp = [];
    for gn_i = 1:length(FASTCORE_gls_tmp0)
        tmp = strsplit(FASTCORE_gls_tmp0{gn_i},trans_sep);
        FASTCORE_gls_tmp = vertcat(FASTCORE_gls_tmp, string(tmp{1,1}));
    end
    %}
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    
	output0 = vertcat(FASTCORE_gls_tmp ,num2cell(zeros(max_gls  - length(FASTCORE_gls_tmp),1)));    
    FASTCORE_gls = [FASTCORE_gls output0];
end
clear max_gls trans_sep output0 FASTCORE_gls_tmp FASTCORE_gls_tmp0 tmp output0 gn_i modl_i

% find each model essential gene list
FASTCORE_gls_esn = [];
max_gls_ens = 500;
for modl_i = 1:size(core_g,2)
    FASTCORE_gls_tmp = core_g{12, modl_i};
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
	output0 = vertcat(string([core_g{5, modl_i} core_g{5, modl_i}]),FASTCORE_gls_tmp ,num2cell(zeros(max_gls_ens  - length(FASTCORE_gls_tmp),2)));
    FASTCORE_gls_esn = [FASTCORE_gls_esn output0];
end
clear max_gls_ens modl_i FASTCORE_gls_tmp output0

fastcore_g = core_g;
clear core_g

% fastcore model gene lists write as csv
FASTCORE_glsT = cell2table((FASTCORE_gls));
writetable(FASTCORE_glsT,'../../MetabolicSLOutput/R3d_fastcore_gls.csv'); clear FASTCORE_glsT;
% fastcore model gmcs lists write as csv
FASTCORE_gmcsT = array2table(FASTCORE_gmcs);
writetable(FASTCORE_gmcsT,'../../MetabolicSLOutput/R3d_fastcore_gmcs.csv'); clear FASTCORE_gmcsT;
% fastcore model mcs lists write as csv
FASTCORE_mcsT = array2table(FASTCORE_mcs);
writetable(FASTCORE_mcsT,'../../MetabolicSLOutput/R3d_fastcore_mcs.csv'); clear FASTCORE_mcsT;
% fastcore model jdl lists write as csv
FASTCORE_fastSLT = array2table(FASTCORE_jdl);
writetable(FASTCORE_fastSLT,'../../MetabolicSLOutput/R3d_FASTCORE_fastSL.csv'); clear FASTCORE_fastSLT;
% fastcore model essential gene lists write as csv
FASTCORE_gls_esn = array2table(FASTCORE_gls_esn);
writetable(FASTCORE_gls_esn,'../../MetabolicSLOutput/R3d_fastcore_essen_genes.csv');
% fastcore model percent of functionalities and portion of contextualized model fba sim to generic model res lists write as csv
FASTCORE_funcSim_ls_esn = array2table(fastcore_g([5,10,11], 1:22));
writetable(FASTCORE_funcSim_ls_esn,'../../MetabolicSLOutput/R3d_fastcore_func_sim.csv'); clear FASTCORE_funcSim_ls_esn;
can_ls = fastcore_g(5, :);
writetable(array2table(can_ls),'../../MetabolicSLOutput/R3d_fastcorecanls.csv');
save ../../MetabolicSLOutput/R3d_fastcore_core_g.mat fastcore_g

clear

%%
% ######### chunk 2 Recon3D iMAT #########

s1 = pwd;
s2 = 'D:\G\thesis\thesis\contextualization\MetabolicSL\code';
if ~strcmp(s1,s2)
    cd(s2);
end
clear s1 s2

Recon3D = readCbModel('..\..\MetabolicSLInput\data\Recon3D.xml');
Recon3D = creategrRulesField(Recon3D);

epsilon = getCobraSolverParams('LP', 'feasTol');

% find modelR204 ATP and biomass rxns to add to contextualized models
atp_idx = find(contains(Recon3D.rxnNames,'ATP demand','IgnoreCase',true));
atp_idx1 = find(contains(Recon3D.rxns,'ATPM','IgnoreCase',true));
bio_idx = find(contains(Recon3D.rxns,'biomass','IgnoreCase',true));
addi_idx = vertcat(atp_idx,bio_idx,atp_idx1);
clear atp_idx bio_idx atp_idx1

T = readtable('..\..\MetabolicSLInput\data\R3D_rnaData_entrz.csv');
gene_ls = T(:,1);
T(:,1) = [];
expression.gene  = table2cell(gene_ls); expressionSp.gene  = table2cell(gene_ls) ;
cantyp = table2cell(readtable('..\..\MetabolicSLInput\data\R3D_cantyp.csv'));
ttl_expressionRxns = [];
ttl_expressionRxns_samples = [];
ttl_expressionRxns_all = [];
for exp_i = 1:size(cantyp,1)
    % for each cancer type aggregated samples
    strt = cantyp{exp_i,2};
    endd = cantyp{exp_i,3};
    expression.value =  table2array(T(:,strt:endd));
    expression.value =  str2double(expression.value);
    [expressionRxns parsedGPR] = mapExpressionToReactions(Recon3D, expression); % defualt: minMax
    ttl_expressionRxns = horzcat(ttl_expressionRxns,expressionRxns);
    disp('*** canty reaction mapping ***') 
    disp('remaining:')
    disp(size(cantyp,1) - exp_i)
    
    % for each sample of cancer type
    for exp_j = 1:size(expression.value,2)
        expressionSp.value = expression.value(:,exp_j);
        [expressionRxnsSP parsedGPRSP] = mapExpressionToReactions(Recon3D, expressionSp); % defualt: minMax
        ttl_expressionRxns_samples = horzcat(ttl_expressionRxns_samples, expressionRxnsSP);
        disp(' sample reaction mapping remaining:') 
        disp(size(expression.value,2) - exp_j)
        disp(' out of:') 
        disp(size(expression.value,2))
    end
    ttl_expressionRxns_all = horzcat(ttl_expressionRxns_all , ttl_expressionRxns_samples);
end 
clear exp_i strt endd expression expressionRxns parsedGPR gene_ls T

[AA,orientation,V] = fastcc(Recon3D,epsilon);

%ttl_expressionRxns_all; cantyp{:,1};
% here medium
Recon3D_cons = RPMI_Medium(Recon3D);             f: 755.0032

optimizeCbModel(Recon3D);
optimizeCbModel(Recon3D_cons)
comp_bounds = RPMI_Medium(Recon3D, Recon3D_cons)

%length(find(res))
core_g = struct([]);
core_g_med = struct([]);

options.tol = epsilon;
clear epsilon

options.core = Recon3D.rxns(addi_idx);
options.solver = 'iMAT';

for i = 1:(size(cantyp,1))
    disp('*********') 
    disp('loop')
    disp(i)
    %disp('sub_loop')
    %disp(1)
    endd = i; 
    levels = ttl_expressionRxns(:,endd);
    core_g{5, i} = cantyp(endd,1);
    
    %disp(3)
    %funcModel = 0;
    
    options.expressionRxns = levels;
    options.threshold_ub = prctile(levels(find(levels > -1)),75);
    options.threshold_lb = prctile(levels(find(levels > -1)),25);
       
    core_g{6, i} = createTissueSpecificModel(Recon3D, options);
    core_g{6, i}.description = char(strcat('RD3_iMAT_',core_g{5, i}));
    core_g{3, i} = createTissueSpecificModel(model, options);
    core_g{3, i}.description = char(strcat('RD3_iMAT_con_',core_g{5, i}));
    %disp('atp demand id:')
    %find(contains(core_g{6, i}.rxnNames,'ATP demand','IgnoreCase',true))
    %disp('biomass id:')
    %idx_bio = find(contains(core_g{6, i}.rxnNames,'biomass','IgnoreCase',true));
    %disp(4)
    tmp_fba = optimizeCbModel(core_g{6, i});
    disp('fba simulation')
    disp(tmp_fba.f) 
end
clear levels endd tmp_fba options i ttl_expressionRxns addi_idx cantyp

for modl_i = 1:size(core_g,2)
    disp(modl_i)
   tmp = core_g{6, modl_i};
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
    core_g{4, modl_i} = results;
    clear results tmp
end
clear tmp modl_i LPproblem results ans

%[A,orientation,V] = fastcc(core_g{6, i},epsilon);
%optimizeCbModel(core_g{6, i})

%{%} test fba
for fbi = 1:size(core_g,2)
    atp_id = find(contains(core_g{6, fbi}.rxnNames,'ATP demand','IgnoreCase',true));
    disp(find(contains(core_g{6, fbi}.rxns,'ATPM','IgnoreCase',true)))
    disp(core_g{6, fbi }.rxns(atp_id))
    disp(core_g{6, fbi }.rxns(find(core_g{6, fbi }.c)))
    %fastcc(core_g{6, fbi},epsilon,printLevel);
    aaa = optimizeCbModel(core_g{6,fbi});
    disp(aaa.full(find(core_g{6, fbi }.c)))
    disp(aaa.full(atp_id))
end
clear fbi atp_id aaa atp_id

% to test models functionalities
[totTest, numTests] = testFunctionality(Recon3D)
gnr_mdl_FBA = optimizeCbModel(Recon3D);
gnr_mdl_FBA.f

[totTest_con, numTests_con] = testFunctionality(model)
gnr_mdl_FBA_con = optimizeCbModel(model);
gnr_mdl_FBA_con.f

for j = 1:size(core_g,2)
    [totTest, numTests] = testFunctionality(core_g{6, j});
    % percent of functionalities
    core_g{10, j} = totTest;
    spc_mdl_FBA = optimizeCbModel(core_g{6, j}); 
    core_g{11, j} = spc_mdl_FBA.f;
    
        [totTest, numTests] = testFunctionality(core_g{3, j});
    % percent of functionalities
    core_g{1, j} = totTest;
    spc_mdl_FBA = optimizeCbModel(core_g{3, j}); 
    core_g{2, j} = spc_mdl_FBA.f;
end
clear gnr_mdl_FBA totTest numTests j spc_mdl_FBA

% thr functionality checking 100 iteration 17 min single core
% I wrote the following lines after a while so load the core strct again,
% if you are doing all sections sequential its not need to load it again
% so do not uncomment
% load('../../MetabolicSLOutput/R3d_iMAT_core_g.mat'); core_g = iMAT_core_g; 
Gmodel = Recon3D;
%outp = {};
ncor = size(core_g,2);
%M = 5;                     % M specifies maximum number of workers
%parfor  (j = 3:ncor,M)
tmp = {strcat('Recon3D','iMAT')};
tl_msg = [];
for  j = 1:ncor
    % 22 coreg_models
    n_iter = 100;
    CS_model = core_g{6, j};
    [CS_totTest, seq_tst, top_1p, top_50p, Iqr, Max ] = functionality_thr(Gmodel, CS_model, n_iter);
     tl_msg = vertcat(tl_msg , horzcat(CS_totTest, top_1p, top_50p, Iqr, seq_tst));
    %outp = horzcat(outp ,msg);
end
% cell to table conversion
tl_msg = array2table((tl_msg));
%tl_msg_n = tl_msg;
tl_msg.Properties.VariableNames([1:4]) = {'NoFuncs' 'top1p' 'top50p' 'IQR' }
writetable(tl_msg,strcat('../../MetabolicSLOutput/func_test_', tmp{1, 1}  ,'.csv'));


% do FastSL on 22 can_specific models of core_g in core_g{7, i}
addpath('.\fastSL_Rxns');
%Eliminate Exhange reactions for lethlaity analysis
for j = 1:size(core_g,2)
	rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList = [rxnsList ,core_g{6, j}.rxns{i_gr}];
        end
    end
    rxnsList = rxnsList';
    %findGenesFromRxns(model,rxnsList(1))
    eliList = rxnsList;

    fastSL(core_g{6, j},0.01,2,eliList); 
end
clear j i_gr rxnsList filePattern eliList
% using func FastSLIO which find output of fastSL from HDD
core_g = FastSLIO(core_g,'RD3_iMAT');
save('../../MetabolicSLOutput/RD3iMAT_core_g.mat', 'core_g', '-v7.3')

%{
%Eliminate Exhange reactions for lethlaity analysis
for j = 1:size(core_g,2)
    eliList = core_g{6, j}.rxns(find(findExcRxns(core_g{6, j}))); 
    %fastSL(core_g{6, j},0.01,2,0); 
    fastSL(core_g{6, j},0.01,2,eliList); 
    %core_g{7, j} = load('fastcoreM_ can name ... _Rxn_lethals.mat');
    %disp(core_g{7, j})
end
%}


% searching for SLs using gMCS 
%parpool(4)
%parfor j = 1:size(core_g,2)
for j = 1:size(core_g,2)
    disp('gmcs:')
    disp(j)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using gMCS 
    model_name = core_g{6, j}.description;
     core_g{6, j} = buildRxnGeneMat(core_g{6, j});
    %model_name = 'Recon2v2';
    n_gMCS = 100;
    max_len_gMCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;

    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    % parpool;
    [gmcs, gMCS_time] = calculateGeneMCS(model_name, core_g{6, j}, n_gMCS, max_len_gMCS, optional_inputs);
    core_g{8, j} = gmcs;
    filePattern    = fullfile('./',strcat('G_',core_g{6, j}.description,'.mat'));
    delete(filePattern);
end
clear gmcs gMCS_time optional_inputs max_len_gMCS n_gMCS model_name j

% searching for SLs using MCS 
%parpool(4)
for j = 1:size(core_g,2)
    disp('mcs:')
    disp(j)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using MCS 
    %model_name = 'Recon2v2';
    n_MCS = 100;
    max_len_MCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;
    
    rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (~isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList  = [rxnsList  ,core_g{6, j}.rxns{i_gr}];
        end
    end
    optional_inputs.rxn_set = rxnsList';
    
    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    %core_g{6, j}.grRules = core_g{6, j}.rules;
    %core_g{6, j} = creategrRulesField(core_g{6, j});
    %core_g{6, j}.rev = zeros(length(core_g{6, j}.rxns),1);
    %core_g{6, j}.rev(find(core_g{6, j}.lb<0)) = 1;
    %core_g{6, j} = buildRxnGeneMat(core_g{6, j});
    %core_g{6, j}.rxnGeneMat = double(core_g{6, j}.rxnGeneMat);
    % parpool;
    [mcs, mcs_time] = calculateMCS(core_g{6, j}, n_MCS, max_len_MCS, optional_inputs);
    core_g{9, j} = mcs;
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
trans_sep = '_AT';
max_gmcs = 100;
iMAT_gmcs = [];
for modl_i = 1:size(core_g,2)
    gmc = gmcs_purifier(core_g{8, modl_i},trans_sep);
    %gmc.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    output0 = vertcat(table2array(gmc) ,zeros(max_gmcs  - size(gmc,1),2));
    iMAT_gmcs  = [iMAT_gmcs  output0];
end
clear trans_sep max_gmcs modl_i gmc output0

% fastSL jdl res purifier
trans_sep = '_AT';
max_jDL = 200;
iMAT_jdl = [];
for modl_i = 1:size(core_g,2)

    output0 = Jsl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jsl,trans_sep);
    output1 = Jdl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jdl,trans_sep);
    output = vertcat(output0,output1);
   [~,idx] = unique(output,'rows','first');
    output = output(sort(idx),:);

    output2 = vertcat(table2array(output) ,zeros(max_jDL - size(output,1),2));
    iMAT_jdl  = [iMAT_jdl  output2];

end
clear max_jDL modl_i JDL output Jdl model output0 trans_sep output output0 output1 output2 idx ans filePattern

% MCS order 2 res purifier
trans_sep = '_AT';
max_mcs = 150;
iMAT_mcs = [];
for modl_i = 1:size(core_g,2)
    mcs = mcs_purifier(core_g{6, modl_i}, core_g{9, modl_i},trans_sep);
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
max_gls = 1500;
% find each model gene list
iMAT_gls = [];
for modl_i = 1:size(core_g,2)
    disp(modl_i)
    iMAT_gls_tmp1 = core_g{6, modl_i}.grRules;
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

iMAT_core_g = core_g;
clear core_g
% fastcore model gene lists write as csv
iMAT_glsT = cell2table(cellstr(iMAT_gls));
writetable(iMAT_glsT,'../../MetabolicSLOutput/R3d_iMAT_gls.csv'); clear iMAT_glsT;
% fastcore model gmcs lists write as csv
iMAT_gmcsT = array2table(iMAT_gmcs);
writetable(iMAT_gmcsT,'../../MetabolicSLOutput/R3d_iMAT_gmcs.csv'); clear iMAT_gmcsT;
% fastcore model mcs lists write as csv
iMAT_mcsT = array2table(iMAT_mcs);
writetable(iMAT_mcsT,'../../MetabolicSLOutput/R3d_iMAT_mcs.csv'); clear iMAT_mcsT;
% fastcore model jdl lists write as csv
iMAT_fastSLT = array2table(iMAT_jdl);
writetable(iMAT_fastSLT,'../../MetabolicSLOutput/R3d_iMAT_fastSL.csv'); clear iMAT_fastSLT;
% fastcore model essential gene lists write as csv
iMAT_gls_esn = array2table(iMAT_gls_esn);
writetable(iMAT_gls_esn,'../../MetabolicSLOutput/R3d_iMAT_essen_genes.csv'); clear iMAT_gls_esn;
% fastcore model percent of functionalities and portion of contextualized model fba sim to generic model res lists write as csv
iMAT_funcSim_ls_esn = array2table(iMAT_core_g([5,10,11], 1:22));
writetable(iMAT_funcSim_ls_esn,'../../MetabolicSLOutput/R3d_iMAT_func_sim.csv'); clear iMAT_funcSim_ls_esn;
ab = iMAT_core_g(5, :);
writetable(array2table(ab),'../../MetabolicSLOutput/R3d_iMATcanls.csv'); clear ab

save('../../MetabolicSLOutput/R3d_iMAT_core_g.mat ', 'iMAT_core_g', '-v7.3')

% percent of functionalities core_g{10, j}
% fraction of wild generic model fba sol core_g{11, j}
% essential genes core_g{12, j}
% '../../MetabolicSLOutput/fastcore_essen_genes.csv'

%%
% ######### chunk 3 Recon2.04 fastcore #########

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

% using fastcc the flux consistent rxns will be determined for next
% sessions
[A,orientation,V] = fastcc(modelR204,epsilon);
clear orientation V

% fastcore models reconstruction
options.solver = 'fastCore';

% reading cancer core gene list for 22 cancer of fastcore 
T = readtable('..\..\MetabolicSLInput\data\cancer_core_genes_entrz.csv');
T(:,1) = [];

%T = readtable('cancer_core_genes_entrz.csv');
%T(:,1) = [];
%T = T1;

% making cell struct which contains all models of fastcore and next step
% will save at lower elements row of the struct
for i = 1:size(T,2)
    disp('loop')
    disp(i)
    disp('sub_loop')
    disp(1)
    %i = 1
    %T(:,i).Properties.VariableNames
    T0 = T{:,i};
    disp(2)
    TF = {T(find(~ismissing(T0,{'NA'})),i)};
    disp(3)
    core_g(1,i) = TF;
    disp(4)
    core_g(2,i) = {ones(size(TF{1,1},1),1)};
    disp(5)
    expression.gene  = string(table2cell(core_g{1, i}));
    disp(6)
    expression.value = core_g{2, i};
    disp(7)
    [expressionRxns parsedGPR] = mapExpressionToReactions(modelR204, expression); % defualt: minMax
    disp(8)
    core_g{3, i} = find(expressionRxns>0);
    disp(9)
    core_g{4, i} = intersect((core_g{3, i}),A);
    % adding biomass and atp demand
    core_g{4, i} = vertcat(core_g{4, i},addi_idx);
    disp(9.1)
    core_g{5, i} = T(:,i).Properties.VariableNames;
    disp(10)
	%options.expressionCol = expressionRxns;
    %expressionCol = expressionRxns;
    disp(11)
    tol = 1e-4;
    %scaling = 1;
    options.epsilon = tol;  
    options.core = core_g{4,i};  
    %options.core = (10:30)';  
    options.printLevel = 2;
    options.solver = 'fastCore';
    core_g{6, i} = createTissueSpecificModel(modelR204, options);
    core_g{6, i}.description = char(strcat('R24_fastcoreM_',core_g{5, i}));
    %{
    core_g{6, i} = creategrRulesField(core_g{6, i});
    core_g{6, i}.rev = zeros(length(core_g{6, i}.rxns),1);
    core_g{6, i}.rev(find(core_g{6, i}.lb<0)) = 1;
    core_g{6, i} = buildRxnGeneMat(core_g{6, i});
    core_g{6, i}.rxnGeneMat = double(core_g{6, i}.rxnGeneMat);
    %}
    disp(12)
    %core_g{7, i} = fastcc(core_g{6, i},epsilon);
    %disp(13)
end
clear i  options tol expressionRxns TF T0 parsedGPR A addi_idx T ans expression epsilon


%{%}
cntr = 0;
for fbi = 1:size(core_g,2)
    %fastcc(core_g{6, fbi},epsilon,printLevel);
    disp('loop:')
    disp(fbi)
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxnNames,'ATP demand','IgnoreCase',true))))
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxns,'ATPM','IgnoreCase',true))))
    disp(core_g{6,fbi}.rxns(find(contains(core_g{6,fbi}.rxns,'biomass','IgnoreCase',true))))
    aaa = optimizeCbModel(core_g{6,fbi});
    disp(core_g{6, fbi }.description)
    disp(core_g{6, fbi }.rxns(find(core_g{6, fbi }.c)))
    disp(aaa.f)
    %disp(aaa.v(find(contains(core_g{6, fbi}.rxnNames,'biomass','IgnoreCase',true))))
    idx = find(contains(core_g{6, fbi}.rxnNames,'ATP demand','IgnoreCase',true));
    disp(core_g{6, fbi }.rxns(idx))
    disp(aaa.v(idx ))
    if aaa.f > 0
        cntr = cntr + 1;
    end
    cntr
end
clear cntr aaa idx ans fbi

% load('../../MetabolicSLOutput/R24_fastcore_core_g.mat');
% core_g = fastcore_g;
% clear fastcore_g

% do FastSL on 22 can_specific models of core_g in core_g{7, i}
addpath('.\fastSL_Rxns');
%Eliminate Exhange reactions for lethlaity analysis
for j = 1:size(core_g,2)
	rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList = [rxnsList ,core_g{6, j}.rxns{i_gr}];
        end
    end
    rxnsList = rxnsList';
    %findGenesFromRxns(model,rxnsList(1))
    eliList = rxnsList;

    fastSL(core_g{6, j},0.01,2,eliList); 
    	
end
clear j rxnsList i_gr eliList

core_g = FastSLIO(core_g,'R24_fastcoreM_');

mem = 'fastcoreM';
% searching for SLs using gMCS 
for j = 1:size(core_g,2)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using gMCS 
    model_name = core_g{6, j}.description;
    %model_name = 'Recon2v2';
    n_gMCS = 100;
    max_len_gMCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;

    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    % parpool;
    [gmcs, gMCS_time] = calculateGeneMCS(model_name, core_g{6, j}, n_gMCS, max_len_gMCS, optional_inputs);
    core_g{8, j} = gmcs;
    filePattern    = fullfile('./',strcat('G_',core_g{6, j}.description,'.mat'));
    delete(filePattern);
end
clear mem gmcs gMCS_time filePattern optional_inputs max_len_gMCS n_gMCS model_name

% searching for SLs using MCS 
for j = 1:size(core_g,2)
    % searching for SLs using gMCS and MCS algorithms
    % searching for SLs using MCS 
    %model_name = 'Recon2v2';
    n_MCS = 100;
    max_len_MCS = 2;
    %optional_inputs.KO = '6240';
    %optional_inputs.separate_transcript = '.';
    optional_inputs.timelimit = 4*60;
    optional_inputs.numWorkers = 6;
    
    rxnsList = {};
    for i_gr = 1:length(core_g{6, j}.rules)
        if (~isempty(core_g{6, j}.rules{i_gr,1}))
            rxnsList  = [rxnsList  ,core_g{6, j}.rxns{i_gr}];
        end
    end
    optional_inputs.rxn_set = rxnsList';
    
    % Finally, we will calculate the gMCSs. Turning the parallel pool on is 
    % recommended if the G matrix has not been calculated yet.
    %core_g{6, j}.grRules = core_g{6, j}.rules;
    %core_g{6, j} = creategrRulesField(core_g{6, j});
    %core_g{6, j}.rev = zeros(length(core_g{6, j}.rxns),1);
    %core_g{6, j}.rev(find(core_g{6, j}.lb<0)) = 1;
    %core_g{6, j} = buildRxnGeneMat(core_g{6, j});
    %core_g{6, j}.rxnGeneMat = double(core_g{6, j}.rxnGeneMat);
    % parpool;
    [mcs, mcs_time] = calculateMCS(core_g{6, j}, n_MCS, max_len_MCS, optional_inputs);
    core_g{9, j} = mcs;
end
clear mcs mcs_time optional_inputs i_gr rxnsList j     n_MCS  max_len_MCS    

% to test each model functionalities
[totTest, numTests] = testFunctionality(modelR204)
gnr_mdl_FBA = optimizeCbModel(modelR204);

for j = 1:size(core_g,2)
    [totTest, numTests] = testFunctionality(core_g{6, j});
    % percent of functionalities
    core_g{10, j} = totTest/numTests;
    spc_mdl_FBA = optimizeCbModel(core_g{6, j}); 
    % portion of contextualized model fba sim to generic model
    core_g{11, j} = spc_mdl_FBA.f/gnr_mdl_FBA.f;
end
clear j spc_mdl_FBA gnr_mdl_FBA totTest numTests

% thr functionality checking 100 iteration 17 min single core
% I wrote the following lines after a while so load the core strct again,
% if you are doing all sections sequential its not need to load it again
% so do not uncomment
% load('../../MetabolicSLOutput/R24_fastcore_core_g.mat'); core_g = fastcore_g;
Gmodel = modelR204;
%outp = {};
ncor = size(core_g,2);
%M = 5;                     % M specifies maximum number of workers
%parfor  (j = 3:ncor,M)
tmp = {strcat('Recon24_','fastcore')};
tl_msg = [];
for  j = 1:ncor
    % 22 coreg_models
    n_iter = 100;
    CS_model = core_g{6, j};
    [CS_totTest, seq_tst, top_1p, top_50p, Iqr, Max] = functionality_thr(Gmodel, CS_model, n_iter);
    tl_msg = vertcat(tl_msg , horzcat(CS_totTest, top_1p, top_50p, Iqr, seq_tst));
    %outp = horzcat(outp ,msg);
end
% cell to table conversion
tl_msg = array2table((tl_msg));
tl_msg.Properties.VariableNames([1:4]) = {'NoFuncs' 'top1p' 'top50p' 'IQR' }
writetable(tl_msg,strcat('../../MetabolicSLOutput/func_test_', tmp{1, 1}  ,'.csv'));


% find essential genes for comparison to CRISPR KO experiments
HasEff = [];
for esnGn_i =1:size(core_g,2)
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(core_g{6, esnGn_i});
    HasEff = horzcat(cellstr(core_g{6, esnGn_i}.genes(find(grRatio<.99))),string(grRatio(find(grRatio<.99))));
    core_g{12, esnGn_i} = HasEff;
end
clear HasEff grRatio grRateKO grRateWT hasEffect delRxns fluxSolution esnGn_i ans


% fastSL jdl res purifier
max_jDL = 7000;
FASTCORE_jdl = [];
%generic model has suffix like modelR204
trans_sep  = '.';
for modl_i = 1:size(core_g,2)
    output0 = Jsl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jsl,trans_sep);
    output1 = Jdl_purifier(core_g{6, modl_i}, core_g{7, modl_i}.Jdl,trans_sep);
    output = vertcat(output0,output1);
    [~,idx] = unique(output,'rows','first');
    output = output(sort(idx),:);

    output2 = vertcat(table2array(output) ,zeros(max_jDL - size(output,1),2));
    FASTCORE_jdl = horzcat(FASTCORE_jdl, output2);
end
clear max_jDL  output0 output JDL modl_i max_jDL modl_i JDL output Jdl model output0 trans_sep output output0 output1 output2 ans idx

% gMCS SL res purifier
max_gmcs = 500;
FASTCORE_gmcs = [];
trans_sep  = '.';
for modl_i = 1:size(core_g,2)
    gmc = gmcs_purifier(core_g{8, modl_i},trans_sep);
    %gmc.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    output0 = vertcat(table2array(gmc) ,zeros(max_gmcs  - size(gmc,1),2));
    FASTCORE_gmcs = [FASTCORE_gmcs output0];
end
clear max_gmcs output0 gmc modl_i gmc trans_sep

% MCS order 2 res purifier
trans_sep  = '.';
max_mcs = 5000;
FASTCORE_mcs = [];
for modl_i = 1:size(core_g,2)
    mcs = mcs_purifier(core_g{6, modl_i}, core_g{9, modl_i},trans_sep);
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    if isempty(mcs)
        output0 = zeros(max_mcs  - size(mcs,1),2);    
    else
        output0 = vertcat(table2array(mcs) ,zeros(max_mcs  - size(mcs,1),2));    
    end
    FASTCORE_mcs = [FASTCORE_mcs output0];
end
clear max_mcs modl_i mcs output0 trans_sep

% find gene list of each model for dd method
% find which model has most genes
max_gls = length(core_g{6, 1}.genes);
for modl_i = 1:size(core_g,2)
    gls_tmp = (length(core_g{6, modl_i}.genes));
    if max_gls < gls_tmp
        max_gls = gls_tmp;
        disp(gls_tmp)
        disp(modl_i)
    end
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
end
clear modl_i gls_tmp 

% find each model gene list
trans_sep = '.';
FASTCORE_gls = [];
for modl_i = 1:size(core_g,2)
    FASTCORE_gls_tmp = core_g{6, modl_i}.genes;
    %{
    FASTCORE_gls_tmp = [];
    for gn_i = 1:length(FASTCORE_gls_tmp0)
        tmp = strsplit(FASTCORE_gls_tmp0{gn_i},trans_sep);
        FASTCORE_gls_tmp = vertcat(FASTCORE_gls_tmp, string(tmp{1,1}));
    end
    %}
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
    
	output0 = vertcat(FASTCORE_gls_tmp ,num2cell(zeros(max_gls  - length(FASTCORE_gls_tmp),1)));    
    FASTCORE_gls = [FASTCORE_gls output0];
end
clear max_gls trans_sep output0 FASTCORE_gls_tmp FASTCORE_gls_tmp0 tmp output0 gn_i modl_i

% find each model essential gene list
FASTCORE_gls_esn = [];
max_gls_ens = 500;
for modl_i = 1:size(core_g,2)
    FASTCORE_gls_tmp = core_g{12, modl_i};
    %mcs.Properties.VariableNames{1} = FASTCORE_core_g{6, modl_i}.description;
	output0 = vertcat(string([core_g{5, modl_i} core_g{5, modl_i}]),FASTCORE_gls_tmp ,num2cell(zeros(max_gls_ens  - length(FASTCORE_gls_tmp),2)));
    FASTCORE_gls_esn = [FASTCORE_gls_esn output0];
end
clear max_gls_ens modl_i FASTCORE_gls_tmp output0

fastcore_g = core_g;
clear core_g

% fastcore model gene lists write as csv
FASTCORE_glsT = cell2table((FASTCORE_gls));
writetable(FASTCORE_glsT,'../../MetabolicSLOutput/R24_fastcore_gls.csv'); clear FASTCORE_glsT;
% fastcore model gmcs lists write as csv
FASTCORE_gmcsT = array2table(FASTCORE_gmcs);
writetable(FASTCORE_gmcsT,'../../MetabolicSLOutput/R24_fastcore_gmcs.csv'); clear FASTCORE_gmcsT;
% fastcore model mcs lists write as csv
FASTCORE_mcsT = array2table(FASTCORE_mcs);
writetable(FASTCORE_mcsT,'../../MetabolicSLOutput/R24_fastcore_mcs.csv'); clear FASTCORE_mcsT;
% fastcore model jdl lists write as csv
FASTCORE_fastSLT = array2table(FASTCORE_jdl);
writetable(FASTCORE_fastSLT,'../../MetabolicSLOutput/R24_FASTCORE_fastSL.csv'); clear FASTCORE_fastSLT;
% fastcore model essential gene lists write as csv
FASTCORE_gls_esn = array2table(FASTCORE_gls_esn);
writetable(FASTCORE_gls_esn,'../../MetabolicSLOutput/R24_fastcore_essen_genes.csv');
% fastcore model percent of functionalities and portion of contextualized model fba sim to generic model res lists write as csv
FASTCORE_funcSim_ls_esn = array2table(fastcore_g([5,10,11], 1:22));
writetable(FASTCORE_funcSim_ls_esn,'../../MetabolicSLOutput/R24_fastcore_func_sim.csv'); clear FASTCORE_funcSim_ls_esn;
can_ls = fastcore_g(5, :);
writetable(array2table(can_ls),'../../MetabolicSLOutput/R24_fastcorecanls.csv');

save('../../MetabolicSLOutput/R24_fastcore_core_g.mat', 'fastcore_g', '-v7.3')

clear

%%
% ######### chunk 4 Recon2.04 iMAT #########

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
% functionality checks mets conversion to metNames based on Recon2.0v4
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


%%
% SL idea Recon2.v4 iMAT
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

% preparing each cell line gene list

for i=1:size(all_SL, 2)
    
end

tb_ko_res = array2table(ko_res);
tb_ko_res.Properties.VariableNames(1:35) = cel_all(1:35, 1);
writetable(tb_ko_res,'../../MetabolicSLOutput/KO_res_bile_celllines_35.csv');

tb_ko_gene_ls = array2table(gene_ls{:,:});
writetable(tb_ko_gene_ls ,'../../MetabolicSLOutput/KO_res_bile_gene_ls.csv');

tb_gen_cell_KO = array2table(gen_cell_KO);
tb_gen_cell_KO.Properties.VariableNames(1:35) = cel_all(1:35, 1);
writetable(tb_gen_cell_KO ,'../../MetabolicSLOutput/cellln_bile_gen_ls.csv');
%%
% SL idea Recon2.v4 iMAT
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
% random biomass for our method
% at first I need to combined two parts of all_SL files cause they contain
% the parameters of ub and lb for imat for model construction without grid
% search
load('../../MetabolicSLOutput/all_SL_save_446_ttestV_cons.mat');
% check if the foler containing results of random biomass cell lines KO
% experiments exist. if not make it
folder_path = '../../MetabolicSLOutput/can_KO'; % Example folder path
if (~exist(folder_path, 'dir')) % Check if folder doesn't exist
    mkdir(folder_path); % Create folder
end
biomass_id = find(contains(modelR204.rxns, 'Biomass_reaction','IgnoreCase',true));
precursor_idx = findReactionPrecursors(modelR204, biomass_id);
precursors_values = full(modelR204.S(precursor_idx, biomass_id));
precursors_values_median = median(full(modelR204.S(precursor_idx, biomass_id)));

biomass_size = (size(precursor_idx,1));

parameter_df = all_SL; expression_table = T1;
for can_i=1:size(cantyp,1)
    disp('*********************************************')
    disp(cantyp{can_i, 2})
    can = (cantyp{can_i, 2});
    starts = cantyp{can_i, 3}-1;
    ends = cantyp{can_i, 4}-1;
    % loop over each cell line of each cancer if its WT has a non zero biomass for the geniun biomass then
    % go for its 100 random biomass experiment
    parameter_df_tmp = parameter_df(1, starts:ends); expression_table_tmp = expression_table(:, starts:ends);
    numcol = size(parameter_df_tmp,2);
    % in first loop lets construct all models then find number of overlaps
    % mets then find a random set of mets in overlap set and continouing
    % the process, do last step 100 times and save each cancer 
    % grid search and put the results of ub and lb in 13 and 14 
    %%% cell lines model construction
    for cel_i =1:numcol
        disp(cel_i);
        clear expression.value levels
        expression.value =  table2array(expression_table_tmp(:,cel_i));
        expression.value =  str2double(expression.value);
        [levels parsedGPRSP] = mapExpressionToReactions(modelR204_cons, expression); % defualt: minMax
        options.expressionRxns = levels;
        res = [];
        for ub_i = 1:10
        for lb_i = 5:19
                %disp(((20 - ub_i) * 5));
                %disp(((20 - lb_i) * 5));
                options.threshold_ub = prctile(levels(find(levels > -1)),(20 - ub_i) * 5);
                options.threshold_lb = prctile(levels(find(levels > -1)),(20 - lb_i) * 5);
                clear tmp_fba tmpl
                tmpl = createTissueSpecificModel(modelR204_cons, options);
                % ckeck for biomass and atp demand rxns activity               
                %id_DM_datp = find(contains(tmpl.rxns, options.core{2},'IgnoreCase',true));
                tmp_fba= optimizeCbModel(tmpl);
                id1 = find(contains(tmpl.rxns, options.core{2},'IgnoreCase',true));
                %precursors = findReactionPrecursors(tmpl, rxn_id);
                %[prerxn, prerxn_id] = findReactionrxnPrecursors(model, reaction_id);
                %res{ub_i,lb_i} = tmp_fba.f;
                %if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.v([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
                if (tmp_fba.f > (gen_cons_FBA.f/10)  & length(tmpl.rxns) < length(modelR204_cons.rxns))
                    break % break lb_i loop
                end
            end
             %if (tmp_fba.f > (gen_cons_FBA.f/3) & tmp_fba.v([id1]) ~= 0 & length(tmpl.rxns) < length(modelR204_cons.rxns))
             if (tmp_fba.f > (gen_cons_FBA.f/10)  & length(tmpl.rxns) < length(modelR204_cons.rxns))
               %all_SL{1, i} = char(strcat( cel_all{i}));
                disp((20 - ub_i) * 5);
                disp((20 - lb_i) * 5);
                parameter_df_tmp{2, cel_i} = tmpl;
                parameter_df_tmp{3, cel_i} = tmp_fba.f;
                parameter_df_tmp{10, cel_i} = ((20 - ub_i) * 5);
                parameter_df_tmp{11, cel_i} = ((20 - lb_i) * 5);
                break % break ub_i loop
            end
        end
    end
    %%% geniune biomass model construction
    ones_matrix = ones(length(gene_ls.Var2),numcol); % Create a matrix of ones
    ko_biomass = num2cell(ones_matrix);
    %%% making ko matrix of genes ko experiments for geniune biomass
    %%% models
    for cel_i = 1:numcol
        % here for the KO expriments of WT biomass models
        if ~eq(parameter_df_tmp{3, cel_i}, 0)
                        res_fba = optimizeCbModel(parameter_df_tmp{2, cel_i}, 'max', 'one');
            % find candid gens for KO
            rxn_id = find(res_fba.v);
            genesList = {};
            for i_gr = 1:length(rxn_id)
                tmp = rxn_id(i_gr);
                if (~isempty(parameter_df_tmp{2, cel_i}.rules{tmp ,1}))
                    genesList = [genesList, parameter_df_tmp{2, cel_i}.genes{find(parameter_df_tmp{2, cel_i}.rxnGeneMat(tmp,:))}];
                end
            end
            clear del_gen
            del_gen = genesList';
            del_gen = unique(del_gen);
            [genels_id, del_geid] = ismember(gene_ls.Var2, del_gen);
            del_gen = gene_ls.Var2(genels_id);
            clear rxn_id tmp genesList genels_id del_geid
            % del_gen list KO expriments
            for d =1:length(del_gen)
                    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(parameter_df_tmp{2, cel_i}, 'FBA', del_gen(d));
                    ko_gid = find(ismember(gene_ls.Var2, del_gen(d)));
                    ko_biomass{ko_gid,cel_i} = grRatio;
            end  
        else
             ko_biomass(:,cel_i) = {NaN};
        end
    end
    num_elements = ~cellfun(@isempty, parameter_df_tmp(3, 1:numcol));
    tb_ko_res0 = array2table(ko_biomass(:,num_elements));
    tb_ko_res0.Properties.VariableNames = parameter_df_tmp(1, num_elements);
    filename = sprintf('../../MetabolicSLOutput/can_KOthr10/ko_thr10geninue_biomass__%s.csv',can); % Create file name using sprintf function
    writetable(tb_ko_res0, filename);  
%{
    %%% random biomass model construction
    threshold_notNat=100;
    n_iteration =1000;
    not_NaN_counter = 0;
    numcol = size(parameter_df_tmp,2);
    num_elements = ~cellfun(@isempty, parameter_df_tmp(3, 1:numcol));
    biomass_gen_val = sum(num_elements);
    list1 = parameter_df_tmp(3, num_elements);
    list1 =  cell2mat(list1');
    % 100 NaN excluded random reps
    for rep_i = 1:n_iteration
        ones_matrix = ones(length(gene_ls.Var2),numcol); % Create a matrix of ones
        ko_rnd_biomass = num2cell(ones_matrix);
        model = parameter_df_tmp{2, 1};
        biomass_id = find(ismember(model.rxns, 'biomass_reaction'));
        biomass_pre_ids = find(model.S(:, biomass_id));
        biomass_coef = model.S(biomass_pre_ids, biomass_id);
        biomass_size = size(biomass_pre_ids,1);
        shuffled_biomass_pre = biomass_pre_ids(randperm(biomass_size));
        shuffled_biomass_pre = model.mets(shuffled_biomass_pre);        
        biomass_coef = biomass_coef(randperm(size(biomass_coef, 1)));
        for cel_i = 1:numcol
            if ~eq(parameter_df_tmp{3, cel_i}, 0)
                    model = parameter_df_tmp{2, cel_i};
                    biomass_id = find(ismember(model.rxns, 'biomass_reaction'));
                    biomass_pre_ids = find(ismember(model.mets, shuffled_biomass_pre));
                    model.S(biomass_pre_ids(1:length(biomass_coef)), biomass_id) = biomass_coef;
                    % store random biomassed models in 12th row
                    parameter_df_tmp{12, cel_i} = model;
                    try
                        res = optimizeCbModel(model,'max');
                    catch ME
                        %rethrow(ME)
                    end
                    % if any of random biomass model has fba sol of 0 then the
                    % loop of the cel_i will break
                    if ~eq(res.f, 0)
                        parameter_df_tmp{(rep_i+12), cel_i} = res.f;
                    else
                        parameter_df_tmp{(rep_i+12), cel_i} = NaN;
                    end
                end
            end
            % if any of random biomass model has fba sol of 0 then the
            % loop of the cel_i will break

        % check for #repeats which are NaN less; if the #random repeat
        % meets the threshold_notNat=100 then update n_iteration value
        num_elements_ls2 = ~cellfun(@isnan, parameter_df_tmp((rep_i +12) , num_elements));        
        if (sum(num_elements_ls2) == sum(num_elements))
            list2 = parameter_df_tmp(rep_i +12, num_elements);
            list2 =  cell2mat(list2');
            [rho, pval] = corr(list1, list2, 'Type', 'Pearson', 'Rows','pairwise');
            parameter_df_tmp{rep_i +12, (numcol+1)} = rho;
            parameter_df_tmp{rep_i +12, (numcol+2)} = pval;
            not_NaN_counter = not_NaN_counter + 1;
            disp(not_NaN_counter);
            %%% making ko matrix of genes ko experiments for random biomass
            %%% models
            for cel_i = 1:numcol
                if ~isempty(parameter_df_tmp{3, cel_i})
                % here for the KO expriments of WT biomass and random 100
                % biomass models
                res_random_fba = optimizeCbModel(parameter_df_tmp{12, cel_i}, 'max', 'one');
                % find candid gens for KO
                rxn_id = find(res_random_fba.v);
                genesList = {};
                for i_gr = 1:length(rxn_id)
                    tmp = rxn_id(i_gr);
                    if (~isempty(parameter_df_tmp{12, cel_i}.rules{tmp ,1}))
                        genesList = [genesList, parameter_df_tmp{12, cel_i}.genes{find(parameter_df_tmp{12, cel_i}.rxnGeneMat(tmp,:))}];
                    end
                end
                clear del_gen
                del_gen = genesList';
                del_gen = unique(del_gen);
                [genels_id, del_geid] = ismember(gene_ls.Var2, del_gen);
                del_gen = gene_ls.Var2(genels_id);
                clear rxn_id tmp genesList genels_id del_geid
                % del_gen list KO expriments
                for d =1:length(del_gen)
                        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(parameter_df_tmp{12, cel_i}, 'FBA', del_gen(d));
                        ko_gid = find(ismember(gene_ls.Var2, del_gen(d)));
                        ko_rnd_biomass{ko_gid,cel_i} = grRatio;
                end  
            end
            end
            tb_ko_res0= [];
            tb_ko_res0 = array2table(ko_rnd_biomass(:,num_elements));
            tb_ko_res0.Properties.VariableNames = parameter_df_tmp(1, num_elements);
            filename = sprintf('../../MetabolicSLOutput/random_biomass_KO/ko_rnd_biomass__%s_repeat%d.csv',can, not_NaN_counter); % Create file name using sprintf function
            writetable(tb_ko_res0, filename);
            if not_NaN_counter == threshold_notNat
                break;
            end        
        end
    end
    %}
end

