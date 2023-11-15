%% report fluxSum
clc, clearvars -except solverOK

load('consistent_model.mat')
model=consistent_model
% FBAsolution=optimizeCbModel(model,'max','zero')

whichData=5 %1:FBA, 2-3:FVA, 4-6:sampling median/mean/std
data=[];
disp('... model & data loading 0% ...')
[NUM,TXT,RAW]=xlsread('stats_model_7_ctrl'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_7_sc1'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_7_sc12'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_7_sc2'); data=[data, NUM(:,whichData)];
disp('... model & data loading 50% ...')
[NUM,TXT,RAW]=xlsread('stats_model_H_ctrl'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_H_sc1'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_H_sc12'); data=[data, NUM(:,whichData)];
[NUM,TXT,RAW]=xlsread('stats_model_H_sc2'); data=[data, NUM(:,whichData)];

glcUptake=data(find(ismember(model.rxns,'EX_glc_D[e]')),:)

disp('... model & data loading done ...')

%% report flux sum per selected metabolites
delete fluxSum.txt
diary fluxSum.txt
diary off

toComps=1:8 %8 %FBAsolutions / Models
nrPathways=1:10 %6


for counterP=1:numel(nrPathways) %for all pathways
    switch nrPathways(counterP)
        case 1
            pathway='Glycolysis'
            metList={'glc_D[c]','g6p[c]','fdp[c]','g3p[c]','dhap[c]','pep[c]','pyr[c]','lac_L[e]'}
        case 2
            pathway='PPP'
            metList={'ru5p_D[c]','s7p[c]','e4p[c]'}
        case 3
            pathway='TCA (cytoplasm only)'
            metList={'cit[c]','icit[c]','akg[c]','succ[c]','fum[c]','oaa[c]'}
        case 4
            pathway='TCA+oxphos (mitochondrial)'
            metList={'cit[m]','icit[m]','akg[m]','succ[m]','succoa[m]','fum[m]','mal_L[m]','oaa[m]','focytC[m]','nadh[m]','fadh2[m]','q10h2[m]','o2[m]'},
        case 5
            pathway='Pyruvate metabolism'
            metList={'mthgxl[e]','mthgxl[c]','g3p[c]','dhap[c]','pyr[c]','aact[c]','acetol[c]','lald_D[c]','lald_L[c]','lgt_S[c]'}
        case 6
            pathway='Glutathione metabolism'
            temp=find(ismember(model.subSystems,pathway));
            metList=findMetsFromRxns(model,model.rxns(temp))'
        case 7
            pathway='Oxidative phosphorylation'
            temp=find(ismember(model.subSystems,pathway));
            metList=findMetsFromRxns(model,model.rxns(temp))'
        case 8
            pathway='ROS detoxification'
            temp=find(ismember(model.subSystems,pathway));
            metList=findMetsFromRxns(model,model.rxns(temp))'
        case 9
            pathway='Fatty acid oxidation'
            metList={'accoa[c]','accoa[m]'}
        case 10
            pathway='Cholesterol (r=ER)'
            metList={'mev_R[c]','5pmev[c]','sql[r]','zymst[r]','chsterol[r]','chsterol[c]'}
    end
    
    resAllMets=nan(numel(toComps),numel(metList));
    for counterC=1:numel(toComps) %for all FBA solutions
        toComp=toComps(counterC)
        
        %% calculate flux sum per metabolite
        %         model=model1;
        v=data(:,toComp);
        temp=repmat(v',size(model.S,1),1);
        fluxes=model.S.*temp;
        fluxSumP=full(sum((fluxes>0).*fluxes,2));
        fluxSumN=full(sum((fluxes<0).*fluxes,2));
        temp=[fluxSumP, fluxSumN];
        
        resAll=[];
        for counterM=1:numel(metList) %for all listed metabolites
            temp2=find(ismember(model.mets,cell2mat(metList(counterM))));
            res=temp(temp2);
            resAll=[resAll; res];
            resAllMets(counterC,counterM)=res;
        end
        resAllt=table(resAll,'RowNames',metList)
        
    end
    
    % report to diary file
    diary on
    pathway
    disp('Flux Sum of the following metabolites:')
    disp(metList)
    disp('Rows:Models; Columns:Metabolites')
    resAllMets
    diary off

end
