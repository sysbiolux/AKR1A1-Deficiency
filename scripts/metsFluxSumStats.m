clearvars -except solverOK, clc, close all
load('../consistent_model.mat')
model_orig=consistent_model;
% model=model_orig;

epsilon=1e-4
nrA=1
nrB=4
cutoffFluxDiff=100
metsOfInterest=model_orig.mets;

metsOfInterest={};
pathway='Glycolysis/gluconeogenesis'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='Pentose phosphate pathway'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='Glutathione metabolism'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='Pyruvate metabolism'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='Oxidative phosphorylation'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='ROS detoxification'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);
pathway='Citric acid cycle'
metList=findMetsFromRxns(model_orig,model_orig.rxns(find(ismember(model_orig.subSystems,pathway))))';
metsOfInterest=union(metsOfInterest,metList);

metsOfInterest=find(ismember(model_orig.mets,metsOfInterest));

%% loading
file=['SamplingResults_medium_1500_model_' num2str(nrA) '.mat'];
load(file)
model1=x.modelSampling
data1=x.samples(:,1:500);
size(data1)
file=['SamplingResults_medium_1500_model_' num2str(nrB) '.mat'];
load(file)
model2=x.modelSampling
data2=x.samples(:,1:500);
size(data2)

disp('... model & data loading done ...')

%% calculate flux sum per metabolite
res1=[];
for counter=1:size(data1,2)
    v=data1(:,counter);
    temp=repmat(v',size(model1.S,1),1);
    fluxes=model1.S.*temp;
    fluxSumP=full(sum((fluxes>0).*fluxes,2));
    fluxSumN=full(sum((fluxes<0).*fluxes,2));
    res1=[res1, fluxSumP];
end
disp('... fluxSum A calculated ...')

res2=[];
for counter=1:size(data2,2)
    v=data2(:,counter);
    temp=repmat(v',size(model2.S,1),1);
    fluxes=model2.S.*temp;
    fluxSumP=full(sum((fluxes>0).*fluxes,2));
    fluxSumN=full(sum((fluxes<0).*fluxes,2));
    res2=[res2, fluxSumP];
end
disp('... fluxSum B calculated ...')

figure
boxplot(res1(1:30,:)','Labels',model1.mets(1:30))
set(gca,'FontSize',10,'XTickLabelRotation',45)

mean(res1(find(ismember(model1.mets,'glc_D[c]')),:))
mean(res2(find(ismember(model2.mets,'glc_D[c]')),:))

%% mapping to original model
[C,IA,IB] = intersect(model_orig.mets,model1.mets,'stable');
resA=zeros(numel(model_orig.mets),500);
resA(IA,:)=res1;

[C,IA,IB] = intersect(model_orig.mets,model2.mets,'stable');
resB=zeros(numel(model_orig.mets),500);
resB(IA,:)=res2;

mean(resA(find(ismember(model_orig.mets,'glc_D[c]')),:))
mean(resB(find(ismember(model_orig.mets,'glc_D[c]')),:))

%% statistical test
stats=[];
for counter=1:size(resA,1)
    A=resA(counter,:);
    B=resB(counter,:);
    P=ranksum(A,B);
    stats=[stats; mean(A), mean(B), mean(B)-mean(A), mean(B)/mean(A), P -log10(P)];
end
stats(1:10,:)

figure
hist(stats(:,6))
title('P values (-log10)')

figure
hist(log10(stats(:,4)))
title('log10 foldchange (mean(B)/mean(A))')

figure
plot(log10(stats(:,4)),stats(:,6),'*')
title('vulcano: log10 foldchange vs -log10(P)')

%% up/down
str = fileread('metsCofactors.txt');
cofactorNames = regexp(str, '\r\n|\r|\n', 'split')'
cofactors=find(ismember(model_orig.metNames,cofactorNames))

up=find((log10(stats(:,4))>0).*(stats(:,6)>0).*(stats(:,3)>cutoffFluxDiff));
up=setdiff(up,cofactors);
up=intersect(up,metsOfInterest);
upt=array2table(stats(up,:),'RowNames',model_orig.mets(up),'VariableNames',{'meanA','meanB','diff','fc','pValue','-log10_p'});
upt.metNames=model_orig.metNames(up);
% upt=[upt(:,end), upt(:,1:(end-1))];
upt=sortrows(upt,6,'descend')

dn=find((log10(stats(:,4))<0).*(stats(:,6)>0).*(stats(:,3)<-cutoffFluxDiff));
dn=setdiff(dn,cofactors);
dn=intersect(dn,metsOfInterest);
dnt=array2table(stats(dn,:),'RowNames',model_orig.mets(dn),'VariableNames',{'meanA','meanB','diff','fc','pValue','-log10_p'});
dnt.metNames=model_orig.metNames(dn);
dnt=sortrows(dnt,6,'descend')

file=['fluxSumStats_' num2str(nrA) '_vs_' num2str(nrB) '.xlsx']
delete(file)
writetable(upt,file,'WriteRowNames',true,'Sheet','Up')
writetable(dnt,file,'WriteRowNames',true,'Sheet','Down')
