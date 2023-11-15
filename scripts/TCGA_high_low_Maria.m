function[Table_ON, Table_OFF]=high_low(LIHC_num, LIHC)
tLIHC_num=LIHC_num';
tLIHC_sorted=sortrows(tLIHC_num,find(contains(LIHC.Genes, 'AKR1A1')));
LIHC_sorted=tLIHC_sorted';
selectON=LIHC_sorted(contains(LIHC.Genes, 'AKR1A1'),:)>=prctile(LIHC_sorted(contains(LIHC.Genes, 'AKR1A1'),:),75);
LIHC_ON=LIHC_sorted(:,selectON);
selectOFF=LIHC_sorted(contains(LIHC.Genes, 'AKR1A1'),:)<=prctile(LIHC_sorted(contains(LIHC.Genes, 'AKR1A1'),:),25);
LIHC_OFF=LIHC_sorted(:,selectOFF);
colnames=LIHC.Properties.VariableNames(LIHC_sorted(end,:));
colnamesON=colnames(selectON);
colnamesOFF=colnames(selectOFF);
Table_OFF=LIHC_OFF(1:end-1,:);
Table_ON=LIHC_ON(1:end-1,:);
Table_OFF=array2table(Table_OFF);
Table_ON=array2table(Table_ON);
Table_ON.Properties.VariableNames=colnamesON;
Table_OFF.Properties.VariableNames=colnamesOFF;


end
