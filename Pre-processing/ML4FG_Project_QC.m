%% Pancreas Data: Load Data 

Pancreas_Data = readtable('pancreas_refseq_rpkms_counts_3514sc.txt');
Metadata_Table = readtable('PancreasMetadata.csv'); 
Pancreas_Samples = string(table2array(readtable('PancreasSampleNames.txt', ...
    'Delimiter','tab', ...
    'ReadVariableNames',false)));

%% Pancreas Data: Segment Data 

PancreasGenelist = string(table2array(Pancreas_Data(1:26179,1)));
X_Pancreas = table2array(Pancreas_Data(1:26179,3:3516));

QC_indices = zeros(1,length(Metadata_Table.CellID));
for i = 1:length(Metadata_Table.CellID)
    QC_indices(i) = find(Pancreas_Samples == Metadata_Table.CellID(i));
end

X_filtered = X_Pancreas(:,QC_indices);
Samples_filtered = Pancreas_Samples(QC_indices);

%% Pancreas Data: t-SNE vs UMAP 
log_X_filtered = log2(X_filtered'+1); 
[X_pca, PcaMapping] = pca(log_X_filtered, 0.6); 

X_umap = sc_umap(X_pca', 3, false, false); 
figure(1)
subplot(1,2,1)
gscatter(X_umap(:,1),X_umap(:,2), Metadata_Table.CellType); 
title('UMAP')
legend('boxoff')

X_tsne = tsne(X_pca, 'NumPCAComponents', 0, 'Perplexity', 50);
subplot(1,2,2)
gscatter(X_tsne(:,1),X_tsne(:,2), Metadata_Table.CellType); 
title('tSNE')

%% Neftel et al GBM Data: Load Data 

[X_GBM,GBM_Genelist,GBM_Barcodelist] = sc_readmtxfile('matrix2.mtx','features2.tsv','barcodes2.tsv',1);
GBM_Metadata = readtable('metadata2.csv'); 

% Cells and genes have already been filtered
% Already logged, normalized and centered?
% Use relative expression as input to UMAP (PCA components = 50)

log_X_GBM = log2((X_GBM/10)+1);
RelativeExpression = log_X_GBM - mean(log_X_GBM,2);

GBM_umap = sc_umap(RelativeExpression, 3, false, false); 
figure(2)
gscatter(GBM_umap(:,1),GBM_umap(:,2), GBM_Metadata.cell_assignment); 
title('UMAP')
legend('boxoff')

%% Export Data Files 
writematrix(log_X_filtered','~/Desktop/Pancreas_LogAndFiltered.txt','Delimiter','tab');
writematrix(RelativeExpression,'~/Desktop/GBM_RelativeExpression.txt','Delimiter','tab');
