% first set the parent folder to be the "experiments_results" locally
drctrRd='.../experiments_results/';
load([drctrRd,'PBMC_5k_DistEuc_FullDim.mat'])

[~,~,QR]=svd(pX_Euc_short_short_sparse_dmEstm);pXQR=pX_Euc_short_short_sparse_dmEstm*QR;
figure,scatter3(pXQR(:,1),pXQR(:,2),pXQR(:,3),[],pY,'filled'),axis equal
title(['Propose method',char(39),'s 3D projection'])
pXumap=run_umap(pXL,'n_components',dmEstm_dmFull+1,'n_neighbors',dmEstm_dmFull+2,'verbose','text');
[~,~,QR]=svd(pXumap);pXQR=pXumap*QR;
figure,scatter3(pXQR(:,1),pXQR(:,2),pXQR(:,3),[],pY,'filled'),axis equal
title('UMAP 3D projection')
[pXQhssn] = compute_mapping(pXL,'tSNE',dmEstm_dmFull+1);
[~,~,QR]=svd(pXQhssn);pXQR=pXQhssn*QR;
figure,scatter3(pXQR(:,1),pXQR(:,2),pXQR(:,3),[],pY,'filled'),axis equal
title('t-SNE Eigenmaps 3D projection')
pXQhssn=laplacian_eigen_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
[~,~,QR]=svd(pXQhssn);pXQR=pXQhssn*QR;
figure,scatter3(pXQR(:,1),pXQR(:,2),pXQR(:,3),[],pY,'filled'),axis equal
title('Laplacian Eigenmaps 3D projection')
pXQhssn=EmbdngFrmGrmn(distances(graph(dstMtrx_dmFull)),dmEstm_dmFull+1);
[~,~,QR]=svd(pXQhssn);pXQR=pXQhssn*QR;
figure,scatter3(pXQR(:,1),pXQR(:,2),pXQR(:,3),[],pY,'filled'),axis equal
title('IsoMap 3D projection')