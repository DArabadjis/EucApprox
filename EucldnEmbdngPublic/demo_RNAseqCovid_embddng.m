%clear all
% first set the parent folder to be the "experiments_results" locally
%drctrRd=;
drctrRd=[drctrRd,'Covid19/2020_Wilk_COVID/'];
load([drctrRd,'560_cell.counts.matrices.mat'])
pX=[sparsevar.intron',sparsevar.exon'];
orioCells=20;
pX=pX(:,sum(pX)>=orioCells);
pX=pX./sum(pX,2);
pXL=log2(1+pX);
%Labels of diagnoses
tabT=readtable([drctrRd,'categories_annotations_2026-04-26.csv'],'Range',[2,1]);
slctOut='cell_type_coarse';
% slctOut='cell_type_fine';
pY=tabT.(slctOut);[catsY,~,pY]=unique(pY);
isCvd=tabT.Status;isVent=tabT.Ventilated;ptntId=tabT.Donor_full;
pY=pY(1:size(pXL,1),:);isCvd=isCvd(1:size(pXL,1),:);ptntId=ptntId(1:size(pXL,1),:);

%distance function
fncD=@(x,y)sqrt(mean((x-y).^2,2));
% fncD=@(x,y)pdist2(y,x,'correlation');
% fncD=@(x,y)1-mean(y.*x,2)./sqrt(mean(y.^2,2).*mean(x.^2,2));

%2 NNs
dstMtrx_dmFull=CmptDstSprs(pXL,fncD,false,3,false);
%Intrinsic dim
dmEstm_dmFull=round(estmDim2NNs(dstMtrx_dmFull,true))
dstMtrx_dmFull=CmptDstSprs(pXL,fncD,false,dmEstm_dmFull+1,false);
save([drctrRd,'2020_Wilk_COVID_graph_constructs.mat'],'dmEstm_dmFull','dstMtrx_dmFull');
[pXEdgeSparse,dstMtrxEdgeSparse]=CmptLclBase(dstMtrx_dmFull);
save([drctrRd,'2020_Wilk_COVID_graph_constructs.mat'],'pXEdgeSparse','dstMtrxEdgeSparse','-append');
%Compute embedding
wVlmFrm=false;wRcrsSlv=true;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];
pX_Euc_short_short_sparse_dmEstm=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxEdgeSparse,pXEdgeSparse(:,1:dmEstm_dmFull+3),dmEstm_dmFull+1,wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);
save([drctrRd,'2020_Wilk_COVID_embeddings.mat'],'pX_Euc_short_short_sparse_dmEstm');
%Evaluate clustering
datOut_dmEstm_Euc_NewCalc_sparse=evalEmbddngClssfr(pX_Euc_short_short_sparse_dmEstm,pY,dmEstm_dmFull+1);
pXQhssn=EmbdngFrmGrmn(distances(graph(dstMtrx_dmFull)),dmEstm_dmFull+1);
datOut_dmEstm_isomap_intrnsc=evalEmbddngClssfr(pXQhssn,pY,dmEstm_dmFull+1)
pXQhssn=lle_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOut_dmEstm_lle_intrnsc=evalEmbddngClssfr(pXQhssn,pY,dmEstm_dmFull+1);
pXQhssn=ltsa_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOut_dmEstm_ltsa_intrnsc=evalEmbddngClssfr(pXQhssn,pY,dmEstm_dmFull+1);
pXQhssn=laplacian_eigen_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOut_dmEstm_laplacian_intrnsc=evalEmbddngClssfr(pXQhssn,pY,dmEstm_dmFull+1);
pXumap=run_umap(pXL,'n_components',dmEstm_dmFull+1,'n_neighbors',dmEstm_dmFull+2,'verbose','text');
datOut_dmEstm_umap=evalEmbddngClssfr(pXumap,pY,dmEstm_dmFull+1);
[pXQhssn] = compute_mapping(pXL,'tSNE',dmEstm_dmFull+1,dmEstm_dmFull+1);
datOut_dmEstm_tsne=evalEmbddngClssfr(pXQhssn,pY,dmEstm_dmFull+1);
save([drctrRd,'2020_Wilk_COVID_classifications_',slctOut,'.mat'],'datOut_dmEstm_Euc_NewCalc_sparse',...
    'datOut_dmEstm_isomap_intrnsc','datOut_dmEstm_lle_intrnsc','datOut_dmEstm_ltsa_intrnsc',...
    'datOut_dmEstm_laplacian_intrnsc','datOut_dmEstm_umap','datOut_dmEstm_tsne');
%Evaluate geometry
datOutQuality_dmEstm_short_short_sparse_dmEstm=evalEmbddng(dstMtrx_dmFull,pX_Euc_short_short_sparse_dmEstm);
pXQhssn=laplacian_eigen_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOutQuality_dmEstm_laplacian_intrnsc=evalEmbddng(dstMtrx_dmFull,pXQhssn)
pXQhssn=lle_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOutQuality_dmEstm_lle_intrnsc=evalEmbddng(dstMtrx_dmFull,pXQhssn)
pXQhssn=ltsa_intrnsc(dstMtrx_dmFull,dmEstm_dmFull+1);
datOutQuality_dmEstm_ltsa_intrnsc=evalEmbddng(dstMtrx_dmFull,pXQhssn)
pXQhssn=EmbdngFrmGrmn(distances(graph(dstMtrx_dmFull)),dmEstm_dmFull+1);
datOutQuality_dmEstm_isomap_intrnsc=evalEmbddng(dstMtrx_dmFull,pXQhssn)
datOutQuality_dmEstm_short_short_sparse_dmEstm=evalEmbddng(dstMtrx_dmFull,pX_Euc_short_short_sparse_dmEstm)
pXumap=run_umap(pXL,'n_components',dmEstm_dmFull+1,'n_neighbors',dmEstm_dmFull+2,'verbose','text');
datOutQuality_dmEstm_umap=evalEmbddng(dstMtrx_dmFull,pXumap)
[pXQhssn] = compute_mapping(pXL,'tSNE',dmEstm_dmFull+1);
datOutQuality_dmEstm_tsne=evalEmbddng(dstMtrx_dmFull,pXQhssn)
save([drctrRd,'2020_Wilk_COVID_geometryEval_',slctOut,'.mat'],'datOutQuality_dmEstm_short_short_sparse_dmEstm',...
    'datOutQuality_dmEstm_laplacian_intrnsc','datOutQuality_dmEstm_lle_intrnsc','datOutQuality_dmEstm_ltsa_intrnsc',...
    'datOutQuality_dmEstm_isomap_intrnsc','datOutQuality_dmEstm_umap','datOutQuality_dmEstm_tsne');