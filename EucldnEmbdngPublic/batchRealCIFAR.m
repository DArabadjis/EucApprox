clear all
%!!!!load CIFAR here!!!!!
drctr='/home/darampatzis/MATLAB Drive/mathimata/prjctsCds/datasets/NeurIPS/CIFAR_train/';
tabY=readtable([drctr,'trainLabels.csv']);
% pY=tabY(:,2).Variables;[catsY,~,pY]=unique(pY);
flss=dir([drctr,'*.png']);
f=rgb2lab(imread([drctr,flss(1).name]));
pXall=nan(length(flss),numel(f));pXall(1,:)=f(:);
for ii=2:length(flss),ii/length(flss),f=rgb2lab(imread([drctr,flss(ii).name]));pXall(ii,:)=f(:);end
pYall=tabY(:,2).Variables;[catsY,~,pYall]=unique(pYall);

N=2500; %the number of images per split
iSlct=floor(rand(size(pXall,1),1)*(size(pXall,1)/N));

fncD=@(x,y)sqrt(mean((x-y).^2,2));
for kSlct=1:max(iSlct)+1
pX=pXall(iSlct==(kSlct-1),:);
dstMtrxSlct(kSlct).dstMtrx=CmptDstSprs(pX,fncD,false,3,false);
dmEstm(kSlct)=round(estmDim2NNs(dstMtrxSlct(kSlct).dstMtrx,true));drawnow
dstMtrxSlct(kSlct).dstMtrx=CmptDstSprs(pX,fncD,false,dmEstm(kSlct)+1,false);
dstMtrxSlct(kSlct).dstMtrx_dmFull=MnmlNghbrhdsODim(false,dmEstm(kSlct)-0.5,pX,fncD);
dstMtrxSlct(kSlct).dstMtrx=mkCnnctd(pX,dstMtrxSlct(kSlct).dstMtrx,fncD);
dstMtrxSlct(kSlct).dstMtrx_dmFull=mkCnnctd(pX,dstMtrxSlct(kSlct).dstMtrx_dmFull,fncD);
close
end

for kSlct=1:max(iSlct)+1
[dstMtrxSlct(kSlct).pXEdge_sparse,dstMtrxSlct(kSlct).dstMtrxEdge_sparse]=CmptLclBase(dstMtrxSlct(kSlct).dstMtrx,inf,false);
[dstMtrxSlct(kSlct).pXEdge_sparse_dmFull,dstMtrxSlct(kSlct).dstMtrxEdge_sparse_dmFull]=CmptLclBase(dstMtrxSlct(kSlct).dstMtrx_dmFull,inf,false);
end
save('CIFAR_graph_constructs.mat','iSlct','dstMtrxSlct','dmEstm','-v7.3')
%!!!parfor!!! if unwilling it should be turned to for
parfor kSlct=1:max(iSlct)+1
posaDims=min(sum(spones(dstMtrxSlct(kSlct).dstMtrx)));
dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_NewCalc=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxSlct(kSlct).dstMtrxEdge_sparse,dstMtrxSlct(kSlct).pXEdge_sparse(:,1:posaDims+2),posaDims,[],false);
posaDims=dmEstm(kSlct)+1;
dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_dmFull_NewCalc=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxSlct(kSlct).dstMtrxEdge_sparse_dmFull,dstMtrxSlct(kSlct).pXEdge_sparse_dmFull(:,1:posaDims+2),posaDims,[],false);
end
save('CIFAR_embedding.mat','dstMtrxEmbddng')
parfor kSlct=1:max(iSlct)+1
datOutClss_short_short_sparse_dmFull_NewCalc(kSlct)=evalEmbddngClssfr(dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_dmFull_NewCalc,pYall(iSlct==(kSlct-1)),size(dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_dmFull_NewCalc,2));
datOut_short_short_sparse_dmFull_NewCalc(kSlct,:)=evalEmbddng(dstMtrxSlct(kSlct).dstMtrx_dmFull,dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_dmFull_NewCalc);
datOutClss_short_short_sparse_NewCalc(kSlct)=evalEmbddngClssfr(dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_NewCalc,pYall(iSlct==(kSlct-1)),size(dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_NewCalc,2));
datOut_short_short_sparse_NewCalc(kSlct,:)=evalEmbddng(dstMtrxSlct(kSlct).dstMtrx,dstMtrxEmbddng(kSlct).pXQ_short_short_sparse_NewCalc);
end
save('CIFAR_embedding.mat','datOutClss_short_short_sparse_dmFull_NewCalc',...
    'datOut_short_short_sparse_dmFull_NewCalc',...
    'datOutClss_short_short_sparse_NewCalc',...
    'datOut_short_short_sparse_NewCalc','-append')
parfor kSlct=1:max(iSlct)+1
pXQhssn=EmbdngFrmGrmn(distances(graph(dstMtrxSlct(kSlct).dstMtrx_dmFull)),dmEstm(kSlct)+1);
datOutClss_isomap_intrnsc_dmFull(kSlct)=evalEmbddngClssfr(pXQhssn,pYall(iSlct==(kSlct-1)),size(pXQhssn,2));
datOut_isomap_intrnsc_dmFull(kSlct,:)=evalEmbddng(dstMtrxSlct(kSlct).dstMtrx_dmFull,pXQhssn);
end
save('CIFAR_embedding.mat','datOutClss_isomap_intrnsc_dmFull',...
    'datOut_isomap_intrnsc_dmFull','-append')
parfor kSlct=1:max(iSlct)+1
pXQhssn=lle_intrnsc(dstMtrxSlct(kSlct).dstMtrx_dmFull,dmEstm(kSlct)+1);
datOutClss_lle_intrnsc_dmFull(kSlct)=evalEmbddngClssfr(pXQhssn,pYall(iSlct==(kSlct-1)),size(pXQhssn,2));
datOut_lle_intrnsc_dmFull(kSlct,:)=evalEmbddng(dstMtrxSlct(kSlct).dstMtrx_dmFull,pXQhssn);
end
save('CIFAR_embedding.mat','datOutClss_lle_intrnsc_dmFull',...
    'datOut_lle_intrnsc_dmFull','-append')
parfor kSlct=1:max(iSlct)+1
pX=pXall(iSlct==(kSlct-1),:);
pXQhssn=compute_mapping(pX,'tSNE',dmEstm(kSlct)+1);
datOutClss_tsne_dmFull(kSlct)=evalEmbddngClssfr(pXQhssn,pYall(iSlct==(kSlct-1)),size(pXQhssn,2));
datOut_tsne_dmFull(kSlct,:)=evalEmbddng(dstMtrxSlct(kSlct).dstMtrx_dmFull,pXQhssn);
end
save('CIFAR_embedding.mat','datOutClss_tsne_dmFull',...
    'datOut_tsne_dmFull','-append')