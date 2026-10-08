clear all
% first set the parent folder to be the "experiments_results" locally
drctrRd='/home/darampatzis/MATLAB Drive/mathimata/prjctsCds/datasets/NeurIPS/';
drctrRd=[drctrRd,'GDC-raw_data/GDC-lung-RNAseq/'];
load([drctrRd,'GDC-lung-RNAseq-tpm.mat'])
pX=pX(:,5:end);pXL=log2(1+pX);
%distance function
fncD=@(x,y)1-mean(y.*x,2)./sqrt(mean(y.^2,2).*mean(x.^2,2));

%2 NNs
dstMtrx_dmFull=CmptDstSprs(pXL,fncD,false,3,false);
%Intrinsic dim
dmEstm_dmFull=round(estmDim2NNs(dstMtrx_dmFull,true))

dstMtrx_dmFull=CmptDstSprs(pXL,fncD,false,dmEstm_dmFull+1,false);
[pXEdgeSparse,dstMtrxEdgeSparse]=CmptLclBase(dstMtrx_dmFull);

wVlmFrm=false;wRcrsSlv=tue;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];
pX_Euc_short_short_sparse_dmEstm=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxEdgeSparse,pXEdgeSparse(:,1:dmEstm_dmFull+3),dmEstm_dmFull+1,wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);