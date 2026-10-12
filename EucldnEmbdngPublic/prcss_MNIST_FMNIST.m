% first set the parent folder to be the "experiments_results" locally
drctrRd='.../experiments_results/';
%either mnist or fmnist
wchOne='mnist'
load([drctrRd,wchOne,'.mat']),

if strcmp(wchOne,'mnist')
    pYall=mnisttest(:,1).Variables;pXall=mnisttest(:,2:end).Variables;
elseif strcmp(wchOne,'fmnist')
    pYall=fashionmnisttest(:,1).Variables;pXall=fashionmnisttest(:,2:end).Variables;
end
pX=pXall;
fncD=@(x,y)sqrt(mean((x-y).^2,2));
dstMtrx=CmptDstSprs(pX,fncD,false,3,false);
dmEstm=round(estmDim2NNs(dstMtrx,true));close
dstMtrx=CmptDstSprs(pX,fncD,false,dmEstm+1,false);
dstMtrx_dmFull=MnmlNghbrhdsODim(false,dmEstm-0.5,pX,fncD);
dstMtrx=mkCnnctd(pX,dstMtrx,fncD);
dstMtrx_dmFull=mkCnnctd(pX,dstMtrx_dmFull,fncD);

[pXEdge_sparse,dstMtrxEdge_sparse]=CmptLclBase(dstMtrx,inf,false);
[pXEdge_sparse_dmFull,dstMtrxEdge_sparse_dmFull]=CmptLclBase(dstMtrx_dmFull,inf,false);
save([drctrRd,wchOne,'_graphConstructs_all.mat'],'dstMtrx','dstMtrx_dmFull','pXEdge_sparse','dstMtrxEdge_sparse','pXEdge_sparse_dmFull','dstMtrxEdge_sparse_dmFull')

wVlmFrm=false;wRcrsSlv=true;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];

posaDims=dmEstm+1;
[pXQ_short_short_sparse_dmFull_NewCalc,errors_dmFull]=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxEdge_sparse_dmFull,pXEdge_sparse_dmFull(:,1:posaDims+2),posaDims,...
	wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);
save([drctrRd,wchOne,'_embeddings_all.mat'],'pXQ_short_short_sparse_dmFull_NewCalc','errors_dmFull')
[pXQ_short_short_sparse_NewCalc,errors]=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxEdge_sparse,pXEdge_sparse(:,1:posaDims+2),posaDims,...
	wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);
save([drctrRd,wchOne,'_embeddings_all.mat'],'pXQ_short_short_sparse_NewCalc','errors',"-append")
