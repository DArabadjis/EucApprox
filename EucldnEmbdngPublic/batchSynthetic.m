clear all
fncD=@(x,y)sqrt(mean((x-y).^2,2));

nPnts=1500;nXprmnts=5;
dSet='swiss';posaDims=3;Nnghbrs=[5:2:15];
%dSet='difficult';posaDims=6;Nnghbrs=[7:2:17];

%dSet='helix';posaDims=2;Nnghbrs=[3:2:13];

wDmFull=false;
wVlmFrm=false;wRcrsSlv=false;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];
if wDmFull
for kSlct=1:nXprmnts
[pX,pL,pT]=generate_data(dSet,nPnts);dstMtrxSlct(kSlct).pT=pT;
dstMtrxSlct(kSlct).dstMtrx=MnmlNghbrhdsODim(false,posaDims-0.5,pX,fncD);
dstMtrxSlct(kSlct).dstMtrx=mkCnnctd(pX,dstMtrxSlct(kSlct).dstMtrx,fncD);
end

for kSlct=1:nXprmnts
[dstMtrxSlct(kSlct).pXEdge_sparse,dstMtrxSlct(kSlct).dstMtrxEdge_sparse]=CmptLclBase(dstMtrxSlct(kSlct).dstMtrx);
end


%wVlmFrm=false;wRcrsSlv=false;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];
parfor ii=1:nXprmnts
[dstMtrxSlctOut(ii).pXQ_short_short_sparse_NewCalc,dstMtrxSlctOut(ii).errorsAll]=EcldnEmbdngDatSgndNewCalc_Edge(...
	dstMtrxSlct(ii).dstMtrxEdge_sparse,dstMtrxSlct(ii).pXEdge_sparse(:,1:posaDims+2),posaDims,...
	wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);
end
else
for qoo=1:length(Nnghbrs)
for kSlct=1:nXprmnts
[pX,pL,pT]=generate_data(dSet,nPnts);dstMtrxSlct(qoo,kSlct).pT=pT;
dstMtrxSlct(qoo,kSlct).dstMtrx=CmptDstSprs(pX,fncD,false,Nnghbrs(qoo),false);
dstMtrxSlct(qoo,kSlct).dstMtrx=mkCnnctd(pX,dstMtrxSlct(qoo,kSlct).dstMtrx,fncD);
end
end

for qoo=1:length(Nnghbrs)
for kSlct=1:nXprmnts
[dstMtrxSlct(qoo,kSlct).pXEdge_sparse,dstMtrxSlct(qoo,kSlct).dstMtrxEdge_sparse]=CmptLclBase(dstMtrxSlct(qoo,kSlct).dstMtrx);
end
end

%wVlmFrm=~false;wRcrsSlv=false;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];
parfor ii=1:(length(Nnghbrs)*nXprmnts)
[qoo,kSlct]=ind2sub([length(Nnghbrs),nXprmnts],ii);
[dstMtrxSlctOut(ii).pXQ_short_short_sparse_NewCalc,dstMtrxSlctOut(ii).errorsAll]=EcldnEmbdngDatSgndNewCalc_Edge(...
	dstMtrxSlct(qoo,kSlct).dstMtrxEdge_sparse,dstMtrxSlct(qoo,kSlct).pXEdge_sparse(:,1:posaDims+2),posaDims,...
	wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut);
end
end
if ~wDmFull
if ~wVlmFrm
save(['batchOut_',dSet,'_nXprmnts_',num2str(nXprmnts),'.mat'],'dstMtrxSlctOut','dstMtrxSlct','Nnghbrs','nXprmnts')
else
save(['batchOut_VlmFrm_',dSet,'_nXprmnts_',num2str(nXprmnts),'.mat'],'dstMtrxSlctOut','dstMtrxSlct','Nnghbrs','nXprmnts')
end
else
if ~wVlmFrm
save(['batchOut_dmFull_',dSet,'_nXprmnts_',num2str(nXprmnts),'.mat'],'dstMtrxSlctOut','dstMtrxSlct','Nnghbrs','nXprmnts')
else
save(['batchOut_dmFull_VlmFrm_',dSet,'_nXprmnts_',num2str(nXprmnts),'.mat'],'dstMtrxSlctOut','dstMtrxSlct','Nnghbrs','nXprmnts')
end
end

delete(gcp('nocreate'))
exit