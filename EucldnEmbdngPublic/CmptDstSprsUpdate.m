function [isrt,indxQ,dsdQ,nNghbrs]=CmptDstSprsUpdate(pX,pXQry,fncD,op,nDims,apst)
% --------------------Inputs---------------------
% nDimsIn - intrinsic dimension that all neighbors should have
% pX - Data Points for supporting the neighboring relations (point,feature)
% pXQry - New Data Points for querying their neighboring relations (point,feature)
% fncD - distance func in the matlab pattern @(x,y)apostash(x,y) 
%       where x: 1 x features and y: points x features
% op - true for probability-like distances (multiplicative aggregation), 
%       false for ordinary distances (additive aggregation)
% nDims - number of neighbors(apst=false)/radius of the neighborhood(apst=true)
% --------------------Outputs---------------------
% nNghbrs - number of neighbors of each point, points x 1
% dsdQ - distances of neighbors from the source point, neighborhoods' points x 1
% isrt - indices of neighbors, neighborhoods' points x 1
% indxQ - indices of the source-points of the neighborhoods, useful for
%       sparse indexing together with isrt

nObs=size(pXQry,1);

dsdQ=inf(nObs,1);isrt=zeros(size(dsdQ));indxQ=[1:nObs]';nNghbrs=isrt;
arx=1;pX=[pX;nan(1,size(pX,2))];
for ii=1:nObs,pX(end,:)=pXQry(ii,:);
    dsdQ0=fncD(pXQry(ii,:),pX);
    if apst
        slct=dsdQ0<=nDims;iSlctOut=find(slct);dsdQslctOut=dsdQ0(iSlctOut);
    else,[dsdQslctOut,iSlctOut]=mink(dsdQ0,nDims);
    end
    % iSlct(iSlct==size(pX,1))=size(pX,1)+ii-1;
    L0=length(iSlctOut);isrt(arx:arx+L0-1)=iSlctOut;dsdQ(arx:arx+L0-1)=dsdQslctOut;
    indxQ(arx:arx+L0-1)=ii;arx=arx+L0;nNghbrs(ii)=L0;ii
end
if op,dsdQ=abs(log(1-dsdQ(:)));end