function [dstMtrx,isrt,indxQ,dsdQ,nNghbrs]=CmptDstSprs(pX,fncD,op,nDims,apst)
% --------------------Inputs---------------------
% pX - Data Points (point,feature)
% fncD - distance func in the matlab pattern @(x,y)apostash(x,y) 
%       where x: 1 x features and y: points x features
% op - true for probability-like distances (multiplicative aggregation), 
%       false for ordinary distances (additive aggregation)
% nDims - number of neighbors(apst=false)/radius of the neighborhood(apst=true)
% --------------------Outputs---------------------
% nNghbrs - number of neighbors of each point  
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% dsdQ - distances of neighbors from the source point 
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% isrt - indices of neighbors 
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% indxQ - indices of the source-points of the neighborhoods, useful for
%       sparse indexing together with isrt
% dstMtrx - points x points sparse matrix having zeros in non-connected
%           points and the distance between conneected points

if apst
dsdQ=inf(size(pX,1),1);isrt=zeros(size(dsdQ));indxQ=[1:size(pX,1)]';
nNghbrs=isrt;arx=1;
for ii=1:size(pX,1),dsdQ0=fncD(pX(ii,:),pX);
    slct=dsdQ0<=nDims;isrt0=find(slct);dsdQ0=dsdQ0(isrt0);L0=length(isrt0);
    % [dsdQ0,isrt1]=sort(dsdQ0);isrt0=isrt0(isrt1);
    isrt(arx:arx+L0-1)=isrt0;dsdQ(arx:arx+L0-1)=dsdQ0;
    indxQ(arx:arx+L0-1)=ii;arx=arx+L0;nNghbrs(ii)=L0;
    % ii
end
isrt=isrt(1:arx-1);dsdQ=dsdQ(1:arx-1);indxQ=indxQ(1:arx-1);
else
nDims=min(nDims,size(pX,1)-1);nDims=nDims+1;nNghbrs=nDims*ones(size(pX,1),1);
dsdQ=inf(size(pX,1),nDims);isrt=zeros(size(pX,1),nDims);
indxQ=repmat([1:size(pX,1)]',1,nDims);
for ii=1:size(pX,1),dsdQ0=fncD(pX(ii,:),pX);
    [dsdQ(ii,:),isrt(ii,:)]=mink(dsdQ0,nDims);
    % ii
end
end
if op
dstMtrx=max(sparse(isrt(:),indxQ(:),abs(log(1-dsdQ(:))),size(pX,1),size(pX,1)), ...
    sparse(indxQ(:),isrt(:),abs(log(1-dsdQ(:))),size(pX,1),size(pX,1)));
else
dstMtrx=max(sparse(isrt(:),indxQ(:),dsdQ(:),size(pX,1),size(pX,1)), ...
    sparse(indxQ(:),isrt(:),dsdQ(:),size(pX,1),size(pX,1)));
end