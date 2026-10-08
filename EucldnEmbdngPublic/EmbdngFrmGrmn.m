function [pXQ,S,V]=EmbdngFrmGrmn(CMmtrx,nDims)
% --------------------Inputs----------------------
% CMmtrx - points x points, all pairwise on-graph distances
% nDims - dimensions of the embedding + 1
% --------------------Outputs---------------------
% pXQ - points x dimensions, the embedding coordinates


if ~issymmetric(CMmtrx),CMmtrx=max(CMmtrx,CMmtrx');end
for ii=1:size(CMmtrx,1),CMmtrx=min(CMmtrx,CMmtrx(:,ii)+CMmtrx(ii,:));end
CMmtrx=CMmtrx.^2;gHat=mean(CMmtrx);gMg=0.5*((gHat-CMmtrx)+gHat');
MxItrr=1000;
[V,S]=svds(gMg,nDims+1,'largest','MaxIterations',MxItrr);
V=V(:,2:end);S=sqrt(S(2:end,2:end));pXQ=V*S;