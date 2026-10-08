function dmsLcl=CmptDms(dstMtrx,isrt,indxQ,apst)
% --------------------Inputs---------------------
% dstMtrx - points x points sparse matrix having zeros in non-connected
%           points and the distance between connected points
% isrt - indices of neighbors 
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% indxQ - indices of the source-points of the neighborhoods, useful for
%       sparse indexing together with isrt
% --------------------Outputs---------------------
% dmsLcl - points x 1, dimensionality of each point's neighborhood

dmsLcl=zeros(size(dstMtrx,1),1);nObs=size(dstMtrx,1);
if apst
for ii=1:nObs,slct=indxQ==ii;
    iNbghr=isrt(slct);CMmtrx=dstMtrx(iNbghr,iNbghr);
    CMmtrx=(1-eye(size(CMmtrx))).*CMmtrx;CMmtrx=distances(graph(CMmtrx));
    CMmtrx=CMmtrx(~any(isinf(CMmtrx)),~any(isinf(CMmtrx)));CMmtrx=CMmtrx.^2;
    if ~isempty(CMmtrx)
    % dstMtrx=pdist(pX(slct,:),fncD);CMmtrx=squareform(dstMtrx);
    CMmtrx=1/2*((CMmtrx(:,1)+CMmtrx(1,:))-CMmtrx);enas=ones(size(CMmtrx,1),1);tT=0;
    % CMmtrx=CMmtrx.*(1-eye(size(CMmtrx)));enas=ones(size(CMmtrx,1),1);tT=1;
    % CMmtrx=[[CMmtrx,enas];[enas',0]];
    % dmsLcl(ii)=rank(CMmtrx);
    s=svd(CMmtrx);s=sqrt(s(2:end-tT));prbs=(2*asin(s./sqrt(cumsum(s.^2)))/pi);
    prbs0=(2*asin(enas(2:end)./sqrt(cumsum(enas(2:end))))/pi);
    dmsLcl(ii)=sum(prbs./prbs0);
    else,dmsLcl(ii)=0;
    end
end
else,enas=ones(size(isrt,2),1);
for ii=1:nObs
    iNbghr=isrt(ii,:);CMmtrx=dstMtrx(iNbghr,iNbghr);
    CMmtrx=(1-eye(size(CMmtrx))).*CMmtrx;CMmtrx=distances(graph(CMmtrx));
    CMmtrx=CMmtrx(~any(isinf(CMmtrx)),~any(isinf(CMmtrx)));CMmtrx=CMmtrx.^2;
    if ~isempty(CMmtrx)
    % gWgI=1/2*(((dstMtrx(iNbghr,ii))+(dstMtrx(iNbghr,ii))')-dstMtrx(iNbghr,iNbghr));
    % CMmtrx=full((1-wDatQQ(iNbghr,ii)).*gWgI);
    % s=svd(CMmtrx);dmsLcl(ii)=sum(s/sqrt(sum(s.^2)));
    % dstMtrx=pdist(pX(isrt(ii,:),:),fncD);CMmtrx=squareform(dstMtrx);
    CMmtrx=1/2*((CMmtrx(:,1)+CMmtrx(1,:))-CMmtrx);tT=0;
    % CMmtrx=CMmtrx.*(1-eye(size(CMmtrx)));tT=1;
    % CMmtrx=[[CMmtrx,enas];[enas',0]];
    % dmsLcl(ii)=rank(CMmtrx);
    s=svd(CMmtrx);s=sqrt(s(2:end-tT));prbs=(2*asin(s./sqrt(cumsum(s.^2)))/pi);
    prbs0=(2*asin(enas(2:end)./sqrt(cumsum(enas(2:end))))/pi);
    dmsLcl(ii)=sum(prbs./prbs0);
    else,dmsLcl(ii)=0;
    end
end
end
dmsLcl(isnan(dmsLcl))=0;
% dmsLcl=dmsLcl-2;