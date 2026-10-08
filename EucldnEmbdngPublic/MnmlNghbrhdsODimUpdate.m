function [isrt,indxQ,dsdQ,nNghbrs]=MnmlNghbrhdsODimUpdate(nDimsIn,pX,pXQry,fncDX,pY,pYQry,fncDY)
% --------------------Inputs---------------------
% nDimsIn - intrinsic dimension that all neighbors should have
% pX - Data Points for supporting the neighboring relations (point,feature)
% pY - Data Points for supporting the local dimension (point,feature)
% pXQry - New Data Points for querying their neighboring relations (point,feature)
% pYQry - New Data Points for querying their neighborhoods' local dimension (point,feature)
% fncDX - distance func for pX 
% fncDY - distance func for pY
%       in the matlab pattern @(x,y)apostash(x,y) 
%       where x: 1 x features and y: points x features
% --------------------Outputs---------------------
% nNghbrs - number of neighbors of each point, points x 1
% dsdQ - distances of neighbors from the source point, neighborhoods' points x 1
% isrt - indices of neighbors, neighborhoods' points x 1
% indxQ - indices of the source-points of the neighborhoods, useful for
%       sparse indexing together with isrt

if ~exist('pY','var')||isempty(pY),pY=pX;pYQry=pXQry;idia=true;else,idia=false;end
if ~exist('fncDY','var')||isempty(fncDY),fncDY=fncDX;end

nObs=size(pXQry,1);

dsdQ=inf(nObs,1);isrt=zeros(size(dsdQ));indxQ=[1:nObs]';nNghbrs=isrt;
arx=1;pX=[pX;nan(1,size(pX,2))];pY=[pY;nan(1,size(pY,2))];
for ii=1:nObs,pX(end,:)=pXQry(ii,:);pY(end,:)=pYQry(ii,:);
    dsdQ0=fncDX(pXQry(ii,:),pX);kNghbrs=1;
    while true,kNghbrs=kNghbrs+1;[dsdQslct,iSlct]=mink(dsdQ0,kNghbrs);CMmtrx=zeros(kNghbrs);
        if ~idia,dsdQslct=fncDY(pYQry(ii,:),pY(iSlct,:));end
        for jj=1:kNghbrs,dsdQ1=fncDX(pX(iSlct(jj),:),pX);
            [~,iSlct1]=mink(dsdQ1,kNghbrs);[iSlct1,iA]=intersect(iSlct,iSlct1);
            if ~idia,CMmtrx(iA,jj)=fncDY(pY(iSlct(jj),:),pY(iSlct1,:));
            else,CMmtrx(iA,jj)=dsdQ1(iSlct1);
            end
        end
        CMmtrx=(1-eye(size(CMmtrx))).*max(CMmtrx,CMmtrx');CMmtrx=distances(graph(CMmtrx));
        CMmtrx=CMmtrx(~any(isinf(CMmtrx)),~any(isinf(CMmtrx)));CMmtrx=CMmtrx.^2;
        if isempty(CMmtrx),dmsLcl=0;
        else
        % dstMtrx=pdist(pX(slct,:),fncD);CMmtrx=squareform(dstMtrx);
        CMmtrx=1/2*((CMmtrx(:,1)+CMmtrx(1,:))-CMmtrx);enas=ones(size(CMmtrx,1),1);tT=0;
        % CMmtrx=CMmtrx.*(1-eye(size(CMmtrx)));enas=ones(size(CMmtrx,1),1);tT=1;
        % CMmtrx=[[CMmtrx,enas];[enas',0]];
        % dmsLcl(ii)=rank(CMmtrx);
        s=svd(CMmtrx);s=sqrt(s(2:end-tT));prbs=(2*asin(s./sqrt(cumsum(s.^2)))/pi);
        prbs0=(2*asin(enas(2:end)./sqrt(cumsum(enas(2:end))))/pi);
        dmsLcl=sum(prbs./prbs0);
        end
        if dmsLcl>=nDimsIn,break;end
    end
    % iSlct(iSlct==size(pX,1))=size(pX,1)+ii-1;
    L0=length(iSlct);isrt(arx:arx+L0-1)=iSlct;dsdQ(arx:arx+L0-1)=dsdQslct;
    indxQ(arx:arx+L0-1)=ii;arx=arx+L0;nNghbrs(ii)=L0;ii
end