function [pX,dstMtrxN]=CmptLclBase(dstMtrx,wSync,dpth)
% Inputs
% dstMtrx : graph adjacency matrix - sparse #points x #points
% dpth : flag for 2nd order neighborhood output (false)
% wSync : flag for spectral synchronization of the transports (false)
% Outputs
% pX : local trivializing transports in the form [edge_row,edge_col,frame] - edges x intrinsic dim + 2 
% dstMtrxN : outpput graph in correspondence with the transports - sparse #points x #points

if ~exist('dpth','var')||isempty(dpth)||~islogical(dpth),dpth=false;end
if ~exist('wSync','var')||isempty(wSync)||~islogical(wSync),wSync=false;end
if ~issparse(dstMtrx),dstMtrx=sparse(dstMtrx);end
dmsLcl=zeros(size(dstMtrx,1),1);nObs=size(dstMtrx,1);
dstMtrx=spdiags(zeros(nObs,1),0,dstMtrx);[isrt,indxQ]=find(dstMtrx);
% nDims=min(nDims,min(sum(spones(dstMtrx))))

% if dpth,nDims=min(sum(spones(spones(dstMtrx)*spones(dstMtrx)))),end,dstMtrx=max(dstMtrx,dstMtrx');
nDims=min(sum(spones(spones(dstMtrx)*spones(dstMtrx))))
pX=zeros(size(isrt,1)+nObs,nDims+2);arx=1;offst=0;dsdQ=inf(size(pX,1),1);
for ii=1:nObs,slct=indxQ==ii;iNbghr=[ii;isrt(slct)];iNbghr0=iNbghr;%iNbghr1=iNbghr0;
    slct=any(indxQ==iNbghr',2);iNbghr=union(iNbghr,isrt(slct));iNbghr1=iNbghr;%if dpth,iNbghr1=iNbghr;else,iNbghr1=iNbghr0;end
    slct=any(indxQ==iNbghr',2);iNbghr=union(iNbghr,isrt(slct));
    CMmtrx=dstMtrx(iNbghr,iNbghr);CMmtrx=(1-eye(size(CMmtrx))).*CMmtrx;
    CMmtrx=distances(graph(CMmtrx));
    % -----------------------------
    % CMmtrx=(iNbghr1==iNbghr')*CMmtrx*(iNbghr==iNbghr1');
    % if ~isempty(CMmtrx),pXQ=EmbdngFrmGrmn(CMmtrx,nDims);
    %     pXQ=(iNbghr0==iNbghr1')*(pXQ-pXQ(iNbghr1==ii,:));
    %     pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),pXQ];
    % end
    % -----------------------------
    for jj=1:size(CMmtrx,1),CMmtrx=min(CMmtrx,CMmtrx(:,jj)+CMmtrx(jj,:));end
    CMmtrx=(iNbghr1==iNbghr')*CMmtrx*(iNbghr==iNbghr1');
    if dpth,dsdQ(arx:arx+length(iNbghr1)-1)=CMmtrx(:,iNbghr1==ii);
    else,dsdQ(arx:arx+length(iNbghr0)-1)=(iNbghr0==iNbghr1')*CMmtrx(:,iNbghr1==ii);
    end
    % CMmtrx=CMmtrx(~any(isinf(CMmtrx)),~any(isinf(CMmtrx)));
    CMmtrx=CMmtrx.^2;
    if ~isempty(CMmtrx)
    % dstMtrx=pdist(pX(slct,:),fncD);CMmtrx=squareform(dstMtrx);
    % CMmtrx=1/2*((CMmtrx(:,1)+CMmtrx(1,:))-CMmtrx);enas=ones(size(CMmtrx,1),1);
    CMmtrx=1/2*((CMmtrx(:,iNbghr1==ii)+CMmtrx(iNbghr1==ii,:))-CMmtrx);
    % CMmtrx=CMmtrx.*(1-eye(size(CMmtrx)));enas=ones(size(CMmtrx,1),1);tT=1;
    % CMmtrx=[[CMmtrx,enas];[enas',0]];
    % dmsLcl(ii)=rank(CMmtrx);
    [Q,s,~]=svd(CMmtrx);%Q=(iNbghr0==iNbghr1')*Q;
    s=sqrt(s(1+offst:nDims+offst,1+offst:nDims+offst));Q=Q(:,1+offst:nDims+offst);
    if dpth
    % pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),Q*s];
    pX(arx:arx+length(iNbghr1)-1,:)=[iNbghr1,ii*ones(size(iNbghr1)),Q*s];
    arx=arx+length(iNbghr1);
    else
    % zOut=(iNbghr0==iNbghr1')*Q*s;
    pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),(iNbghr0==iNbghr1')*Q*s];
    arx=arx+length(iNbghr0);
    end
    end
    % -----------------------------
    % arx=arx+length(iNbghr0);
    % arx=arx+length(iNbghr1);
end
% dstMtrxN=max(sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2)),...
%     sparse(pX(:,2),pX(:,1),dsdQ,size(dstMtrx,1),size(dstMtrx,2)));
dstMtrxN=sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2));
if wSync,pX=SpectrlSync(pX);end
% indsCvrd=zeros(nObs,1);indsAll=true(nObs,1);prsm=ones(1,nDims);
% ii=1;
% while any(indsAll)
% indsCvrd(ii)=find(indsAll,1,'first');tel=ii;indsAll(indsCvrd(ii))=false;
% while (ii<=tel),arx=indsCvrd(ii);
%     slct=pX(:,2)==arx;iNbghr=pX(slct,1);pXC=pX(slct,3:end);
%     iNbghr0=setdiff(iNbghr,indsCvrd(indsCvrd>0));
%     for jj=1:length(iNbghr0),slctIn=pX(:,2)==iNbghr0(jj);iNbghrIn=pX(slctIn,1);
%         pXCin=pX(slctIn,3:end);pXCin=pXCin-pXCin(iNbghrIn==arx,3:end);
%         antstxa=(iNbghr==iNbghrIn');pXCin=antstxa*pXCin;pXC0=antstxa'*pXC;
%         pXCin=pXCin(any(antstxa,2),:);pXC0=pXC0(any(antstxa',2),:);
%         [qL,~,qR]=svd(pXC0'*pXCin);prsm(end)=sign(det(qL*qR'));
%         pX(slctIn,3:end)=pX(slctIn,3:end)*(qR.*prsm)*qL';
%     end
%     indsAll(iNbghr0)=false;
%     indsCvrd(tel+1:tel+length(iNbghr0))=iNbghr0;tel=tel+length(iNbghr0);
%     ii=ii+1;
% end
% end

function pX=SpectrlSync(pX)
nObs=max(pX(:,2));nDims=size(pX,2)-2;dkts=(1:nDims^2);
prsm=ones(1,nDims);pQR=nan(size(pX,1),nDims^2+2);
arx=1;
for ii=1:nObs
    slct=pX(:,2)==ii;iNbghr0=pX(slct,1);iNbghr=unique(pX(any(pX(:,1)==iNbghr0',2),1));
    pXC=pX(slct,3:end);
    for jj=1:length(iNbghr),slctIn=pX(:,2)==iNbghr(jj);iNbghrIn=pX(slctIn,1);
        pXCin=pX(slctIn,3:end);antstxa=(iNbghr==iNbghrIn');pXCin=antstxa*pXCin;%pXC0=antstxa'*pXC;
        pXCin=pXCin(any(antstxa,2),:);pXC0=pXC(any(antstxa,2),:);
        pXCin=pXCin-mean(pXCin);pXC0=pXC0-mean(pXC0);
        [qL,~,qR]=svd(pXC0'*pXCin);%prsm(end)=sign(det(qL*qR'));
        RotMat=(qR.*prsm)*qL';pQR(arx,:)=[iNbghr(jj),ii,RotMat(:)'];
        arx=arx+1;
    end
end
hMat0=((pQR(:,2)-1).*dkts(end)+dkts)';hMat=[mod(hMat0(:)-1,nDims)+1];
hMat0=((pQR(:,1)-1).*dkts(end)+dkts)';hMat=[hMat,floor((hMat0(:)-1)/nDims)+1];
hMat0=pQR(:,3:end)';hMat=[hMat,hMat0(:)];
hMat=sparse(hMat(:,1),hMat(:,2),hMat(:,3),nObs*nDims,nObs*nDims);
% opts.isreal = true;opts.issym = true;
% [Qin,~] = eigs(hMat, nDims, 'largestreal', opts);
% opts.isreal = true;opts.issym = true;
[Qin,~] = svds(hMat, nDims);
for ii=1:nObs
    [qL,~,qR]=svd(Qin((ii-1)*nDims+(1:nDims),:));%prsm(end)=sign(det(qR*qL'));
    RotMat=(qL.*prsm)*qR';
    slct=pX(:,2)==ii;pX(slct,3:end)=pX(slct,3:end)*RotMat;
end