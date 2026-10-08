function [pX,frmDat,dstMtrxN]=LclFrmsE(dstMtrx,nDims,dpth,wVlmFrm,wBtch,pY)
if ~exist('dpth','var')||isempty(dpth)||~islogical(dpth),dpth=false;end
if ~exist('wVlmFrm','var')||isempty(wVlmFrm)||~islogical(wVlmFrm),wVlmFrm=false;end
if ~exist('wBtch','var')||isempty(wBtch)||~islogical(wBtch),wBtch=true;end
if ~exist('pY','var')||isempty(pY),wDfrm=false;else,wDfrm=true;end
if ~issparse(dstMtrx),dstMtrx=sparse(dstMtrx);end
dmsLcl=zeros(size(dstMtrx,1),1);nObs=size(dstMtrx,1);
dstMtrx=spdiags(zeros(nObs,1),0,dstMtrx);
nDims=min(nDims,min(sum(spones(dstMtrx))))
if wBtch
    geitBin=spones(max(dstMtrx,speye(size(dstMtrx,1),size(dstMtrx,2))));nNeighbors=sum(geitBin);
    [isrt,indxQ]=find(geitBin);
    geitBin1=geitBin*geitBin;[isrt1,indxQ1]=find(geitBin1);nNeighbors1=sum(spones(geitBin1));
    geitBin2=geitBin1*geitBin;[isrt2,indxQ2]=find(geitBin2);nNeighbors2=sum(spones(geitBin2));
else,[isrt,indxQ]=find(dstMtrx);
end
if dpth,nDims=min(sum(spones(spones(dstMtrx)*spones(dstMtrx)))),end,dstMtrx=max(dstMtrx,dstMtrx');
pX=zeros(size(isrt,1)+(1-wBtch)*nObs,nDims+2);arx=1;offst=0;dsdQ=inf(size(pX,1),1);
arx0=1;arx1=arx0;arx2=arx0;
for ii=1:nObs
    if wBtch,%iNbghr0=find(geitBin(:,ii));
        % if dpth,iNbghr1=find(geitBin1(:,ii));else,iNbghr1=iNbghr0;end
        % iNbghr=find(geitBin2(:,ii));
        
        iNbghr0=isrt(arx0:arx0+nNeighbors(ii)-1);arx0=arx0+nNeighbors(ii);
        if dpth,iNbghr1=isrt1(arx1:arx1+nNeighbors1(ii)-1);arx1=arx1+nNeighbors1(ii);
        else,iNbghr1=iNbghr0;
        end
        iNbghr=isrt2(arx2:arx2+nNeighbors2(ii)-1);arx2=arx2+nNeighbors2(ii);

        % slct=indxQ==ii;iNbghrA=[ii;isrt(slct)];
        % slct=any(indxQ==iNbghrA',2);iNbghrB=union(iNbghrA,isrt(slct));
        % if dpth,iNbghrQ=iNbghrB;else,iNbghrQ=iNbghrA;end
        % slct=any(indxQ==iNbghrB',2);iNbghrC=union(iNbghrB,isrt(slct));
        % if ~isempty(setdiff(iNbghr0,iNbghrA))||~isempty(setdiff(iNbghr1,iNbghrQ))||~isempty(setdiff(iNbghr,iNbghrC))
        %     setdiff(iNbghr0,iNbghrA)'
        %     setdiff(iNbghr1,iNbghrQ)'
        %     setdiff(iNbghr,iNbghrC)'
        %     pause,
        % end
    else
    slct=indxQ==ii;iNbghr=[ii;isrt(slct)];iNbghr0=iNbghr;%iNbghr1=iNbghr0;
    slct=any(indxQ==iNbghr',2);iNbghr=union(iNbghr,isrt(slct));if dpth,iNbghr1=iNbghr;else,iNbghr1=iNbghr0;end
    slct=any(indxQ==iNbghr',2);iNbghr=union(iNbghr,isrt(slct));
    end
    CMmtrx=dstMtrx(iNbghr,iNbghr);CMmtrx=(1-eye(size(CMmtrx))).*CMmtrx;
    CMmtrx=distances(graph(CMmtrx));
    % -----------------------------
    % CMmtrx=(iNbghr1==iNbghr')*CMmtrx*(iNbghr==iNbghr1');
    % if ~isempty(CMmtrx),pXQ=EmbdngFrmGrmn(CMmtrx,nDims);
    %     pXQ=(iNbghr0==iNbghr1')*(pXQ-pXQ(iNbghr1==ii,:));
    %     pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),pXQ];
    % end
    % -----------------------------
    % for jj=1:size(CMmtrx,1),CMmtrx=min(CMmtrx,CMmtrx(:,jj)+CMmtrx(jj,:));end
    CMmtrx=(iNbghr1==iNbghr')*CMmtrx*(iNbghr==iNbghr1');
    for jj=1:size(CMmtrx,1),CMmtrx=min(CMmtrx,CMmtrx(:,jj)+CMmtrx(jj,:));end
    dsdQ(arx:arx+length(iNbghr1)-1)=CMmtrx(:,iNbghr1==ii);
    CMmtrx=CMmtrx.^2;
    if ~isempty(CMmtrx)
    % dstMtrx=pdist(pX(slct,:),fncD);CMmtrx=squareform(dstMtrx);
    % CMmtrx=1/2*((CMmtrx(:,1)+CMmtrx(1,:))-CMmtrx);enas=ones(size(CMmtrx,1),1);
    CMmtrx=1/2*((CMmtrx(:,iNbghr1==ii)+CMmtrx(iNbghr1==ii,:))-CMmtrx);
    % CMmtrx=CMmtrx.*(1-eye(size(CMmtrx)));enas=ones(size(CMmtrx,1),1);tT=1;
    % CMmtrx=[[CMmtrx,enas];[enas',0]];
    % dmsLcl(ii)=rank(CMmtrx);
    if wDfrm,dY=pY(iNbghr1)-pY(ii);CMmtrx=CMmtrx-(((CMmtrx*dY)*(CMmtrx*dY)')./(dY'*CMmtrx*dY));end
    [Q,s,~]=svd(CMmtrx);%Q=(iNbghr0==iNbghr1')*Q;
    s=diag(s);
    if wVlmFrm
    volC=ConvHullApprx(Q(:,1+offst:nDims+offst).*sqrt(s(1+offst:nDims+offst)'),7,1e3,1e2,1e-2,1e-2);
    s=(s./sqrt(volC));
    end
    s0=1./sqrt(s(1+offst:nDims+offst))';s0(s(1+offst:nDims+offst)==0)=0;Q=Q(:,1+offst:nDims+offst);pXz=Q.*s0;%max(abs(pXz))
    % pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),Q*s];
    pX(arx:arx+length(iNbghr1)-1,:)=[iNbghr1,ii*ones(size(iNbghr1)),pXz-0*mean(pXz)];
    end
    % -----------------------------
    % arx=arx+length(iNbghr0);
    arx=arx+length(iNbghr1);ii/nObs
end

% indsCvrd=zeros(nObs,1);indsAll=true(nObs,1);prsm=ones(1,nDims);
% ii=1;
% while any(indsAll)
% indsCvrd(ii)=find(indsAll,1,'first');tel=ii;indsAll(indsCvrd(ii))=false;
% while (ii<=tel),arx=indsCvrd(ii);
%     slct=pX(:,2)==arx;iNbghr=pX(slct,1);pXC=pX(slct,3:end);
%     iNbghr0=setdiff(iNbghr,indsCvrd(indsCvrd>0));
%     for jj=1:length(iNbghr0),slctIn=pX(:,2)==iNbghr0(jj);iNbghrIn=pX(slctIn,1);
%         pXCin=pX(slctIn,3:end);antstxa=(iNbghr==iNbghrIn');pXCin=antstxa*pXCin;%pXC0=antstxa'*pXC;
%         pXCin=pXCin(any(antstxa,2),:);pXC0=pXC(any(antstxa,2),:);
%         pXCin=pXCin-mean(pXCin);pXC0=pXC0-mean(pXC0);
%         [qL,~,qR]=svd(pXC0'*pXCin);prsm(end)=sign(det(qL*qR'));
%         pX(slctIn,3:end)=pX(slctIn,3:end)*(qR.*prsm)*qL';
%     end
%     indsAll(iNbghr0)=false;
%     indsCvrd(tel+1:tel+length(iNbghr0))=iNbghr0;tel=tel+length(iNbghr0);
%     ii=ii+1;
% end
% end

% pX=SpectrlSync(pX);

for ii=1:nDims
    frmDat(ii).E=sparse(pX(:,1),pX(:,2),pX(:,2+ii),size(dstMtrx,1),size(dstMtrx,2));
end
% dstMtrxN=max(sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2)),...
%     sparse(pX(:,2),pX(:,1),dsdQ,size(dstMtrx,1),size(dstMtrx,2)));
dstMtrxN=sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2));

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
