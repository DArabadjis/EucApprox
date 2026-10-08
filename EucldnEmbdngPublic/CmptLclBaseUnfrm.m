function [pX,dstMtrxN]=CmptLclBaseUnfrm(dstMtrx,nDims)
if ~issparse(dstMtrx),dstMtrx=sparse(dstMtrx);end
dmsLcl=zeros(size(dstMtrx,1),1);nObs=size(dstMtrx,1);
dstMtrx=spdiags(zeros(nObs,1),0,dstMtrx);[isrt,indxQ]=find(dstMtrx);
nDims=min(nDims,min(sum(spones(dstMtrx))));

nDims=min(sum(spones(spones(dstMtrx)*spones(dstMtrx))))
pX=zeros(size(isrt,1)+nObs,nDims+2);arx=1;offst=0;dsdQ=inf(size(pX,1),1);
for ii=1:nObs,slct=indxQ==ii;iNbghr=[ii;isrt(slct)];iNbghr0=iNbghr;iNbghr1=iNbghr0;
    slct=any(indxQ==iNbghr',2);iNbghr=union(iNbghr,isrt(slct));iNbghr1=iNbghr;
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
    dsdQ(arx:arx+length(iNbghr1)-1)=CMmtrx(:,iNbghr1==ii);
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
    EcldCrds=Q*s;%EcldCrds0=EcldCrds(iNbghr1~=ii,1:end-1);
    % nNghbrs=size(CMmtrx,1)-1;jCmb=nchoosek(1:nNghbrs,nDims);EcldCrds1=zeros(nNghbrs);
    % for iCmb=1:size(jCmb,1)
    %     EcldCrds1(jCmb(iCmb,:),:)=EcldCrds1(jCmb(iCmb,:),:)+...
    %         EcldCrds0(jCmb(iCmb,:),:)./det(EcldCrds0(jCmb(iCmb,:),:));
    % end
    % EcldCrds1=EcldCrds1/size(jCmb,1);EcldCrds(iNbghr1~=ii,1:end-1)=EcldCrds1;

    % sDet=diag(s);sDet=mean(log(sDet(1:end-1)));EcldCrds=EcldCrds.*exp(-sDet);
    volC=ConvHullApprx(EcldCrds,7,1e3,1e2,1e-2,1e-2);volC=(volC)^(1/min(size(EcldCrds,2),7));
    EcldCrds=EcldCrds./volC;
    % pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),Q*s];
    pX(arx:arx+length(iNbghr1)-1,:)=[iNbghr1,ii*ones(size(iNbghr1)),EcldCrds];
    
    end
    % -----------------------------
    % arx=arx+length(iNbghr0);
    arx=arx+length(iNbghr1);ii/nObs
end
% dstMtrxN=max(sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2)),...
%     sparse(pX(:,2),pX(:,1),dsdQ,size(dstMtrx,1),size(dstMtrx,2)));
dstMtrxN=sparse(pX(:,1),pX(:,2),dsdQ,size(dstMtrx,1),size(dstMtrx,2));
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