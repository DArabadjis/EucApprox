function errsLcl=cmptMtrc(pXQ,dstMtrx)

if ~issparse(dstMtrx),dstMtrx=sparse(dstMtrx);end
nObs=size(dstMtrx,1);
dstMtrx=spdiags(zeros(nObs,1),0,dstMtrx);[isrt,indxQ]=find(dstMtrx);

nDims=min(sum(spones(spones(dstMtrx)*spones(dstMtrx))))
offst=0;
errGrm=inf(nObs,1);errAng=errGrm;errDst=errGrm;rcll=errGrm;cntnt=errGrm;
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
    CMmtrx=(iNbghr1==iNbghr')*CMmtrx*(iNbghr==iNbghr1');CMmtrxD=CMmtrx(:,iNbghr1==ii);
    
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
    % pX(arx:arx+length(iNbghr0)-1,:)=[iNbghr0,ii*ones(size(iNbghr0)),Q*s];
    
    CMmtrxEstm=pXQ(iNbghr1,:)-pXQ(ii,:);CMmtrxEstmD=sqrt(sum(CMmtrxEstm.^2,2));
    CMmtrxEstmAng=acos((CMmtrxEstm./CMmtrxEstmD)*(CMmtrxEstm./CMmtrxEstmD)');
    CMmtrxEstmInv=CMmtrxEstm'*Q/s;CMmtrxEstmInv=CMmtrxEstmInv*CMmtrxEstmInv';
    errGrm(ii)=norm(CMmtrxEstmInv-eye(size(CMmtrxEstmInv)),"fro");

    CMmtrxAng=acos(1./CMmtrxD.*(CMmtrx./CMmtrxD'));
    CMmtrxAng=CMmtrxAng(iNbghr1~=ii,iNbghr1~=ii);
    CMmtrxEstmAng=CMmtrxEstmAng(iNbghr1~=ii,iNbghr1~=ii);
    errAng(ii)=sum(abs(CMmtrxEstmAng-CMmtrxAng),'all')/(size(CMmtrxAng,1)*(size(CMmtrxAng,1)-1));
    errDst(ii)=mean(abs(CMmtrxEstmD(iNbghr1~=ii)./CMmtrxD(iNbghr1~=ii)-1));
    end
    dstAll=sqrt(sum((pXQ-pXQ(ii,:)).^2,2));[~,ipXsrt]=sort(dstAll);[~,ipXsrt]=sort(ipXsrt);
    cntnt(ii)=2/(2*nObs-3*length(iNbghr1)-1)*...
        mean(max(0,ipXsrt(iNbghr1)-length(iNbghr1)));
    cmmn=ipXsrt(iNbghr1)<=length(iNbghr1);rcll(ii)=mean(cmmn);
end
errsLcl.errDst=errDst;errsLcl.errAng=errAng;errsLcl.errGrm=errGrm/(size(pXQ,2)^2);
errsLcl.cntnt=cntnt;errsLcl.rcll=rcll;