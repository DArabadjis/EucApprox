function errsGlbl=cmptMtrcGlbl(pXQ,dstMtrx)

orioN=1e4;nObs=size(dstMtrx,1);
if ~issparse(dstMtrx),dstMtrx=sparse(dstMtrx);end
dstMtrx=spdiags(zeros(nObs,1),0,dstMtrx);
nNghbrs=(sum(spones(dstMtrx)));
dstMtrxG=graph(dstMtrx,'omitselfloops');
trst=inf(nObs,1);prcsn=trst;prsnCorr=trst;sprmnCorr=trst;kndllCorr=trst;
if nObs<orioN
    dsts=distances(dstMtrxG);dsts=min(dsts,dsts');
    [~,iNghbrs0]=sort(dsts);[~,iNghbrs0]=sort(iNghbrs0);
    for ii=1:nObs,dpXQ=sqrt(sum((pXQ-pXQ(ii,:)).^2,2));
        [~,iNghbrs]=mink(dpXQ,nNghbrs(ii)+1);
        trst(ii)=2/(2*nObs-3*(nNghbrs(ii)+1)-1)*...
            mean(max(0,iNghbrs0(iNghbrs,ii)-(nNghbrs(ii)+1)));
        cmmn=iNghbrs0(iNghbrs,ii)<=(nNghbrs(ii)+1);prcsn(ii)=mean(cmmn);
        prsnCorr(ii)=corr(dsts(:,ii),dpXQ);sprmnCorr(ii)=corr(dsts(:,ii),dpXQ,'type','Spearman');
        kndllCorr(ii)=corr(dsts(:,ii),dpXQ,'type','Kendall');
    end
    dsts=squareform(dsts);strss=(pdist(pXQ)-dsts).^2./mean(dsts.^2);
else,strss=trst;
    for ii=1:nObs
        dsts=distances(dstMtrxG,ii);
        [~,iNghbrs0]=sort(dsts);[~,iNghbrs0]=sort(iNghbrs0);
        dstEmbd=sqrt(sum((pXQ-pXQ(ii,:)).^2,2));[~,iNghbrs]=mink(dstEmbd,nNghbrs(ii)+1);
        trst(ii)=2/(2*nObs-3*(nNghbrs(ii)+1)-1)*...
            mean(max(0,iNghbrs0(iNghbrs)-(nNghbrs(ii)+1)));
        cmmn=iNghbrs0(iNghbrs)<=(nNghbrs(ii)+1);prcsn(ii)=mean(cmmn);
        strss(ii)=(dstEmbd-dsts).^2./mean(dsts.^2);
        prsnCorr(ii)=corr(dsts(:,ii),dstEmbd);
        sprmnCorr(ii)=corr(dsts(:,ii),dstEmbd,'type','Spearman');
        kndllCorr(ii)=corr(dsts(:,ii),dstEmbd,'type','Kendall');
    end
end
errsGlbl.strss=strss;errsGlbl.trst=trst;errsGlbl.prcsn=prcsn;
errsGlbl.prsnCorr=prsnCorr;errsGlbl.sprmnCorr=sprmnCorr;errsGlbl.kndllCorr=kndllCorr;