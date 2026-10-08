function dstMtrxOut=mkCnnctd(pX,dstMtrx,fncD)
[isrt,indxQ,dsdQ]=find(dstMtrx-spdiags(spdiags(dstMtrx,0),0,size(dstMtrx,1),size(dstMtrx,2)));
dkts0=(1:size(dstMtrx,1))';
isrt=[dkts0;isrt];indxQ=[dkts0;indxQ];dsdQ=[zeros(size(dstMtrx,1),1);dsdQ];
cnnCmp=CmptCnnctCmps(dsdQ,isrt,indxQ,[],true);
cnnCmp=cnnCmp(:,end); dstMtrxOut=dstMtrx;
while true,length(unique(cnnCmp))
if any(cnnCmp~=cnnCmp(1))
    [cmpLbls0,~,cnnCmp]=unique(cnnCmp);cmpLbls=length(cmpLbls0);
    [~,iDtx]=sort(sum(cnnCmp==[1:cmpLbls]));antstxa=[1:cmpLbls];
    dstMinAll=inf(cmpLbls,1);srcs=zeros(cmpLbls,1);trgts=zeros(cmpLbls,1);
    for ii=1:cmpLbls-1,ind0=find(cnnCmp==iDtx(ii));sma0=pX(ind0,:);
        for jj=1:cmpLbls,ind1=find(cnnCmp==iDtx(jj));sma1=pX(ind1,:);dstMin=inf;
            for iPnts=1:size(sma1,1)
                [dst01,indMin01]=min(fncD(sma1(iPnts,:),sma0));
                if dst01<dstMin,dstMin=dst01;ind0Min=ind0(indMin01);ind1Min=ind1(iPnts);end
            end
            if dstMin<dstMinAll(ii),dstMinAll(ii)=dstMin;srcs(ii)=ind0Min;trgts(ii)=ind1Min;antstxa(iDtx(ii))=iDtx(jj);end
            if dstMin<dstMinAll(jj),dstMinAll(jj)=dstMin;srcs(jj)=ind1Min;trgts(jj)=ind0Min;antstxa(iDtx(jj))=iDtx(ii);end
        end
    end
    for ii=1:cmpLbls,dstMtrxOut(trgts(ii),srcs(ii))=dstMinAll(ii);
        dstMtrxOut(srcs(ii),trgts(ii))=dstMinAll(ii);
        cmpLbls0(iDtx(end-ii+1))=cmpLbls0(antstxa(iDtx(end-ii+1)));
    end,cnnCmp=cmpLbls0(cnnCmp);
else,break;
end
end