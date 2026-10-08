function W=LnrQdrtrMtrx(dstMtrx,nDims,pY)
if ~exist('pY','var')||isempty(pY),wDfrm=false;else,wDfrm=true;end
nDimsMx=7;
uncvrdM=sparse([],[],[],size(dstMtrx,1),size(dstMtrx,2),nnz(dstMtrx));
W=sparse([],[],[],size(dstMtrx,1),size(dstMtrx,2),nnz(dstMtrx));
[isrtIn,indxQIn]=find(max(spones(dstMtrx),speye(size(dstMtrx,1),size(dstMtrx,2))));
arx=0;
for arx=1:size(dstMtrx,1),arx
    slct=indxQIn==arx;iNbghrOutArx=isrtIn(slct);slct=any(isrtIn==iNbghrOutArx',2);
    iNbghrOut=unique(indxQIn(slct));iNbghrOut0=iNbghrOut(uncvrdM(iNbghrOut,arx)==0);
    iNbghr=iNbghrOut;iNbghr0=iNbghr;
    slct=any(indxQIn==iNbghr',2);iNbghr=union(iNbghr,isrtIn(slct));iNbghr1=iNbghr;
    slct=any(indxQIn==iNbghr',2);iNbghr=union(iNbghr,isrtIn(slct));iNbghrAll=iNbghr;
    CMmtrx=full(dstMtrx(iNbghr,iNbghr));
    CMmtrx=(1-eye(size(CMmtrx))).*max(CMmtrx,CMmtrx');
    CMmtrx=distances(graph(CMmtrx));CMmtrx=min(CMmtrx,CMmtrx');
    CMmtrx=CMmtrx.^2;
    for qq=1:size(CMmtrx,2),CMmtrx=min(CMmtrx,CMmtrx(:,qq)+CMmtrx(qq,:));end
    % CMmtrx=(iNbghr0==iNbghr')*CMmtrx*(iNbghr==iNbghr0');iNbghr=iNbghr0;
    gWgI=1/2*((CMmtrx(:,iNbghr==arx)+CMmtrx(iNbghr==arx,:))-CMmtrx);%(rank(gWgI))
    if wDfrm,dY=pY(iNbghr)-pY(arx);gWgI=gWgI-(((gWgI*dY)*(gWgI*dY)')./(dY'*gWgI*dY));end,%(rank(gWgI)),pause
    [Q,s,~]=svd(gWgI);Q=Q(:,1:nDims);s=sqrt(diag(s));s=diag(s(1:nDims));
    xMat=Q*s;xMat=(iNbghrOutArx==iNbghr')*xMat;
    % --------------------------/*test--------------------------
    % figure,hold on,axis equal
    % --------------------------test*/--------------------------
    for ii=1:length(iNbghrOut0)
        slct=indxQIn==iNbghrOut0(ii);iNbghr=isrtIn(slct)==iNbghrOutArx';iNbghr=iNbghr(any(iNbghr,2),:);
        % slct=indxQIn==iNbghrOut0(ii);iNbghr=intersect(isrtIn(slct),iNbghrOut);iNbghr0=iNbghr;
        % slct=any(indxQIn==iNbghr',2);iNbghr=union(iNbghr,isrtIn(slct));iNbghr1=iNbghr;
        % slct=any(indxQIn==iNbghr',2);iNbghr=union(iNbghr,isrtIn(slct));iNbghrAll=iNbghr;
        % CMmtrx=full(dstMtrx(iNbghr,iNbghr));
        % CMmtrx=(1-eye(size(CMmtrx))).*max(CMmtrx,CMmtrx');
        % CMmtrx=distances(graph(CMmtrx));CMmtrx=min(CMmtrx,CMmtrx');
        % CMmtrx=CMmtrx.^2;
        % for qq=1:size(CMmtrx,2),CMmtrx=min(CMmtrx,CMmtrx(:,qq)+CMmtrx(qq,:));end
        % CMmtrx=(iNbghr0==iNbghr')*CMmtrx*(iNbghr==iNbghr0');iNbghr=iNbghr0;
        % --------------------------/*test--------------------------
        % [iNbghrOut0(ii),arx],(iNbghr*iNbghrOutArx)',plot(xMat(:,1),xMat(:,2),'.')
        % --------------------------test*/--------------------------
        if size(iNbghr,1)<=nDims,volC=0;
        else
            % gWgI=1/2*((CMmtrx(:,iNbghr==iNbghrOut0(ii))+CMmtrx(iNbghr==iNbghrOut0(ii),:))-CMmtrx);%(rank(gWgI))
            % if wDfrm,dY=pY(iNbghr)-pY(iNbghrOut0(ii));gWgI=gWgI-(((gWgI*dY)*(gWgI*dY)')./(dY'*gWgI*dY));end,%(rank(gWgI)),pause
            % [Q,s,~]=svd(gWgI);Q=Q(:,1:nDims);s=sqrt(diag(s));s=diag(s(1:nDims));xMat=Q*s;
            volC=ConvHullApprx(iNbghr*xMat,nDimsMx-1,1e3,1e2,0.01,0.01);
        end
        % --------------------------/*test--------------------------
        % [k,v]=convhulln(iNbghr*xMat);(iNbghr(k,:)*iNbghrOutArx)',v
        % plot(iNbghr(k,:)*xMat(:,1),iNbghr(k,:)*xMat(:,2)),pause
        % --------------------------test*/--------------------------
        if iNbghrOut0(ii)==arx,W(arx,arx)=2*volC;uncvrdM(arx,arx)=1;
        else,W(iNbghrOut0(ii),arx)=volC;%W(arx,iNbghrOut0(ii))=volC;
            uncvrdM(iNbghrOut0(ii),arx)=1;%uncvrdM(arx,iNbghrOut0(ii))=1;
        end
    end
    % --------------------------/*test--------------------------
    % close
    % --------------------------test*/--------------------------
end