function [fOut,Qout,gDatOut,lplcOpD]=EcldnMppngGrphDSprsSgndQnewCalc_Edge(fDat, ...
    wDat,nDims,Qmat,gDatIn,lplcOpD,LclFrm,wQdrtr,DivMtrx,isrt,wRcrsSlv,wAMG,prcndM,wPrllInv)
% fDat -> [target,source,edge-features]
nObsAll=size(wDat,1);
if ~exist('isrt','var')||isempty(isrt)
    wDat0=wDat;[isrt(:,1),isrt(:,2)]=find(wDat0);
end

if ~exist('wRcrsSlv','var')||isempty(wRcrsSlv),wRcrsSlv=nObsAll>2000;wRcrsSlv=false;end
if ~exist('wAMG','var')||isempty(wAMG),wAMG=false;end
if ~exist('wPrllInv','var')||isempty(wPrllInv),wPrllInv=false;end

if wRcrsSlv&&(~exist('prcndM','var')||isempty(prcndM))
    opts.type= 'ict';opts.droptol=1e-3;opts.diagcomp=1e-3;
    prcndM=ichol(lplcOpD,opts);
end

%gDat=nan(nDims,nDims,nObsAll);
if ~exist('Qmat','var')||isempty(Qmat),Qmat=repmat(eye(nDims),1,1,nObsAll);
    % lplcOpD=lsqminnorm(lplcOpD,speye(size(lplcOpD,1),size(lplcOpD,2)));
    % geitBin=spones(max(wDat,speye(size(wDat,1),size(wDat,2))));wQdrtrN=(wQdrtr*geitBin)';
    gDat0.m=zeros(nObsAll,length(LclFrm),length(fDat));%gDat0.mg=gDat0.m;
    for ii=1:size(gDat0.m,2),u=LclFrm(ii).E;%ug=LclFrm(ii).E.*wQdrtrN;
        for jj=1:size(gDat0.m,3)
            gDat0.m(:,ii,jj)=sum(u.*fDat(jj).E);
            % gDat0.mg(:,ii,jj)=sum(ug.*fDat(jj).E);
            % gDat0.m(:,ii,jj)=Qmat(ii,jj,:);
        end
    end
else
    gDat0=gDatIn;%gDat=gDatIn;
end
gOut=zeros(size(gDat0.m,1),size(gDat0.m,3));
for iDims=1:length(DivMtrx)
    % gOut=gOut+DivMtrx(iDims).E*squeeze(gDat0.m(:,iDims,:));
    gOut=gOut+DivMtrx(iDims).E*squeeze(gDat0.m(:,iDims,:))-sum(DivMtrx(iDims).E)'.*squeeze(gDat0.m(:,iDims,:));
    % gOut=gOut+0.5*(DivMtrx(iDims).E*squeeze(gDat0.m(:,iDims,:))+LclFrm(iDims).E*squeeze(gDat0.mg(:,iDims,:)));
end
if wRcrsSlv,fOut=nan(size(gOut));mxItr=100;gOut=gOut-mean(gOut);
    tolInv=1e-15*std(gOut);
    if wPrllInv
        parfor kFout=1:size(fOut,2)
        if wAMG,[fOut(:,kFout),flag,relres]=pcg(lplcOpD,gOut(:,kFout),tolInv(kFout),mxItr,prcndM(kFout).M);
        else,[fOut(:,kFout),flag,relres]=pcg(lplcOpD,gOut(:,kFout),tolInv(kFout),2000,prcndM,prcndM');
	end
        %relres
        end
    else
    for kFout=1:size(fOut,2)
        if wAMG,[fOut(:,kFout),flag,relres]=pcg(lplcOpD,gOut(:,kFout),tolInv(kFout),mxItr,prcndM);
        else,[fOut(:,kFout),flag,relres]=pcg(lplcOpD,gOut(:,kFout),tolInv(kFout),2000,prcndM,prcndM');
        end
        %relres
    end
    end
else,fOut=lsqminnorm(lplcOpD,gOut);
end
% if sPwDat,fOut=lsqminnorm(lplcOpD,-gOut);else,fOut=lplcOpD*(-gOut);end
% fOut=lplcOpD*(gOut);

fOut=fOut-mean(fOut);
dfDat0=zeros(nObsAll,nDims,nDims);
for ii=1:length(LclFrm),ug=DivMtrx(ii).E'*fOut-sum(DivMtrx(ii).E)'.*fOut;%u=LclFrm(ii).E'*fOut;
    for jj=1:size(gDat0.m,3)
        dfDat0(:,:,jj)=dfDat0(:,:,jj)+ug.*(gDat0.m(:,ii,jj));
        % dfDat0(:,:,jj)=dfDat0(:,:,jj)+0.5*(ug.*(gDat0.m(:,ii,jj))+u.*(gDat0.mg(:,ii,jj)));
    end
end
Qout=nan(size(Qmat));sAllRot=0;
% for ii=1:size(fOut,1)
%     slct=isrt(:,2)==ii;iNbghr=isrt(slct,1);
%     dfDat0=(fOut(iNbghr,:)-fOut(ii,:)).*(wDatQQ(iNbghr,ii));
%     dfDat=gDat0(ii).mg'*dfDat0;sAllRot=sAllRot+dfDat;
%     [vVcs,~,uUVcs]=svd(dfDat);
%     dDgns=[ones(1,nDims-1),(sign(det(Qmat(:,:,ii)'*vVcs*uUVcs')))];
%     % Qout(:,:,ii)=(Qmat(:,:,ii)'*vVcs*diag(dDgns)*uUVcs')';
%     Qout(:,:,ii)=uUVcs*diag(dDgns)*vVcs';
% end
for ii=1:size(fOut,1)
    dfDat=squeeze(dfDat0(ii,:,:))';sAllRot=sAllRot+dfDat;
    [vVcs,~,uUVcs]=svd(dfDat);
    dDgns=[ones(1,nDims-1),(sign(det(Qmat(:,:,ii)'*vVcs*uUVcs')))];
    % Qout(:,:,ii)=(Qmat(:,:,ii)'*vVcs*diag(dDgns)*uUVcs')';
    Qout(:,:,ii)=uUVcs*diag(dDgns)*vVcs';
    % sAllRot=sAllRot+Qout(:,:,ii)'*(dfDat'*dfDat);
end
[vVcs,~,uUVcs]=svd(sAllRot);dDgns=[ones(1,nDims-1),(sign(det(vVcs*uUVcs')))];
QoutAll=uUVcs*diag(dDgns)*vVcs';Qout=pagemtimes(Qout,QoutAll');
Qout=pagemtimes(Qout,Qmat);
% -----------------------------------------
% iNan=find(squeeze(any(any(isnan(Qout)|isinf(Qout)))));
% wDat0=wDat;wDat0(iNan,iNan)=-inf;[~,isrt0]=sort(1-wDat0);
% Qout0=Qout;nghbrs=nDims;
% for ii=1:length(iNan),iNbghr=isrt0(1:nghbrs,iNan(ii));
%     Qslct=0;
%     for jj=1:nghbrs,Qslct=Qslct+1/nghbrs*logm(Qout(:,:,iNbghr(jj)));end
%     Qout0(:,:,iNan(ii))=expm(Qslct);Qout0(:,:,iNan(ii))*Qout0(:,:,iNan(ii))'
% end,Qout=Qout0;%gDatOut=gDat;
% % for ii=1:nObsAll
% %     gDatOut(:,:,ii)=gDat(:,:,ii)*Qmat(:,:,ii)*Qout(:,:,ii)';
% % end
% -----------------------------------------
u=permute(gDat0.m,[2,3,1]);
gDatOut.m=permute(pagemtimes(u,pagemtimes(Qmat,pagetranspose(Qout))),[3,1,2]);
% ug=permute(gDat0.mg,[2,3,1]);
% gDatOut.mg=permute(pagemtimes(ug,pagemtimes(Qmat,pagetranspose(Qout))),[3,1,2]);
% -----------------------------------------
% gDatOut=gDat0;
% for ii=1:nObsAll
%     %     dQW=tensorprod(Qmat(:,:,wNbghbrhd(:,ii)),Qout(:,:,wNbghbrhd(:,ii)),2,2);
%     %     dQWslct=repmat(logical(eye(nDims0(ii))),[nDims,nDims,1,1]);
%     %     dQW=reshape(dQW(dQWslct),[nDims,nDims,nDims0(ii)]);
%     dQW=pagemtimes(Qmat(:,:,logical(wNbghbrhd(:,ii))),'none',Qout(:,:,logical(wNbghbrhd(:,ii))),'transpose');
%     gDatOut(ii).aQmQ=permute(pagemtimes(permute(gDat0(ii).aQmQ,[3,2,1]),dQW),[3,2,1]);
%     %     dQWa=Qmat(:,:,wNbghbrhd(:,ii));dQWb=Qout(:,:,wNbghbrhd(:,ii));
%     %     dQWa=reshape(permute(dQWa,[2,1,3]),[size(dQWa,2),size(dQWa,1)*size(dQWa,3)]);
%     %     dQWb=reshape(permute(dQWb,[2,1,3]),[size(dQWb,2),size(dQWb,1)*size(dQWb,3)]);
%     %     dQWa=reshape(dQWa'*dQWb);
%     gDatOut(ii).m=gDat0(ii).m*Qmat(:,:,ii)*Qout(:,:,ii)';
%     gDatOut(ii).mg=gDat0(ii).mg*Qmat(:,:,ii)*Qout(:,:,ii)';
% end
% % for ii=1:nObsAll,iNbghr=find(wNbghbrhd(:,ii));
% %     aQmQ=nan(nDims0(ii),nDims);
% %     for jj=1:nDims0(ii)%,size(wDatQQ),size(wNbghbrhd(:,iNbghr(jj)))
% %         lWdat=wDatQQ(wNbghbrhd(:,iNbghr(jj)),iNbghr(jj));
% %         aQmQ(jj,:)=gDatOut(iNbghr(jj)).mg(isrt0(ii,iNbghr(jj)),:)*...
% %             (lWdat(isrt0(ii,iNbghr(jj))));
% %     end
% %     norm(gDatOut(ii).aQmQ-aQmQ)
% % end,%pause
% % % for ii=1:size(Qout,1),for jj=1:size(Qout,2)
% % %         Qout(ii,jj,:)=gDiort(squeeze(Qout(ii,jj,:)),wDat,nDims);
% % % end,end


function fOut=gDiort(fDat,wDat,nghbrs)
fOut=fDat;slct=isinf(fOut)|isnan(fOut);wDat(slct,slct)=-inf;
[~,isrt]=sort(1-wDat);
fOut(slct)=mean(fOut(isrt(1:nghbrs,slct)))';