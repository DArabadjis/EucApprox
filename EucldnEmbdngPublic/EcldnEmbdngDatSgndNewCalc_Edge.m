function [pXN,errsLcl,QBst,gDatBst,lplcOpD,LclFrm,DivMtrx]=EcldnEmbdngDatSgndNewCalc_Edge(wDat,pX,nDims,wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut,wPrllQ)
% -------------------------------------------------------
% Inputs
% wDat : graph adjacency matrix - sparse #points x #points
% pX : local transports in the form [edge_row,edge_col,frame] - edges x intrinsic dim + 2 
% nDims : embedding dimension
% wVlmFrm : flag for volume normalization (false)
% wRcrsSlv : flag for iterative sparse solver (false)
% wAMG : flag for multigrid (false)
% wPlot : flag for plotting 3D projections of the embedding
% tlrnc : convergence tolerancre (1e-15)
% wErrOut : flag for keeping per iteration strains (true)

% Outputs
% pXN : the embedding - #points x nDims
% errsLcl : local strain - [strains of iterations,strain of output]
% QBst : orthogonal alignment matrix - nDims x nDims x #points 
% lplcOpD : Laplacian op
% LclFrm : local metric-trivializing frames - nDims tuple of sparse #points x #points, LclFrm(k).E 
% DivMtrx : Divergence op

if ~exist('wPlot','var')||isempty(wPlot),wPlot=true;end
if ~exist('wVlmFrm','var')||isempty(wVlmFrm)||~islogical(wVlmFrm),wVlmFrm=false;end
if ~exist('wRcrsSlv','var')||isempty(wRcrsSlv),wRcrsSlv=size(wDat,1)>2000;wRcrsSlv=false;end
if ~exist('wAMG','var')||isempty(wAMG),wAMG=false;end
if ~exist('tlrnc','var')||isempty(tlrnc)||tlrnc<1e-15,tlrnc=1e-15;end
if ~exist('wErrOut','var')||isempty(wErrOut),wErrOut=~false;end
if ~exist('wPrllQ','var')||isempty(wPrllQ),wPrllQ=true;end

tlrnc=10^floor(log10(tlrnc));
orioConv=0.1;orioConvL=5*tlrnc;nFailsIn=inf;orioConvL00=tlrnc;
orioConvL11=5*sqrt(tlrnc);nFailsIn11=100;

% orioConv=0.1;orioConvL=5e-15;nFailsIn=inf;orioConvL00=1e-15;
% orioConvL11=5e-7;nFailsIn11=100;
% -------------------------------------------------------

wDat0=wDat;nObsAll=size(wDat,1);
for ii=1:size(wDat0,2),wDat0(ii,ii)=0;end
[isrt(:,1),isrt(:,2)]=find(wDat0);nDims0=sum((spones(wDat0)));
isrt0=nan(size(isrt));isrt0(:,2)=isrt(:,2);
for ii=1:size(wDat,1),slct=isrt(:,2)==ii;isrt0(slct,1)=(1:nDims0(ii));end
if size(pX,2)>(nDims+2),pXAA=pX(:,1:nDims+2);
    for ii=1:nObsAll,ii/nObsAll,
        slct=pX(:,2)==ii;dfDat=pX(slct,3:end);
        [vVcs,sVls]=svd(full(dfDat));diafSz=size(sVls,2)-nDims;
        if diafSz<0,vVcs=[vVcs,zeros(size(vVcs,1),-diafSz)];
            sVls=[sVls,zeros(size(sVls,1),-diafSz)];sVls=[sVls;zeros(-diafSz,size(sVls,2))];
        end
        % [vVcs,sVls]=svds(dfDat,nDimSing,'largest','MaxIterations',MxItrr);
        vVcs=vVcs(:,1:nDims);sVls=diag(sVls);sVls=diag(sVls(1:nDims));
        pXAA(slct,3:end)=vVcs*((sVls));
    end
    pX=pXAA;
end
for ii=1:size(pX,2)-2
    pXDat(ii).E=sparse(pX(:,1),pX(:,2),pX(:,ii+2),size(wDat,1),size(wDat,2));
end

[lplcOpD,LclFrm,wQdrtr,DivMtrx]=DvrgncOpEstmWdfrm_epekGeit_fst(max(spones(wDat),speye(size(wDat,1),size(wDat,2))),wDat,nDims,[],wVlmFrm);
wPrllInv=false;pPrllInv=[];
if wRcrsSlv,%lplcOpD=lplcOpD-(eigs(lplcOpD,1,'smallestabs'))*speye(size(lplcOpD,1));(eigs(lplcOpD,1,'smallestabs'))
    warning('off','MATLAB:pcg:tooSmallTolerance');warning('off','MATLAB:minres:tooSmallTolerance');
    if wAMG
        %lplcOpDMdfd=lplcOpD-(eigs(lplcOpD,1,'smallestabs'))*speye(size(lplcOpD,1));
        amgClss=hypre_amg_mex('create',lplcOpD);prcndM=@(v)hypre_amg_mex('apply',amgClss,v);
        % lplcOpD=lplcOpD-mean(lplcOpD);lplcOpD=lplcOpD-mean(lplcOpD,2);
        % prcndM=@(v)[hypre_amg_mex('apply',amgClss,v(1:nObsAll));v(nObsAll+1)/nObsAll];
        % enas=ones(size(lplcOpD,1),1);lplcOpD=[lplcOpD,enas;enas',0];
    else,
        lplcOpDMdfd=lplcOpD-(eigs(lplcOpD,1,'smallestabs'))*speye(size(lplcOpD,1));(eigs(lplcOpDMdfd,1,'smallestabs'))
        try,opts.type='ict';opts.droptol=1e-3;opts.diagcomp=1e-2;prcndM=ichol(lplcOpDMdfd,opts);
		if wPrllQ
        tstTm=2*rand(nObsAll,1)-1;tolInv=1e-15*std(tstTm);
        tic;pcg(lplcOpD,tstTm,tolInv,2000,prcndM,prcndM');tstTm=toc;
        if (nDims*tstTm)>1,wPrllInv=true;if isempty(gcp("nocreate")),parpool("threads",nDims);end,end
		end
        catch lthosMhn, lthosMhn.identifier
	%tstTm=2*rand(nObsAll,nDims)-1;tic;lsqminnorm(lplcOpD,tstTm);tstTm=toc
        %if (tstTm)>2.5, 'hey AMG!'
	%    wPrllInv=true;if isempty(gcp("nocreate")),delete(gcp("nocreate"));end,parpool("Processes",nDims);
	%    wAMG=true;parfor kPrcss=1:nDims,amgClss=hypre_amg_mex('create',lplcOpD);prcndM(kPrcss).M=@(v)hypre_amg_mex('apply',amgClss,v);end
	%else,wRcrsSlv=false;prcndM=[];
	%end
	errTrm=sum(abs(lplcOpD),2);errTrm=(max(errTrm-diag(lplcOpD),0))+1.5*(2.2e-16)*max(errTrm);
	lplcOpDMdfd=lplcOpD+errTrm.*speye(size(lplcOpD,1));(eigs(lplcOpDMdfd,1,'smallestabs'))
	opts.type='ict';opts.droptol=1e-3;opts.diagcomp=1e-2;prcndM=ichol(lplcOpDMdfd,opts);
		if wPrllQ
        tstTm=2*rand(nObsAll,1)-1;tolInv=1e-15*std(tstTm);
        tic;pcg(lplcOpD,tstTm,tolInv,2000,prcndM,prcndM');tstTm=toc;
        if (nDims*tstTm)>1,wPrllInv=true;if isempty(gcp("nocreate")),parpool("threads",nDims);end,end
		end
	%wRcrsSlv=false;prcndM=[];
        end
    end
else,prcndM=[];
end
[fOut,Qout,gDatOut,lplcOpDOut]=EcldnMppngGrphDSprsSgndQnewCalc_Edge(pXDat,wDat,nDims,[],[],lplcOpD,LclFrm,wQdrtr,DivMtrx,isrt,wRcrsSlv,wAMG,prcndM,wPrllInv);
% [fOut,Qout,gDatOut,lplcOpDout]=EcldnMppngGrphDSprsSgndQ_Edge(pX,wDat,nDims,[],[],[],szNbghbrhd,isrt,rdcDms,wNrmlzd);

% [fOut,Qout,gDatOut]=EcldnMppngGrphCrr((dstMtrx),wDat,nDims);
errMin=inf;dErrMin=inf;
if wPlot
if nDims==2,figure,plot(fOut(:,1),fOut(:,2),'.'),axis equal,drawnow
elseif nDims==3,figure,plot3(fOut(:,1),fOut(:,2),fOut(:,3),'.'),axis equal,drawnow
else [~,~,QRf]=svd(fOut);fPlot=fOut*QRf;plot3(fPlot(:,1),fPlot(:,2),fPlot(:,3),'.'),axis equal,drawnow
end
end,iters=0;errsLcl=[];
while true,iters=iters+1;
tic,[fOutN,QoutN,gDatOutN]=EcldnMppngGrphDSprsSgndQnewCalc_Edge(fOut,wDat,nDims,Qout,gDatOut,lplcOpDOut,LclFrm,wQdrtr,DivMtrx,isrt,wRcrsSlv,wAMG,prcndM,wPrllInv);toc
% if wPrll,
%     % if isempty(gcp('nocreate')),parpool('threads',24);end
%     tic,[fOutN,QoutN,gDatOutN]=EcldnMppngGrphDSprsSgndQfst_prll(fOut,wDat,nDims,Qout,gDatOut,lplcOpDout,szNbghbrhd,isrt,isrt0,nDims0);toc
% else,tic,[fOutN,QoutN,gDatOutN]=EcldnMppngGrphDSprsSgndQfst(fOut,wDat,nDims,Qout,gDatOut,lplcOpDout,szNbghbrhd,isrt,isrt0,nDims0,rdcDms,wNrmlzd);toc
% end
% % [fOutN,QoutN,gDatOutN]=EcldnMppngGrphDSprsSgndQ(fOut,wDat,nDims,Qout,gDatOut,lplcOpDout,szNbghbrhd,isrt);
% % [fOutN,QoutN,gDatOutN]=EcldnMppngGrphCrr(fOut,wDat,nDims,Qout,gDatOut);
errAll=0;%for ii=1:size(QoutN,3),errAll=errAll+1/size(QoutN,3)*(QoutN(:,:,ii)'*Qout(:,:,ii)-eye(size(Qout(:,:,ii)))).^2;end
for ii=1:size(QoutN,3),errAll=errAll+1/size(QoutN,3)*...
        (QoutN(:,:,ii)/Qout(:,:,ii)-eye(size(Qout(:,:,ii)))).^2;
end
errAll=sqrt(mean(errAll,'all')),fOut=fOutN;Qout=QoutN;gDatOut=gDatOutN;

% errMin-errAll
% if ((errMin-errAll)>orioConvL)||(errMin>orioConv)||(iters<4)
%     fBst=fOut;errMin=errAll;flDIn=0;
% else,flDIn=flDIn+1;
%     if (flDIn>(nFailsIn*nDims))||(abs(errMin-errAll)<orioConvL00),break;end
% end
dErrAll=errMin-errAll
if (dErrAll>orioConvL)||(errMin>orioConv)||(iters<4)
    fBst=fOut;flDIn=0;if true||dErrAll>0,errMin=errAll;end
    QBst=Qout;gDatBst=gDatOut;
else,flDIn=flDIn+1;
    if (flDIn>(nFailsIn*nDims))||(abs(dErrAll)<orioConvL00),break;end
end
% % if (dErrAll>0)&&(errMin<orioConv)&&(iters>4)
% % if (abs(dErrMin-dErrAll)>orioConvL00),dErrMin=dErrAll;flDIn0=0;
% % else,flDIn0=flDIn0+1;if flDIn0>(nFailsIn*nDims),break;end
% % end
% % end
if (dErrAll>0)&&(errMin<orioConv)&&(iters>4)
if (dErrMin>dErrAll),dErrMin=dErrAll;flDIn0=0;
else,flDIn0=flDIn0+((dErrAll-dErrMin)<(orioConvL)),if flDIn0>(nFailsIn*nDims^2),break;end
end
else,flDIn0=0;
end
if (dErrAll<0)&&(errMin<orioConvL11)&&(iters>(nFailsIn11*nDims))
if (abs(dErrAll)<errMin),flDIn1=0;
else,flDIn1=flDIn1+1,if flDIn1>(nFailsIn11*nDims^2),break;end
end
else,flDIn1=0;
end

if wPlot
if nDims==2,plot(fOut(:,1),fOut(:,2),'.'),axis equal,drawnow
elseif nDims==3,plot3(fOut(:,1),fOut(:,2),fOut(:,3),'.'),axis equal,drawnow
else [~,~,QRf]=svd(fOut);fPlot=fOut*QRf;plot3(fPlot(:,1),fPlot(:,2),fPlot(:,3),'.'),axis equal,drawnow
end
end
if wErrOut,errsLcl(:,iters)=[LclStrainEval(fOutN,DivMtrx,gDatOutN);errAll];end
iters
end
errsLcl=[errsLcl,[LclStrainEval(fBst,DivMtrx,gDatOutN);errAll]];
pXN=fBst;%if wPrll,delete(gcp('nocreat]e'));end
if wRcrsSlv&&wAMG,%hypre_amg_mex('destroy', amgClss);
    warning('on','MATLAB:pcg:tooSmallTolerance');
    warning('on','MATLAB:minres:tooSmallTolerance');
end

function errLclStrain=LclStrainEval(fOut,DivMtrx,gDat0)
[nObsAll,nDims]=size(fOut);
dfDat0=zeros(nObsAll,nDims,nDims);
for ii=1:nDims,ug=DivMtrx(ii).E'*fOut-sum(DivMtrx(ii).E)'.*fOut;
    for jj=1:nDims
        dfDat0(:,:,jj)=dfDat0(:,:,jj)+(ug.*ug(:,jj))+(squeeze(gDat0.m(:,ii,:)).*gDat0.m(:,ii,jj))-2*(ug.*gDat0.m(:,ii,jj));
    end
end
errLclStrain=1/2*mean(diag(squeeze(mean(dfDat0))));