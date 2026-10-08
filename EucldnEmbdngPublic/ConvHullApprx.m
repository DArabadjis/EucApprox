function volC=ConvHullApprx(pX,nDimsIn,nDrctMx,nDrctMn,eps_rel,cnfdnc)

if ~exist('nDimsIn','var')||isempty(nDimsIn),nDimsIn=5;end
if ~exist('nDrctMx','var')||isempty(nDrctMx),nDrctMx=5e3;end
if ~exist('nDrctMn','var')||isempty(nDrctMn),nDrctMn=5e2;end
if ~exist('eps_rel','var')||isempty(eps_rel),eps_rel=0.05;end
if ~exist('cnfdnc','var')||isempty(cnfdnc),cnfdnc=0.05;end

pX=pX(:,1:min(end,size(pX,1)-1));
[nElms,nDims]=size(pX);
if nDims<=nDimsIn,try,[~,volC]=convhulln(pX,{'Qt','Qs','Qx'});catch,volC=0;end
else,nDrct=0;vlsAll=[];pX=pX-mean(pX);
    while true
        if nDrct>=nDrctMx,warning('Reached mxIters before meeting tolerance');break;end
        nDrct=nDrct+1;
        [Qp,~]=qr(randn(nDims,nDimsIn),'econ');
        try,[~,vlsAll(nDrct)]=convhulln(pX*Qp,{'Qt','Qs','Qx'});catch,vlsAll(nDrct)=0;end
        volC=mean(vlsAll);
        nDrctU=ceil((var(vlsAll,1)/nDrct/(eps_rel*volC)^2)*log(2/cnfdnc));
        if nDrct>=nDrctU && nDrct>=0*nDrctMn,break;end
    end
end


