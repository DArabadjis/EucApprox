function [lplcOpD,pX,wQdrtr,DivMtrx]=DvrgncOpEstmWdfrm_epekGeit_fst(wQdrtr,dstMtrx,nDims,dpth,wVlmFrm,wBtch,pY) 
% dstMtrx : distance graph matrix #points x #points sparse
% wQdrtr : quadrature graph matrix #points x #points sparse
% pX : struct with the local frames' edge functions, 
%      1 x dimensions stuct -> vcFld(#).E : #points(target) x #points(source) sparse
% DivMtrx : struct with the basic divergence operator pX*wQdrtr, 
%           1 x dimensions stuct -> vcFld(#).E : #points(target) x #points(source) sparse
% lplcOpD : struct with the basic Laplacian operator {DivMtrx_i*pX_i'}_i, 
%           1 x dimensions stuct -> vcFld(#).E : #points(target) x #points(source) sparse

if ~exist('wVlmFrm','var')||isempty(wVlmFrm)||~islogical(wVlmFrm),wVlmFrm=false;end
if ~exist('wBtch','var')||isempty(wBtch)||~islogical(wBtch),wBtch=true;end
if ~exist('dpth','var')||isempty(dpth)||~islogical(dpth),dpth=false;end
if ~exist('pY','var')||isempty(pY),wDfrm=false;else,wDfrm=true;end

if wDfrm,[~,pX]=LclFrmsE(dstMtrx,nDims,dpth,wVlmFrm,wBtch,pY);
else,[~,pX]=LclFrmsE(dstMtrx,nDims,dpth,wVlmFrm,wBtch);
end
if isempty(wQdrtr)
if wDfrm,wQdrtr=LnrQdrtrMtrx(dstMtrx,nDims,pY);else,wQdrtr=LnrQdrtrMtrx(dstMtrx,nDims);end
geitBin=spones(max(dstMtrx,speye(size(dstMtrx,1),size(dstMtrx,2))));wQdrtr=(min(wQdrtr,wQdrtr'));
wQdrtrN=(wQdrtr*geitBin)';
wQdrtrS=sum(wQdrtr,2)';
else,wQdrtrN=wQdrtr;
end
% wQdrtr=speye(size(wQdrtr,1),size(wQdrtr,2)).*wQdrtr;
nDimsVcFld=length(pX);
for ii=1:nDimsVcFld,DivMtrx(ii).E=pX(ii).E.*wQdrtrN;end
% for ii=1:nDimsVcFld,DivMtrx(ii).E=DivMtrx(ii).E.*spones(dstMtrx);end
lplcOpD=sparse([],[],[],size(dstMtrx,1),size(dstMtrx,2),nnz(dstMtrx));
for ii=1:nDimsVcFld
    % lplcOpD=lplcOpD+(DivMtrx(ii).E)*(pX(ii).E)';
    e0=sum(DivMtrx(ii).E);e1=DivMtrx(ii).E.*e0;
    lplcOpD=lplcOpD+(DivMtrx(ii).E)*(DivMtrx(ii).E)'-e1-e1'+spdiags((e0.^2)',0,size(dstMtrx,1),size(dstMtrx,2));
end
% lplcOpD=0.5*(lplcOpD+lplcOpD');
% lplcOpD=lplcOpD.*geitBin;
