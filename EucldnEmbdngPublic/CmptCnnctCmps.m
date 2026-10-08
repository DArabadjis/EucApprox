function [CnnCmp]=CmptCnnctCmps(dsdQ,isrt,indxQ,nNghbrs,apst)
% --------------------Inputs---------------------
% dsdQ - distances of neighbors from the source point 
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% isrt - indices of neighbors 
%       points x nDims (apst=false) / neighborhoods' points x 1 (apst=true)
% indxQ - indices of the source-points of the neighborhoods, useful for
%       sparse indexing together with isrt
% --------------------Outputs---------------------
% CnnCmp - points x merge depth, index of the connected cover of each point
%       at each level of covers' merging

if apst,CnnCmp=zeros(max(indxQ),1);dsdQMx=inf(size(CnnCmp));
    for ii=1:size(dsdQMx,1),dsdQMx(ii)=max(dsdQ(indxQ==ii));end
    % dsdQMx=splitapply(@max,dsdQ,indxQ);dsdQMx=dsdQMx./nNghbrs;
    % dsdQMx=1./nNghbrs;
else,dsdQMx=max(dsdQ,[],2);CnnCmp=zeros(size(isrt,1),1);
end
indsAll=true(size(CnnCmp));indCvrd=zeros(size(CnnCmp));
ii=0;
while any(indsAll),ii=ii+1;iFnd=find(indsAll);
[~,iArx]=min(dsdQMx(iFnd));iArx=iFnd(iArx(1));
[indCvrdN]=CrtCmpnnt(iArx,isrt,indxQ,indCvrd,apst);
slctQQ=(indCvrd==0)&(indCvrdN>0);indCvrd=indCvrdN;
CnnCmp(indCvrd(slctQQ))=ii;indsAll(indCvrd(slctQQ))=false;mean(indsAll)
end

if ~apst,isrt=isrt';isrt=isrt(:);indxQ=indxQ';indxQ=indxQ(:);end
CnnCmpN=CnnCmp;
for ii=1:size(CnnCmp,1),slct=isrt==ii;slct=any(indxQ==(indxQ(slct)'),2);
    CnnCmpN(ii)=min(CnnCmp(isrt(indxQ==ii)));
end
while any(CnnCmpN~=CnnCmp(:,end))
    CnnCmp=[CnnCmp,CnnCmpN];CnnCmpN0=CnnCmp(:,end);
    for ii=1:size(CnnCmpN0,1),slct=isrt==ii;slct=any(indxQ==(indxQ(slct)'),2);
        CnnCmpN(ii)=min(CnnCmpN0(isrt(slct)));
    end
end

function indCvrd=CrtCmpnnt(iArx,isrt,indxQ,indCvrdIn,apst)
if ~exist('indCvrdIn','var')||isempty(indCvrdIn),prthEpan=true;
else,prthEpan=false;
end
if prthEpan,indCvrd=zeros(max(indxQ),1);tel=1;ii=1;indCvrd(1)=iArx;
else,indCvrd=indCvrdIn;ii=find(indCvrd==0,1,'first');indCvrd(ii)=iArx;tel=ii;
end
if apst
while ii<=tel,iArx=indCvrd(ii);dInd=setdiff(isrt(indxQ==iArx),indCvrd(1:tel));
indCvrd(tel+1:tel+length(dInd))=dInd;ii=ii+1;tel=tel+length(dInd);tel-ii
end
else
while ii<=tel,iArx=indCvrd(ii);dInd=setdiff(isrt(iArx,:),indCvrd(1:tel));
indCvrd(tel+1:tel+length(dInd))=dInd;ii=ii+1;tel=tel+length(dInd);tel-ii
end
end