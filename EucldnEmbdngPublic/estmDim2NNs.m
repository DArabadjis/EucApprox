function [dmEstm,uDim,vDim]=estmDim2NNs(dstMtrx,wplot)
if ~exist('wplot','var')||isempty(wplot),wplot=false;end
nObs=size(dstMtrx,2);
rdNN2=inf(nObs,2);
for ii=1:size(dstMtrx,2),rdNN2(ii,:)=mink(dstMtrx(dstMtrx(:,ii)>0,ii),2);end

rdNN2=rdNN2(:,2)./rdNN2(:,1);[~,HrdNN2]=sort(rdNN2);[~,HrdNN2]=sort(HrdNN2);
uDim=log(rdNN2);vDim=-log(1-HrdNN2./nObs);uDim=uDim(~isinf(vDim));vDim=vDim(~isinf(vDim));
tAng=1/2*atan2(2*mean(uDim.*vDim),mean(uDim.^2)-mean(vDim.^2));
dmEstm=tan(tAng);
if wplot
    figure,plot(uDim,vDim,'.'),hold on
    plot([0,max(uDim)],[0,dmEstm*max(uDim)]);axis equal;
end