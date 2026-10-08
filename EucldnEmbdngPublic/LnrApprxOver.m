function [fApprx,kSBst,errBst]=LnrApprxOver(pX,pY,tStRto)
vNrm=sqrt(max(sum(pX.^2,2)));
[V00,S00,W00]=svd(pX','econ');kSmx=max(find(diag(S00)>=vNrm,1,'last'),2);
Nx=size(pX,1);lngthTst=floor(tStRto*Nx)+1;kSmx=min(lngthTst,kSmx);
% for kS=2:kSmx
errApprx=inf(kSmx-1,ceil(1/tStRto));errBst=inf(kSmx-1,1);
for kRps=1:ceil(1/tStRto)
[~,indsRnd]=sort(rand(Nx,1));[~,indsRnd0]=sort(indsRnd);
trnSt=indsRnd0>lngthTst;[V0,S0,W0]=svd(pX(trnSt,:)','econ');
vNrmIn=sqrt(max(sum(pX(trnSt,:).^2,2)));
kSmxIn=min(kSmx,max(2,find(diag(S0)>=vNrmIn,1,'last')));
for kS=2:kSmxIn
S=diag(S0);orioS=S(kS);
V=V0(:,S>=orioS);W=W0(:,S>=orioS);S=diag(S(S>=orioS));
fApprx=pY(trnSt,:)'*W*diag(1./diag(S))*V';
fApprx0=pX(~trnSt,:)*fApprx';
errApprx(kS-1,kRps)=mean((fApprx0-pY(~trnSt,:)).^2,'all');
end,%[kS,errApprx/ceil(1/tStRto)]
% if errApprx<errBst,errBst=errApprx;kSBst=kS;end
end
% [errBst,kSBst]=min(mean(errApprx,2));
for kS=1:size(errApprx,1)
    errBst(kS)=mean(errApprx(kS,~isinf(errApprx(kS,:))));
end
[errBst,kSBst]=min(errBst);
S=diag(S00);orioS=S(kSBst);
V=V00(:,S>=orioS);W=W00(:,S>=orioS);S=diag(S(S>=orioS));
fApprx=pY'*W*diag(1./diag(S))*V';
errBst=errBst/ceil(1/tStRto);