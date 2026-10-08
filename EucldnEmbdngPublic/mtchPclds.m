function [pXInUNQry0,err]=mtchPclds(pXInUNQry,pXInUN)
[~,~,qUnQ]=svd(pXInUN);[~,~,qUnQQry]=svd(pXInUNQry);
pXInUNQry0=pXInUNQry*qUnQQry*qUnQ';cIndcs0=zeros(size(pXInUNQry0,1),1);
cIndcs=zeros(size(pXInUN,1),1);
err=inf;
while true
for ii=1:size(pXInUNQry0,1)
    [~,cIndcs0(ii)]=min(sum((pXInUNQry0(ii,:)-pXInUN).^2,2));
end
for ii=1:size(pXInUN,1)
    [~,cIndcs(ii)]=min(sum((pXInUNQry0-pXInUN(ii,:)).^2,2));
end
[qUnQ,~,qUnQQry]=svd(pXInUN(cIndcs0,:)'*pXInUNQry0+pXInUN'*pXInUNQry0(cIndcs,:));
sgVls=[ones(1,size(pXInUN,2)-1),(sign(det(qUnQ))*sign(det(qUnQQry)))];
pXInUNQry0N=pXInUNQry0*qUnQQry*diag(sgVls)*qUnQ';
errN=sqrt(mean(sum((pXInUNQry0N-pXInUN(cIndcs0,:)).^2,2))+ ...
    mean(sum((pXInUNQry0N(cIndcs,:)-pXInUN).^2,2)))
if errN<err,pXInUNQry0=pXInUNQry0N;err=errN;else,break;end
end