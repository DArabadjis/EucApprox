function datOut=evalEmbddngClssfr(pXQhssn,pY,Nnghbrs,fncD,fncPrb)

if ~exist('fncD','var')||isempty(fncD),fncD=@(x,y)sqrt(mean((y-x).^2,2));end
if ~exist('Nnghbrs','var')||isempty(Nnghbrs),Nnghbrs=size(pXQhssn,2)+1;end
if ~exist('fncPrb','var')||isempty(Nnghbrs),fncPrb=@(x)1-1./(1+exp(-x));end

sclX=sqrt(mean(fncD(mean(pXQhssn),pXQhssn).^2));pXQhssn=(pXQhssn-mean(pXQhssn))./sclX;
[pYun,~,pYlbls]=unique(pY);indsRnk=[1:Nnghbrs-1]';
nObs=size(pXQhssn,1);prbs=zeros(nObs,length(pYun));mnRnk=inf(nObs,length(pYun));
frqnt=prbs;
for iObs=1:nObs,[dsdMn,iMn]=mink(fncD(pXQhssn(iObs,:),pXQhssn),Nnghbrs);
    lclLbls=pY(iMn(iMn~=iObs));lclLblsBin=lclLbls==pYun';slct=any(lclLblsBin);
    prbsIn=fncPrb(dsdMn(iMn~=iObs).*lclLblsBin(:,slct));
    prbsIn0=prbsIn;prbsIn0(prbsIn==0)=1;
    prbs(iObs,slct)=exp(sum(log(prbsIn))./sum(lclLblsBin(:,slct)));
    mnRnk(iObs,slct)=sum(indsRnk.*prbsIn)./sum(prbsIn);
    frqnt(iObs,slct)=sum(prbsIn);
end
datOut.prbs=prbs;datOut.mnRnk=mnRnk;datOut.frqnt=frqnt./mean(pY==pYun');

[~,datOut.prbsY]=max(datOut.prbs,[],2);datOut.prbsAcc=mean(pYlbls==datOut.prbsY);
[~,datOut.frqntY]=max(datOut.frqnt,[],2);datOut.frqntAcc=mean(pYlbls==datOut.frqntY);
[~,datOut.mnRnkY]=min(datOut.mnRnk,[],2);datOut.mnRnkAcc=mean(pYlbls==datOut.mnRnkY);

for iLbls=1:length(pYun),mesaY=pYlbls==iLbls;
    [recall,precision,thr]=perfcurve(mesaY,datOut.prbs(:,iLbls),true,'xCrit','reca','yCrit','prec');
    [datOut.prbsF1(iLbls,1),iMx]=max(sqrt(recall.*precision));
    datOut.prbsF1(iLbls,2)=thr(iMx);datOut.prbsPR=[recall,precision,thr];

    [recall,precision,thr]=perfcurve(mesaY,1./datOut.mnRnk(:,iLbls),true,'xCrit','reca','yCrit','prec');
    [datOut.mnRnkF1(iLbls,1),iMx]=max(sqrt(recall.*precision));
    datOut.mnRnkF1(iLbls,2)=thr(iMx);datOut.mnRnkPR=[recall,precision,thr];

    [recall,precision,thr]=perfcurve(mesaY,datOut.frqnt(:,iLbls),true,'xCrit','reca','yCrit','prec');
    [datOut.frqntF1(iLbls,1),iMx]=max(sqrt(recall.*precision));
    datOut.frqntF1(iLbls,2)=thr(iMx);datOut.frqntPR=[recall,precision,thr];
end

datOut.prbs_ARI=adjusted_rand_index(pY,datOut.prbsY);
datOut.frqnt_ARI=adjusted_rand_index(pY,datOut.frqntY);
datOut.mnRnk_ARI=adjusted_rand_index(pY,datOut.frqntY);
datOut.prbs_NMI=normalized_mutual_info(pY,datOut.prbsY,'avg');
datOut.frqnt_NMI=normalized_mutual_info(pY,datOut.frqntY,'avg');
datOut.mnRnk_NMI=normalized_mutual_info(pY,datOut.frqntY,'avg');

function ari = adjusted_rand_index(y_true, y_pred)
[~,~,u] = unique(y_true);[~,~,v] = unique(y_pred);
n = numel(u);
% Contingency table n_ij
C = accumarray([u,v], 1, [max(u), max(v)]);
% Combinatorial helpers
comb2 = @(x) x.*(x-1)/2;

sum_nij_2 = sum(comb2(C), 'all');
ai = sum(C,2);bj = sum(C,1);
sum_ai_2 = sum(comb2(ai));sum_bj_2 = sum(comb2(bj));
n2 = comb2(n);

expected_index = (sum_ai_2 * sum_bj_2) / n2;
max_index = 0.5*(sum_ai_2 + sum_bj_2);

ari = (sum_nij_2 - expected_index) / (max_index - expected_index);
if isnan(ari), ari = 0; end


function nmi = normalized_mutual_info(y_true, y_pred, method)

if ~exist('method','var')||isempty(method)||~any(strcmp(method,{'sqrt','avg'})),method = 'sqrt';end

[~,~,u] = unique(y_true);[~,~,v] = unique(y_pred);
n = numel(u);

% Contingency and probabilities
C = accumarray([u,v], 1, [max(u), max(v)]);
Pij =C/n;Pi=sum(Pij,2);Pj=sum(Pij,1)';

% Mutual information (natural logs)

% vectorized MI: sum Pij log(Pij / (Pi*Pj))
[rr,cc,vals] = find(Pij); %#ok<ASGLU> rr,cc indices; vals = Pij>0 entries
mi=(vals.*log(vals./(Pi(rr).*Pj(cc))));mi(isnan(mi))=0;mi=sum(mi);

% Entropies
H1=-sum(Pi(Pi>0).*log(Pi(Pi>0)));H2 = -sum(Pj(Pj>0).*log(Pj(Pj>0)));

switch lower(method)
    case 'sqrt'
        denom = sqrt(H1*H2);
        if denom==0, nmi = 0; else, nmi = mi / denom; end
    case 'avg'
        denom = 0.5*(H1+H2);
        if denom==0, nmi = 0; else, nmi = 2*mi / (H1+H2); end
    otherwise
        error('Unknown NMI normalization "%s". Use "sqrt" or "avg".', method);
end
if isnan(nmi), nmi = 0; end

