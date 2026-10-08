function datOut=evalEmbddngClssfrKfold(pXQhssn,pY,Nnghbrs,kFld,fncD,fncPrb)

if ~exist('fncD','var')||isempty(fncD),fncD=@(x,y)sqrt(mean((y-x).^2,2));end
if ~exist('Nnghbrs','var')||isempty(Nnghbrs),Nnghbrs=size(pXQhssn,2)+1;end
if ~exist('fncPrb','var')||isempty(Nnghbrs),fncPrb=@(x)1-1./(1+exp(-x));end

sclX=sqrt(mean(fncD(mean(pXQhssn),pXQhssn).^2));pXQhssn=(pXQhssn-mean(pXQhssn))./sclX;
[pYun,~,pYlbls]=unique(pY);indsRnk=[1:Nnghbrs-1]';
lblsKfld=floor(rand(size(pY,1),1)*kFld);
datOut.prbsF1=zeros(length(pYun),1);
datOut.prbs_ARI=0;datOut.prbs_NMI=0;datOut.prbsAcc=0;
for iFld=0:kFld-1,slctKfld=lblsKfld==iFld;pXQhssnIn=pXQhssn(slctKfld,:);
nObs=size(pXQhssnIn,1);prbs=zeros(nObs,length(pYun));
for iObs=1:nObs,[dsdMn,iMn]=mink(fncD(pXQhssnIn(iObs,:),pXQhssn(~slctKfld,:)),Nnghbrs);
    lclLbls=pY(iMn);lclLblsBin=lclLbls==pYun';slct=any(lclLblsBin);
    prbsIn=fncPrb(dsdMn.*lclLblsBin(:,slct));
    prbsIn0=prbsIn;prbsIn0(prbsIn==0)=1;
    prbs(iObs,slct)=exp(sum(log(prbsIn))./sum(lclLblsBin(:,slct)));
end

[~,prbsY]=max(prbs,[],2);datOut.prbsAcc=datOut.prbsAcc+log(mean(pYlbls(slctKfld)==prbsY))/kFld;%(mean(pYlbls(slctKfld)==prbsY))
thrAll=[];
for iLbls=1:length(pYun),mesaY=pYlbls(slctKfld)==iLbls;
    [recall,precision,thr]=perfcurve(mesaY,prbs(:,iLbls),true,'xCrit','reca','yCrit','prec');
    thrAll0=unique([thrAll;thr]);
    if ~isempty(thrAll)
        rcAll=interp1(thrAll,rcAll,thrAll0);prcAll=interp1(thrAll,prcAll,thrAll0);
        rcAll=[rcAll,interp1(thr(2:end),recall(2:end),thrAll0)];
        prcAll=[prcAll,interp1(thr(2:end),precision(2:end),thrAll0)];thrAll=thrAll0;
    else,thrAll=thrAll0;prcAll=precision(2:end);rcAll=recall(2:end);
    end
end
[mxF1,iMx]=max(sqrt(rcAll.*prcAll));datOut.prbsF1(:,1)=datOut.prbsF1(:,1)+log(mxF1')/kFld;%mxF1
datOut.prbs_ARI=datOut.prbs_ARI+adjusted_rand_index(pY(slctKfld),prbsY)/kFld;%adjusted_rand_index(pY(slctKfld),prbsY)
datOut.prbs_NMI=datOut.prbs_NMI+normalized_mutual_info(pY(slctKfld),prbsY,'avg')/kFld;%normalized_mutual_info(pY(slctKfld),prbsY,'avg'),pause
end
datOut.prbsAcc=exp(datOut.prbsAcc);datOut.prbsF1=exp(datOut.prbsF1);
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

