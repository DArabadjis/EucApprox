function datOut=evalEmbddng(dstMtrx,pXQhssn,pT)
% Inputs
% dstMtrx : graph adjacency matrix - sparse #points x #points
% pXQhssn : the embedding - #points x nDims
% pT : manifold generating coordinates (empty [])

% Outputs
% datOut : table of measures of local & global embedding consistency

if ~exist('pT','var'),pT=[];end
if ~isempty(pT)
    Amap=@(pX,pT)(pX-mean(pX))*lsqminnorm(pX-mean(pX),pT-mean(pT));
    pXQproc=Amap(pXQhssn,pT);[errhsssn,pXQproc]=procrustes(pT,pXQproc);
    figure,corrplot([pT,pXQproc])
else,errhsssn=nan;
end
errsLcl=cmptMtrc(pXQhssn,dstMtrx);errsGlbl=cmptMtrcGlbl(pXQhssn,dstMtrx);

errDst=mean(errsLcl.errDst);errAng=mean(errsLcl.errAng)*180/pi;
errGrm=sqrt(mean(errsLcl.errGrm(~isinf(errsLcl.errGrm)&~isnan(errsLcl.errGrm))));
cntnt=exp(mean(log(1-errsLcl.cntnt)));
rcll=exp(mean(log(errsLcl.rcll)));

strss=sqrt(mean(errsGlbl.strss));trst=exp(mean(log(1-errsGlbl.trst)));
prcsn=exp(mean(log(errsGlbl.prcsn)));
prsnCorr=exp(mean(log(abs(errsGlbl.prsnCorr))));
sprmnCorr=exp(mean(log(abs(errsGlbl.sprmnCorr))));
kndllCorr=exp(mean(log(abs(errsGlbl.kndllCorr))));

varNams=["Embeddings' affine-Procrustes error","Local distances' absolute relative error","Local angles' absolute error(degrees)",...
    "Local metric reconstruction error (local strain)","Neighborhoods' continuity",...
    "Neighborhoods' trustworthiness","Neighborhoods' recall","Neighborhoods' precision","Neighborhoods' F1",...
    "Global stress","Global distances' Pearson","Global distances' Spearman","Global distances' Kendall"];
datOut=table(errhsssn,errDst,errAng,errGrm,cntnt,trst,rcll,prcsn,sqrt(rcll*prcsn),...
    strss,prsnCorr,sprmnCorr,kndllCorr,'VariableNames',varNams);