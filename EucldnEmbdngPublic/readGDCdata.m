function [pX,genes,tOut,Diagnwseis]=readGDCdata(drctRd,varRd)
osTos=computer('arch');if strcmp(osTos(1:min(3,end)),'win'),flSp='\';else,flSp='/';end
%-----------------data files----------------
drctData=dir(drctRd);drctData=drctData([drctData(:).isdir]);drctData=drctData(3:end);
nFldrs=length(drctData);fldrs=string({drctData(:).name});iSlctd=true(nFldrs,1);flsTSV=fldrs;
iPlc=1;
for ii=1:nFldrs,arxeio=dir([drctRd,drctData(ii).name,flSp,'*.tsv']);
    if isempty(arxeio),iSlctd(iPlc)=false;iPlc=iPlc+1;
    else,nEntrs=length(arxeio);iSlctd(iPlc:iPlc+nEntrs-1)=true;
        flsTSV(iPlc:iPlc+nEntrs-1)=string({arxeio(:).name});
        fldrs(iPlc:iPlc+nEntrs-1)=string(drctData(ii).name);
        iPlc=iPlc+nEntrs;
    end
end
fldrs=fldrs(iSlctd);flsTSV=flsTSV(iSlctd);
%-----------------cases----------------
mtdata=dir([drctRd,'*.json']);txtJSON=fileread([drctRd,flSp,mtdata(1).name]);
txtJSON=jsondecode(txtJSON);
fls_cases=strings(size(flsTSV));iSlctd=false(length(flsTSV),1);
for ii=1:length(txtJSON),a=txtJSON{ii};
    if isfield(a,'associated_entities')&&~isempty(a.associated_entities)&&...
            isfield(a.associated_entities,'case_id')&&~isempty(a.associated_entities.case_id)
        iFnd=[];
        if isfield(a,'file_name')&&~isempty(a.file_name)
            iFnd=find(flsTSV==string(a.file_name),1);
        end
        if isempty(iFnd)&&isfield(a,'file_id')&&~isempty(a.file_id)
            iFnd=find(fldrs==string(a.file_id));
        end
        if ~isempty(iFnd)
            fls_cases(iFnd)=string(a.associated_entities.case_id);
            iSlctd(iFnd)=true;
        end
    end
end,fls_cases=fls_cases(iSlctd);fldrs=fldrs(iSlctd);flsTSV=flsTSV(iSlctd);
%-----------------clinical----------------
clnclFl=[drctRd,'clinical.tsv'];
opts=detectImportOptions(clnclFl,'FileType','text','Delimiter','\t');
opts=setvartype(opts,repmat("string",1,length(opts.VariableNames)));
clnclFl=readtable(clnclFl,opts);
VrsRd={'cases_case_id','project_project_id','diagnoses_primary_diagnosis',...
    'diagnoses_ajcc_pathologic_stage','diagnoses_ajcc_clinical_stage',...
    'demographic_gender','cases_primary_site','diagnoses_diagnosis_is_primary_disease'};
clnclFl=clnclFl(:,VrsRd);clnclFl=unique(clnclFl,'rows');
isPrmr0=(clnclFl(:,'diagnoses_diagnosis_is_primary_disease').Variables=="false");
isPrmr1=(clnclFl(:,'diagnoses_diagnosis_is_primary_disease').Variables=="true");
[Diagnwseis,~,dLbls]=unique(clnclFl(:,'diagnoses_primary_diagnosis').Variables);
diagns=false(length(flsTSV),length(Diagnwseis));stgLvl=strings(size(diagns));
sitePrm=strings(size(diagns));iSlctd=true(length(flsTSV),1);gender=strings(length(flsTSV),1);
for ii=1:length(flsTSV)
    iFnd=find(clnclFl(:,'cases_case_id').Variables==fls_cases(ii));
    if ~isempty(iFnd),iFnd1=find(isPrmr1(iFnd));
        if isempty(iFnd1),iFnd1=find(isPrmr0(iFnd));end
        if isempty(iFnd1),iSlctd(ii)=false;
        else
            diagns(ii,dLbls(iFnd(iFnd1)))=true;
            stgLvl(ii,dLbls(iFnd(iFnd1)))=clnclFl(iFnd(iFnd1),'diagnoses_ajcc_pathologic_stage').Variables;
            gender(ii)=clnclFl(iFnd(iFnd1),'demographic_gender').Variables;
        end
    end
end
fls_cases=fls_cases(iSlctd);fldrs=fldrs(iSlctd);flsTSV=flsTSV(iSlctd);
diagns=diagns(iSlctd,:);stgLvl=stgLvl(iSlctd,:);gender=gender(iSlctd);
% size(fls_cases),size(fldrs),size(flsTSV)
tOut=table(fls_cases',fldrs',flsTSV',gender,diagns,stgLvl);
genes="";
for ii=1:length(flsTSV),ii
    tbl=readtable([drctRd,char(fldrs(ii)),flSp,char(flsTSV(ii))],...
        'FileType','text','Delimiter','\t','CommentStyle','#');
    tbl=tbl(5:end,:);
    % genes=union(genes,unique(tbl(:,'gene_id').Variables));
    if ii==1,genes=tbl(:,'gene_id').Variables;slct=startsWith(genes, "ENSG");
        genes=[genes(slct),tbl(slct,'gene_name').Variables];
        pX=nan(length(flsTSV),size(genes,1));pX(ii,:)=tbl(slct,varRd).Variables;
    else,genes=tbl(:,'gene_id').Variables;
        [~,ind0,ind1]=intersect(genes(:,1),tbl(:,'gene_id').Variables, 'stable');
        pX(ii,ind0)=(tbl(ind1,varRd).Variables);pX=pX(:,ind0);genes=genes(ind0,:);
    end
end