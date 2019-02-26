%import cell array --> ptCellArray
%having trouble saving the snp data so might not work--> running on indel

deleterow= false(size(ptIndels(1,:)));

I= []
J= []
K= []

for n= 1:size(ptIndels(:,1))
    if ptIndels{n, 10}==1 & ptIndels{n, 17}==1
        %ptIndels{n, 'clinvar_status'}==1 & ptIndels{n, 'HGMD_status'}==1
        vec1= ptIndels(n,:)
        i(n,:)= vec1
    elseif ptIndels{n, 10}==1 & ptIndels{n, 17}==0
        vec2= ptIndels(n,:)
        j(n,:)= vec2
    elseif ptIndels{n, 10}==0 & ptIndels{n, 17}==1
        vec3= ptIndels(n,:)
        k(n,:)= vec3
    end
    
end

I=i(~cellfun(@isempty, i(:,1)), :);
J=j(~cellfun(@isempty, j(:,1)), :);
K=k(~cellfun(@isempty, k(:,1)), :);


T1 = cell2table(I,'VariableNames',{'chr','pos','ref','alt','chr_hg19','pos_hg19','ref_hg19','alt_hg19','ref_hg19ref_hg38','clinvar_status','clinvar_rs','clinvar_clnsig','clinvar_trait','clinvar_review','clinvar_hgvs','clinvar_var_source','HGMD_status','HGMD_ACC_NUM','HGMD_HGVS_cdna','HGMD_disease','HGMD_pmid','HGMD_Variant_class'});
writetable(T1,'filteredIndels_cm.csv');

T2 = cell2table(J,'VariableNames',{'chr','pos','ref','alt','chr_hg19','pos_hg19','ref_hg19','alt_hg19','ref_hg19ref_hg38','clinvar_status','clinvar_rs','clinvar_clnsig','clinvar_trait','clinvar_review','clinvar_hgvs','clinvar_var_source','HGMD_status','HGMD_ACC_NUM','HGMD_HGVS_cdna','HGMD_disease','HGMD_pmid','HGMD_Variant_class'});
writetable(T2,'filteredIndels_c.csv');

T3 = cell2table(K,'VariableNames',{'chr','pos','ref','alt','chr_hg19','pos_hg19','ref_hg19','alt_hg19','ref_hg19ref_hg38','clinvar_status','clinvar_rs','clinvar_clnsig','clinvar_trait','clinvar_review','clinvar_hgvs','clinvar_var_source','HGMD_status','HGMD_ACC_NUM','HGMD_HGVS_cdna','HGMD_disease','HGMD_pmid','HGMD_Variant_class'});
writetable(T3,'filteredIndels_m.csv');



% %% now for the snp's
% 
% deleterow= false(size(PtSNP(1,:)));
% 
% S= []
% 
% for m= 1:size(PtSNP(:,1))
%     if PtSNP{m, 9}==1 && PtSNP{m, 16}==1
%         %ptIndels{n, 'clinvar_status'}==1 & ptIndels{n, 'HGMD_status'}==1
%         vecc= PtSNP(m,:)
%         s(m,:)= vecc
%     end
%     
% end
% 
% S=s(~cellfun(@isempty, s(:,1)), :);
% 
% Ts = cell2table(S,'VariableNames',{});
% writetable(Ts,'filteredSNPs.csv');
% 









