# A pipeline to quickly prioritize candidate disease-causing WES variants 

# Remove dict items
def entries_to_remove(entries, the_dict):
    for key in entries:
        if key in the_dict:
            del the_dict[key]
# Parents' variants
parents_var = {}
def parent_var(file):
    for line in file:
        line = line.strip('\n')
        line = line.split('\t')
        if line[0]!='chr':
            key = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]
            parents_var[key] = 0
with open('anno_snp1','r') as f6:
    p_snp_file = f6.readlines()   
with open('anno_indel1','r') as f7: 
    p_indel_file = f7.readlines()  
with open('anno_snp2','r') as f8:
    m_snp_file = f8.readlines()   
with open('anno_indel2','r') as f9: 
    m_indel_file = f9.readlines()
with open('anno_dbnsfp1','r',encoding='latin-1') as f10:
    p_dbnsfp_file = f10.readlines()   
with open('anno_dbnsfp2','r',encoding='latin-1') as f11: 
    m_dbnsfp_file = f11.readlines()  
parent_var(p_snp_file)
parent_var(p_indel_file)
parent_var(m_snp_file)
parent_var(m_indel_file)
parent_var(p_dbnsfp_file)
parent_var(m_dbnsfp_file)
f6.close()
f7.close()
f8.close()
f9.close()
f10.close()
f11.close()
# Patient's variants
with open('anno_snp3','r') as f0:
    snp_file = f0.readlines()
with open('anno_indel3','r') as f3:     #WGSA indel file
    indel_file = f3.readlines()
clinvar = {}
hgmd = {}
indel_gene = {}
anno_filter = ['stopgain','stoploss','splicing','nonsynonymous','startgain','startloss','frameshift']
snpeff_filter = ['missense','splice','stop_gain','stop_lost','start_gain','start_lost','frameshift_variant']
vep_filter  = snpeff_filter
anno_conseq = {}
snpeff_conseq = {}
vep_conseq = {}
def indel(file):
    for line in file:
        line = line.strip('\n')
        line = line.split('\t')
        key = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]
        if line[129]!='.':
            clinvar_value = line[129:135]
            clinvar[key] = clinvar_value
        if line[135]!='.':
            hgmd_value = line[135:140] 
            hgmd[key] = hgmd_value
        indel_gene[key]=[line[11],line[39],line[85]]
        anno_result = line[9].split('|')
        for i in anno_result:
            if i in anno_filter:
                anno_test = 1
            else:
                anno_test = 0
        snpeff_result = line[33].split('|')
        for i in snpeff_result:
            if i in snpeff_filter:
                snpeff_test = 1
            else:
                snpeff_test = 0
        vep_result = line[83].split('|')
        for i in vep_result:
            if i in vep_filter:
                vep_test = 1
            else:
                vep_test = 0
        anno_conseq[key] = anno_test
        snpeff_conseq[key] = snpeff_test
        vep_conseq[key] = vep_test
def snp(file):
    for line in file:
        line = line.strip('\n')
        line = line.split('\t')
        key = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]
        if line[13]!='.':
            clinvar_value = line[13:19] 
            clinvar[key] = clinvar_value
        if line[19]!='.':
            hgmd_value = line[19:25] 
            hgmd[key] = hgmd_value
        
        
snp(snp_file)
f0.close()
indel(indel_file)
f3.close()
# User input gene list
r_gene = {}
exp_gene = {}
with open('related_genes','r') as f4: # genelist from literature search/prior knowledge
    r_gene_file = f4.readlines()
with open('final_genes','r') as f5: # genelist from literature search/prior knowledge
    exp_gene_file = f5.readlines() 
for line in r_gene_file:
    line = line.strip('\n')
    line = line.split('\t')
    r_gene[line[1]] = 0
for line in exp_gene_file:
    line = line.strip('\n')
    line = line.split('\t')
    exp_gene[line[1]] = 0
f4.close()
f5.close()
# User input variant file
with open('candidate_var','r') as f1:    #input list of patient's variants
    variant_file = f1.readlines()
var_list = {}
for line in variant_file:
    line = line.strip('\n')
    line = line.split('\t')
    key = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]
    var_list[key] = 0
f1.close()
#other files
phenotype = {}
with open('pheno','r') as pheno:
    pheno = pheno.readlines()
for line in pheno:
    line = line.strip('\n')
    line = line.split('\t')
    line[0]= line[0].lower
    phenotype[line[0]] = 0
    
with open('anno_dbnsfp3','r',encoding = 'latin-1') as f2:    #dbnsfp out
    dbnsfp_file = f2.readlines()
    
tier1 = open('tier1','w+')    # report variants
tier2 = open('tier2','w+')
tier3 = open('tier3','w+')
tier4 = open('tier4','w+')
dbnsfp_gene = {}
cur_clinvar = {}
cur_hgmd = {}
no_var_list = {}
yes_var_list = {}
dbnsfp_damaging = {}
for line in dbnsfp_file:
    line = line.strip('\n')
    line = line.split('\t')
    key = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]
    dbnsfp_gene[key] = line[12]
    if line[4] not in ['X','.'] and line[5] not in ['X','.'] and line[71]=='D' and line[64]=='D':
        dbnsfp_damaging[key] = 0
# start pipeline tier 1 & 2    
for key,value in var_list.items():
    if key in clinvar:
        cur_clinvar[key] = clinvar[key] 
    if key in hgmd:
        cur_hgmd[key] = hgmd[key] 
    if key not in clinvar and key not in hgmd:
        no_var_list[key] = 0
for key,value in cur_clinvar.items():
    if value[2].lower() in phenotype:
        yes_var_list[key] = 0
    else:
        no_var_list[key] = 0
for key,value in cur_hgmd.items():
    if value[3].lower() in phenotype:
        yes_var_list[key] = 0
    else:
        no_var_list[key] = 0
        
if len(yes_var_list)>0:
    for key,value in yes_var_list.items():
        if key in parents_var:
            tier1.write(key+'\n')
        else:
            tier2.write(key+'\n')
tier1.close()
tier2.close()
# damaging

def damaging(position_identifier):
    key = position_identifier
    if key in dbnsfp_damaging:
        return True
    if (key in anno_conseq and anno_conseq[key]==1) or (key in snpeff_conseq and snpeff_conseq[key]==1) or (key in vep_conseq and vep_conseq[key]==1):
        return True
# tier 3

for key,value in no_var_list.items():
    if key in dbnsfp_gene and dbnsfp_gene[key] in r_gene:
        if damaging(key):
            tier3.write(key+'\n')
            entries_to_remove(key,no_var_list)
    elif key in indel_gene:
        anno_genes = indel_gene[key][0].split('|')
        snpeff_genes = indel_gene[key][1].split('|')
        vep_genes = indel_gene[key][2].split('|')
        for genes in r_gene:
            if (genes in anno_genes or genes in snpeff_genes or genes in vep_genes):
                if damaging(key):
                    tier3.write(key+'\n')
                    entries_to_remove(key,no_var_list)
tier3.close()
# tier 4
for key,value in no_var_list.items():
    if key in dbnsfp_gene and dbnsfp_gene[key] in exp_gene:
        vvv += 1
        if damaging(key):
            #tier4.write(key+'\n')
            tier4.write(dbnsfp_gene[key]+'\n')
    if key in indel_gene and damaging(key):
        anno_genes = indel_gene[key][0].split('|')
        snpeff_genes = indel_gene[key][1].split('|')
        vep_genes = indel_gene[key][2].split('|')
        for genes in exp_gene:
            if genes in anno_genes or genes in snpeff_genes or genes in vep_genes:
                #tier4.write(key+'\n')
                tier4.write(genes+'\n')
                break
tier4.close()
