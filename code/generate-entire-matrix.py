import os
import pandas as pd 
import numpy as np
from multiprocess import Pool, cpu_count

cwd0=os.path.dirname(os.getcwd())

cwd_results = os.path.join(cwd0, 'results/all-muts')
cwd_refdata=os.path.join(cwd0,'data/reference_data_files/')

if not os.path.exists(cwd_results):
    os.makedirs(cwd_results)

setlen=lambda x:len(set(x)) 


# Find Studies
cwd=os.path.join(cwd0,'data/cbioportal_raw_EXOME139') 

list_study_dirs=['all_phase2_target_2018_pub','acbc_mskcc_2015','acc_tcga_pan_can_atlas_2018','acyc_mda_2015','acyc_mskcc_2013','acyc_sanger_2013','all_stjude_2013','all_stjude_2015','all_stjude_2016','aml_target_2018_pub','angs_project_painter_2018','bfn_duke_nus_2015','blca_bgi','blca_cornell_2016','blca_dfarber_mskcc_2014','blca_tcga_pub_2017','brca_bccrc','brca_broad','brca_igr_2015','brca_mbcproject_wagle_2017','brca_sanger','brca_tcga_pan_can_atlas_2018','ccrcc_irc_2014','ccrcc_utokyo_2013','cesc_tcga_pan_can_atlas_2018','chol_jhu_2013','chol_nccs_2013','chol_tcga_pan_can_atlas_2018','cll_iuopa_2015','cllsll_icgc_2011','coadread_dfci_2016','coadread_genentech','coadread_tcga_pan_can_atlas_2018','ctcl_columbia_2015','desm_broad_2015','dlbc_broad_2012','dlbc_tcga_pan_can_atlas_2018','dlbcl_dfci_2018','dlbcl_duke_2017','egc_tmucih_2015','es_dfarber_broad_2014','es_iocurie_2014','esca_tcga_pan_can_atlas_2018','escc_icgc','escc_ucla_2014','gbc_shanghai_2014','gbm_tcga_pan_can_atlas_2018','hcc_inserm_fr_2015','hnsc_broad','hnsc_jhu','hnsc_mdanderson_2013','hnsc_tcga_pan_can_atlas_2018','kich_tcga_pan_can_atlas_2018','kirc_bgi','kirc_tcga_pan_can_atlas_2018','kirp_tcga_pan_can_atlas_2018','laml_tcga_pan_can_atlas_2018','lcll_broad_2013','lgg_ucsf_2014','lgggbm_tcga_pub','lihc_amc_prv','lihc_riken','lihc_tcga_pan_can_atlas_2018','luad_broad','luad_mskcc_2015','luad_tcga_pan_can_atlas_2018','lusc_tcga_pan_can_atlas_2018','mbl_broad_2012','mbl_icgc','mbl_pcgp','mbl_sickkids_2016','mcl_idibips_2013','mds_tokyo_2011','mel_tsam_liang_2017','meso_tcga_pan_can_atlas_2018','mm_broad','mpnst_mskcc','mrt_bcgsc_2016','nbl_amc_2012','nbl_broad_2013','nbl_target_2018_pub','nbl_ucologne_2015','nccrcc_genentech_2014','nepc_wcm_2016','nhl_bcgsc_2011','nhl_bcgsc_2013','npc_nusingapore','nsclc_tcga_broad_2016','ov_tcga_pan_can_atlas_2018','paac_jhu_2014','paad_icgc','paad_qcmg_uq_2016','paad_tcga_pan_can_atlas_2018','paad_utsw_2015','panet_arcnet_2017','panet_jhu_2011','panet_shanghai_2013','past_dkfz_heidelberg_2013','pcnsl_mayo_2015','pcpg_tcga_pan_can_atlas_2018','plmeso_nyu_2015','prad_broad','prad_broad_2013','prad_cpcg_2017','prad_eururol_2017','prad_fhcrc','prad_mich','prad_p1000','prad_su2c_2015','prad_tcga_pan_can_atlas_2018','rms_nih_2014','rt_target_2018_pub','sarc_tcga_pan_can_atlas_2018','sclc_cancercell_gardner_2017','sclc_clcgp','sclc_jhu','sclc_ucologne_2015','skcm_broad','skcm_broad_brafresist_2012','skcm_broad_dfarber','skcm_tcga_pan_can_atlas_2018','skcm_ucla_2016','skcm_yale','stad_pfizer_uhongkong','stad_tcga_pan_can_atlas_2018','stad_uhongkong','stad_utokyo','stes_tcga_pub','tet_nci_2014','tgct_tcga_pan_can_atlas_2018','thca_tcga_pan_can_atlas_2018','thym_tcga_pan_can_atlas_2018','uccc_nih_2017','ucec_tcga_pan_can_atlas_2018','ucs_jhu_2014','ucs_tcga_pan_can_atlas_2018','um_qimr_2016','uvm_tcga_pan_can_atlas_2018','vsc_cuk_2018','wt_target_2018_pub']

print('Number of Studies included:',len(list_study_dirs),'\n The following studies are not present in cbioportal_raw_EXOME139/ folder. Please download by running script provided in the folder: ',[istudy for istudy in list_study_dirs if not os.path.isdir(os.path.join(cwd,istudy))]) # Second value should be null vector to ensure all studies are downloaded and unzipped in current directory.

# Import Clinical Data
# when clinical folder defined in cwd0 directory
cldirlist0=[os.path.join(cwd0,'data/clinical',idir) for idir in os.listdir(os.path.join(cwd0,"data/clinical")) if os.path.isdir(os.path.join(cwd0,'data/clinical',idir))]
cldirlist1=[os.path.join(idir,idir1) for idir in cldirlist0 for idir1 in os.listdir(idir)]

# limit to included studies
clinical_file_list=[os.path.join(idir,'data_clinical_sample_v2.txt') for idir in cldirlist1 if any([row in idir for row in list_study_dirs])] 

N_LongiSamples={}
list_df_clinical=[]
for istudypath in clinical_file_list:
    Study_ID=istudypath.split('/')[-2]
    df_clinical_istudy=pd.read_csv(istudypath,sep='\t',comment='#',dtype=str)
    # Remove longitudinal samples
    N_LongiSamples[Study_ID]=(len(df_clinical_istudy.PATIENT_ID)-setlen(df_clinical_istudy.PATIENT_ID)) # Record how many longitudinal samples are removed within every study
    df_clinical_istudy=df_clinical_istudy.drop_duplicates(subset='PATIENT_ID',keep='first')
    df_clinical_istudy['STUDY_ID']=Study_ID
    keepcols=['STUDY_ID','PATIENT_ID','SAMPLE_ID','CODE'] # These columns are ALWAYS present in the curated clinical file.
    list_df_clinical=list_df_clinical+[df_clinical_istudy[keepcols]]
    del df_clinical_istudy
df_clinical=pd.concat(list_df_clinical).reset_index(drop=True) # A unified clinical matrix with ALL the clinical data.
print('Total Number of Studies:',len(list_study_dirs),'\nClinical Data Avaliable for Samples(=Patients or cases)):',len(df_clinical),'\nNumber of Longitudinal Samples Dropped: ',sum(N_LongiSamples.values()))


# Note list of histological codes
histcodes=list(set(df_clinical.CODE.sort_values()))


# Import cbioportal studies and sort reverse chronologically
df_cbiostudies=pd.read_excel(os.path.join(cwd_refdata,'Table2_v6.xlsx'),dtype=str)
df_cbiostudies['Year']=df_cbiostudies['Year'].astype(int)

df_clinical['Year']=[df_cbiostudies.Year[df_cbiostudies.Study_ID==stid].values[0] for stid in df_clinical.STUDY_ID.values]


# Filter repeated samples acrross studies within a given histological code. Note that this filter assumes that dataframes preseve indices upon copying and slicing. If that functinality changes, the codes need to change.
list_redunPID=[]
df_clinical_NR=df_clinical[:].copy() # initialize non redundant clinical samples
for icode in histcodes:
    # Find all patient IDs within a given hist code
    df_patient=df_clinical[df_clinical.CODE==icode]
    # Find repeating patient IDs
    for pid in set(df_patient.PATIENT_ID.values):
        if sum(df_patient.PATIENT_ID==pid)>1:
            df_tmp_remove=df_patient[df_patient.PATIENT_ID==pid].sort_values('Year',ascending=False)
            for idx_drop in df_tmp_remove.index[1:]:
                df_clinical_NR.drop(idx_drop,inplace=True)
            del df_tmp_remove
    df_patient=df_patient[[sum(df_patient.PATIENT_ID==pid)>1 for pid in df_patient.PATIENT_ID.values]].sort_values('SAMPLE_ID')
    list_redunPID=list_redunPID+[df_patient]    
    del df_patient

# %%
print('Clinical Data Processed. Number of Redundant Samples:',sum([len(dfi) for dfi in list_redunPID]))

# %%
df_redundants=pd.concat(list_redunPID).reset_index(drop=True)

# %%
print('Set of ROSETTA codes with redundant samples and corresponding studies:\n',[[icode]+list(set(df_redundants.STUDY_ID[df_redundants.CODE==icode])) for icode in histcodes if len(set(df_redundants.STUDY_ID[df_redundants.CODE==icode]))>0])

# %%
print("Number of Patient IDs which have been counted more than once:", len([pid for pid in set(df_redundants.PATIENT_ID) if sum(df_clinical.PATIENT_ID==pid)>1]))
print("Validation: Number of redundant samples in post-filtration :",len([pid for pid in set(df_redundants.PATIENT_ID) if sum(df_clinical_NR.PATIENT_ID==pid)>1]))

# %%
df_redundants.to_excel(os.path.join(cwd_results,'Redundant_Patient_ID.xlsx'))


# List All genes
list_mutfile=[]
iter1=0
for istudyRaw in list_study_dirs:
    os.chdir(os.path.join(cwd,istudyRaw))
    # import the mutations extended file
    if os.path.exists(os.path.join(cwd,istudyRaw,'data_mutations_extended.txt')):
        df_mut_istudy=pd.read_csv('data_mutations_extended.txt',sep='\t',dtype=str,comment='#',encoding='cp1252',quoting=3,keep_default_na=False)
    else:
        df_mut_istudy=pd.read_csv('data_mutations.txt',sep='\t',dtype=str,comment='#',encoding='cp1252',quoting=3,keep_default_na=False)
    keepcols=['Hugo_Symbol','Chromosome','Variant_Classification','Tumor_Sample_Barcode', 'HGVSp'] # Only these columns are needed
    df_mut_istudy=df_mut_istudy[keepcols]
    df_mut_istudy.dropna(subset=['Hugo_Symbol'],inplace=True,how='any')

    var_list=['translation_start_site','nonsense_mutation', 'frame_shift_del', 'frame_shift', 'frame_shift_ins', 'missense_mutation', 'missense', 'nonsense', 'in_frame_del',  'in_frame_ins','nonstop_mutation']
    df_mut_istudy.Variant_Classification=df_mut_istudy.Variant_Classification.apply(lambda x:x.lower())
    df_mut_istudy=df_mut_istudy[df_mut_istudy.Variant_Classification.isin(var_list)]
    
    # limit to clinically curated samples
    df_clinical_istudy=df_clinical_NR[df_clinical_NR.STUDY_ID==istudyRaw]
    df_mut_istudy=df_mut_istudy[df_mut_istudy.Tumor_Sample_Barcode.isin(df_clinical_istudy.SAMPLE_ID)]
    
    # Include only Unique genes within each sample
    df_mut_istudy.drop_duplicates(subset=['Tumor_Sample_Barcode','Hugo_Symbol'],keep='first',inplace=True)

    # Add ICD-code information to the mutation data
    df_mut_istudy['CODE']=[df_clinical_istudy.CODE[df_clinical_istudy.SAMPLE_ID==sid].values[0] for sid in df_mut_istudy.Tumor_Sample_Barcode]
    
    # store mutations file in memory
    df_mut_istudy['STUDY_ID']=istudyRaw
    list_mutfile=list_mutfile+[df_mut_istudy]
    
    del df_mut_istudy
    del df_clinical_istudy
    
    iter1+=1
    #clear_output()
    print('Study Number ',iter1,' done:', istudyRaw)# track progress in case of errors.


df_mutfiles=pd.concat(list_mutfile).reset_index(drop=True)

df_mutfiles.Hugo_Symbol=df_mutfiles.Hugo_Symbol.astype(str).apply(lambda x:x.upper())

# List All genes
geneset=list(set(df_mutfiles.Hugo_Symbol.values))
remlist=[' ','Unknown',np.nan,'na','NA','NaN','NAN']
for igene in remlist:
    if igene in geneset:
        geneset.remove(igene)

geneset=sorted(geneset)
print(len(geneset)) # Number of genes included in our study


# remove NaN values
df_mutfiles.drop(df_mutfiles[[not(row) for row in df_mutfiles.Hugo_Symbol.isin(geneset)]].index,inplace=True)


df_mutfiles.to_csv(os.path.join(cwd_results,'df_mutfiles.csv'),sep='\t',index=False)


# Define Null Output Matrix: genes vs hist as Integer
df_Out=pd.DataFrame(0,columns=['All']+histcodes,index=geneset+['Total'],dtype=int)
# Construct the output matrix by iteratating over rows of reduced dataframe thereby covering each gene mutated for each sample
def fun_Outdf(inpvec):
    df_Out1=inpvec[1]
    inpdf=inpvec[0]
    for idx,row in inpdf.iterrows():
        igene=row.Hugo_Symbol
        icode=row.CODE
        df_Out1.loc[igene,icode]=df_Out1.loc[igene,icode]+1
        df_Out1.loc[igene,'All']=df_Out1.loc[igene,'All']+1
    return df_Out1

ncores=cpu_count()-2# number of cores the task can be split into
imarkers=[i*int(len(df_mutfiles)/ncores) for i in range(ncores+1)]
imarkers[-1]=len(df_mutfiles)
datalist=[[df_mutfiles[imarkers[i]:imarkers[i+1]],df_Out[:].copy(deep=True)] for i in range(ncores)]
if __name__ == '__main__': #necessary acc multiprocessing docs
    po=Pool(ncores) 
    list_resdf=list(po.map(fun_Outdf,datalist)) 
    po.close() 
    po.join()

for idf in list_resdf:
    df_Out=df_Out+idf


#Check correct counts: All these numbers are same if all mutations were counted once
df_Out.drop(columns='All').sum().sum(),df_Out['All'].sum(),len(df_mutfiles)


# This is an important point. We could choose to count the total number of samples however, then if the same sample showed up in two different studies, it will be counted as one. Instead, we consider the possibility that since some studies name samples simply by numbers, it is possible that two studies have similar name but different samples. So we count unique Study_ID+Sample_ID.
# This is consistent since we already removed redundant Patient IDs with the same histology in a previous step.


df_mutfiles['StidSid']=[df_mutfiles.loc[idx,'STUDY_ID']+'_'+df_mutfiles.loc[idx,'Tumor_Sample_Barcode'] for idx in df_mutfiles.index]
# Count number of cases 'Total' within each histology : Perform this action after filtering for curated and sequenced samples
TotalRow=[setlen(df_mutfiles.StidSid.values)]+[setlen(df_mutfiles.StidSid[df_mutfiles.CODE==ihist].values) for ihist in df_Out.columns[1:]]
df_Out.loc['Total']=TotalRow
print(TotalRow)


for item in ['n/a','NAN','NA','na','nan']:
    if item in df_Out.index:
        df_Out=df_Out.drop(index=item)

df_Out=df_Out.drop(columns=df_Out.columns[df_Out.loc['Total']==0])


# To test genomics pipeline pre-gene renaming, the interim file df_out can be output by un-commenting this command.
#df_Out.to_excel(os.path.join(cwd_interim,'df_Out_preGeneRename.xlsx'),index_label='Hugo_Symbol')



df_chDegen=pd.read_csv(os.path.join(cwd_refdata,'Genelist_ManyChromosomes.xlsx'))
# Only MARCH1 MARCH2 and SEPT15 have any real issues. This is very very likely due to excel errors someone made in the past by copying data incorrectly without realizing.
degenchlist=df_chDegen[['Hugo_Symbol','Chromosome Locations','Entries_per_Location']].applymap(lambda x:str(x).upper()).values


#Import alternate gene nomenclature file from cbioportal: https://docs.cbioportal.org/3.-cbioportal-maintenance/updating-gene-and-gene_alias-tables
#Homo_sapien.gene_info.gz ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
dfGeneNames=pd.read_csv(os.path.join(cwd_refdata,'Homo_sapiens_gene_info_GM.txt'),sep='\t',dtype=str)
dfGeneNames=dfGeneNames.astype(str).applymap(lambda x:x.upper())
dfGeneNames.Synonyms=[str(row).split('|') for row in dfGeneNames.Synonyms.values]
listSynonyms=[elem for row in dfGeneNames.Synonyms.values for elem in row]
def FindGeneName(igene):
    retgene=np.nan
    if (igene in set(dfGeneNames.Symbol)) or (igene not in set(listSynonyms)) or (igene in [row[0] for row in degenchlist]):
        return igene
    else:
        retgene=dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0]
    # proceed to rename if the chromosome no is same.    
    chno_retgene=str(dfGeneNames[dfGeneNames.Symbol==retgene].chromosome.values[0])
    chno_igene=str(df_mutfiles[df_mutfiles.Hugo_Symbol==igene].Chromosome.values[0])# to cover simple renaming situations
    return retgene if ((chno_retgene==chno_igene) and (igene in geneset)) else igene


dfC=df_Out[:].copy(deep=True)
genesrenamed=[]
genesadded=[]
for igene in dfC.index:
    newgene=FindGeneName(igene)
    if (newgene != igene):
        if (newgene in dfC.index.values):
            dfC.loc[newgene]=dfC.loc[newgene].copy()+dfC.loc[igene].copy()
            dfC.drop(igene,inplace=True)
            genesadded=genesadded+[igene]
        else:
            dfC.loc[newgene]=dfC.loc[igene].copy()
            dfC.drop(igene,inplace=True)
            genesrenamed=genesrenamed+[igene]

dfCadd=df_Out.loc[genesadded]
dfCadd['Hugo_Symbol_parent']=[FindGeneName(igene) for igene in dfCadd.index]
dfCadd=dfCadd[[dfCadd.columns[-1]]+list(dfCadd.columns[:-1])]
dfCadd.to_excel(os.path.join(cwd_results,'Genes_Added.xlsx'))

dfCren=df_Out.loc[genesrenamed]
dfCren['Hugo_Symbol_new']=[FindGeneName(igene) for igene in dfCren.index]
dfCren=dfCren[[dfCren.columns[-1]]+list(dfCren.columns[:-1])]
dfCren.to_excel(os.path.join(cwd_results,'Genes_Renamed.xlsx'))

print('len(genesrenamed),len(genesadded): ',len(genesrenamed),len(genesadded))


idx1=list(dfC.index)
idx1=sorted(idx1)
idx1.remove('Total')
idx1=idx1+['Total']
dfC=dfC.loc[idx1]


colist=dfC.columns.sort_values()
colist=[colist[-1]]+list(colist[:-1])
dfC=dfC[colist]


dfC.to_csv(os.path.join(cwd_results,'all-counts-matrix.txt'),header=True,sep='\t',index_label='Hugo_Symbol')

print('Done. Results stored in all-counts-matrix.txt')



