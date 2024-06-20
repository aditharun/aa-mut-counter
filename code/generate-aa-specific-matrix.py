import os
import pandas as pd 
import numpy as np
from multiprocess import Pool, cpu_count
import argparse
import sys

#read in from command line
parser = argparse.ArgumentParser()

parser.add_argument('--aa', type=str, help='3 letter amino acid code with first letter capitalized')

args = parser.parse_args()

amino_acid = args.aa

cwd0=os.path.dirname(os.getcwd())
cwd_results = os.path.join(cwd0, 'results')
cwd_refdata=os.path.join(cwd0,'data/reference_data_files/')

outpath = os.path.join(cwd0, "results", "counts-" + amino_acid + ".txt")

mutfiles_path = os.path.join(cwd_results, "all-muts", "df_mutfiles.csv")

if not os.path.isfile(mutfiles_path):
    print(f"File does not exist: {mutfiles_path}")
    sys.exit()

setlen=lambda x:len(set(x)) # Calculate length of set of a list.

df_mutfiles=pd.read_csv(mutfiles_path,sep='\t',dtype=str)

df_mutfiles_aa=df_mutfiles[df_mutfiles["HGVSp"].str.endswith(amino_acid) & df_mutfiles["HGVSp"].notna()].reset_index()

geneset=list(set(df_mutfiles.Hugo_Symbol.values))
histcodes=list(set(df_mutfiles.CODE.sort_values()))
geneset_aa=list(set(df_mutfiles_aa.Hugo_Symbol.values))

geneset_specific_aa=list(set(df_mutfiles_aa["Hugo_Symbol"]+df_mutfiles_aa["HGVSp"]))

geneset_aa_list=[]
for i in geneset_specific_aa:
    geneset_aa_list.append(i)
    
df_Out_aa=pd.DataFrame(0,columns=['All']+histcodes,index=geneset_aa_list+['Total'],dtype=int)


def fun_Outdf(inpvec):
    df_Out1=inpvec[1]
    inpdf=inpvec[0]
    for idx,row in inpdf.iterrows():
        igene=str(row.Hugo_Symbol)
        imutsite=str(row.HGVSp)
        icode=row.CODE
        """
        print("Values:")
        print("igene: "+str(igene))
        print("imut: "+str(imutsite))
        print("icode: "+str(icode))
        print([igene+'_'+imutsite,icode])
        print(df_Out1)
        """
        df_Out1.loc[igene+imutsite,icode]=df_Out1.loc[igene+imutsite,icode]+1
        df_Out1.loc[igene+imutsite,'All']=df_Out1.loc[igene+imutsite,'All']+1
    return df_Out1


df_Out_aa=fun_Outdf([df_mutfiles_aa, df_Out_aa])


#Check correct counts: All these numbers are same if all mutations were counted once
#in order to just get cystine, carry out all steps up until this one, but create two copies of df out, one with just cystine and one normal
#then we continue with  everything, performing operations on the cystine filtered result, but once done we take the non cystine
#filter out indices not in cystine result, and take total and all values from normal to maintain counts. 

df_Out_aa.drop(columns='All').sum().sum(),df_Out_aa['All'].sum(),len(df_mutfiles_aa)

# %%
# This is an important point. We could choose to count the total number of samples however, then if the same sample showed up in two different studies, it will be counted as one. Instead, we consider the possibility that since some studies name samples simply by numbers, it is possible that two studies have similar name but different samples. So we count unique Study_ID+Sample_ID.
# This is consistent since we already removed redundant Patient IDs with the same histology in a previous step.

df_mutfiles_aa['StidSid']=[df_mutfiles_aa.loc[idx,'STUDY_ID']+'_'+df_mutfiles_aa.loc[idx,'Tumor_Sample_Barcode'] for idx in df_mutfiles_aa.index]
# Count number of cases 'Total' within each histology : Perform this action after filtering for curated and sequenced samples
TotalRow=[setlen(df_mutfiles_aa.StidSid.values)]+[setlen(df_mutfiles_aa.StidSid[df_mutfiles_aa.CODE==ihist].values) for ihist in df_Out_aa.columns[1:]]

df_Out_aa.loc['Total']=TotalRow

for item in ['n/a','NAN','NA','na','nan']:
    if item in df_Out_aa.index:
        df_Out_aa=df_Out_aa.drop(index=item)

df_Out_aa=df_Out_aa.drop(columns=df_Out_aa.columns[df_Out_aa.loc['Total']==0])


# To test genomics pipeline pre-gene renaming, the interim file df_out can be output by un-commenting this command.
#df_Out.to_excel(os.path.join(cwd_interim,'df_Out_preGeneRename.xlsx'),index_label='Hugo_Symbol')
#df_Out_cystine.to_csv(os.path.join(cwd_interim,'df_Out_cystine_preGeneRename.csv'),index_label='Hugo_Symbol')

df_chDegen=pd.read_csv(os.path.join(cwd_refdata,'Genelist_ManyChromosomes.xlsx'))
# Only MARCH1 MARCH2 and SEPT15 have any real issues. This is very very likely due to excel errors someone made in the past by copying data incorrectly without realizing.
degenchlist=df_chDegen[['Hugo_Symbol','Chromosome Locations','Entries_per_Location']].applymap(lambda x:str(x).upper()).values


# Identifying chromosomal degeneracies takes half an hour using 6 cores at >2Ghz. If not updating raw data, uncomment this step to skip next block.
def fun_chromosome(inplist):
    setlen=lambda x:len(set(x))
    df_mutfiles_aa=inplist[0]
    geneset=inplist[1]
    return [igene for igene in geneset if setlen(df_mutfiles_aa[df_mutfiles_aa.Hugo_Symbol==igene].Chromosome)>1]

ncores=cpu_count()-2# number of cores the task can be split into
imarkers=[i*int(len(geneset)/ncores) for i in range(ncores+1)]
imarkers[-1]=len(geneset)
gslist=[[df_mutfiles_aa,geneset[imarkers[i]:imarkers[i+1]]] for i in range(ncores)]
if __name__ == '__main__': #necessary for multiprocessing docs
    po=Pool(ncores) # invoke 6 pooled threads/processes. 
    list_degen=list(po.map(fun_chromosome,gslist)) 
    po.close() 
    po.join()

degenchlist=[elem for row in list_degen for elem in row]
degenchlist=[[igene,set(df_mutfiles_aa[df_mutfiles_aa.Hugo_Symbol==igene].Chromosome)] for igene in degenchlist]
degenchlist=[[igene,ichlist,[sum(df_mutfiles_aa[df_mutfiles_aa.Hugo_Symbol==igene].Chromosome==ich) for ich in ichlist]] for igene,ichlist in degenchlist]

# Check for chromosome assignment issues in the genomic dataset. > THey are present but minimal.
df_chDegen=pd.DataFrame(degenchlist, columns=['Hugo_Symbol','Chromosome Locations','Entries_per_Location'])

#df_chDegen.to_csv(os.path.join(cwd_interim,'Genelist_ManyChromosomes.xlsx'))# Only MARCH1 MARCH2 and SEPT15 have any real issues. This is very very likely due to excel errors someone made in the past by copying data incorrectly without realizing.


#Import alternate gene nomenclature file from cbioportal: https://docs.cbioportal.org/3.-cbioportal-maintenance/updating-gene-and-gene_alias-tables
#Homo_sapien.gene_info.gz ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

dfGeneNames=pd.read_csv(os.path.join(cwd_refdata,'Homo_sapiens_gene_info_GM.txt'),sep='\t',dtype=str)
dfGeneNames=dfGeneNames.astype(str).applymap(lambda x:x.upper())
dfGeneNames.Synonyms=[str(row).split('|') for row in dfGeneNames.Synonyms.values]
listSynonyms=[elem for row in dfGeneNames.Synonyms.values for elem in row]
#print(dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0])
def FindGeneName(igene):
    retgene=np.nan
    if (igene in set(dfGeneNames.Symbol)) or (igene not in set(listSynonyms)) or (igene in [row[0] for row in degenchlist]):
        return igene
    else:
        retgene=dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0]
    # proceed to rename if the chromosome no is same.    
    chno_retgene=str(dfGeneNames[dfGeneNames.Symbol==retgene].chromosome.values[0])
    chno_igene=str(df_mutfiles_aa[df_mutfiles_aa.Hugo_Symbol==igene].Chromosome.values[0])# to cover simple renaming situations
    return retgene if ((chno_retgene==chno_igene) and (igene in geneset)) else igene

dfC=df_Out_aa[:].copy(deep=True)
genesrenamed=[]
genesadded=[]
#dfGeneNames = dfGeneNames[~dfGeneNames['Symbol'].isin(cystineGeneList)]
#retgene=dfGeneNames[[igene in row for row in dfGeneNames.Synonyms]].Symbol.values[0]
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


dfCadd=df_Out_aa.loc[genesadded]
dfCadd['Hugo_Symbol_parent']=[FindGeneName(igene) for igene in dfCadd.index]
dfCadd=dfCadd[[dfCadd.columns[-1]]+list(dfCadd.columns[:-1])]
#dfCadd.to_excel(os.path.join(cwd_interim,'Genes_Added.xlsx'))

dfCren=df_Out_aa.loc[genesrenamed]
dfCren['Hugo_Symbol_new']=[FindGeneName(igene) for igene in dfCren.index]
dfCren=dfCren[[dfCren.columns[-1]]+list(dfCren.columns[:-1])]
#dfCren.to_excel(os.path.join(cwd_interim,'Genes_Renamed.xlsx'))

print('len(genesrenamed),len(genesadded): ',len(genesrenamed),len(genesadded))

idx1=list(dfC.index)
idx1=sorted(idx1)
idx1.remove('Total')
idx1=idx1+['Total']
dfC=dfC.loc[idx1]

colist=dfC.columns.sort_values()
colist=[colist[-1]]+list(colist[:-1])
dfC=dfC[colist]

dfC.to_csv(os.path.join(cwd0, outpath), header=True, sep='\t', index_label='Hugo_Symbol')

print('Done. Results stored in results/ folder for ' + amino_acid)



