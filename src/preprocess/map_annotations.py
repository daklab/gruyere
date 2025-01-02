import os, sys
import numpy as np
import pickle
from tqdm import tqdm
import pandas as pd
import re
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
import warnings
warnings.filterwarnings('ignore')

variant_type = sys.argv[1]  # "coding" or "noncoding"
chromosome = int(sys.argv[2])
output_dir = sys.argv[3]

VARIANTS_EXTRACT = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract/{variant_type}/ADSP.chr{chromosome}.txt"
ANNOTATION_FILE_PATH = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/WGSA_annotations/ADSP_annotated_chr{chromosome}.annotated.snp.gz"
OUTPUT_DIR = f"/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/wgsa_annotated/{variant_type}"

print(variant_type, chromosome, VARIANTS_EXTRACT, ANNOTATION_FILE_PATH)

coding_annotations = ['chr','pos','MAP20','phyloP17way_primate','phyloP30way_mammalian',
                     'phastCons30way_mammalian','phastCons17way_primate_rankscore', 'integrated_fitCons_score',
                     'H1-hESC_fitCons_score', 'GenoCanyon_score',
                     'RegulomeDB_score', 'CADD_raw', 'CADD_phred', 'DANN_score', 'fathmm-MKL_coding_score',
                     'fathmm-MKL_coding_group', 'fathmm-XF_score', 'Eigen-raw','Eigen-PC-raw', 
                     'FANTOM5_CAGE_peak_robust', 
                     'GERP_RS', 'bStatistic', 'SpliceAI_DS_AG',
                     'SpliceAI_DS_AL','SpliceAI_DS_DG','SpliceAI_DS_DL',
                     'gnomAD_genomes_POPMAX_AF', 'gnomAD_genomes_AFR_AF','gnomAD_genomes_AMR_AF', 
                     'gnomAD_genomes_NFE_AF', 'BayesDel_addAF_score', 'Eigen_coding_or_noncoding', 
                     'fathmm-XF_coding_or_noncoding', 'ref','alt', 'LIST-S2_score', 'MVP_score']


noncoding_annotations = ['chr','pos','MAP20','phyloP17way_primate','phyloP30way_mammalian',
                     'phastCons30way_mammalian','phastCons17way_primate_rankscore', 'integrated_fitCons_score',
                     'H1-hESC_fitCons_score', 'GenoCanyon_score',
                     'RegulomeDB_score', 'CADD_raw', 'CADD_phred', 'DANN_score','Eigen-raw','Eigen-PC-raw', 'EnhancerFinder_brain_enhancer',
                     'FANTOM5_CAGE_peak_robust', 'Roadmap_E074_GenoSkyline_Plus_score', 'Roadmap_E068_GenoSkyline_Plus_score',
                     'Roadmap_E069_GenoSkyline_Plus_score', 'Roadmap_E072_GenoSkyline_Plus_score', 
                     'Roadmap_E067_GenoSkyline_Plus_score', 'Roadmap_E073_GenoSkyline_Plus_score',
                     'Roadmap_E070_GenoSkyline_Plus_score', 'Roadmap_E030_GenoSkyline_Plus_score',
                     'Roadmap_E050_GenoSkyline_Plus_score', 'Roadmap_E051_GenoSkyline_Plus_score', 
                     'Roadmap_E124_GenoSkyline_Plus_score',
                     'LINSIGHT', 'funseq2_noncoding_score',
                     'fathmm-MKL_non-coding_score', 'fathmm-MKL_non-coding_group', 'fathmm-XF_score','GERP_RS', 'bStatistic', 
                     'gnomAD_genomes_POPMAX_AF', 'gnomAD_genomes_AFR_AF','gnomAD_genomes_AMR_AF', 
                     'gnomAD_genomes_NFE_AF', 'Eigen_coding_or_noncoding', 
                     'fathmm-XF_coding_or_noncoding', 'ref','alt', 'SpliceAI_DS_AG',
                     'SpliceAI_DS_AL','SpliceAI_DS_DG','SpliceAI_DS_DL']

def get_extracted_df(df, to_extract):
    anno = pd.DataFrame()
    to_extract.columns = ['chr','start','end','gene']
    for index,row in to_extract.iterrows():
        start = row['start'] 
        end = row['end'] 
        if len(df[(df['pos']>=start) & (df['pos']<=end)]) > 0:
            temp = df[(df['pos']>=start) & (df['pos']<=end)]
            temp['Gene'] = row['gene']
            anno = pd.concat([anno, temp])
    return anno

def get_nonsyn(value):
    out = re.findall(r"[-+]?(?:\d*\.*\d+)", value)
    if out == []:
        return 0
    else:
        try:
            return float(max(out))
        except:
            return 0

def process_spliceai_column(df, columns):
    for column in columns:
        df[column] = df[column].replace('.', np.NaN)
        df[column] = df[column].apply(lambda x: max(map(float, x.split(';'))) if isinstance(x, str) else x)
        df[column] = df[column].fillna(0)
    return df
        
def make_coding_numerical(df):
    # go through each column and make the dataframe numerical
    df['LIST-S2_nonsyn'] = df['LIST-S2_score'].apply(get_nonsyn)
    df['MVP_missense'] = df['MVP_score'].apply(get_nonsyn)
    df = df.replace(".", np.NaN)
    # impute gnomad scores
    df['gnomAD_genomes_POPMAX_AF'] = df['gnomAD_genomes_POPMAX_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['gnomAD_genomes_AFR_AF'] = df['gnomAD_genomes_AFR_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)#mins of other populations are 0
    df['gnomAD_genomes_AMR_AF'] = df['gnomAD_genomes_AMR_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['gnomAD_genomes_NFE_AF'] = df['gnomAD_genomes_NFE_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['Eigen_coding_or_noncoding'] = np.where(df['Eigen_coding_or_noncoding']=='c', 1, 0)
    df['fathmm-XF_coding_or_noncoding'] = np.where(df['fathmm-XF_coding_or_noncoding']=='coding', 1, 0)
    df = process_spliceai_column(df, ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL'])
    # make regulomedb quantitative
    for val in df['RegulomeDB_score'].unique():
        new_val = float(re.search(r'\d+', str(val)).group())
        df.loc[df["RegulomeDB_score"] == val, "RegulomeDB_score"] = new_val
    df.loc[df['FANTOM5_CAGE_peak_robust']=='Y','FANTOM5_CAGE_peak_robust'] = 1
    df.loc[df['FANTOM5_CAGE_peak_robust']!=1,'FANTOM5_CAGE_peak_robust'] = 0
    df['BayesDel_addAF_score'] = df['BayesDel_addAF_score'].str.split("|").str[0] # keeping first val if many
    df['BayesDel_addAF_score'] = df['BayesDel_addAF_score'].replace(".", np.NaN)    
    df.loc[df['BayesDel_addAF_score'].isna(),'BayesDel_addAF_score'] = 0 
    # fathmm coding group 
    df.loc[df['fathmm-MKL_coding_group'].isna(),'fathmm-MKL_coding_group'] = ""
    not_all = df[df['fathmm-MKL_coding_group']!="ALL"][['fathmm-MKL_coding_group']]
    df.loc[df['fathmm-MKL_coding_group'] =="ALL",'fathmm-MKL_coding_group'] = "".join(set("".join(list(not_all['fathmm-MKL_coding_group'].unique()))))
    for i in set("".join(list(df['fathmm-MKL_coding_group'].unique()))):
        df['fathmm-MKL_coding_group_' + i] = df['fathmm-MKL_coding_group'].str.contains(i).mul(1)
    df.loc[df["RegulomeDB_score"].isna(), "RegulomeDB_score"] = 0 
    df['variant_id'] = "chr" + df['chr'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['ref'] + "/" + df['alt']
    df = df.drop(['fathmm-MKL_coding_group', 'ref','alt', 'fathmm-XF_coding_or_noncoding','LIST-S2_score','MVP_score'], axis = 1)
    df = df.set_index(["variant_id"])
    for column in df.columns:
        try: df[column] = df[column].astype(float)
        except:
            df = df.drop(column, axis = 1)
            print(column, "COULD NOT BE REPRESENTED AS FLOAT AND IS BEING REMOVED")
    print("before removing NA", df.shape)
    df = df.dropna(how = 'any') # remove NA values
    print("after removing NA", df.shape)
    return df

def make_noncoding_numerical(df):
    # go through each column and make the dataframe numerical
    df = df.replace(".", np.NaN)
   # impute gnomad scores
    df['gnomAD_genomes_POPMAX_AF'] = df['gnomAD_genomes_POPMAX_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['gnomAD_genomes_AFR_AF'] = df['gnomAD_genomes_AFR_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)#mins of other populations are 0
    df['gnomAD_genomes_AMR_AF'] = df['gnomAD_genomes_AMR_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['gnomAD_genomes_NFE_AF'] = df['gnomAD_genomes_NFE_AF'].fillna(np.min(df['gnomAD_genomes_POPMAX_AF'].astype(float))/2)
    df['Eigen_coding_or_noncoding'] = np.where(df['Eigen_coding_or_noncoding']=='c', 1, 0)
    df['fathmm-XF_coding_or_noncoding'] = np.where(df['fathmm-XF_coding_or_noncoding']=='coding', 1, 0)
    df = process_spliceai_column(df, ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL'])
    
    df.loc[df["RegulomeDB_score"].isna(), "RegulomeDB_score"] = 0 # make regulomedb quantitative
    df.loc[df["GERP_RS"].isna(), "GERP_RS"] = 0 
    for val in df['RegulomeDB_score'].unique():
        new_val = float(re.search(r'\d+', str(val)).group())
        df.loc[df["RegulomeDB_score"] == val, "RegulomeDB_score"] = new_val
    df.loc[df['EnhancerFinder_brain_enhancer']=='Y','EnhancerFinder_brain_enhancer'] = 1
    df.loc[df['EnhancerFinder_brain_enhancer']!=1,'EnhancerFinder_brain_enhancer'] = 0
    df.loc[df['FANTOM5_CAGE_peak_robust']=='Y','FANTOM5_CAGE_peak_robust'] = 1
    df.loc[df['FANTOM5_CAGE_peak_robust']!=1,'FANTOM5_CAGE_peak_robust'] = 0
    # fathmm noncoding group
    df.loc[df['fathmm-MKL_non-coding_group'].isna(),'fathmm-MKL_non-coding_group'] = ""
    not_all = df[df['fathmm-MKL_non-coding_group']!="ALL"][['fathmm-MKL_non-coding_group']]
    df.loc[df['fathmm-MKL_non-coding_group'] =="ALL",'fathmm-MKL_non-coding_group'] = "".join(set("".join(list(not_all['fathmm-MKL_non-coding_group'].unique()))))
    for i in set("".join(list(df['fathmm-MKL_non-coding_group'].unique()))):
        df['fathmm-MKL_non-coding_group_' + i] = df['fathmm-MKL_non-coding_group'].str.contains(i).mul(1)
    df['variant_id'] = "chr" + df['chr'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['ref'] + "/" + df['alt']
    df = df.drop(['fathmm-MKL_non-coding_group', 'ref','alt'], axis = 1)
    df = df.set_index(["variant_id"])
    for column in df.columns:
        try: df[column] = df[column].astype(float)
        except:
            df = df.drop(column, axis = 1)
            print(column, "COULD NOT BE REPRESENTED AS FLOAT AND IS BEING REMOVED")
    print("before removing NA", df.shape)
    df = df.dropna(how = 'any') # remove NA values
    print("after removing NA", df.shape)
    return df


def main():
    to_extract = pd.read_csv(VARIANTS_EXTRACT, sep = "\t", header = None)
    if variant_type == "coding": annotations = coding_annotations
    if variant_type == "noncoding": annotations = noncoding_annotations

    output = os.path.join(OUTPUT_DIR, variant_type, f"WGSA_chr{chromosome}.csv")
    if os.path.isfile(output):
        print(output, "already exists, continue on")
        return None
    if not os.path.exists(os.path.join(OUTPUT_DIR)):
        os.mkdir(os.path.join(OUTPUT_DIR))
    print("working on ", output)
    df = pd.read_csv(ANNOTATION_FILE_PATH, sep="\t", usecols = annotations)
    print("Finished reading annotation file.")
    if variant_type == "coding":
        df = df[(df['Eigen_coding_or_noncoding'] =='c') | (df['fathmm-XF_coding_or_noncoding']=='coding')]
        coding_anno = get_extracted_df(df, to_extract)
        coding_anno = make_coding_numerical(coding_anno)
        coding_anno.to_csv(output)
        print("done with ", output)
    else:
        df = df[(df['Eigen_coding_or_noncoding'] =='n') | (df['fathmm-XF_coding_or_noncoding']=='noncoding')]
        noncoding_anno = get_extracted_df(df, to_extract)
        noncoding_anno = make_noncoding_numerical(noncoding_anno)
        noncoding_anno.to_csv(output)
        print("done with ", output)
    return

if __name__ == "__main__":
    main()



