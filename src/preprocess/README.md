## Here, we detail all steps to get from VCF to per-gene genotype and functional annotation tensors for modelling.

Input data:
- VCF file: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated/ADSP.chr{CHR}.vcf.gz`
    - Individuals to keep: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract2/keep_qc.txt`
        - File will be used in --keep by plink. Formatted FAM_ID \t IID
        - Criteria: Controls are 65+, >90% genotyping rate, unrelated.
    - Variant regions to keep: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/to_extract2/{CELL}/ADSP.chr{CHR}.txt`
        - File will be used in --extract range by plink. Formatted CHR \t START \t END \t GENE
        - Script to generate file: `/gpfs/commons/home/adas/gruyere/src/preprocess/make_files_for_plink.py`
- Functional annotations: 
    - WGSA: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/WGSA_annotations/ADSP_annotated_chr{CHR}.annotated.snp.gz`
    - VEP
- Covariates: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/phenotypes/36K_QC_filtered_final.csv`
- ABC Scores by cell type: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_all/glass_lab_fastq/processed_files/ABC_data/ABC_results_{CELL}/{CELL}/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.tsv`
    - Scripts to run ABC scores: `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/code/ABC/abc.sh`
    - My scripts are outdated and result in hg19. Using Chirag's runs from hg38.
    
    
    
Preprocessing scripts and steps:

1) Extract genotypes of interest: `/gpfs/commons/home/adas/gruyere/src/preprocess/extract_genotypes.sh`.
    - This script extracts genotypes for unrelated individuals. 
    - We select coding and non-coding variants of interest that pass QC (genotyping rate > 90%, individual missingness > 90%). 
    - Genotypes are saved `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/{coding/noncoding}/ADSP.chr{CHR}.vcf`
2) Map variants to functional annotations: 
    - Mapping VEP/missense annotations: Run `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/VEP/ensembl-vep/VEP_full.sh`
    - Mapping WGSA annotations: Run `/gpfs/commons/home/adas/gruyere/src/preprocess/map_annotations.sh`