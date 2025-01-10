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

1) Map variants to functional annotations: 
    - Mapping WGSA annotations: Run `/gpfs/commons/home/adas/gruyere/src/preprocess/map_annotations.sh`

2) Extract genotypes of interest:
    - First create input files for PLINK: Run `/gpfs/commons/home/adas/gruyere/src/preprocess/make_input_plink.sh`
    - Then extract genotypes for unrelated individuals: `/gpfs/commons/home/adas/gruyere/src/preprocess/extract_genotypes.sh`
    - We select coding and non-coding variants of interest that pass QC (genotyping rate > 90%, individual missingness > 90%). 
        - Genotypes are saved `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/genotypes/{coding/noncoding}/ADSP.chr{CHR}.vcf`
        - Compress genotypes `/gpfs/commons/home/adas/gruyere/src/preprocess/vcf_comp.sh`

3) Efficiently get VEP and missense annotations:
    - Run `/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/VEP/ensembl-vep/VEP_full.sh`
    
4) Calculate LD matrices with GWAS significant loci if conditional analysis
    - Run `/gpfs/commons/home/adas/gruyere/src/preprocess/LD_gwas.sh`
    
4) Generate per-gene tensors:
    - Selecting only variants that are present both in genotype and annotation datasets (steps 1 & 2), generate efficient tensors for loading
    - Genotype tensors: per chromosome individual by genotype (separate for coding, non-coding)
    - Annotation tensors: per chromosome variant by annotation (separate for coding, non-coding)
    - Variant-Gene index tensors: per chromosome variant-gene index mapping (separate for coding, non-coding & includes ABC annotations for non-coding)
    - Run `/gpfs/commons/home/adas/gruyere/src/preprocess/gene_tensors.sh`
