# gruyere: Genome-wide Rare Variant EnRichment Evaluation <img src="figures/logo3.png" align="right" height="100"/>

We develop a genome-wide rare variant association test designed for identifying trait-associated loci and functional annotations. This repository accompanies our recent preprint: **Leveraging functional annotations to map rare variants associated with Alzheimer’s disease with gruyere**.


![gruyere model](figures/overview.png)

## Installation

gruyere is written in Python. You can load gruyere along with required dependencies with the following:

``` r
git clone https://github.com/daklab/gruyere.git
cd gruyere
pip install -r requirements.txt OR conda create --name gruyere --file requirements.txt
```

## Overview: Inputs and Outputs
**Model Inputs** 
- G: Genotypes for N individuals and P variants [P x N]. Index should contain gene name that variant maps to. Can optionally include variant id "gene_variantID"
- Z: Functional annotations for P variants and Q annotations [P x Q]. Index should contain gene name that variant maps to. Can optionally include variant id "gene_variantID"
- XY: Individual-level covariates for N individuals and C covariates [NxC] and "Diagnosis" column for binary or continuous phenotypes

**Model Outputs (Joint analysis)**
- Joint analysis: 


## Example:
1) Adjust `example_data/inputs.yaml`accordingly

``` r
sbatch example_bash/joint_analysis.sh
sbatch example_bash/per_gene.sh
```

