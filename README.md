# joint-model-analysis

> [!WARNING]  
> This project is still under development.

## todo

- [ ] Introdution
- [ ] Workflow Steps
- [X] Inputs
- [X] Outputs
- [X] Components

## Introdution

## Workflow Steps

## Inputs

- `analysis_script`: Path to the R script for analysis.
- `serum_csv`: Path to the serum data CSV file.
- `demographics_csv`: Path to the demographics data CSV file.
- `protein_name`: Name of the protein to analyze.
- `docker`: Docker image to use for the analysis.
- `memory_gb`: Memory allocation for the task (default: `24` GB).
- `cpu`: CPU allocation for the task (default: `16`).

## Outputs

- `summary`: CSV file containing the joint model results.
- `rds`: RDS file containing the joint model object.

## Components

- Docker image
  - `rocker/tidyverse`
- R packages
  - survival
  - JMbayes2
  - nlme
