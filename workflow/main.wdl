version 1.0

import "./tasks/analysis.wdl" as process

workflow main {
    input {
        File serum_csv
        File demographics_csv
        String protein_name
        Int memory_gb = 24
        Int cpu = 16
    }
    String docker = "ghcr.io/anand-imcm/joint-model-analysis:1.0.0"
    call process.analysis{
        input:
            serum=serum_csv,
            demog=demographics_csv,
            protein=protein_name,
            memory=memory_gb,
            cpu=cpu,
            docker=docker
    }
    output{
        File summary = analysis.csv
        File rds = analysis.rds
    }
}
