version 1.0

task analysis {
    input {
        File serum
        File demog
        String protein
        String docker
        Int memory
        Int cpu
        Int disk_size = 500
        String disk_type = "SSD"
    }
    command <<<
        Rscript /scripts/analysis.R ~{serum} ~{demog} ~{protein}
    >>>
    output {
        File csv = protein + "_joint_results.csv"
        File rds = protein + "_joint_results.RDS"
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_size} ~{disk_type}"
    }
}