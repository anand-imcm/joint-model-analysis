version 1.0

task analysis {
    input {
        File script
        File serum
        File demog
        String protein
        String docker
        Int memory
        Int cpu
    }
    Array[File] all_data = flatten([[serum], [demog]])
    Int disk_size_gb = ceil(size(all_data, "GB")) + 5
    command <<<
        Rscript ~{script} ~{serum} ~{demog} ~{protein}
    >>>
    output {
        File csv = protein + "_joint_results.csv"
        File rds = protein + "_joint_results.RDS"
    }
    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}