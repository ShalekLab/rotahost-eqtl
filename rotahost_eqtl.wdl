version 1.0

# RotaHost sc-eQTL + sc-pQTL WDL workflow
# Published on Dockstore: https://dockstore.org/organizations/ShalekLab
#
# Implements sc-eQTLGen WG3 approach (Kaptijn et al. 2026)
# Supports 6 analysis modes:
#   global | longitudinal | iga_interaction | temporal_igay | temporal_igan | adt
#
# Toy test parameters (chr22, B_L1, 10 perms):
#   cell_types = ["B"]
#   chromosomes = ["22"]
#   n_permutations = 10
#   analysis_mode = "global"
#
# Full run parameters:
#   cell_types = ["B","CD4_T","CD8_T","NK","Mono","DC","gdT"]
#   chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
#   n_permutations = 1000
#   analysis_mode = "global" (run separately for each mode)
#
# Cost monitoring: estimated $0.03 per task (n1-standard-4 preemptible, ~40 min)
# Set GCP budget alert at $2 for safety; abort if exceeded

workflow rotahost_eqtl {
    input {
        # Analysis configuration
        String  analysis_mode      = "global"   # global|longitudinal|iga_interaction|temporal_igay|temporal_igan|adt
        Array[String] cell_types   = ["B"]      # L1 cell types to analyze
        Array[String] chromosomes  = ["22"]     # chromosomes to test
        Int     n_permutations     = 10         # 10 for toy test; 1000 for production
        Float   maf_filter         = 0.05
        Float   hwe_filter         = 0.0001
        Int     cis_window         = 1000000    # ±1 Mb

        # GCS paths to input files (upload from local before submission)
        String  gcs_bgen_prefix               # gs://.../ directory containing chr{N}.bgen files
        String  gcs_kinship_tsv               # gs://.../kinship_matrix.tsv
        File    gene_annotation_tsv           # local or GCS path to gene annotation
        File    sample_mapping_tsv            # local or GCS path (no header)

        # Per-cell-type input directories (GCS prefix)
        String  gcs_phenotype_dir             # gs://.../phenotypes/{analysis_mode}/
        String  gcs_covariate_dir             # gs://.../covariates/{analysis_mode}/

        # Optional: interaction term for analyses 3-5
        String? interaction_term              # "iga_status" or "timepoint_group"

        # Docker image (GHCR → ShalekLab)
        String  docker_image = "ghcr.io/shaleklab/rotahost-eqtl:latest"

        # Compute resources
        Int     cpu     = 2
        Int     mem_gb  = 8
        Int     disk_gb = 30
        Int     preemptible_tries = 3
    }

    # Scatter over all (cell_type, chromosome) pairs
    scatter (ct in cell_types) {
        scatter (chrom in chromosomes) {
            call run_eqtl_task {
                input:
                    cell_type        = ct,
                    chromosome       = chrom,
                    analysis_mode    = analysis_mode,
                    n_permutations   = n_permutations,
                    maf_filter       = maf_filter,
                    hwe_filter       = hwe_filter,
                    cis_window       = cis_window,
                    bgen_prefix      = "~{gcs_bgen_prefix}chr~{chrom}",
                    kinship_tsv      = gcs_kinship_tsv,
                    annotation_tsv   = gene_annotation_tsv,
                    sample_map_tsv   = sample_mapping_tsv,
                    phenotype_dir    = gcs_phenotype_dir,
                    covariate_dir    = gcs_covariate_dir,
                    interaction_term = interaction_term,
                    docker_image     = docker_image,
                    cpu              = cpu,
                    mem_gb           = mem_gb,
                    disk_gb          = disk_gb,
                    preemptible      = preemptible_tries,
            }
        }
    }

    output {
        Array[Array[File]] h5_results = run_eqtl_task.h5_output
        Array[Array[File]] log_files  = run_eqtl_task.log_file
    }

    meta {
        author:      "Sergio Triana, Shalek Lab"
        email:       "strianas@mit.edu"
        description: "RotaHost sc-eQTL + sc-pQTL pipeline (sc-eQTLGen WG3 approach)"
    }
}

task run_eqtl_task {
    input {
        String  cell_type
        String  chromosome
        String  analysis_mode
        Int     n_permutations
        Float   maf_filter
        Float   hwe_filter
        Int     cis_window
        String  bgen_prefix       # GCS prefix without .bgen extension
        String  kinship_tsv       # GCS path or local path
        File    annotation_tsv
        File    sample_map_tsv
        String  phenotype_dir     # GCS directory prefix
        String  covariate_dir     # GCS directory prefix
        String? interaction_term
        String  docker_image
        Int     cpu
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String limix_script = if (defined(interaction_term))
        then "/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py"
        else "/limix_qtl/Limix_QTL/run_QTL_analysis.py"

    String inter_flag = if (defined(interaction_term))
        then "-int ~{interaction_term}"
        else ""

    # Phenotype file: denoised_adt for adt mode, expression for others
    String pheno_suffix = if (analysis_mode == "adt")
        then "_denoised_adt.csv"
        else "_expression.tsv"

    command <<<
        set -euo pipefail

        # Download BGEN files from GCS
        gsutil cp "~{bgen_prefix}.bgen"     ./chr~{chromosome}.bgen
        gsutil cp "~{bgen_prefix}.bgen.bgi" ./chr~{chromosome}.bgen.bgi

        # Download kinship matrix
        gsutil cp "~{kinship_tsv}" ./kinship_matrix.tsv

        # Download phenotype and covariate files for this cell type
        gsutil cp "~{phenotype_dir}/~{cell_type}~{pheno_suffix}" \
                  ./phenotype_~{cell_type}.tsv
        gsutil cp "~{covariate_dir}/~{cell_type}_covariates.tsv" \
                  ./covariate_~{cell_type}.tsv

        mkdir -p output/

        python ~{limix_script} \
            --bgen   chr~{chromosome} \
            -af      ~{annotation_tsv} \
            -cf      ./covariate_~{cell_type}.tsv \
            -pf      ./phenotype_~{cell_type}.tsv \
            -rf      ./kinship_matrix.tsv \
            -smf     ~{sample_map_tsv} \
            -od      output/ \
            -gr      ~{chromosome} \
            -np      ~{n_permutations} \
            -maf     ~{maf_filter} \
            -hwe     ~{hwe_filter} \
            -c -gm gaussnorm \
            -w       ~{cis_window} \
            ~{inter_flag} \
            > output/run.log 2>&1

        echo "Completed: ~{cell_type} chr~{chromosome}" >> output/run.log
    >>>

    output {
        File h5_output = "output/qtl_results_~{chromosome}.h5"
        File log_file  = "output/run.log"
    }

    runtime {
        docker:       docker_image
        cpu:          cpu
        memory:       "~{mem_gb} GB"
        disks:        "local-disk ~{disk_gb} SSD"
        preemptible:  preemptible
        maxRetries:   preemptible
    }

    meta {
        description: "Run limix_qtl cis-eQTL or interaction-eQTL for one cell type × chromosome"
    }
}
