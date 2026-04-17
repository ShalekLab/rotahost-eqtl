version 1.0

# RotaHost sc-eQTL + sc-pQTL WDL workflow
# Published on Dockstore: https://dockstore.org/organizations/ShalekLab
#
# Implements sc-eQTLGen WG3 approach (Kaptijn et al. 2026)
# Two-step design:
#   1. compute_windows — splits gene annotation into genomic windows (~215 genes each)
#   2. run_eqtl_task   — runs limix-qtl per (cell_type × window)
#
# Bgen input modes:
#   Mode 1 (genome-wide):      genome_bgen_path = "gs://.../all_chroms.bgen"
#   Mode 2 (per-chromosome):   chr_bgen_paths = ["gs://.../chr1.bgen", ..., "gs://.../chr22.bgen"]
#
# Toy test:
#   cell_types = ["B"], genes_per_window = 9999 (one window per chr)
#   chr_bgen_paths = ["gs://.../chr22.bgen"], chromosomes = ["22"]
#
# Full run: 5 analyses × 7 CT × ~53 windows ≈ 1,855 tasks
# Cost: ~$400 | Wall time: ~8h (all tasks parallel)

workflow rotahost_eqtl {
    input {
        String        analysis_mode    = "longitudinal"
        Array[String] cell_types       = ["B"]
        Int           n_permutations   = 100
        Float         maf_filter       = 0.05
        Float         hwe_filter       = 0.0001
        Int           cis_window       = 1000000

        # --- Bgen input: provide ONE mode ---
        # Mode 1: genome-wide bgen (single file, simplest UI)
        String        genome_bgen_path     = ""
        String        genome_bgi_path      = ""
        # Mode 2: per-chromosome bgens (22 entries each)
        Array[String] chromosomes          = []
        Array[String] chr_bgen_paths       = []
        Array[String] chr_bgi_paths        = []

        # Windowing: target ~100 genes/window ≈ 3.3h at 120s/gene (long-read pipeline pattern)
        # Lower = shorter tasks = lower preemption rate; higher = fewer tasks = lower overhead
        # 100 → ~3,255 tasks, ~20% preemption rate (recommended for full run)
        # 215 → ~1,855 tasks, ~40% preemption rate
        Int           genes_per_window     = 100

        # Shared input files
        File          donor_re_tsv
        File          gene_annotation_tsv
        File          sample_mapping_tsv

        # Per-cell-type input files (parallel to cell_types)
        Array[File]   phenotype_files
        Array[File]   covariate_files

        # Interaction term (empty = no interaction; "iga_status" | "pre_post")
        String interaction_term = ""

        String  docker_image       = "ghcr.io/shts2123/rotahost-eqtl:latest"
        Int     cpu                = 2
        Int     mem_gb             = 8
        Int     disk_gb            = 30
        # preemptible_tries=1: try once preemptible, then fall back to non-preemptible
        # (long-read pipeline / spectronaut pattern for tasks >1h)
        # 0 = always non-preemptible (most reliable); 3 = most preemptible savings
        Int     preemptible_tries  = 1

        # Optional: GCS folder to copy merged h5 results after merge_results
        # e.g. "gs://fc-secure-.../results/eqtl"
        # Leave empty to keep outputs only in the submission path
        String  results_gcs_dir    = ""
        # Disk for merge task — genome-wide run: ~53 windows × ~70MB h5 + merged output
        Int     merge_disk_gb      = 50
    }

    # Step 1: compute genomic windows from gene annotation
    call compute_windows {
        input:
            gene_annotation_tsv = gene_annotation_tsv,
            genes_per_window    = genes_per_window,
            cis_window          = cis_window,
            genome_bgen_path    = genome_bgen_path,
            genome_bgi_path     = genome_bgi_path,
            chromosomes         = chromosomes,
            chr_bgen_paths      = chr_bgen_paths,
            chr_bgi_paths       = chr_bgi_paths,
    }

    # Step 2: scatter over (cell_type × window)
    scatter (ct_idx in range(length(cell_types))) {
        scatter (wi in range(length(compute_windows.genomic_windows))) {
            call run_eqtl_task {
                input:
                    cell_type        = cell_types[ct_idx],
                    genomic_window   = compute_windows.genomic_windows[wi],
                    analysis_mode    = analysis_mode,
                    n_permutations   = n_permutations,
                    maf_filter       = maf_filter,
                    hwe_filter       = hwe_filter,
                    cis_window       = cis_window,
                    # String → File coercion: Cromwell localizes these GCS URIs
                    bgen_file        = compute_windows.bgen_per_window[wi],
                    bgen_bgi_file    = compute_windows.bgi_per_window[wi],
                    donor_re_tsv     = donor_re_tsv,
                    annotation_tsv   = gene_annotation_tsv,
                    feature_filter_file = compute_windows.feature_filter_files[wi],
                    sample_map_tsv   = sample_mapping_tsv,
                    phenotype_file   = phenotype_files[ct_idx],
                    covariate_file   = covariate_files[ct_idx],
                    interaction_term = interaction_term,
                    docker_image     = docker_image,
                    cpu              = cpu,
                    mem_gb           = mem_gb,
                    disk_gb          = disk_gb,
                    preemptible      = preemptible_tries,
            }
        }

        # Gather: merge all window h5s for this cell type, apply genome-wide FDR
        call merge_results {
            input:
                h5_files            = run_eqtl_task.h5_output,
                cell_type           = cell_types[ct_idx],
                analysis_mode       = analysis_mode,
                interaction_term    = interaction_term,
                docker_image        = docker_image,
                results_gcs_dir     = results_gcs_dir,
                disk_gb             = merge_disk_gb,
        }
    }

    output {
        Array[String]      windows_used = compute_windows.genomic_windows
        Array[Array[File]] h5_results   = run_eqtl_task.h5_output
        Array[Array[File]] log_files    = run_eqtl_task.log_file
        Array[File]        top_results  = merge_results.top_results
        Array[File]        merged_h5    = merge_results.merged_h5
        Array[File]        metrics      = merge_results.metrics
    }

    meta {
        author:      "Sergio Triana, Shalek Lab"
        email:       "strianas@mit.edu"
        description: "RotaHost sc-eQTL + sc-pQTL (sc-eQTLGen WG3, Kaptijn et al. 2026). Auto-windowing, flexible bgen input."
    }
}

# ============================================================================
# Task 1: Compute genomic windows from gene annotation
# ============================================================================
# Splits genes into chunks of ~genes_per_window, padded by ±cis_window.
# Accepts genome-wide OR per-chromosome bgen paths (as Strings, no localization).
# Outputs parallel arrays of windows + bgen GCS URIs for String→File coercion.

task compute_windows {
    input {
        File          gene_annotation_tsv
        Int           genes_per_window   = 100
        Int           cis_window         = 1000000

        # Mode 1: genome-wide
        String        genome_bgen_path   = ""
        String        genome_bgi_path    = ""
        # Mode 2: per-chromosome
        Array[String] chromosomes        = []
        Array[String] chr_bgen_paths     = []
        Array[String] chr_bgi_paths      = []
    }

    command <<<
        python3 - << 'PYEOF'
        import csv, json, math, os

        genome_bgen = "~{genome_bgen_path}".strip()
        genome_bgi  = "~{genome_bgi_path}".strip()
        raw_chroms  = "~{sep=',' chromosomes}"
        raw_bgens   = "~{sep=',' chr_bgen_paths}"
        raw_bgis    = "~{sep=',' chr_bgi_paths}"
        cis         = ~{cis_window}
        gpw         = ~{genes_per_window}

        chroms = [x for x in raw_chroms.split(',') if x]
        bgens  = [x for x in raw_bgens.split(',') if x]
        bgis   = [x for x in raw_bgis.split(',') if x]

        bgen_map = dict(zip(chroms, bgens)) if not genome_bgen else {}
        bgi_map  = dict(zip(chroms, bgis))  if not genome_bgen else {}

        # Read gene annotation (autosomal only) — keep (start, feature_id) tuples
        genes = {}
        with open("~{gene_annotation_tsv}") as f:
            for row in csv.DictReader(f, delimiter='\t'):
                c = row['chromosome']
                if not c.isdigit():
                    continue
                if not genome_bgen and c not in bgen_map:
                    continue
                genes.setdefault(int(c), []).append((int(row['start']), row['feature_id']))

        os.makedirs("feature_filters", exist_ok=True)
        windows, win_bgens, win_bgis = [], [], []
        # Use global zero-padded sequence numbers so glob() returns files in window order
        global_idx = 0
        for chrom in sorted(genes):
            pairs = sorted(genes[chrom])  # sort by (start, feature_id)
            n_win = max(1, math.ceil(len(pairs) / gpw))
            chunk = math.ceil(len(pairs) / n_win)
            for i in range(n_win):
                sl = pairs[i * chunk : (i + 1) * chunk]
                starts = [p[0] for p in sl]
                fids   = [p[1] for p in sl]
                w_start = max(1, starts[0] - cis)
                w_end   = starts[-1] + cis
                windows.append(f"{chrom}:{w_start}-{w_end}")
                win_bgens.append(genome_bgen if genome_bgen else bgen_map[str(chrom)])
                win_bgis.append(genome_bgi  if genome_bgen else bgi_map[str(chrom)])
                # Write per-window feature filter TSV (limix-qtl reads index_col=0)
                # 6-digit zero-pad supports up to 999,999 windows; sorts alphabetically = window order
                ff_path = f"feature_filters/win{global_idx:06d}.tsv"
                with open(ff_path, 'w') as out:
                    out.write("feature_id\n")
                    for fid in fids:
                        out.write(fid + "\n")
                global_idx += 1

        print(json.dumps({
            "windows":  windows,
            "bgens":    win_bgens,
            "bgis":     win_bgis,
        }))
        PYEOF
    >>>

    output {
        Array[String] genomic_windows      = read_json(stdout())["windows"]
        Array[String] bgen_per_window      = read_json(stdout())["bgens"]
        Array[String] bgi_per_window       = read_json(stdout())["bgis"]
        # glob returns files in alphabetical order — matches window order via zero-padded names
        Array[File]   feature_filter_files = glob("feature_filters/win*.tsv")
    }

    runtime {
        docker:      "python:3.11-slim"
        cpu:         1
        memory:      "2 GB"
        disks:       "local-disk 10 SSD"
        preemptible: 0
    }

    meta {
        description: "Splits gene annotation into genomic windows for parallelized eQTL mapping."
    }
}

# ============================================================================
# Task 2: Run limix-qtl per cell type × genomic window
# ============================================================================

task run_eqtl_task {
    input {
        String  cell_type
        String  genomic_window    # "CHR:START-END"
        String  analysis_mode
        Int     n_permutations
        Float   maf_filter
        Float   hwe_filter
        Int     cis_window
        File    bgen_file         # localized by Terra/Cromwell (String→File coercion)
        File    bgen_bgi_file     # matching index file
        File    donor_re_tsv
        File    annotation_tsv
        File    feature_filter_file  # per-window TSV: only test these features (no boundary duplicates)
        File    sample_map_tsv
        File    phenotype_file
        File    covariate_file
        String  interaction_term
        String  docker_image
        Int     cpu
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    # Parse window string "CHR:START-END" for output filename
    # limix-qtl writes {prefix}_{chr}_{start}_{end}.h5 when -gr is specified
    String chromosome   = sub(genomic_window, ":.*", "")
    String window_range = sub(genomic_window, ".*:", "")
    String window_start = sub(window_range, "-.*", "")
    String window_end   = sub(genomic_window, ".*-", "")

    # NOTE: bgen_prefix must be computed in bash (not WDL String expression)
    # because WDL File expressions resolve to GCS URIs in String context,
    # not the localized local path.

    String limix_script = if (interaction_term != "")
        then "/limix_qtl/Limix_QTL/run_interaction_QTL_analysis.py"
        else "/limix_qtl/Limix_QTL/run_QTL_analysis.py"

    String inter_flag = if (interaction_term != "")
        then "-it ~{interaction_term}"
        else ""

    # Output h5 prefix differs by analysis type
    String h5_prefix = if (interaction_term != "") then "iqtl_results" else "qtl_results"

    command <<<
        set -euo pipefail

        echo "=== RotaHost eQTL: ~{cell_type} window ~{genomic_window} ==="
        echo "Analysis: ~{analysis_mode} | Perms: ~{n_permutations} | Donor RE: chip-indexed identity"

        mkdir -p output/

        # Compute bgen prefix from the localized file path (strip .bgen extension)
        BGEN_LOCAL="~{bgen_file}"
        BGEN_PREFIX="${BGEN_LOCAL%.bgen}"

        # Filter sample_mapping to only samples present in:
        #   (a) this CT's covariate file
        #   (b) this CT's phenotype file
        #   (c) the donor random effect matrix (kinship file)
        # Without (c), interaction script gets NaN kinship rows for unmatched chips.
        python << 'FILTERPY'
        import pandas as pd
        smap = pd.read_csv("~{sample_map_tsv}", sep="\t", header=None, names=["chip","sample"])
        cov  = pd.read_csv("~{covariate_file}", sep="\t", index_col=0)
        ph   = pd.read_csv("~{phenotype_file}", sep="\t", index_col=0, nrows=1)
        kin  = pd.read_csv("~{donor_re_tsv}", sep="\t", index_col=0, nrows=1)  # just need col names
        valid_samples = set(cov.index) & set(ph.columns)
        valid_chips   = set(kin.columns)
        keep = smap[smap["sample"].isin(valid_samples) & smap["chip"].isin(valid_chips)]
        keep.to_csv("filtered_sample_mapping.tsv", sep="\t", header=False, index=False)
        print(f"Sample mapping filtered: {len(smap)} -> {len(keep)} "
              f"(unique chips: {smap['chip'].nunique()} -> {keep['chip'].nunique()})")
        FILTERPY

        python -u ~{limix_script} \
            --bgen   "$BGEN_PREFIX" \
            -af      ~{annotation_tsv} \
            -ff      ~{feature_filter_file} \
            -cf      ~{covariate_file} \
            -pf      ~{phenotype_file} \
            -rf      ~{donor_re_tsv} \
            -smf     filtered_sample_mapping.tsv \
            -od      output/ \
            -gr      ~{genomic_window} \
            -np      ~{n_permutations} \
            -maf     ~{maf_filter} \
            -hwe     ~{hwe_filter} \
            -c -gm gaussnorm \
            -w       ~{cis_window} \
            ~{inter_flag} \
            2>&1 | tee output/run.log

        echo "Completed: ~{cell_type} ~{genomic_window}" | tee -a output/run.log

    >>>

    output {
        File h5_output = "output/~{h5_prefix}_~{chromosome}_~{window_start}_~{window_end}.h5"
        File log_file  = "output/run.log"
    }

    runtime {
        docker:       docker_image
        cpu:          cpu
        memory:       "~{mem_gb} GB"
        disks:        "local-disk ~{disk_gb} SSD"
        preemptible:  preemptible
        maxRetries:   preemptible
        # NOTE: checkpointFile is NOT supported in GCPBatch mode (Terra backend since June 2025).
        # The Docker image includes h5 append-mode patches for future use when/if
        # Terra adds GCPBatch checkpoint support.
    }

    meta {
        description: "limix-qtl cis-eQTL per cell type × genomic window. Checkpoint-enabled for preemption resilience."
    }
}

# ============================================================================
# Task 3: Merge per-window h5 results and apply genome-wide FDR
# ============================================================================
# Runs after all scatter tasks for one cell type complete.
# Reads all window h5 files, extracts best SNP per gene, applies BH-FDR.

task merge_results {
    input {
        Array[File] h5_files
        String      cell_type
        String      analysis_mode
        String      interaction_term
        String      docker_image
        # Optional GCS folder to copy all 3 output files
        # e.g. "gs://fc-secure-.../eqtl_results"
        # All 35 runs (7 CT × 5 analyses) can share one folder — filenames are unique
        String      results_gcs_dir = ""
        Int         disk_gb         = 50
    }

    String h5_prefix       = if (interaction_term != "") then "iqtl" else "eqtl"
    String merged_h5_name  = "~{h5_prefix}_~{cell_type}_~{analysis_mode}_merged.h5"
    String top_results_name = "~{h5_prefix}_~{cell_type}_~{analysis_mode}_top_results.tsv"
    String metrics_name     = "~{h5_prefix}_~{cell_type}_~{analysis_mode}_metrics.tsv"

    command <<<
        set -euo pipefail

        python3 -u - << 'EOF'
        import h5py, numpy as np, os, sys

        h5_paths = "~{sep=' ' h5_files}".split()
        cell_type = "~{cell_type}"
        analysis = "~{analysis_mode}"
        merged_name     = "~{merged_h5_name}"
        top_results_name = "~{top_results_name}"
        metrics_name    = "~{metrics_name}"
        results_gcs_dir = "~{results_gcs_dir}".strip()

        print(f"Merging {len(h5_paths)} h5 files for {cell_type} ({analysis})")

        header = "feature_id\tsnp_id\tp_value\tbeta\tbeta_se\tempirical_feature_p_value\tn_snps_tested"
        rows = []
        seen_genes = set()  # genes appear in multiple windows due to ±cis_window overlap

        # Write merged h5 (all variants for all genes — full summary stats for plotting/colocalization)
        with h5py.File(merged_name, 'w') as out_f:
            for h5_path in h5_paths:
                with h5py.File(h5_path, 'r') as f:
                    for gene in f.keys():
                        data = f[gene][:]
                        if len(data) == 0 or gene in seen_genes:
                            continue
                        seen_genes.add(gene)
                        out_f.create_dataset(gene, data=data)
                        best_idx = int(np.argmin(data['empirical_feature_p_value']))
                        rows.append((
                            gene,
                            data['snp_id'][best_idx].decode(),
                            float(data['p_value'][best_idx]),
                            float(data['beta'][best_idx]),
                            float(data['beta_se'][best_idx]),
                            float(data['empirical_feature_p_value'][best_idx]),
                            len(data),
                        ))

        print(f"  Total genes: {len(rows)}")
        print(f"  Merged h5 written: {merged_name} ({os.path.getsize(merged_name) / 1e6:.1f} MB)")

        if len(rows) == 0:
            print("WARNING: No genes found in any h5 file")
            with open(top_results_name, 'w') as f:
                f.write(header + "\tfdr_gene\n")
            with open(metrics_name, 'w') as f:
                f.write("cell_type\tanalysis_mode\tn_genes\tn_egenes_05\tn_egenes_20\ttop_gene\ttop_emp_p\n")
                f.write(f"{cell_type}\t{analysis}\t0\t0\t0\tNA\tNA\n")
            sys.exit(0)

        # Sort by empirical p-value
        rows.sort(key=lambda x: x[5])

        # BH-FDR correction across all genes
        emp_pvals = np.array([r[5] for r in rows])
        n = len(emp_pvals)
        ranked = np.argsort(emp_pvals)
        fdr = np.zeros(n)
        cummin = 1.0
        for i in range(n - 1, -1, -1):
            idx = ranked[i]
            bh = emp_pvals[idx] * n / (i + 1)
            cummin = min(cummin, bh)
            fdr[idx] = min(cummin, 1.0)

        n_egenes_05 = int(np.sum(fdr < 0.05))
        n_egenes_20 = int(np.sum(fdr < 0.20))
        print(f"  eGenes (FDR<0.05): {n_egenes_05} | eGenes (FDR<0.20): {n_egenes_20}")

        # top_results tsv — best SNP per gene + FDR, named uniquely
        with open(top_results_name, 'w') as f:
            f.write(header + "\tfdr_gene\n")
            for i, row in enumerate(rows):
                f.write("\t".join(str(x) for x in row) + f"\t{fdr[i]:.6g}\n")
        print(f"  Written to {top_results_name}")

        # metrics tsv — one-line summary for quick QC
        top_gene = rows[0][0]
        top_emp_p = rows[0][5]
        with open(metrics_name, 'w') as f:
            f.write("cell_type\tanalysis_mode\tn_genes\tn_egenes_05\tn_egenes_20\ttop_gene\ttop_emp_p\n")
            f.write(f"{cell_type}\t{analysis}\t{n}\t{n_egenes_05}\t{n_egenes_20}\t{top_gene}\t{top_emp_p:.6g}\n")
        print(f"  Written to {metrics_name}")

        # Upload all 3 files to results_gcs_dir if provided
        if results_gcs_dir:
            from google.cloud import storage
            client = storage.Client()
            base = results_gcs_dir.rstrip("/")
            for local_file in [merged_name, top_results_name, metrics_name]:
                dest = base + "/" + local_file
                dest_no_prefix = dest[len("gs://"):]
                bucket_name, blob_path = dest_no_prefix.split("/", 1)
                bucket = client.bucket(bucket_name)
                blob = bucket.blob(blob_path)
                blob.upload_from_filename(local_file)
                print(f"  Uploaded {local_file} → {dest}")

        EOF
    >>>

    output {
        File top_results = top_results_name
        File merged_h5   = merged_h5_name
        File metrics     = metrics_name
    }

    runtime {
        docker:      docker_image
        cpu:         1
        memory:      "4 GB"
        disks:       "local-disk ~{disk_gb} SSD"
        preemptible: 0
    }

    meta {
        description: "Merge per-window eQTL h5 results (full variant stats), compute genome-wide BH-FDR, write metrics summary."
    }
}
