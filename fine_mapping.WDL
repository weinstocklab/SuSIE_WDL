version 1.0


workflow susie_finemapping {

    input {

        ## Summary statistics
        ## Concatenated REGENIE output parquet with all phenotypes included.
        File            summary_stats_parquet

        ## Genotype files for LD computation.
        ## Accepts the REGENIE step2 chunk manifest or step2 chromosome manifest.
       
        File            step2_chunk_manifest

        ## GWAS covariates file used to define the analysis sample.
        ## The workflow extracts IID/EID and builds a PLINK keep file internally.
        File            covariates_file

        ## Fine-mapping parameters 
        Float           p_threshold     = 5e-8   # genome-wide significance
        Float           min_maf         = 0.001  # minimum minor allele frequency for locus discovery
        Int             window_kb       = 500     # +/- kb around index SNP
        Int             max_causal      = 10      # L parameter in SuSiE
        Float           min_abs_corr    = 0.5     # min |r| for credible sets
        Float           coverage        = 0.95    # credible set coverage
        Int             max_iter        = 500     # SuSiE max iterations
        Float           tol             = 1e-3    # SuSiE convergence tol
        Boolean         refine          = true    # use SuSiE refinement
        Boolean         estimate_residual_variance = true

        ## Plink2 binary 
        File            plink2_binary

        ## Runtime 
        Int             threads     = 4
        String          r_docker    = "rocker/r-ver:4.3.0"
    }

    call extract_sample_keep {
        input:
            covariates_file = covariates_file
    }

    ## Step 1: Identify significant loci per phenotype.
    call identify_loci {
        input:
            summary_stats_parquet = summary_stats_parquet,
            step2_chunk_manifest  = step2_chunk_manifest,
            p_threshold           = p_threshold,
            min_maf               = min_maf,
            window_kb             = window_kb
    }

 
    ## Steps 2-4: For each significant locus, run full LD computation + SuSiE pipeline

    Array[Array[String]] loci   = read_tsv(identify_loci.loci_manifest)
    Int                  n_loci = length(loci) - 1   # remove header row

    scatter (i in range(n_loci)) {

        String locus_phenotype = loci[i + 1][0]
        String locus_chrom     = loci[i + 1][1]
        String locus_index_snp = loci[i + 1][3]
        String locus_start     = loci[i + 1][4]
        String locus_end       = loci[i + 1][5]
        String locus_pvar_path = loci[i + 1][6]
        String locus_psam_path = loci[i + 1][7]
        String locus_pgen_path = loci[i + 1][8]

        ## 2. Compute LD matrix
        call compute_ld {
            input:
                pvar          = "dx://" + locus_pvar_path,
                psam          = "dx://" + locus_psam_path,
                pgen          = "dx://" + locus_pgen_path,
                plink2_prefix = basename(locus_pvar_path, ".pvar"),
                plink2_binary = plink2_binary,
                chrom         = locus_chrom,
                start         = locus_start,
                end           = locus_end,
                sample_keep   = extract_sample_keep.sample_keep_file,
                out_prefix    = "ld_~{locus_chrom}_~{locus_start}_~{locus_end}_~{locus_phenotype}",
                threads       = threads
        }

        ## 3. Prepare per-locus summary stats
        call prepare_locus_sumstats {
            input:
                summary_stats_parquet = summary_stats_parquet,
                phenotype             = locus_phenotype,
                chrom                 = locus_chrom,
                start                 = locus_start,
                end                   = locus_end,
                variant_ids           = compute_ld.variant_ids
        }

        ## 4. Run SuSiE
        call run_susie {
            input:
                z_scores                   = prepare_locus_sumstats.z_scores_tsv,
                ld_matrix                  = compute_ld.ld_matrix,
                variant_ids                = compute_ld.variant_ids,
                phenotype                  = locus_phenotype,
                chrom                      = locus_chrom,
                start                      = locus_start,
                end                        = locus_end,
                index_snp                  = locus_index_snp,
                n_samples                  = compute_ld.n_samples,
                max_causal                 = max_causal,
                min_abs_corr               = min_abs_corr,
                coverage                   = coverage,
                max_iter                   = max_iter,
                tol                        = tol,
                refine                     = refine,
                estimate_residual_variance = estimate_residual_variance,
                docker                     = r_docker
        }

        ## 5. Create locuszoom-style plot with SuSiE overlays
        call plot_locus {
            input:
                z_scores          = prepare_locus_sumstats.z_scores_tsv,
                pip_tsv           = run_susie.pip_tsv,
                credible_sets_tsv = run_susie.credible_sets_tsv,
                ld_matrix         = compute_ld.ld_matrix,
                variant_ids       = compute_ld.variant_ids,
                phenotype         = locus_phenotype,
                chrom             = locus_chrom,
                start             = locus_start,
                end               = locus_end,
                index_snp         = locus_index_snp,
                p_threshold       = p_threshold,
                docker            = r_docker
        }
    }

    ## Step 6: Combine all loci results
    
    call combine_results {
        input:
            pip_tsvs          = run_susie.pip_tsv,
            credible_set_tsvs = run_susie.credible_sets_tsv
    }

    output {
        File          merged_pip           = combine_results.merged_pip
        File          merged_credible_sets = combine_results.merged_credible_sets
        Array[File]   susie_rds_all        = run_susie.susie_rds
        Array[File]   susie_logs           = run_susie.log
        Array[File]   locuszoom_pdfs       = plot_locus.plot_pdf
        Array[File]   locuszoom_pngs       = plot_locus.plot_png
        Array[File]   locuszoom_plot_data  = plot_locus.plot_data_tsv
        File          loci_manifest        = identify_loci.loci_manifest
    }
}


## TASKS

task extract_sample_keep {

    input {
        File covariates_file
    }

    command <<<
        set -euo pipefail

        python3 -m venv .venv
        .venv/bin/pip install pandas --quiet

        .venv/bin/python <<'PYEOF'
        import pandas as pd

        df = pd.read_csv("~{covariates_file}", sep=None, engine="python", dtype=str)

        normalized = {}
        for col in df.columns:
            key = col.strip().lstrip("#").upper()
            normalized[key] = col

        sample_col = None
        for candidate in ("IID", "EID"):
            if candidate in normalized:
                sample_col = normalized[candidate]
                break

        if sample_col is None:
            raise ValueError(
                "covariates_file must contain an IID or EID column in the header"
            )

        sample_ids = (
            df[sample_col]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        sample_ids = sample_ids[sample_ids != ""].drop_duplicates()

        keep_df = pd.DataFrame({"#IID": sample_ids})
        keep_df.to_csv("sample_keep.tsv", sep="\t", index=False)

        print(f"Derived sample keep file from column: {sample_col}")
        print(f"Samples retained: {len(keep_df)}")
        PYEOF
    >>>

    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        File sample_keep_file = "sample_keep.tsv"
    }
}


task identify_loci {

    input {
        File    summary_stats_parquet
        File    step2_chunk_manifest
        Float   p_threshold
        Float   min_maf
        Int     window_kb
    }

    command <<<
        set -euo pipefail

        python3 -m venv .venv
        .venv/bin/pip install duckdb pandas numpy --quiet

        .venv/bin/python <<'PYEOF'
        import duckdb
        import pandas as pd
        import numpy as np

        p_thresh  = ~{p_threshold}
        min_maf   = ~{min_maf}
        window_kb = ~{window_kb}
        log10p_thresh = -np.log10(p_thresh)

        con = duckdb.connect()
        sig = con.execute(
            """
            SELECT
                phenotype,
                CHROM,
                CAST(POS AS BIGINT) AS POS,
                ID,
                CAST(A1FREQ AS DOUBLE) AS A1FREQ,
                CAST(LOG10P AS DOUBLE) AS LOG10P
            FROM read_parquet(?)
            WHERE LOG10P >= ?
              AND least(CAST(A1FREQ AS DOUBLE), 1.0 - CAST(A1FREQ AS DOUBLE)) >= ?
            ORDER BY phenotype, LOG10P DESC, POS ASC
            """,
            ["~{summary_stats_parquet}", float(log10p_thresh), float(min_maf)],
        ).df()

        sig["CHROM"] = (
            sig["CHROM"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.upper()
            .replace({"23": "X"})
        )
        sig["POS"] = sig["POS"].astype(int)

        loci = []
        for pheno, pheno_df in sig.groupby("phenotype", sort=True):
            used = set()
            for _, row in pheno_df.iterrows():
                snp_key = (row["CHROM"], int(row["POS"]))
                if snp_key in used:
                    continue

                chrom  = row["CHROM"]
                pos    = int(row["POS"])
                win_st = max(1, pos - window_kb * 1000)
                win_en = pos + window_kb * 1000
                snp_id = row["ID"] if pd.notna(row["ID"]) else f"{chrom}:{pos}"

                loci.append({
                    "phenotype": pheno,
                    "chrom": chrom,
                    "pos": pos,
                    "index_snp": snp_id,
                    "window_start": win_st,
                    "window_end": win_en,
                })

                in_window = (
                    (pheno_df["CHROM"] == chrom)
                    & (pheno_df["POS"] >= win_st)
                    & (pheno_df["POS"] <= win_en)
                )
                for hit in pheno_df.loc[in_window, ["CHROM", "POS"]].itertuples(index=False):
                    used.add((hit.CHROM, int(hit.POS)))

        loci_df = pd.DataFrame(
            loci,
            columns=[
                "phenotype",
                "chrom",
                "pos",
                "index_snp",
                "window_start",
                "window_end",
            ],
        )

        manifest = pd.read_csv("~{step2_chunk_manifest}", sep="\t", dtype=str)
        required_cols = [
            "chrom",
            "plink2_chrom_pvar",
            "plink2_chrom_psam",
            "plink2_chrom_pgen",
        ]
        has_required_cols = set(required_cols).issubset(manifest.columns)

        if not has_required_cols:
            probe = pd.read_csv("~{step2_chunk_manifest}", sep="\t", header=None, dtype=str)

            if probe.shape[1] == 10:
                manifest = pd.read_csv(
                    "~{step2_chunk_manifest}",
                    sep="\t",
                    header=None,
                    names=[
                        "chrom",
                        "start",
                        "end",
                        "range",
                        "plink2_chrom_pvar_name",
                        "plink2_chrom_psam_name",
                        "plink2_chrom_pgen_name",
                        "plink2_chrom_pvar",
                        "plink2_chrom_psam",
                        "plink2_chrom_pgen",
                    ],
                    dtype=str,
                )
            elif probe.shape[1] == 7:
                manifest = pd.read_csv(
                    "~{step2_chunk_manifest}",
                    sep="\t",
                    header=None,
                    names=[
                        "chrom",
                        "plink2_chrom_pvar_name",
                        "plink2_chrom_psam_name",
                        "plink2_chrom_pgen_name",
                        "plink2_chrom_pvar",
                        "plink2_chrom_psam",
                        "plink2_chrom_pgen",
                    ],
                    dtype=str,
                )
            else:
                raise ValueError(
                    "step2_chunk_manifest must be either a 10-column chunk manifest "
                    "or a 7-column chromosome manifest"
                )

        if not manifest.empty and str(manifest.iloc[0]["chrom"]).strip().lower() == "chrom":
            manifest = manifest.iloc[1:].copy()

        missing_cols = [c for c in required_cols if c not in manifest.columns]
        if missing_cols:
            raise ValueError(f"step2_chunk_manifest is missing columns: {missing_cols}")

        manifest = manifest.copy()
        manifest["chrom"] = (
            manifest["chrom"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.upper()
            .replace({"23": "X"})
        )
        has_windows = {"start", "end"}.issubset(manifest.columns)
        if has_windows:
            manifest["start"] = manifest["start"].astype(int)
            manifest["end"] = manifest["end"].astype(int)

        if loci_df.empty:
            merged = pd.DataFrame(columns=[
                "phenotype",
                "chrom",
                "pos",
                "index_snp",
                "window_start",
                "window_end",
                "pvar_dx_path",
                "psam_dx_path",
                "pgen_dx_path",
            ])
        else:
            rows = []
            for locus in loci_df.to_dict(orient="records"):
                if has_windows:
                    matches = manifest[
                        (manifest["chrom"] == locus["chrom"])
                        & (manifest["start"] <= locus["pos"])
                        & (manifest["end"] >= locus["pos"])
                    ]
                else:
                    matches = manifest[manifest["chrom"] == locus["chrom"]]

                if matches.empty:
                    raise ValueError(
                        f"No step2 manifest row covers {locus['chrom']}:{locus['pos']} "
                        f"for phenotype {locus['phenotype']}"
                    )
                if len(matches) > 1:
                    raise ValueError(
                        f"Multiple step2 manifest rows cover {locus['chrom']}:{locus['pos']} "
                        f"for phenotype {locus['phenotype']}"
                    )

                match = matches.iloc[0]
                locus["pvar_dx_path"] = match["plink2_chrom_pvar"]
                locus["psam_dx_path"] = match["plink2_chrom_psam"]
                locus["pgen_dx_path"] = match["plink2_chrom_pgen"]
                rows.append(locus)

            merged = pd.DataFrame(rows)

        merged.to_csv("merged_loci_manifest.tsv", sep="\t", index=False)
        print(f"Total loci to fine-map: {len(merged)}")
        PYEOF
    >>>

    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        File loci_manifest = "merged_loci_manifest.tsv"
    }
}


task compute_ld {

    input {
        File    pvar
        File    psam
        File    pgen
        String  plink2_prefix
        File    plink2_binary
        String  chrom
        String  start
        String  end
        File    sample_keep
        String  out_prefix
        Int     threads
    }

    command <<<
        set -euo pipefail

        echo "$(date -u) Computing LD for ~{chrom}:~{start}-~{end}"

        ln -sf ~{psam} ~{plink2_prefix}.psam
        ln -sf ~{pvar} ~{plink2_prefix}.pvar
        ln -sf ~{pgen} ~{plink2_prefix}.pgen

        apt-get update -qq && apt-get install -y unzip --quiet
        unzip ~{plink2_binary}

        ## Extract variants in the locus window
        ./plink2 \
            --pfile ~{plink2_prefix} \
            --no-pheno \
            --chr ~{chrom} \
            --from-bp ~{start} \
            --to-bp ~{end} \
            --keep ~{sample_keep} \
            --maf 0.001 \
            --make-pgen \
            --threads ~{threads} \
            --out ~{out_prefix}.extracted

        ## Count samples
        tail -n +2 ~{out_prefix}.extracted.psam | wc -l > n_samples.txt

        ## Compute full symmetric LD matrix (r, not r^2) expected by susie_rss()
        ./plink2 \
            --pfile ~{out_prefix}.extracted \
            --no-pheno \
            --r-unphased square \
            --threads ~{threads} \
            --out ~{out_prefix}

        ## Ordered variant IDs matching LD matrix rows/cols
        awk 'NR>1 {print $3}' ~{out_prefix}.extracted.pvar > ~{out_prefix}.variant_ids.txt

        echo "$(date -u) Done computing LD"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
        dx_instance_type: "mem2_ssd1_v2_x8"
    }

    parameter_meta {
        pgen: { description: "input pgen for LD computation", patterns: ["*.pgen"], stream: true }
        pvar: { description: "input pvar",                    patterns: ["*.pvar"], stream: true }
        psam: { description: "input psam",                    patterns: ["*.psam"], stream: true }
    }

    output {
        File   ld_matrix   = "~{out_prefix}.unphased.vcor"
        File   variant_ids = "~{out_prefix}.variant_ids.txt"
        String n_samples   = read_string("n_samples.txt")
    }
}


task prepare_locus_sumstats {

    input {
        File    summary_stats_parquet
        String  phenotype
        String  chrom
        String  start
        String  end
        File    variant_ids
    }

    command <<<
        set -euo pipefail

        python3 -m venv .venv
        .venv/bin/pip install duckdb pandas numpy --quiet

        .venv/bin/python <<'PYEOF'
        import duckdb
        import pandas as pd
        import numpy as np

        pheno  = "~{phenotype}"
        chrom  = "~{chrom}".upper()
        win_st = int("~{start}")
        win_en = int("~{end}")

        ld_snps = pd.read_csv("~{variant_ids}", header=None, names=["ID"])

        con = duckdb.connect()
        window = con.execute(
            """
            SELECT
                ID,
                CHROM,
                CAST(POS AS BIGINT) AS POS,
                ALLELE0,
                ALLELE1,
                CAST(BETA AS DOUBLE) AS BETA,
                CAST(SE AS DOUBLE) AS SE,
                CAST(A1FREQ AS DOUBLE) AS A1FREQ,
                CAST(N AS DOUBLE) AS N,
                CAST(LOG10P AS DOUBLE) AS LOG10P
            FROM read_parquet(?)
            WHERE phenotype = ?
              AND replace(upper(CHROM), 'CHR', '') IN (?, ?)
              AND CAST(POS AS BIGINT) BETWEEN ? AND ?
            """,
            ["~{summary_stats_parquet}", pheno, chrom, "23" if chrom == "X" else chrom, win_st, win_en],
        ).df()

        window["CHROM"] = (
            window["CHROM"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.upper()
            .replace({"23": "X"})
        )
        window["POS"]   = window["POS"].astype(int)
        window["Z"] = window["BETA"].astype(float) / window["SE"].astype(float)
        window["P"] = 10 ** (-window["LOG10P"].astype(float))

        merged = ld_snps.merge(
            window[["ID", "CHROM", "POS", "ALLELE0", "ALLELE1",
                    "BETA", "SE", "Z", "A1FREQ", "N", "LOG10P", "P"]],
            on="ID",
            how="left"
        )

        n_missing = merged["Z"].isna().sum()
        if n_missing > 0:
            print(f"WARNING: {n_missing} LD variants have no z-score (set to 0)")
            merged["Z"] = merged["Z"].fillna(0)

        merged.to_csv("zscore_~{phenotype}_~{chrom}_~{start}_~{end}.tsv", sep="\t", index=False)
        print(f"Prepared {len(merged)} variants for SuSiE")
        PYEOF
    >>>

    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        File z_scores_tsv = "zscore_~{phenotype}_~{chrom}_~{start}_~{end}.tsv"
    }
}


task run_susie {

    input {
        File    z_scores
        File    ld_matrix
        File    variant_ids
        String  phenotype
        String  chrom
        String  start
        String  end
        String  index_snp
        String  n_samples
        Int     max_causal
        Float   min_abs_corr
        Float   coverage
        Int     max_iter
        Float   tol
        Boolean refine
        Boolean estimate_residual_variance
        String  docker
    }

    command <<<
        set -euo pipefail

        Rscript - <<'REOF'

        if (!requireNamespace("susieR",     quietly = TRUE))
            install.packages("susieR",     repos = "https://cloud.r-project.org", quiet = TRUE)
        if (!requireNamespace("data.table", quietly = TRUE))
            install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)

        library(susieR)
        library(data.table)

        phenotype   <- "~{phenotype}"
        chrom       <- "~{chrom}"
        locus_start <- "~{start}"
        locus_end   <- "~{end}"
        index_snp   <- "~{index_snp}"
        n           <- as.integer("~{n_samples}")
        L           <- as.integer("~{max_causal}")
        min_abs_corr          <- ~{min_abs_corr}
        coverage              <- ~{coverage}
        max_iter              <- as.integer("~{max_iter}")
        tol                   <- ~{tol}
        refine                <- "~{refine}" == "true"
        estimate_residual_var <- "~{estimate_residual_variance}" == "true"
        out_prefix <- paste0("susie_", phenotype, "_", chrom, "_", locus_start, "_", locus_end)

        ## Load z-scores
        cat(format(Sys.time()), "Loading z-scores\n")
        zdf  <- fread("~{z_scores}")
        zvec <- setNames(as.numeric(zdf$Z), zdf$ID)

        ## Load LD matrix
        ## plink2 --r-unphased square outputs a full symmetric r matrix.
        cat(format(Sys.time()), "Loading LD matrix\n")
        R_mat <- as.matrix(fread("~{ld_matrix}", header = FALSE))
        R_mat <- (R_mat + t(R_mat)) / 2   # enforce exact symmetry
        diag(R_mat) <- 1.0

        snp_ids <- readLines("~{variant_ids}")
        rownames(R_mat) <- snp_ids
        colnames(R_mat) <- snp_ids
        stopifnot(length(zvec) == nrow(R_mat))

        ## Check LD / z-score consistency (Zou et al. 2022)
        cat(format(Sys.time()), "Checking z-score / LD consistency\n")
        tryCatch({
            check        <- kriging_rss(z = zvec, R = R_mat, n = n)
            suspect_snps <- rownames(check$conditional_dist)[
                abs(check$conditional_dist[, "z"]) > 2
            ]
            if (length(suspect_snps) > 0)
                warning(paste("Potentially inconsistent SNPs:", paste(suspect_snps, collapse = ", ")))
            writeLines(suspect_snps, paste0(out_prefix, ".suspect_snps.txt"))
        }, error = function(e) {
            cat("kriging_rss check failed (non-fatal):", conditionMessage(e), "\n")
            writeLines("check_failed", paste0(out_prefix, ".suspect_snps.txt"))
        })

        ## Run SuSiE-RSS
        cat(format(Sys.time()), "Running susie_rss\n")
        susie_fit <- tryCatch(
            susie_rss(
                z                          = zvec,
                R                          = R_mat,
                n                          = n,
                L                          = L,
                coverage                   = coverage,
                min_abs_corr               = min_abs_corr,
                max_iter                   = max_iter,
                tol                        = tol,
                refine                     = refine,
                estimate_residual_variance = estimate_residual_var,
                verbose                    = TRUE
            ),
            error = function(e) { cat("susie_rss failed:", conditionMessage(e), "\n"); NULL }
        )

        saveRDS(susie_fit, file = paste0(out_prefix, ".rds"))

        ## Extract PIPs
        pip_vec <- if (!is.null(susie_fit)) susie_fit$pip else rep(NA_real_, length(zvec))

        pip_df <- data.frame(
            phenotype   = phenotype,
            chrom       = chrom,
            locus_start = locus_start,
            locus_end   = locus_end,
            index_snp   = index_snp,
            ID          = names(zvec),
            Z           = as.numeric(zvec),
            PIP         = pip_vec
        )
        fwrite(pip_df, paste0(out_prefix, ".pip.tsv"), sep = "\t")

        ## Extract credible sets
        cs_rows <- list()
        if (!is.null(susie_fit) && !is.null(susie_fit$sets$cs)) {
            for (cs_name in names(susie_fit$sets$cs)) {
                idx <- susie_fit$sets$cs[[cs_name]]
                cs_rows[[cs_name]] <- data.frame(
                    phenotype    = phenotype,
                    chrom        = chrom,
                    locus_start  = locus_start,
                    locus_end    = locus_end,
                    index_snp    = index_snp,
                    credible_set = cs_name,
                    cs_size      = length(idx),
                    min_abs_corr = susie_fit$sets$purity[cs_name, "min.abs.corr"],
                    ID           = names(zvec)[idx],
                    Z            = zvec[idx],
                    PIP          = pip_vec[idx]
                )
            }
            cs_df <- rbindlist(cs_rows, fill = TRUE)
        } else {
            cs_df <- data.frame(
                phenotype = character(), chrom = character(),
                locus_start = character(), locus_end = character(),
                index_snp = character(), credible_set = character(),
                cs_size = integer(), min_abs_corr = numeric(),
                ID = character(), Z = numeric(), PIP = numeric()
            )
        }
        fwrite(cs_df, paste0(out_prefix, ".credible_sets.tsv"), sep = "\t")

        cat(format(Sys.time()), "Done\n")
        REOF
    >>>

    runtime {
        docker: "~{docker}"
        dx_instance_type: "mem2_ssd1_v2_x4"
    }

    output {
        File pip_tsv           = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.pip.tsv"
        File credible_sets_tsv = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.credible_sets.tsv"
        File susie_rds         = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.rds"
        File log               = stdout()
    }
}


task plot_locus {

    input {
        File   z_scores
        File   pip_tsv
        File   credible_sets_tsv
        File   ld_matrix
        File   variant_ids
        String phenotype
        String chrom
        String start
        String end
        String index_snp
        Float  p_threshold
        String docker
    }

    command <<<
        set -euo pipefail

        Rscript - <<'REOF'
        if (!requireNamespace("data.table", quietly = TRUE))
            install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
        if (!requireNamespace("ggplot2", quietly = TRUE))
            install.packages("ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE)
        if (!requireNamespace("cowplot", quietly = TRUE))
            install.packages("cowplot", repos = "https://cloud.r-project.org", quiet = TRUE)

        library(data.table)
        library(ggplot2)
        library(cowplot)

        phenotype <- "~{phenotype}"
        chrom <- "~{chrom}"
        locus_chr <- if (grepl("^chr", chrom, ignore.case = TRUE)) chrom else paste0("chr", chrom)
        locus_start <- as.numeric("~{start}")
        locus_end <- as.numeric("~{end}")
        index_snp <- "~{index_snp}"
        p_threshold <- ~{p_threshold}
        out_prefix <- paste0("locuszoom_", phenotype, "_", chrom, "_", locus_start, "_", locus_end)

        zdf <- fread("~{z_scores}")
        pip_df <- fread("~{pip_tsv}")

        cs_df <- if (file.exists("~{credible_sets_tsv}") && file.info("~{credible_sets_tsv}")$size > 0) {
            fread("~{credible_sets_tsv}")
        } else {
            data.table(ID = character(), credible_set = character())
        }

        variant_ids <- readLines("~{variant_ids}")
        R_mat <- as.matrix(fread("~{ld_matrix}", header = FALSE))
        rownames(R_mat) <- variant_ids
        colnames(R_mat) <- variant_ids

        lead_id <- NA_character_
        if (index_snp %in% variant_ids) {
            lead_id <- index_snp
        } else if (nrow(pip_df[!is.na(PIP)]) > 0) {
            lead_id <- pip_df[order(-PIP)][["ID"]][1]
        } else if (nrow(zdf[!is.na(P)]) > 0) {
            lead_id <- zdf[order(P)][["ID"]][1]
        }

        r2_df <- data.table(ID = variant_ids, r2 = NA_real_)
        if (!is.na(lead_id) && lead_id %in% variant_ids) {
            r2_df[, r2 := as.numeric(R_mat[lead_id, ])^2]
        }

        cs_annot <- cs_df[, .(
            credible_sets = paste(unique(credible_set), collapse = ","),
            in_credible_set = .N > 0
        ), by = ID]

        plot_df <- merge(zdf, pip_df[, .(ID, PIP)], by = "ID", all.x = TRUE)
        plot_df <- merge(plot_df, cs_annot, by = "ID", all.x = TRUE)
        plot_df <- merge(plot_df, r2_df, by = "ID", all.x = TRUE)

        plot_df[is.na(PIP), PIP := 0]
        plot_df[is.na(in_credible_set), in_credible_set := FALSE]
        plot_df[is.na(credible_sets), credible_sets := ""]
        plot_df[, neg_log10_p := -log10(pmax(P, 1e-300))]

        assign_gene_tracks <- function(gene_df) {
            if (nrow(gene_df) == 0) return(gene_df)

            gene_df <- gene_df[order(start, end)]
            tracks <- integer(nrow(gene_df))
            tracks[1] <- 1L
            buffer <- 1000

            if (nrow(gene_df) > 1) {
                for (i in 2:nrow(gene_df)) {
                    for (track in seq_len(i)) {
                        overlaps <- FALSE
                        for (j in seq_len(i - 1)) {
                            if (tracks[j] == track && gene_df$start[i] <= gene_df$end[j] + buffer) {
                                overlaps <- TRUE
                                break
                            }
                        }
                        if (!overlaps) {
                            tracks[i] <- track
                            break
                        }
                    }
                }
            }

            gene_df[, track := tracks]
            gene_df
        }

        gencode_url <- "https://raw.githubusercontent.com/weinstocklab/gwasplot/main/data/gencode.rda"
        utils::download.file(gencode_url, destfile = "gencode.rda", mode = "wb", quiet = TRUE)
        load("gencode.rda")
        gencode_dt <- as.data.table(gencode)

        genes <- gencode_dt[
            type == "gene" &
            gene_biotype == "protein_coding" &
            chrom == locus_chr &
            start <= locus_end &
            end >= locus_start,
            .(gene_id, gene_name, start, end)
        ]

        exons <- gencode_dt[
            type == "CDS" &
            chrom == locus_chr &
            start <= locus_end &
            end >= locus_start,
            .(gene_id, start, end)
        ]

        plot_df[, ld_bin := fifelse(
            is.na(r2), "Unknown",
            fifelse(r2 >= 0.8, "0.8-1.0",
            fifelse(r2 >= 0.6, "0.6-0.8",
            fifelse(r2 >= 0.4, "0.4-0.6",
            fifelse(r2 >= 0.2, "0.2-0.4",
            fifelse(r2 > 0, "0-0.2", "0"))))))
        ]
        plot_df[, ld_bin := factor(ld_bin, levels = c("Unknown", "0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))]

        label_df <- unique(rbindlist(list(
            plot_df[ID == lead_id],
            plot_df[order(-PIP)][1:min(5, .N)],
            plot_df[in_credible_set == TRUE][1:min(5, .N)]
        ), fill = TRUE), by = "ID")

        p_assoc <- ggplot(plot_df, aes(x = POS, y = neg_log10_p)) +
            geom_point(aes(color = ld_bin, size = pmax(PIP, 0.02), shape = in_credible_set), alpha = 0.85) +
            geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "firebrick") +
            geom_text(
                data = label_df,
                aes(label = ID),
                size = 2.5,
                vjust = -0.7,
                check_overlap = TRUE,
                show.legend = FALSE
            ) +
            scale_color_manual(
                values = c(
                    "Unknown" = "grey70",
                    "0" = "steelblue4",
                    "0-0.2" = "deepskyblue3",
                    "0.2-0.4" = "darkolivegreen3",
                    "0.4-0.6" = "goldenrod2",
                    "0.6-0.8" = "darkorange2",
                    "0.8-1.0" = "firebrick3"
                ),
                drop = FALSE
            ) +
            scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
            scale_size_continuous(range = c(1.2, 5), name = "SuSiE PIP") +
            scale_x_continuous(
                limits = c(locus_start, locus_end),
                expand = c(0, 0)
            ) +
            labs(
                title = paste0("Locus plot: ", phenotype, " chr", chrom, ":", locus_start, "-", locus_end),
                subtitle = paste0("Lead/index SNP: ", lead_id),
                x = NULL,
                y = expression(-log[10](P)),
                color = expression(LD~(r^2)),
                shape = "Credible set"
            ) +
            theme_bw(base_size = 11) +
            theme(
                plot.title = element_text(face = "bold"),
                legend.position = "right",
                panel.grid.minor = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                plot.margin = margin(t = 5, r = 10, b = 0, l = 5, unit = "pt")
            )

        if (nrow(genes) > 0) {
            genes <- assign_gene_tracks(copy(genes))
            max_track <- max(genes$track)
            gene_body_height <- 0.15
            exon_height <- 0.4
            track_spacing <- 0.6

            genes[, `:=`(
                y_center = (track - 1) * track_spacing,
                y_min = y_center - gene_body_height / 2,
                y_max = y_center + gene_body_height / 2,
                label_y = y_center + exon_height / 2 + 0.1
            )]

            if (nrow(exons) > 0) {
                exons <- merge(exons, genes[, .(gene_id, track, y_center)], by = "gene_id", all.x = TRUE)
                exons <- exons[!is.na(track)]
                exons[, `:=`(
                    y_min = y_center - exon_height / 2,
                    y_max = y_center + exon_height / 2
                )]
            }

            p_genes <- ggplot() +
                geom_rect(
                    data = genes,
                    aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
                    fill = "#cccccc",
                    color = "#999999",
                    alpha = 0.7
                )

            if (nrow(exons) > 0) {
                p_genes <- p_genes +
                    geom_rect(
                        data = exons,
                        aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max),
                        fill = "#2c3e50",
                        color = "#34495e",
                        alpha = 0.9
                    )
            }

            p_genes <- p_genes +
                geom_text(
                    data = genes,
                    aes(x = (start + end) / 2, y = label_y, label = gene_name),
                    size = 3.8,
                    hjust = 0.5,
                    vjust = 0,
                    fontface = "italic",
                    color = "#2c3e50"
                ) +
                scale_y_continuous(
                    limits = c(-0.2, max_track * track_spacing + 0.1),
                    expand = expansion(mult = c(0, 0))
                ) +
                scale_x_continuous(
                    limits = c(locus_start, locus_end),
                    expand = c(0, 0)
                ) +
                labs(x = paste0("Position on chr", chrom)) +
                theme_void() +
                theme(
                    axis.title.x = element_text(size = 10, margin = margin(t = 4)),
                    plot.margin = margin(t = 0, r = 10, b = 5, l = 5, unit = "pt")
                )
        } else {
            p_genes <- ggplot() +
                geom_blank() +
                scale_x_continuous(
                    limits = c(locus_start, locus_end),
                    expand = c(0, 0)
                ) +
                labs(x = paste0("Position on chr", chrom)) +
                theme_void() +
                theme(
                    axis.title.x = element_text(size = 10, margin = margin(t = 4)),
                    plot.margin = margin(t = 0, r = 10, b = 5, l = 5, unit = "pt")
                )
        }

        p <- cowplot::plot_grid(p_assoc, p_genes, ncol = 1, align = "v", rel_heights = c(5, 1.2))

        ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = 7)
        ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 7, dpi = 200)
        fwrite(plot_df, paste0(out_prefix, ".plot_data.tsv"), sep = "\t")
        REOF
    >>>

    runtime {
        docker: "~{docker}"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }

    output {
        File plot_pdf      = "locuszoom_~{phenotype}_~{chrom}_~{start}_~{end}.pdf"
        File plot_png      = "locuszoom_~{phenotype}_~{chrom}_~{start}_~{end}.png"
        File plot_data_tsv = "locuszoom_~{phenotype}_~{chrom}_~{start}_~{end}.plot_data.tsv"
    }
}


task combine_results {

    input {
        Array[File] pip_tsvs
        Array[File] credible_set_tsvs
    }

    command <<<
        set -euo pipefail

        python3 -m venv .venv
        .venv/bin/pip install pandas --quiet

        .venv/bin/python <<'PYEOF'
        import pandas as pd

        pip_files = ["~{sep='", "' pip_tsvs}"]
        cs_files  = ["~{sep='", "' credible_set_tsvs}"]

        merged_pip = pd.concat(
            [pd.read_csv(f, sep="\t") for f in pip_files],
            ignore_index=True
        ).sort_values(["phenotype", "chrom", "locus_start", "PIP"],
                      ascending=[True, True, True, False])
        merged_pip.to_csv("all_loci_pip.tsv", sep="\t", index=False)
        print(f"Total PIPs written: {len(merged_pip)}")

        cs_dfs    = [d for d in (pd.read_csv(f, sep="\t") for f in cs_files) if len(d) > 0]
        merged_cs = pd.concat(cs_dfs, ignore_index=True) if cs_dfs else pd.DataFrame()
        merged_cs.to_csv("all_loci_credible_sets.tsv", sep="\t", index=False)
        print(f"Total credible set entries written: {len(merged_cs)}")
        PYEOF
    >>>

    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem2_ssd1_v2_x2"
    }

    output {
        File merged_pip           = "all_loci_pip.tsv"
        File merged_credible_sets = "all_loci_credible_sets.tsv"
    }
}
