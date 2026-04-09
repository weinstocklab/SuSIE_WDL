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
        Float           max_geno_missing = 0.05  # mirror GWAS-side PLINK missingness filter
        Boolean         pass_only_var_filter = true
        String?         plink_var_filter          # optional PLINK --var-filter exception names
        Int             window_kb       = 1000     # +/- kb around index SNP
        Int             max_causal      = 10      # L parameter in SuSiE
        Float           min_abs_corr    = 0.5     # min |r| for credible sets
        Float           coverage        = 0.95    # credible set coverage
        Int             max_iter        = 70      # SuSiE max iterations
        Float           tol             = 0.002    # SuSiE convergence tol
        Boolean         refine          = true    # use SuSiE refinement
        Boolean         estimate_residual_variance = true

        ## Plink2 binary 
        File            plink2_binary

        ## Runtime 
        Int             threads     = 4
        String          r_docker    = "jweinstk/fine_mapping:4.5.3"
        Int             variant_count_threshold = 4000  # threshold for dynamic instance selection (>= threshold = large)
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
                threads       = threads,
                max_geno_missing = max_geno_missing,
                min_maf       = min_maf,
                pass_only_var_filter = pass_only_var_filter,
                plink_var_filter = plink_var_filter
        }

        ## Determine instance type based on variant count
        call determine_instance_type {
            input:
                variant_ids = compute_ld.variant_ids,
                threshold = variant_count_threshold
        }

        ## 3. Prepare per-locus summary stats
        call prepare_locus_sumstats {
            input:
                summary_stats_parquet = summary_stats_parquet,
                phenotype             = locus_phenotype,
                chrom                 = locus_chrom,
                start                 = locus_start,
                end                   = locus_end,
                variant_ids           = compute_ld.variant_ids,
                ld_variants_pvar      = compute_ld.extracted_pvar,
                min_maf               = min_maf
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
                min_maf                    = min_maf,
                refine                     = refine,
                estimate_residual_variance = estimate_residual_variance,
                docker                     = r_docker,
                instance_type              = determine_instance_type.instance_type
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

    call check_indel_pass {
        input:
            loci_manifest = identify_loci.loci_manifest,
            pvar_files = compute_ld.extracted_pvar
    }

    output {
        File          merged_pip           = combine_results.merged_pip
        File          merged_credible_sets = combine_results.merged_credible_sets
        Array[File]   susie_rds_all        = run_susie.susie_rds
        Array[File]   susie_logs           = run_susie.log
        Array[File]   susie_excluded_snps  = run_susie.excluded_snps_tsv
        Array[File]   susie_diagnostics    = run_susie.diagnostics_txt
        File          indel_pass_summary   = check_indel_pass.indel_summary
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

        def clean_series(series):
            out = series.fillna("").astype(str).str.strip()
            out = out[out != ""]
            return out

        iid_series = clean_series(df[normalized["IID"]]).drop_duplicates() if "IID" in normalized else pd.Series(dtype=str)
        eid_series = clean_series(df[normalized["EID"]]).drop_duplicates() if "EID" in normalized else pd.Series(dtype=str)

        if "FID" in normalized and "IID" in normalized:
            fid_iid = pd.DataFrame({
                "FID": clean_series(df[normalized["FID"]]),
                "IID": clean_series(df[normalized["IID"]]),
            })
            fid_iid = fid_iid[(fid_iid["FID"] != "") & (fid_iid["IID"] != "")]
            fid_iid = fid_iid.drop_duplicates()
            if fid_iid.empty:
                raise ValueError("covariates_file has FID/IID headers but no non-empty rows")
            fid_iid.to_csv("sample_keep.tsv", sep="\t", index=False, header=False)
            strategy = "FID+IID"
        elif not iid_series.empty:
            iid_series.to_frame(name="IID").to_csv("sample_keep.tsv", sep="\t", index=False, header=False)
            strategy = "IID"
        elif not eid_series.empty:
            eid_series.to_frame(name="IID").to_csv("sample_keep.tsv", sep="\t", index=False, header=False)
            strategy = "EID"
        else:
            raise ValueError("covariates_file must contain at least one non-empty IID or EID value")

        print(f"Derived PLINK keep file using strategy: {strategy}")
        print(f"Detected IID column: {'IID' in normalized}")
        print(f"Detected EID column: {'EID' in normalized}")
        print(f"Detected FID column: {'FID' in normalized}")
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

        # Merge overlapping intervals within each phenotype and chromosome
        if not loci_df.empty:
            merged_loci = []
            for pheno in loci_df["phenotype"].unique():
                pheno_df = loci_df[loci_df["phenotype"] == pheno]
                for chrom in pheno_df["chrom"].unique():
                    chrom_loci = pheno_df[pheno_df["chrom"] == chrom].sort_values("window_start").copy()
                    
                    if len(chrom_loci) == 0:
                        continue
                    
                    # Merge overlapping intervals
                    intervals = []
                    current = chrom_loci.iloc[0].to_dict()
                    
                    for idx in range(1, len(chrom_loci)):
                        row = chrom_loci.iloc[idx]
                        # Check if current interval overlaps with next interval
                        if row["window_start"] <= current["window_end"]:
                            # Merge: extend the end, keep the most significant (first) index_snp
                            current["window_end"] = max(current["window_end"], row["window_end"])
                        else:
                            # No overlap, save current and start new interval
                            intervals.append(current)
                            current = row.to_dict()
                    
                    # Don't forget the last interval
                    intervals.append(current)
                    merged_loci.extend(intervals)
            
            loci_df = pd.DataFrame(merged_loci)

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
        Float   max_geno_missing
        Float   min_maf
        Boolean pass_only_var_filter
        String? plink_var_filter
        String  instance_type = "mem2_ssd1_v2_x4"  # small by default
    }

    command <<<
        set -euo pipefail

        echo "$(date -u) Computing LD for ~{chrom}:~{start}-~{end}"
        printf "not_started\n" > ld_status.txt

        ln -sf ~{psam} ~{plink2_prefix}.psam
        ln -sf ~{pvar} ~{plink2_prefix}.pvar
        ln -sf ~{pgen} ~{plink2_prefix}.pgen

        apt-get update -qq && apt-get install -y unzip --quiet
        unzip ~{plink2_binary}

        var_filter_args=()
        if [[ -n "~{default="" plink_var_filter}" ]]; then
            var_filter_args=(--var-filter "~{default="" plink_var_filter}")
            echo "$(date -u) Applying PLINK var filter with exceptions: ~{default="" plink_var_filter}"
        elif [[ "~{pass_only_var_filter}" == "true" ]]; then
            var_filter_args=(--var-filter)
            echo "$(date -u) Applying PLINK PASS-only var filter"
        fi

        ## Extract variants in the locus window
        ./plink2 \
            --pfile ~{plink2_prefix} \
            --no-pheno \
            --chr ~{chrom} \
            --from-bp ~{start} \
            --to-bp ~{end} \
            --keep ~{sample_keep} \
            --maf ~{min_maf} \
            --geno ~{max_geno_missing} \
            "${var_filter_args[@]}" \
            --make-pgen \
            --memory 12000 \
            --threads ~{threads} \
            --out ~{out_prefix}.extracted

        sample_count=$(tail -n +2 ~{out_prefix}.extracted.psam | wc -l)
        if [[ "${sample_count}" -eq 0 ]]; then
            echo "No samples remaining after sample/MAF/missingness/window filters" >&2
            exit 13
        fi
        printf "%s\n" "${sample_count}" > n_samples.txt

        ## Ordered variant IDs matching LD matrix rows/cols
        awk 'NR>1 {print $3}' ~{out_prefix}.extracted.pvar > ~{out_prefix}.variant_ids.txt

        variant_count=$(wc -l < ~{out_prefix}.variant_ids.txt)
        echo "$(date -u) Variants retained after filtering: ${variant_count}"

        if [[ "${variant_count}" -eq 0 ]]; then
            echo "No variants remaining in locus after sample/MAF/missingness/variant-filter/window filters" >&2
            exit 14
        elif [[ "${variant_count}" -eq 1 ]]; then
            printf "1\n" > ~{out_prefix}.unphased.vcor
            printf "single_variant_identity\n" > ld_status.txt
        else
            ## Compute full symmetric LD matrix (r, not r^2) expected by susie_rss()
            ./plink2 \
                --pfile ~{out_prefix}.extracted \
                --no-pheno \
                --r-unphased square allow-ambiguous-allele \
                --memory 12000 \
                --threads ~{threads} \
                --out ~{out_prefix}

            if [[ ! -e "~{out_prefix}.unphased.vcor" && -e "~{out_prefix}.unphased.vcor1" ]]; then
                mv "~{out_prefix}.unphased.vcor1" "~{out_prefix}.unphased.vcor"
            fi

            if [[ -s "~{out_prefix}.unphased.vcor" ]]; then
                printf "plink_generated\n" > ld_status.txt
            else
                echo "WARNING: PLINK did not emit ~{out_prefix}.unphased.vcor; using identity fallback for ${variant_count} variants" >&2
                awk -v n="${variant_count}" 'BEGIN {
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= n; ++j) {
                            printf("%s%d", (j == 1 ? "" : " "), (i == j ? 1 : 0))
                        }
                        printf("\n")
                    }
                }' > ~{out_prefix}.unphased.vcor
                printf "plink_missing_identity_fallback\n" > ld_status.txt
            fi
        fi

        echo "$(date -u) Done computing LD"
    >>>

    runtime {
        docker: "ghcr.io/rgcgithub/regenie/regenie:v4.1.gz"
        dx_instance_type: "~{instance_type}"
    }

    parameter_meta {
        pgen: { description: "input pgen for LD computation", patterns: ["*.pgen"], stream: true }
        pvar: { description: "input pvar",                    patterns: ["*.pvar"], stream: true }
        psam: { description: "input psam",                    patterns: ["*.psam"], stream: true }
    }

    output {
        File   ld_matrix   = "~{out_prefix}.unphased.vcor"
        File   variant_ids = "~{out_prefix}.variant_ids.txt"
        File   extracted_pvar = "~{out_prefix}.extracted.pvar"
        String n_samples   = read_string("n_samples.txt")
        String ld_status   = read_string("ld_status.txt")
    }
}


task determine_instance_type {

    input {
        File    variant_ids
        Int     threshold
    }

    command <<<
        variant_count=$(wc -l < ~{variant_ids})
        if [[ ${variant_count} -ge ~{threshold} ]]; then
            echo "mem2_ssd1_v2_x8"
        else
            echo "mem2_ssd1_v2_x4"
        fi
    >>>

    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }

    output {
        String instance_type = read_string(stdout())
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
        File    ld_variants_pvar
        Float   min_maf
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
        min_maf = ~{min_maf}

        ld_snps = pd.read_csv("~{variant_ids}", header=None, names=["ID"], dtype=str)
        ld_pvar = pd.read_csv("~{ld_variants_pvar}", sep="\t", dtype=str)
        if "#CHROM" in ld_pvar.columns:
            ld_pvar = ld_pvar.rename(columns={"#CHROM": "CHROM"})
        required_pvar_cols = {"CHROM", "POS", "ID", "REF", "ALT"}
        missing_pvar_cols = required_pvar_cols - set(ld_pvar.columns)
        if missing_pvar_cols:
            raise ValueError(f"ld_variants_pvar is missing columns: {sorted(missing_pvar_cols)}")
        ld_pvar = ld_pvar[["CHROM", "POS", "ID", "REF", "ALT"]].copy()

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
              AND least(CAST(A1FREQ AS DOUBLE), 1.0 - CAST(A1FREQ AS DOUBLE)) >= ?
            """,
            ["~{summary_stats_parquet}", pheno, chrom, "23" if chrom == "X" else chrom, win_st, win_en, min_maf],
        ).df()

        window["CHROM"] = (
            window["CHROM"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.upper()
            .replace({"23": "X"})
        )
        window["POS"] = window["POS"].astype(int)
        window["ID"] = window["ID"].astype(str)
        window["ALLELE0"] = window["ALLELE0"].astype(str).str.upper()
        window["ALLELE1"] = window["ALLELE1"].astype(str).str.upper()
        with np.errstate(divide="ignore", invalid="ignore"):
            window["Z"] = window["BETA"].astype(float) / window["SE"].astype(float)
        window.loc[~np.isfinite(window["Z"]), "Z"] = np.nan
        window["P"] = 10 ** (-window["LOG10P"].astype(float))

        ld_pvar["CHROM"] = (
            ld_pvar["CHROM"]
            .astype(str)
            .str.replace("^chr", "", regex=True)
            .str.upper()
            .replace({"23": "X"})
        )
        ld_pvar["POS"] = ld_pvar["POS"].astype(int)
        ld_pvar["ID"] = ld_pvar["ID"].astype(str)
        ld_pvar["REF"] = ld_pvar["REF"].astype(str).str.upper()
        ld_pvar["ALT"] = ld_pvar["ALT"].astype(str).str.upper()

        def normalize_variant_id(series):
            return (
                series.fillna("")
                .astype(str)
                .str.strip()
                .str.upper()
                .str.replace("^CHR", "", regex=True)
                .str.replace(r"[-_/]", ":", regex=True)
            )

        def unordered_allele_key(a, b):
            a = a.fillna("").astype(str).str.upper()
            b = b.fillna("").astype(str).str.upper()
            left = a.where(a <= b, b)
            right = b.where(a <= b, a)
            return left + ":" + right

        ld_pvar["norm_id"] = normalize_variant_id(ld_pvar["ID"])
        ld_pvar["coord_alleles"] = (
            ld_pvar["CHROM"].astype(str)
            + ":"
            + ld_pvar["POS"].astype(str)
            + ":"
            + unordered_allele_key(ld_pvar["REF"], ld_pvar["ALT"])
        )

        window["norm_id"] = normalize_variant_id(window["ID"])
        window["coord_alleles"] = (
            window["CHROM"].astype(str)
            + ":"
            + window["POS"].astype(str)
            + ":"
            + unordered_allele_key(window["ALLELE0"], window["ALLELE1"])
        )
        print(f"Window summary-stat rows in locus: {len(window)}")

        ss_cols = [
            "ID", "norm_id", "coord_alleles", "CHROM", "POS", "ALLELE0", "ALLELE1",
            "BETA", "SE", "Z", "A1FREQ", "N", "LOG10P", "P"
        ]
        by_id = ld_pvar.merge(
            window[ss_cols],
            on="norm_id",
            how="left",
            suffixes=("", "_ss"),
        )
        n_id_matches = by_id["Z"].notna().sum()

        unmatched = by_id["Z"].isna()
        by_coord = ld_pvar.loc[unmatched].merge(
            window[ss_cols],
            on="coord_alleles",
            how="left",
            suffixes=("", "_ss"),
        )
        if not by_coord.empty:
            by_coord = by_coord.drop(columns=["norm_id_ss"], errors="ignore")
            for col in by_coord.columns:
                by_id.loc[unmatched, col] = by_coord[col].values
        n_coord_matches = by_coord["Z"].notna().sum() if not by_coord.empty else 0

        merged = ld_snps.merge(
            by_id[[
                "ID", "CHROM_ss", "POS_ss", "ALLELE0", "ALLELE1",
                "BETA", "SE", "Z", "A1FREQ", "N", "LOG10P", "P", "REF", "ALT"
            ]].rename(columns={"CHROM_ss": "CHROM", "POS_ss": "POS"}),
            on="ID",
            how="left",
        )

        same_orientation = (merged["ALLELE0"] == merged["REF"]) & (merged["ALLELE1"] == merged["ALT"])
        flipped_orientation = (merged["ALLELE0"] == merged["ALT"]) & (merged["ALLELE1"] == merged["REF"])
        needs_flip = flipped_orientation & ~same_orientation
        allele_mismatch = ~(same_orientation | flipped_orientation) & merged["Z"].notna()
        merged.loc[needs_flip, "Z"] = -merged.loc[needs_flip, "Z"]
        merged.loc[allele_mismatch, "Z"] = np.nan
        merged.loc[allele_mismatch, "BETA"] = np.nan
        merged = merged.drop(columns=["REF", "ALT"])

        n_flipped = needs_flip.sum()
        n_allele_mismatch = allele_mismatch.sum()
        print(
            f"Summary-stat matching: "
            f"id_matches={int(n_id_matches)} "
            f"coord_allele_matches={int(n_coord_matches)} "
            f"flipped={int(n_flipped)} "
            f"allele_mismatch={int(n_allele_mismatch)}"
        )

        n_missing = (~np.isfinite(merged["Z"])).sum()
        if n_missing > 0:
            print(f"WARNING: {n_missing} LD variants have missing/non-finite z-scores")

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
        Float   min_maf
        Boolean refine
        Boolean estimate_residual_variance
        String  docker
        String  instance_type = "mem2_ssd1_v2_x8"  # small by default for SuSiE
    }

    command <<<
        set -euo pipefail

        Rscript - <<'REOF'

        library(susieR)
        library(data.table)
        library(Rfast)

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
        diagnostics_file <- paste0(out_prefix, ".diagnostics.txt")
        log_diag <- function(...) {
            msg <- paste(..., collapse = " ")
            cat(format(Sys.time()), msg, "\n")
            cat(format(Sys.time()), msg, "\n", file = diagnostics_file, append = TRUE)
        }

        ## Load z-scores
        log_diag("Loading z-scores")
        zdf  <- fread("~{z_scores}")
        z_all <- setNames(as.numeric(zdf$Z), zdf$ID)
        log_diag("Loaded", length(z_all), "z-scores")

        ## Load LD matrix
        ## plink2 --r-unphased square outputs a full symmetric r matrix.
        snp_ids <- readLines("~{variant_ids}")
        log_diag("Loading LD matrix for", length(snp_ids), "variant IDs")
        raw_ld <- fread("~{ld_matrix}", header = FALSE, data.table = FALSE, fill = TRUE)
        log_diag("Raw LD table dimensions:", nrow(raw_ld), "rows x", ncol(raw_ld), "cols")

        if (ncol(raw_ld) == length(snp_ids) + 1) {
            first_col <- as.character(raw_ld[[1]])
            first_num <- suppressWarnings(as.integer(first_col))
            if (all(!is.na(first_num)) && identical(first_num, seq_len(length(first_num)))) {
                log_diag("Dropping leading numeric index column from LD matrix")
                raw_ld <- raw_ld[, -1, drop = FALSE]
            } else if (length(first_col) == length(snp_ids) && all(first_col == snp_ids)) {
                log_diag("Dropping leading variant ID column from LD matrix")
                raw_ld <- raw_ld[, -1, drop = FALSE]
            }
        }

        if (nrow(raw_ld) != length(snp_ids) || ncol(raw_ld) != length(snp_ids)) {
            log_diag(
                "LD matrix dimension mismatch after preprocessing:",
                "rows =", nrow(raw_ld),
                "cols =", ncol(raw_ld),
                "expected =", length(snp_ids)
            )
            stop(
                sprintf(
                    "LD matrix shape mismatch: got %d x %d for %d variant IDs",
                    nrow(raw_ld), ncol(raw_ld), length(snp_ids)
                )
            )
        }

        R_mat <- as.matrix(data.frame(lapply(raw_ld, as.numeric), check.names = FALSE))
        storage.mode(R_mat) <- "double"
        rownames(R_mat) <- snp_ids
        colnames(R_mat) <- snp_ids
        if (length(z_all) != nrow(R_mat)) {
            log_diag(
                "z-score / LD dimension mismatch:",
                "z =", length(z_all),
                "ld_rows =", nrow(R_mat)
            )
            stop(
                sprintf(
                    "z-score / LD dimension mismatch: %d z-scores vs %d LD rows",
                    length(z_all), nrow(R_mat)
                )
            )
        }

        ## First drop variants without usable z-scores; unmatched GWAS variants
        ## should not be allowed to poison the LD QC for the matched subset.
        bad_z <- !is.finite(z_all)
        keep_finite_z <- !bad_z
        log_diag(
            "Initial z-score counts:",
            "total =", length(z_all),
            "nonfinite_z =", sum(bad_z),
            "finite_z =", sum(keep_finite_z)
        )

        excluded_df <- data.frame(
            ID = names(z_all),
            exclude_nonfinite_z = bad_z,
            exclude_nonfinite_ld = FALSE,
            kept_for_susie = FALSE
        )

        if (!any(keep_finite_z)) {
            fwrite(excluded_df, paste0(out_prefix, ".excluded_snps.tsv"), sep = "\t")
            writeLines(head(names(z_all), 20), paste0(out_prefix, ".bad_z_ids.txt"))
            stop("No SNPs have finite z-scores for SuSiE")
        }

        zvec <- z_all[keep_finite_z]
        R_mat <- R_mat[keep_finite_z, keep_finite_z, drop = FALSE]

        ## Enforce exact symmetry / valid diagonal before LD QC.
        R_mat <- (R_mat + t(R_mat)) / 2
        diag(R_mat) <- 1.0

        ld_oob <- abs(R_mat) > (1 + 1e-8)
        ld_oob[is.na(ld_oob)] <- FALSE
        if (any(ld_oob)) {
            log_diag(
                "Clipping", sum(ld_oob),
                "LD entries outside [-1,1]; max abs before clip =",
                max(abs(R_mat[ld_oob]))
            )
            R_mat[R_mat > 1] <- 1
            R_mat[R_mat < -1] <- -1
        }

        bad_ld_entry <- !is.finite(R_mat)
        bad_ld_snp <- apply(bad_ld_entry, 1, any) | apply(bad_ld_entry, 2, any)
        log_diag(
            "Post-z-filter LD non-finite counts:",
            "subset_total =", length(zvec),
            "ld_entries =", sum(bad_ld_entry),
            "ld_rows =", sum(bad_ld_snp)
        )

        keep_snps <- !bad_ld_snp
        n_bad_z <- sum(bad_z)
        n_bad_ld_snp <- sum(bad_ld_snp)
        n_kept <- sum(keep_snps)
        log_diag(
            "SuSiE SNP filtering:",
            "total =", length(z_all),
            "nonfinite_z =", n_bad_z,
            "nonfinite_ld_rows =", n_bad_ld_snp,
            "kept =", n_kept
        )

        excluded_df$exclude_nonfinite_ld[keep_finite_z] <- bad_ld_snp
        excluded_df$kept_for_susie[keep_finite_z] <- keep_snps
        fwrite(excluded_df, paste0(out_prefix, ".excluded_snps.tsv"), sep = "\t")

        if (n_kept == 0) {
            bad_ld_ids <- names(zvec)[bad_ld_snp]
            bad_z_ids <- names(z_all)[bad_z]
            writeLines(head(bad_ld_ids, 20), paste0(out_prefix, ".bad_ld_ids.txt"))
            writeLines(head(bad_z_ids, 20), paste0(out_prefix, ".bad_z_ids.txt"))
            stop("No SNPs remain for SuSiE after excluding non-finite z-scores / LD rows")
        }

        R_mat <- R_mat[keep_snps, keep_snps, drop = FALSE]
        zvec <- zvec[keep_snps]

        ## Check LD / z-score consistency (Zou et al. 2022)
        log_diag("Checking z-score / LD consistency")
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
        log_diag("Running susie_rss with", length(zvec), "SNPs")
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
            error = function(e) { log_diag("susie_rss failed:", conditionMessage(e)); NULL }
        )

        saveRDS(susie_fit, file = paste0(out_prefix, ".rds"))

        ## Extract PIPs
        pip_vec <- rep(NA_real_, length(z_all))
        names(pip_vec) <- names(z_all)
        if (!is.null(susie_fit)) {
            pip_vec[names(zvec)] <- susie_fit$pip
        }

        pip_df <- data.frame(
            phenotype   = phenotype,
            chrom       = chrom,
            locus_start = locus_start,
            locus_end   = locus_end,
            index_snp   = index_snp,
            ID          = names(z_all),
            Z           = as.numeric(z_all),
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
                    PIP          = pip_vec[names(zvec)[idx]]
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

        log_diag("Done")
        REOF
    >>>

    runtime {
        docker: "~{docker}"
        dx_instance_type: "~{instance_type}"
    }

    output {
        File pip_tsv           = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.pip.tsv"
        File credible_sets_tsv = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.credible_sets.tsv"
        File susie_rds         = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.rds"
        File excluded_snps_tsv = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.excluded_snps.tsv"
        File diagnostics_txt   = "susie_~{phenotype}_~{chrom}_~{start}_~{end}.diagnostics.txt"
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

        library(data.table)
        library(ggplot2)
        library(cowplot)

        # Function to clean variant IDs: remove DRAGEN: prefix and replace colons with underscores
        clean_variant_id <- function(id) {
            id <- gsub("^DRAGEN:", "", id)
            id <- gsub(":", "_", id)
            return(id)
        }

        phenotype <- "~{phenotype}"
        chrom <- "~{chrom}"
        locus_chr <- if (grepl("^chr", chrom, ignore.case = TRUE)) chrom else paste0("chr", chrom)
        locus_start <- as.numeric("~{start}")
        locus_end <- as.numeric("~{end}")
        index_snp <- "~{index_snp}"
        p_threshold <- ~{p_threshold}
        out_prefix <- paste0("locuszoom_", phenotype, "_", chrom, "_", locus_start, "_", locus_end)

        zdf <- as.data.table(fread("~{z_scores}"))
        pip_df <- as.data.table(fread("~{pip_tsv}"))

        cs_df <- if (file.exists("~{credible_sets_tsv}") && file.info("~{credible_sets_tsv}")$size > 0) {
            as.data.table(fread("~{credible_sets_tsv}"))
        } else {
            data.table(ID = character(), credible_set = character())
        }

        if (!"ID" %in% names(zdf)) stop("z_scores is missing ID column")
        if (!"ID" %in% names(pip_df)) stop("pip_tsv is missing ID column")
        if (!"PIP" %in% names(pip_df)) pip_df[, PIP := NA_real_]
        if (!"ID" %in% names(cs_df)) cs_df[, ID := character()]
        if (!"credible_set" %in% names(cs_df)) cs_df[, credible_set := character()]

        zdf[, ID := as.character(ID)]
        pip_df[, `:=`(ID = as.character(ID), PIP = as.numeric(PIP))]
        cs_df[, `:=`(ID = as.character(ID), credible_set = as.character(credible_set))]

        variant_ids <- as.character(readLines("~{variant_ids}"))
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

        load("/usr/local/lib/R/gencode.rda")
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

        # Create display version of variant IDs (cleaned: remove DRAGEN: prefix, replace : with _)
        plot_df[, display_id := sapply(ID, clean_variant_id)]
        display_lead_id <- if (!is.na(lead_id)) clean_variant_id(lead_id) else NA_character_

        label_df <- unique(rbindlist(list(
            plot_df[ID == lead_id],
            plot_df[order(-PIP)][1:min(5, .N)],
            plot_df[in_credible_set == TRUE][1:min(5, .N)]
        ), fill = TRUE), by = "ID")

        p_assoc <- ggplot(plot_df, aes(x = POS, y = PIP)) +
            geom_point(aes(color = ld_bin, size = pmax(PIP, 0.02), shape = in_credible_set), alpha = 0.85) +
            geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "firebrick") +
            geom_text(
                data = label_df,
                aes(label = display_id),
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
                title = paste0("Locus plot: ", phenotype, " chr", chrom, "_", locus_start, "_", locus_end),
                subtitle = paste0("Lead/index SNP: ", display_lead_id),
                x = NULL,
                y = "PIP",
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

            genes[, y_center := (track - 1) * track_spacing]
            genes[, `:=`(
                y_min   = y_center - gene_body_height / 2,
                y_max   = y_center + gene_body_height / 2,
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
                        fill = "white",
                        color = "#4682b4",
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
                    color = "black"
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
        dx_instance_type: "mem2_ssd1_v2_x8"
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

task check_indel_pass {
    input {
        File loci_manifest
        Array[File] pvar_files
    }
    command <<<
        set -euo pipefail
        python3 -m venv .venv
        .venv/bin/pip install pandas --quiet

        .venv/bin/python3 <<'PYEOF'
import pandas as pd
import os
import subprocess

loci = pd.read_csv("~{loci_manifest}", sep="\t")
os.makedirs("indel_pass_results", exist_ok=True)

# Get pvar files from WDL input
pvar_files_found = ["~{sep='", "' pvar_files}"]
print(f"Found {len(pvar_files_found)} pvar files", flush=True)
for f in pvar_files_found:
    print(f"  {f}", flush=True)

summary_rows = []
pvar_regions = []

# Build a list of (chrom, start, end, filename) from pvar files
for pvar_file in pvar_files_found:
    # Extract chrom and positions from filename (e.g., ld_12_13343355_15343355_residual.extracted.pvar)
    basename = os.path.basename(pvar_file)
    if basename.startswith("ld_") and basename.endswith(".extracted.pvar"):
        parts = basename.split("_")
        if len(parts) >= 4:
            try:
                chrom = parts[1]
                start_pos = int(parts[2])
                end_pos = int(parts[3])
                pvar_regions.append((chrom, start_pos, end_pos, pvar_file))
                print(f"Parsed: chr{chrom}:{start_pos}-{end_pos} -> {basename}", flush=True)
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse {basename}: {e}", flush=True)

for _, row in loci.iterrows():
    chrom    = str(row["chrom"]).replace("chr", "")
    pos      = int(row["pos"])
    locus_id = f"chr{chrom}_{pos}"
    window   = 10000
    start    = max(0, pos - window)
    end      = pos + window

    # Find the pvar file whose range contains this locus position
    local_pvar = None
    for pvar_chrom, pvar_start, pvar_end, pvar_file in pvar_regions:
        if pvar_chrom == chrom and pvar_start <= pos <= pvar_end:
            local_pvar = pvar_file
            break

    if not local_pvar or not os.path.exists(local_pvar):
        print(f"  {locus_id}: NO_FILE (searched {len(pvar_regions)} regions)", flush=True)
        summary_rows.append({"locus": locus_id, "chrom": chrom,
                             "pos": pos, "verdict_flags": "NO_FILE"})
        continue

    result = subprocess.run(
        f"grep -v '^#' {local_pvar} | awk -v s={start} -v e={end} '$2+0>=s && $2+0<=e'",
        shell=True, capture_output=True, text=True
    )
    lines = [l.split("\t") for l in result.stdout.strip().split("\n") if l]

    flags = []
    multi = sum(1 for l in lines if len(l) > 4 and "," in l[4])
    indels    = sum(1 for l in lines if len(l) > 4
                    and abs(int(l[1]) - pos) <= 1000 and len(l[4]) > 50)
    exact     = [l for l in lines if len(l) > 1 and int(l[1]) == pos]
    n_alleles = len(exact)

    if multi     > 0: flags.append("MULTIALLELIC_NEARBY")
    if indels    > 0: flags.append("LONG_INDEL_NEARBY")
    if n_alleles > 1: flags.append("MULTIPLE_ROWS_AT_POS")

    summary_rows.append({
        "locus":                   locus_id,
        "chrom":                   chrom,
        "pos":                     pos,
        "n_alleles_at_pos":        n_alleles,
        "multiallelic_count_10kb": multi,
        "long_indel_nearby":       indels > 0,
        "variant_density_10kb":    len(lines),
        "verdict_flags":           ";".join(flags) or "NONE"
    })
    print(f"  {locus_id} -> {';'.join(flags) or 'NONE'}", flush=True)

pd.DataFrame(summary_rows).to_csv(
    "indel_pass_results/indel_pass_summary.tsv", sep="\t", index=False
)
print(f"Done. Wrote {len(summary_rows)} loci.", flush=True)
PYEOF
    >>>
    output {
        File indel_summary = "indel_pass_results/indel_pass_summary.tsv"
    }
    runtime {
        docker: "python:3.11-slim"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}
