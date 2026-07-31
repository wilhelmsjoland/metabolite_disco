import os, shutil, glob

onstart:
    os.makedirs(f"{config['output']}/logs", exist_ok=True)

onsuccess:
    logs = sorted(glob.glob(".snakemake/log/*.snakemake.log"))
    if logs:
        shutil.copy(logs[-1], f"{config['output']}/logs/{os.path.basename(logs[-1])}")

onerror:
    logs = sorted(glob.glob(".snakemake/log/*.snakemake.log"))
    if logs:
        shutil.copy(logs[-1], f"{config['output']}/logs/{os.path.basename(logs[-1])}")

rule met_disco:
    input:
        f"{config['output']}/snakemake_objects/11_pca.rds",
        f"{config['output']}/snakemake_objects/12_volcano.rds",
        f"{config['output']}/snakemake_objects/14_intersecting_features.rds",
        f"{config['output']}/snakemake_objects/16_mz_predictions.rds",
        f"{config['output']}/snakemake_objects/19_filter_matches.rds"

rule setup:
    output:
        f"{config['output']}/snakemake_objects/01_setup.rds"
    params:
        output = config["output"],
        data_path = config["data_path"],
        meta_file = config["meta_file"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/01_setup.R"

rule bpc:
    input:
        f"{config['output']}/snakemake_objects/01_setup.rds"
    output:
        f"{config['output']}/snakemake_objects/02_bpc.rds"
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/02_bpc.R"

rule internal_standard:
    input:
        setup = f"{config['output']}/snakemake_objects/01_setup.rds",
        bpc = f"{config['output']}/snakemake_objects/02_bpc.rds"
    output:
        f"{config['output']}/snakemake_objects/03_internal_standard.rds"
    params:
        output = config["output"],
        internal_standard = config["internal_standard"],
        is_adduct = config["is_adduct"],
        ppm_global = config["ppm_global"],
        sn_threshold = config["sn_threshold"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/03_internal_standard.R"

rule peak_calling:
    input:
        setup = f"{config['output']}/snakemake_objects/01_setup.rds",
        internal_std = f"{config['output']}/snakemake_objects/03_internal_standard.rds"
    output:
        f"{config['output']}/snakemake_objects/04_peak_calling.rds"
    params:
        output = config["output"],
        ppm_global = config["ppm_global"],
        sn_threshold = config["sn_threshold"],
        mzdiff = config["mzdiff"],
        beta_cor_threshold = config["beta_cor_threshold"],
        beta_snr_threshold = config["beta_snr_threshold"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/04_peak_calling.R"

rule alignment:
    input:
        setup = f"{config['output']}/snakemake_objects/01_setup.rds",
        bpc = f"{config['output']}/snakemake_objects/02_bpc.rds",
        internal_std = f"{config['output']}/snakemake_objects/03_internal_standard.rds",
        peak_calling = f"{config['output']}/snakemake_objects/04_peak_calling.rds"
    output:
        f"{config['output']}/snakemake_objects/05_alignment.rds"
    params:
        output = config["output"],
        ppm_global = config["ppm_global"],
        bw_first_grouping = config["bw_first_grouping"],
        peak_anchor_sd = config["peak_anchor_sd"],
        min_fraction_align = config["min_fraction_align"],
        extra_peaks = config["extra_peaks"],
        span = config["span"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/05_alignment.R"

rule correspondence:
    input:
        alignment = f"{config['output']}/snakemake_objects/05_alignment.rds",
        internal_std = f"{config['output']}/snakemake_objects/03_internal_standard.rds"
    output:
        f"{config['output']}/snakemake_objects/06_correspondence.rds"
    params:
        output = config["output"],
        ppm_global = config["ppm_global"],
        bw_first_grouping = config["bw_first_grouping"],
        bw_second_grouping = config["bw_second_grouping"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/06_correspondence.R"

rule gap_filling:
    input:
        correspondence = f"{config['output']}/snakemake_objects/06_correspondence.rds"
    output:
        f"{config['output']}/snakemake_objects/07_gap_filling.rds"
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/07_gap_filling.R"

rule filter_features:
    input:
        gap_filling = f"{config['output']}/snakemake_objects/07_gap_filling.rds"
    output:
        f"{config['output']}/snakemake_objects/08_filter_features.rds"
    params:
        output = config["output"],
        missingness = config["missingness"],
        sn_threshold = config["sn_threshold"],
        beta_cor_threshold = config["beta_cor_threshold"],
        beta_snr_threshold = config["beta_snr_threshold"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/08_filter_features.R"

rule scaling:
    input:
        filter_features = f"{config['output']}/snakemake_objects/08_filter_features.rds"
    output:
        f"{config['output']}/snakemake_objects/09_scaling.rds"
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/09_scaling.R"

rule limma:
    input:
        scaling = f"{config['output']}/snakemake_objects/09_scaling.rds",
        setup = f"{config['output']}/snakemake_objects/01_setup.rds"
    output:
        f"{config['output']}/snakemake_objects/10_limma.rds"
    params:
        output = config["output"],
        gap_filling = config["gap_filling"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/10_limma.R"

rule pca:
    input:
        limma = f"{config['output']}/snakemake_objects/10_limma.rds",
        setup = f"{config['output']}/snakemake_objects/01_setup.rds"
    output:
        f"{config['output']}/snakemake_objects/11_pca.rds"
    threads: 1
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/11_pca.R"

rule volcano:
    input:
        limma = f"{config['output']}/snakemake_objects/10_limma.rds"
    output:
        f"{config['output']}/snakemake_objects/12_volcano.rds"
    threads: 1
    params:
        output = config["output"],
        qvalue = config["qvalue"],
        gap_filling = config["gap_filling"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/12_volcano.R"

rule upset:
    input:
        limma = rules.limma.output[0]
    output:
        f"{config['output']}/snakemake_objects/13_upset.rds"
    threads: 1
    params:
        output = config["output"],
        qvalue = config["qvalue"],
        gap_filling = config["gap_filling"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/13_upset.R"

rule intersecting_features:
    input:
        upset = rules.upset.output[0]
    output:
        f"{config['output']}/snakemake_objects/14_intersecting_features.rds"
    threads: 1
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/14_intersecting_features.R"

rule prep_annotation_biotransformation:
    input:
        filter_features = rules.filter_features.output[0],
        limma = rules.limma.output[0]
    output:
        f"{config['output']}/snakemake_objects/15_prep_annotation_biotransformation.rds"
    threads: 1
    params:
        output = config["output"],
        data_path = config["data_path"],
        polarity = config["polarity"],
        mass_shift_path = config["mass_shift_path"],
        rpairs_path = config["rpairs_path"],
        qvalue = config["qvalue"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/15_prep_annotation_biotransformation.R"

rule mz_predictions:
    input:
        limma = rules.limma.output[0],
        prep_biot = rules.prep_annotation_biotransformation.output[0]
    output:
        f"{config['output']}/snakemake_objects/16_mz_predictions.rds"
    threads: config["cores"]
    params:
        output = config["output"],
        metabolite_search = config["metabolite_search"],
        polarity = config["polarity"],
        mass_shift_path = config["mass_shift_path"],
        ppm_match = config["ppm_match"],
        all_vs_all = config["all_vs_all"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/16_mz_predictions.R"

rule annotation:
    input:
        limma = rules.limma.output[0],
        prep_biot = rules.prep_annotation_biotransformation.output[0],
        filter_features = rules.filter_features.output[0]
    output:
        f"{config['output']}/snakemake_objects/17_annotation.rds"
    threads: config["cores"]
    params:
        output = config["output"],
        polarity = config["polarity"],
        ppm_match = config["ppm_match"],
        seed = config["seed"],
        cores = config["cores"]
    script:
        "scripts/17_annotation.R"

rule biotransformer:
    input:
        filter_features = rules.filter_features.output[0]
    output:
        f"{config['output']}/snakemake_objects/18_biotransformer.rds"
    # One process per molecule (see scripts/18_biotransformer.R) - never
    # needs more threads than there are SMILES strings to predict from
    threads: min(config["cores"], len(config["smiles"].split(";")))
    params:
        output = config["output"],
        smiles = config["smiles"],
        polarity = config["polarity"],
        ppm_match = config["ppm_match"],
        seed = config["seed"],
        cores = config["cores"],
        annotate_path = config.get("annotate_path")
    script:
        "scripts/18_biotransformer.R"

rule filter_matches:
    input:
        annotation = rules.annotation.output[0],
        biotransformer = rules.biotransformer.output[0]
    output:
        f"{config['output']}/snakemake_objects/19_filter_matches.rds"
    threads: 1
    params:
        output = config["output"],
        seed = config["seed"],
        cores = config["cores"],
        annotate_path = config.get("annotate_path")
    script:
        "scripts/19_filter_matches.R"
