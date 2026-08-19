#!/usr/bin/env python3
import subprocess
import pandas as pd

def run_subprocess(x, cwd=None, check=True, shell=True):
    if cwd is None:
        cwd = repo_dir

    subprocess.run(
		x,
		cwd=cwd,
		check=check,
		shell=shell
)

# Two different roots: the xlsx inputs live under the project folder, the
# scripts under the repo inside it - the subprocess calls below use paths
# relative to repo_dir, everything else hangs off work_dir.
work_dir = "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm"
repo_dir = f"{work_dir}/metabolite_disco"
output_dir = "/Volumes/bluecub/aglycone_release_100um_24h/output"
mol_info_dir = f"{work_dir}/molekyler/information/aglycone_release_100um_24h"
exp_mtx = pd.read_excel(
    io=f"{mol_info_dir}/experiment_matrix.xlsx",
    sheet_name="Sheet1"
)
clean_nms_glycone_list = pd.read_excel(
    io=f"{mol_info_dir}/clean_names_glycone_list.xlsx",
    sheet_name="Sheet1"
)
clean_glycone_metadata = pd.read_excel(
    io=f"{mol_info_dir}/clean_glycone_metadata.xlsx",
    sheet_name="Sheet1"
)

run_df = (
    exp_mtx[["experiment", "condition", "group"]]
    .groupby("experiment", as_index=False)
    .agg({
        "condition": "first",
        "group": lambda g: sorted(set(g)),
    })
    .merge(
        clean_nms_glycone_list,
        left_on="condition",
        right_on="glycoside",
        how="left"
    )
    .drop(columns="condition")
    .dropna(subset="aglycone")
    .drop_duplicates(subset="experiment")
    .merge(
        clean_glycone_metadata[["molecule", "SMILES"]],
        left_on="aglycone",
        right_on="molecule",
        how="left"
    )
    .rename(columns={"SMILES": "aglycone_smiles"})
    .drop(columns="molecule")
    .merge(
        clean_glycone_metadata[["molecule", "SMILES"]],
        left_on="glycoside",
        right_on="molecule",
        how="left"
    )
    .drop(columns="molecule")
    .rename(columns={"SMILES": "glycoside_smiles"})
    .sort_values("experiment")
)

for idx, row in run_df.iterrows():
    input_dir = f"{output_dir}/{row['experiment']}"
    aglycone = row["aglycone"]
    aglycone_smiles = row["aglycone_smiles"]
    glycoside = row["glycoside"]
    glycoside_smiles = row["glycoside_smiles"]
    exp_groups = list(filter(lambda x: "ycfa_glucose" not in x, row["group"]))
    # groups = ",".join(exp_groups)
    groups = ",".join(f"{g}" for g in exp_groups)

    run_subprocess(
        f"python script_tools/calc_similarity.py \
        --input '{input_dir}' \
        --smiles '{aglycone},{aglycone_smiles};{glycoside},{glycoside_smiles}' \
        --radius 2 \
        --fpsize 4096"
    )

    run_subprocess(
        f"Rscript script_tools/filter_features.R \
        --input {input_dir} \
        --grouping {groups} \
        --similarity_filter 0.1 \
        --lfc 2 \
        --qval 0.05 \
        --beta_cor 0.65"
    )

    run_subprocess(
        f"Rscript script_tools/chromatogram_plots.R \
        --input {input_dir} \
        --features {input_dir}/report/retained_features.csv"
    )

    run_subprocess(
        f"python script_tools/create_report.py \
        --input {input_dir} \
        --smiles '{aglycone},{aglycone_smiles};{glycoside},{glycoside_smiles}' \
        --chromatogram {input_dir}/report/features.parquet \
        --similarity_cutoff 0.1 \
        --mcs_cutoff 0"
    )