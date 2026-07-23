#!/usr/bin/env python3
import subprocess
import pandas as pd

def run_subprocess(x, cwd=None, check=True, shell=True):
    if cwd is None:
        cwd = work_dir

    subprocess.run(
		x,
		cwd=cwd,
		check=check,
		shell=shell
)

work_dir = "/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm"
output_dir = "/Volumes/bluecub/aglycone_release_100um_24h/output"
mol_info_dir = f"{work_dir}/molekyler/information"
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
    exp_mtx[["experiment", "condition"]]
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

    run_subprocess(
        f"Rscript script_tools/filter_intensity.R \
        -i {input_dir} \
        -s 0.2 \
        -f 5"
    )

    run_subprocess(
        f"Rscript script_tools/chromatogram_plots.R \
        -i {input_dir} \
        -f {input_dir}/report/retained_features.csv"
    )

    run_subprocess(
        f"conda run -n psm_chem python script_tools/create_report.py \
        -i {input_dir} \
        -s '{aglycone_smiles}' \
        -c {input_dir}/report/features.parquet \
        -t 0.2 \
        -m 0"
    )