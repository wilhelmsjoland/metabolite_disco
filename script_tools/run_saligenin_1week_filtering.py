#!/usr/bin/env python3
# %% 
import subprocess
import os
from pathlib import Path
import polars as pl

def run_subprocess(x, cwd=None, check=True, shell=True):
    if cwd is None:
        cwd = repo_dir

    subprocess.run(
		x,
		cwd=cwd,
		check=check,
		shell=shell
)

work_dir = Path(
    f"/cfs/klemming/projects/supr/sjoland_naiss/data/"
    f"260805_saligenin_1week/output/260805_saligenin_1week_ms1"
)
repo_dir = Path(
    "/cfs/klemming/projects/supr/sjoland_naiss/project/psm/metabolite_disco"
)
metadata_dir = Path(
    f"/cfs/klemming/projects/supr/sjoland_naiss/data/"
    f"260805_saligenin_1week/metadata"
    )
unfinished_runs = ["Ser"]

exp_list = sorted(os.listdir(work_dir))
for i in exp_list:
    if i not in unfinished_runs:
        # extract grouping metadata
        meta_f_path = metadata_dir.joinpath(f"{i}.csv")
        meta_f = pl.read_csv(meta_f_path)
        groups = (
            meta_f
            .filter(pl.col("group").str.ends_with("_P"))
            .unique("group")
            ["group"][0] 
        )
        input_dir = Path(work_dir).joinpath(f"{i}")
        aglycone = "saligenin"
        aglycone_smiles = "C1=CC=C(C(=C1)CO)O"

        run_subprocess(
            f"python script_tools/calc_similarity.py \
            --input '{input_dir}' \
            --smiles '{aglycone},{aglycone_smiles}' \
            --radius 2 \
            --fpsize 4096"
        )

        run_subprocess(
            f"Rscript script_tools/filter_features.R \
            --input {input_dir} \
            --grouping {groups} \
            --similarity_filter 0.05 \
            --lfc 1 \
            --qval 0.1 \
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
            --smiles '{aglycone},{aglycone_smiles}' \
            --chromatogram {input_dir}/report/features.parquet \
            --similarity_cutoff 0.1 \
            --mcs_cutoff 0"
        )
