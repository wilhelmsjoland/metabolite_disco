# %%
import os
from pathlib import Path
import polars as pl

work_dir = Path(
    f"/Users/wilhelm/Proton/01_juniper/01_arbete/01_projekt/03_psm/"
    f"metabolite_disco"
)
exp_dir = Path("/Volumes/bluecub/aglycone_release_100um_24h/output")
out_path = Path("/Volumes/bluecub/aglycone_release_100um_24h/combined_results")
out_save = Path(out_path, "aglycone_release_100um_24h_hits.parquet")

exp_names = (
    [i for i in os.listdir(exp_dir) if os.path.isdir(os.path.join(exp_dir, i))]
)

exp_names.sort()

exp_path = Path(
    exp_dir,
    exp_names[0],
    "report",
    "similar_hits.parquet"
)

data = pl.read_parquet(exp_path)

all_hits = pl.concat([
    pl.read_parquet(Path(exp_dir, name, "report", "similar_hits.parquet")) 
    for name in exp_names
])

all_hits.write_parquet(file=out_save)


# %%
