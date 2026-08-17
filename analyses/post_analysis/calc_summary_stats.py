# %%
data_with_info = (
    data
    .join(
        other=exp_info,
        on="experiment",
        how="left"
    )
    .filter(
        # pl.col("source") == "biotransformer",
        pl.col("adduct") == "[M-H]-",
    )
)
# %%
unique_feats_per_exp = (
    data_with_info
    .group_by("experiment")
    .agg(
        pl.col("feature").n_unique().alias("feats"),
        pl.col("condition").first(),
        pl.col("strain").first()
    )
)

# Features per experiment
all_hits_points = alt.Chart(unique_feats_per_exp).mark_circle(size = 100).encode(
    x=alt.X(
        "condition:N",
        sort=alt.EncodingSortField(
            field="feats",
            order="descending"
        )
    ),
    y=alt.Y("feats:Q"),
    fill=alt.Fill("strain:N"),
)

all_hits_points.save(fp=chart_path.joinpath("all_hits.pdf"), format="pdf")

# %%
all_exps = data_with_info["condition"].unique().sort().to_list()
chart_path = Path(
    "/Volumes/bluecub/aglycone_release_100um_24h/combined_results/hit_heatmaps"
)

for i in all_exps[:]:
    strains_in_condition = (
        data_with_info
        .filter(pl.col("condition") == i)
        .get_column("strain")
        .unique()
        .sort()
        .to_list()
    )

    for s in strains_in_condition:
        tink = (
            data_with_info
            .filter(
                pl.col("condition") == i,
                pl.col("strain") == s,
                # pl.col("sim") > 0.1
            )
            .with_columns(
                pl.col("skeleton_key").n_unique().over("feature").alias("total_hits")
            )
            .filter(
                pl.col("sim").rank(method="ordinal", descending=True).over("feature") <= 5
            )
            .unique(subset=["feature", "skeleton_key"])
        )

        tink_img = tink.with_columns(
            pl.col("smiles")
            .map_elements(smiles_to_data_uri, return_dtype=pl.String)
            .alias("structure")
        )

        tink_img = tink_img.with_columns(
            (pl.col("feature") + " (" + pl.col("total_hits").cast(pl.String) + ")").alias("feature_label")
        )

        hmp = alt.Chart(tink_img).mark_rect().encode(
            x=alt.X(
                "skeleton_key:N",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                axis=None
            ),
            y=alt.Y("feature_label:N", title=None),
            color=alt.Color(
                "sim:Q",
                scale=alt.Scale(domain=[0, 1]),
                legend=alt.Legend(title="Tanimoto similarity")
            ),
            tooltip=["skeleton_key", "smiles"],
        )

        structures = hmp.mark_image(width=50, height=50).encode(
            url="structure:N",
            x=alt.X(
                "skeleton_key:N",
                sort=alt.EncodingSortField(field="sim", order="descending"),
                axis=None
            ),
            y=alt.value(20),
        ).properties(height=50, width=alt.Step(50))

        full_hmp = (
            alt.vconcat(hmp, structures, spacing=0)
            .resolve_scale(x="shared")
            .configure_view(strokeWidth=0)
        )

        full_hmp.save(fp=chart_path.joinpath(f"{i}_{s}.pdf"), format="pdf")

# %%
