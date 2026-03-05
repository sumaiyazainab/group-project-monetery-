#!/usr/bin/env python
# coding: utf-8

# In[1]:


# =========================
# VISUALISATION MODULE
# =========================
# In this file, I generate the required visuals for the Directed Evolution portal:
# 1) A Top 10 performers table (ranked by Activity Score)
# 2) A per-generation distribution plot (violin plot)
# 3) A 3D activity landscape (sequence diversity on X/Y, activity on Z)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsRegressor


# ============================================================
# 1) TOP PERFORMERS TABLE (Top 10 leaderboard)
# ============================================================
# Here I build a clean “Top 10 variants” table.
# This is useful because it quickly shows the best variants according to Activity_Score.
# I include essential experiment fields AND mutation summary fields so the user can
# see both performance and how many mutations were needed to achieve it.

ESSENTIAL_FIELDS = [
    "Plasmid_Variant_Index",
    "Parent_Plasmid_Variant",
    "Directed_Evolution_Generation",
    "DNA_Quantification_fg",
    "Protein_Quantification_pg",
    "Control",
]

MUTATION_SUMMARY_FIELDS = [
    "mutation_count",
    "synonymous",
    "nonsynonymous",
    "truncating",
]

def top_performers_table(df, n=10, exclude_controls=True):
    """
    I create a Top N table ranked by Activity_Score.

    What this shows:
    - Which variants are performing best overall (highest Activity_Score)
    - Which generation they came from
    - How many mutations they contain (mutation_count)
    - What type of mutations they have (synonymous vs nonsynonymous, truncation)

    Why I do it this way:
    - A table is the fastest “at a glance” summary for scientists.
    - Excluding Control rows avoids WT appearing in the Top 10 (mutants are the goal).
    """

    df = df.copy()

    # I check the minimum required columns exist before I try to rank anything.
    required = ["Activity_Score"] + ESSENTIAL_FIELDS[:-1]  # everything except Control is essential
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for Top Performers table: {missing}")

    # I usually remove WT/control rows from the leaderboard because we want top mutants.
    if exclude_controls and "Control" in df.columns:
        df = df[df["Control"] == False].copy()

    # I rank highest Activity_Score at the top.
    df = df.sort_values("Activity_Score", ascending=False).copy()

    # I build the final columns for display. If some columns are missing, I skip them safely.
    cols = ESSENTIAL_FIELDS + ["Activity_Score"] + MUTATION_SUMMARY_FIELDS
    cols = [c for c in cols if c in df.columns]

    top = df[cols].head(n).reset_index(drop=True)

    # I add a rank column so it looks like a “leaderboard”.
    top.insert(0, "Rank", range(1, len(top) + 1))

    # I round the score so it looks neat in the portal UI.
    top["Activity_Score"] = top["Activity_Score"].astype(float).round(4)

    return top


# --- Testing example ---
# top10 = top_performers_table(df_scored, n=10)
# top10



# ============================================================
# 2) ACTIVITY DISTRIBUTION BY GENERATION (Violin Plot)
# ============================================================
# Here I show how the Activity Score distribution changes across generations.
# This is one of the most important plots because it tells us:
# “Is the directed evolution experiment improving activity over time?”
#
# Why violin plot?
# - It shows the full distribution shape (not just mean/median).
# - Wide violin = lots of diversity; narrow = population converging.
# - If violins shift upward generation by generation, selection is working.

def plot_activity_distribution_by_generation(
    df,
    score_col="Activity_Score",
    gen_col="Directed_Evolution_Generation",
    exclude_controls=True,
    control_col="Control",
    title="Activity Score Distribution by Generation",
    save_path=None,
    dpi=200,
):
    """
    I plot the distribution of Activity_Score for each generation (violin plot).

    What this shows:
    - Each generation = one violin
    - The vertical axis = Activity Score
    - If the violins gradually rise, it suggests that later generations are better
      (the experiment is heading in the right direction).
    """

    plot_df = df.copy()

    # I remove controls so the plot reflects mutant progress, not WT baselines.
    if exclude_controls and control_col in plot_df.columns:
        plot_df = plot_df[plot_df[control_col] == False]

    # I remove any rows missing generation or activity score.
    plot_df = plot_df.dropna(subset=[gen_col, score_col])

    # I sort generations so they appear in the correct order on the x-axis.
    gens = sorted(plot_df[gen_col].unique())

    # I build one list of Activity Score values per generation (required by violinplot).
    data = []
    kept_gens = []
    for g in gens:
        vals = plot_df.loc[plot_df[gen_col] == g, score_col].astype(float).values
        vals = vals[~np.isnan(vals)]
        if len(vals) > 0:
            data.append(vals)
            kept_gens.append(g)

    if len(data) == 0:
        raise ValueError("No data available to plot after filtering.")

    # Plot
    plt.figure(figsize=(10, 5))

    # I show mean + median so the user can quickly see the central trend.
    plt.violinplot(data, showmeans=True, showmedians=True, showextrema=True)

    plt.xticks(np.arange(1, len(kept_gens) + 1), [str(g) for g in kept_gens])
    plt.xlabel("Generation")
    plt.ylabel("WT-normalised Activity Score")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # I save a PNG if needed (useful for Flask/Django static files or PDF reports).
    if save_path:
        if not save_path.lower().endswith(".png"):
            save_path += ".png"
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    # I always show it in Jupyter so I can preview the figure easily.
    plt.show()

    return save_path


# --- Testing example ---
# plot_activity_distribution_by_generation(df_scored)



# ============================================================
# 3) 3D ACTIVITY LANDSCAPE (PCA + KNN Surface)
# ============================================================
# This is the “bonus” high-impact visual.
#
# The biological story:
# - Similar protein sequences should cluster together.
# - Some clusters will have high Activity Score (peaks).
# - This gives a “fitness landscape” view of sequence space.
#
# How I build it:
# 1) I convert protein sequences into numeric features using k-mer counts.
# 2) I reduce those features down to 2D using PCA (Dim1, Dim2).
# 3) I use Activity Score as the height (Z axis).
# 4) I smooth the surface using KNN regression to create a mountain-like topography.
#
# Why PCA?
# - It’s fast, and explainable (good for assessment + reproducibility).
# Why KNN smoothing?
# - The raw points are sparse; KNN makes a continuous “terrain” that looks like a landscape.

def make_3d_activity_landscape(
    df,
    seq_col="protein_seq",
    score_col="Activity_Score",
    exclude_controls=True,
    control_col="Control",
    kmer_k=3,
    grid_size=70,
    knn_neighbors=7,
    knn_weights="distance",
    title="3D Activity Landscape (Sequence diversity → PCA → Activity Score)",
    save_html=None,
):
    """
    I create an interactive 3D activity landscape.

    Axes meaning:
    - X = PCA Dim 1 (sequence diversity)
    - Y = PCA Dim 2 (sequence diversity)
    - Z = Activity Score (performance)

    What this shows:
    - Peaks = high activity variants
    - Valleys = low activity variants
    - Clusters = similar sequences
    """

    plot_df = df.copy()

    # I remove controls so the peaks represent mutant improvements.
    if exclude_controls and control_col in plot_df.columns:
        plot_df = plot_df[plot_df[control_col] == False].copy()

    # I keep only rows where sequence + Activity Score exist.
    plot_df = plot_df.dropna(subset=[seq_col, score_col]).copy()
    plot_df[seq_col] = plot_df[seq_col].astype(str)

    if len(plot_df) < 5:
        raise ValueError("Not enough rows with sequence + Activity_Score to build landscape.")

    # 1) Turn sequences into k-mer feature vectors
    vectorizer = CountVectorizer(analyzer="char", ngram_range=(kmer_k, kmer_k))
    X = vectorizer.fit_transform(plot_df[seq_col])

    # 2) PCA reduces high-dimensional k-mer space into 2D coordinates
    coords = PCA(n_components=2, random_state=42).fit_transform(X.toarray())
    plot_df["Dim1"] = coords[:, 0]
    plot_df["Dim2"] = coords[:, 1]

    # 3) Create a grid across the PCA space so I can draw a continuous surface
    x_lin = np.linspace(plot_df["Dim1"].min(), plot_df["Dim1"].max(), grid_size)
    y_lin = np.linspace(plot_df["Dim2"].min(), plot_df["Dim2"].max(), grid_size)
    xx, yy = np.meshgrid(x_lin, y_lin)
    grid_xy = np.c_[xx.ravel(), yy.ravel()]

    # 4) KNN regression predicts activity score at each grid point (this is the smoothing)
    xy = plot_df[["Dim1", "Dim2"]].values
    z = plot_df[score_col].astype(float).values

    knn = KNeighborsRegressor(n_neighbors=knn_neighbors, weights=knn_weights)
    knn.fit(xy, z)
    zz = knn.predict(grid_xy).reshape(grid_size, grid_size)

    # 5) Build interactive 3D plot (Plotly)
    fig = go.Figure()

    # Surface = smoothed landscape
    fig.add_trace(go.Surface(
        x=xx, y=yy, z=zz,
        opacity=0.85,
        colorbar=dict(title="Activity Score"),
        name="Landscape"
    ))

    # Points = real variants
    hover_text = []
    for _, r in plot_df.iterrows():
        hover_text.append(
            f"Variant: {r.get('Plasmid_Variant_Index', '')}<br>"
            f"Gen: {r.get('Directed_Evolution_Generation', '')}<br>"
            f"Score: {float(r[score_col]):.3f}<br>"
            f"Mutations: {r.get('mutation_count', '')}"
        )

    fig.add_trace(go.Scatter3d(
        x=plot_df["Dim1"],
        y=plot_df["Dim2"],
        z=plot_df[score_col],
        mode="markers",
        marker=dict(size=3, opacity=0.7),
        text=hover_text,
        hoverinfo="text",
        name="Variants"
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="Sequence diversity (PCA Dim 1)",
            yaxis_title="Sequence diversity (PCA Dim 2)",
            zaxis_title="Activity Score",
        ),
        width=950,
        height=700,
        margin=dict(l=0, r=0, b=0, t=50),
    )

    # Saving as HTML is ideal for a web portal (interactive + easy to embed/download)
    if save_html:
        if not save_html.lower().endswith(".html"):
            save_html += ".html"
        fig.write_html(save_html, include_plotlyjs=True, full_html=True)

    return fig


# --- Testing example ---
# fig = make_3d_activity_landscape(df_scored, save_html="activity_landscape.html")
# fig.show()


# In[ ]:




