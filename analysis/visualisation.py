import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsRegressor


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
    Top N variants ranked by Activity_Score, with essential fields + mutation summary.
    """
    df = df.copy()

    # Safety: ensure required columns exist
    required = ["Activity_Score"] + ESSENTIAL_FIELDS[:-1]  # everything except Control is absolutely required
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for Top Performers table: {missing}")

    # Exclude WT controls from leaderboard
    if exclude_controls and "Control" in df.columns:
        df = df[df["Control"] == False].copy()

    # Sort by score
    df = df.sort_values("Activity_Score", ascending=False).copy()

    # Build final column list (only keep those that exist)
    cols = ESSENTIAL_FIELDS + ["Activity_Score"] + MUTATION_SUMMARY_FIELDS
    cols = [c for c in cols if c in df.columns]

    top = df[cols].head(n).reset_index(drop=True)

    # Cosmetic formatting
    top.insert(0, "Rank", range(1, len(top) + 1))
    top["Activity_Score"] = top["Activity_Score"].astype(float).round(4)

    return top

# Usage:
#top10 = top_performers_table(df_scored, n=10)
#top10


# In[26]:


import numpy as np
import matplotlib.pyplot as plt

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
    plot_df = df.copy()

    if exclude_controls and control_col in plot_df.columns:
        plot_df = plot_df[plot_df[control_col] == False]

    plot_df = plot_df.dropna(subset=[gen_col, score_col])
    gens = sorted(plot_df[gen_col].unique())

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

    plt.figure(figsize=(10, 5))
    plt.violinplot(data, showmeans=True, showmedians=True, showextrema=True)

    plt.xticks(np.arange(1, len(kept_gens) + 1), [str(g) for g in kept_gens])
    plt.xlabel("Generation")
    plt.ylabel("WT-normalised Activity Score")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save if requested
    if save_path:
        if not save_path.lower().endswith(".png"):
            save_path += ".png"
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    # Always show in Jupyter
    plt.close()

    return save_path

    #plot_activity_distribution_by_generation(df_scored)


# In[27]:


#plot_activity_distribution_by_generation(df_scored)


# In[28]:


import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsRegressor

def make_3d_activity_landscape(
    df,
    seq_col="protein_seq",
    score_col="Activity_Score",
    exclude_controls=True,
    control_col="Control",
    kmer_k=3,                 # 3-mers for proteins is a good default
    grid_size=70,             # higher = smoother but slower
    knn_neighbors=7,          # smoothing strength
    knn_weights="distance",
    title="3D Activity Landscape (Sequence diversity → PCA → Activity Score)",
    save_html=None,           # e.g. "activity_landscape.html"
):
    plot_df = df.copy()

    # Optional: drop controls so peaks reflect mutants
    if exclude_controls and control_col in plot_df.columns:
        plot_df = plot_df[plot_df[control_col] == False].copy()

    # Keep only valid rows
    plot_df = plot_df.dropna(subset=[seq_col, score_col]).copy()
    plot_df[seq_col] = plot_df[seq_col].astype(str)

    if len(plot_df) < 5:
        raise ValueError("Not enough rows with sequence + Activity_Score to build landscape.")

    # 1) Convert sequences → k-mer count vectors
    vectorizer = CountVectorizer(analyzer="char", ngram_range=(kmer_k, kmer_k))
    X = vectorizer.fit_transform(plot_df[seq_col])

    # 2) PCA → 2D coordinates
    coords = PCA(n_components=2, random_state=42).fit_transform(X.toarray())
    plot_df["Dim1"] = coords[:, 0]
    plot_df["Dim2"] = coords[:, 1]

    # 3) Build grid over PCA space
    x_lin = np.linspace(plot_df["Dim1"].min(), plot_df["Dim1"].max(), grid_size)
    y_lin = np.linspace(plot_df["Dim2"].min(), plot_df["Dim2"].max(), grid_size)
    xx, yy = np.meshgrid(x_lin, y_lin)
    grid_xy = np.c_[xx.ravel(), yy.ravel()]

    # 4) Smooth surface via KNN regression: (Dim1, Dim2) -> Activity_Score
    xy = plot_df[["Dim1", "Dim2"]].values
    z = plot_df[score_col].astype(float).values

    knn = KNeighborsRegressor(n_neighbors=knn_neighbors, weights=knn_weights)
    knn.fit(xy, z)
    zz = knn.predict(grid_xy).reshape(grid_size, grid_size)

    # 5) Build interactive 3D figure
    fig = go.Figure()

    # Surface = mountain range (smoothed activity)
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

    # Save for web portal
    if save_html:
        if not save_html.lower().endswith(".html"):
            save_html += ".html"
        fig.write_html(save_html, include_plotlyjs=True, full_html=True)

    return fig


# In[29]:


#fig = make_3d_activity_landscape(df_scored, save_html="activity_landscape.html")
#fig.show()

