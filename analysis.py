"""Compute the data behind every view of the atlas website.

Views are every single microbiome plus every pair of microbiomes listed in
microbiomes.json. Each view gets two maps, both made with the interactive
UMAP + HDBSCAN tool:

  * DNA-embedding map: UMAP + HDBSCAN on DNABERT-S genome embeddings.
  * Functional map:    UMAP + HDBSCAN on eggNOG relative-abundance profiles
                       (root-level orthologous groups; Bray-Curtis distance),
                       with the top COG groups per functional cluster.

Pairs also get one Mann-Whitney U test per gene family between the two
microbiomes (Benjamini-Hochberg FDR).

Outputs go to docs/views/<view id>/; build_site.py turns them into HTML pages.

Usage:
    python analysis.py              # all views
    python analysis.py gut--vaginal # only the named view(s)
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import false_discovery_control, mannwhitneyu
from sklearn.metrics import adjusted_rand_score

from atlas_config import SITE, WORK, load_config, views

TOOL = SITE.parents[1] / "Interactive UMAP + HDBSCAN Tool" / "interactive-umap-hdbscan" / "interactive_umap_hdbscan.py"

MIN_PREVALENCE_TESTED = 0.01   # gene family must be present in >= 1% of the view's genomes ...
MIN_GENOMES_TESTED = 5         # ... and in at least this many genomes
FUNCTIONAL_METRIC = "braycurtis"
# HDBSCAN min_samples per map, chosen by seed-stability sweeps: DNA-embedding maps collapse into
# 1-2 clusters at >= 15, so they use 5; functional maps were most stable at ~10.
MIN_SAMPLES = {"dna": 5, "functional": 10}
MIN_CLUSTER_PREVALENCE = 0.10  # top COG terms must be present in >= 10% of the cluster's genomes
TOP_N_SAVED = 50               # per-cluster ranked COG terms kept in the CSV
PSEUDOCOUNT = 1e-6
ALPHA = 0.05


class EggnogStore:
    """Root-level eggNOG OG counts for all genomes (built by Random Forest/build_feature_matrix.py)."""

    def __init__(self, features_dir):
        d = Path(features_dir)
        self.meta = pd.read_csv(d / "sample_metadata.csv")
        self.meta["sample_id"] = self.meta["sample_id"].astype(str)
        self.meta = self.meta.set_index("sample_id")
        self.X = sp.load_npz(d / "eggnog_og_matrix.npz").tocsr()
        self.feats = np.array(json.load(open(d / "eggnog_og_features.json")))
        self.desc = json.load(open(d / "eggnog_og_descriptions.json"))

    def gene_counts(self):
        return self.meta["n_genes"].rename("Gene_count")

    def relative_abundance(self, genome_ids):
        rows = self.meta.index.get_indexer(genome_ids)
        if (rows < 0).any():
            missing = list(np.asarray(genome_ids)[rows < 0][:5])
            raise SystemExit(f"Error: {(rows < 0).sum()} genomes have no eggNOG annotation, e.g. {missing}")
        X = self.X[rows]
        n_present = np.asarray((X > 0).sum(axis=0)).ravel()
        keep = n_present >= max(MIN_GENOMES_TESTED, np.ceil(MIN_PREVALENCE_TESTED * len(rows)))
        X = X[:, keep].tocsc()
        ra = (sp.diags(1.0 / self.meta["n_genes"].values[rows]) @ X).tocsc()
        return ra, (X > 0).tocsc(), self.feats[keep]


def tool(*args):
    subprocess.run([sys.executable, str(TOOL), *map(str, args)], check=True)


def cluster_with_tool(view, kind, inputs, extra_args, annotations):
    """Run UMAP + HDBSCAN via the tool. inputs: [(label, csv)]. Returns (coords, report)."""
    out = WORK / view["id"] / kind
    args = ["--out", out, "--annotations", annotations, "--plotly-js", "cdn",
            "--min-samples", MIN_SAMPLES[kind], *extra_args]
    for label, path in inputs:
        args += ["--input", label, path]
    tool(*args)
    coords = pd.read_csv(out / "coordinates_and_clusters.csv", index_col=0)
    return coords, (out / "report.txt").read_text(encoding="utf-8")


def map_subtitle(coords_csv, metric, min_samples):
    n = len(pd.read_csv(coords_csv, usecols=[0]))
    return (f"{n:,} genomes; UMAP {metric.replace('braycurtis', 'Bray–Curtis')}, n_neighbors=30, min_dist=0.1; "
            f"HDBSCAN min_cluster_size={max(2, round(0.05 * n))}, min_samples={min_samples}")


def make_maps(view, kind, coords_csv, annotations, title, other_cluster_col):
    """Interactive maps for one clustering, recoloured several ways, written to docs/views/<id>/maps/."""
    names = [m["name"] for m in view["members"]]
    subtitle = map_subtitle(coords_csv, FUNCTIONAL_METRIC if kind == "functional" else "euclidean", MIN_SAMPLES[kind])
    colors = [x for m in view["members"] for x in ("--colors", m["name"], m["light"])]
    maps = SITE / "views" / view["id"] / "maps"
    maps.mkdir(parents=True, exist_ok=True)
    variants = [("Cluster", "cluster", "coloured by cluster")]
    if len(names) > 1:
        variants.append(("Dataset", "microbiome", "coloured by microbiome"))
    variants.append((other_cluster_col, "other", f"coloured by {other_cluster_col.replace('_', ' ').lower()}"))
    for col, suffix, what in variants:
        tmp = WORK / view["id"] / f"map_{kind}_{suffix}"
        tool("--precomputed", coords_csv, "--annotations", annotations, "--color-by", col, *colors,
             "--plotly-js", "cdn", "--stability-runs", "0", "--title", f"{title}, {what}",
             "--subtitle", subtitle, "--out", tmp)
        (tmp / "interactive_umap.html").replace(maps / f"{kind}_{suffix}.html")


def pair_test(ra, present, feats, desc, is_b, a, b):
    """Mann-Whitney U per OG; rank-biserial > 0 means higher in microbiome b."""
    n_b, n_a = is_b.sum(), (~is_b).sum()
    u_all, p_all = [], []
    for start in range(0, ra.shape[1], 1000):
        block = ra[:, start:start + 1000].toarray()
        u, p = mannwhitneyu(block[is_b], block[~is_b], axis=0, alternative="two-sided", method="asymptotic")
        u_all.append(u)
        p_all.append(p)
    u, p = np.concatenate(u_all), np.nan_to_num(np.concatenate(p_all), nan=1.0)

    mean_a = np.asarray(ra[~is_b].mean(axis=0)).ravel()
    mean_b = np.asarray(ra[is_b].mean(axis=0)).ravel()
    res = pd.DataFrame({
        "OG": feats,
        "Description": [desc.get(f, "") for f in feats],
        f"Prevalence_{a}": np.asarray(present[~is_b].mean(axis=0)).ravel(),
        f"Prevalence_{b}": np.asarray(present[is_b].mean(axis=0)).ravel(),
        f"Mean_per_1000_genes_{a}": mean_a * 1000,
        f"Mean_per_1000_genes_{b}": mean_b * 1000,
        f"log2FC_{b}_vs_{a}": np.log2((mean_b + PSEUDOCOUNT) / (mean_a + PSEUDOCOUNT)),
        f"Rank_biserial_{b}_vs_{a}": 2 * u / (n_b * n_a) - 1,
        "P_value": p,
    })
    res["Q_value_BH"] = false_discovery_control(res["P_value"].values, method="bh")
    res["Significant_q<0.05"] = res["Q_value_BH"] < ALPHA
    res["Higher_in"] = np.where(res[f"Rank_biserial_{b}_vs_{a}"] > 0, b, a)
    return res.sort_values("Q_value_BH")


def cluster_top_cogs(ra, present, feats, desc, labels):
    is_cog = np.char.startswith(feats.astype(str), "COG")
    ra, present, feats = ra[:, is_cog], present[:, is_cog], feats[is_cog]
    clustered = labels >= 0
    out = []
    for cl in sorted(set(labels[clustered])):
        inn, rest = labels == cl, clustered & (labels != cl)
        m_in = np.asarray(ra[inn].mean(axis=0)).ravel()
        m_rest = np.asarray(ra[rest].mean(axis=0)).ravel() if rest.any() else np.zeros(len(feats))
        p_in = np.asarray(present[inn].mean(axis=0)).ravel()
        p_rest = np.asarray(present[rest].mean(axis=0)).ravel() if rest.any() else np.zeros(len(feats))
        fc = np.log2((m_in + PSEUDOCOUNT) / (m_rest + PSEUDOCOUNT))
        eligible = np.where(p_in >= MIN_CLUSTER_PREVALENCE)[0]
        order = eligible[np.argsort(-fc[eligible])][:TOP_N_SAVED]
        for rank, i in enumerate(order, 1):
            out.append({
                "Cluster": cl, "Rank": rank, "COG": feats[i],
                "Description": desc.get(feats[i], ""),
                "log2_fold_enrichment": fc[i],
                "Prevalence_in_cluster": p_in[i], "Prevalence_in_other_clusters": p_rest[i],
                "Mean_per_1000_genes_in_cluster": m_in[i] * 1000,
                "Mean_per_1000_genes_in_other_clusters": m_rest[i] * 1000,
            })
    return pd.DataFrame(out)


def cluster_summary(labels, microbiome, names):
    s = pd.crosstab(labels, microbiome).reindex(columns=names, fill_value=0)
    s.index.name = "Cluster"
    s["Total"] = s.sum(axis=1)
    return s


def mini_maps(coords, img_dir, prefix):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    img_dir.mkdir(parents=True, exist_ok=True)
    themes = {"light": ("#c3c2b7", "#0b0b0b"), "dark": ("#55544f", "#ffffff")}
    x, y, lab = coords["UMAP_1"].values, coords["UMAP_2"].values, coords["Cluster"].values
    for cl in sorted(set(lab)):
        for theme, (bg, fg) in themes.items():
            fig, ax = plt.subplots(figsize=(3, 3), dpi=160)
            ax.scatter(x[lab != cl], y[lab != cl], s=1.2, c=bg, linewidths=0)
            ax.scatter(x[lab == cl], y[lab == cl], s=2.2, c=fg, linewidths=0)
            ax.set_axis_off()
            ax.set_aspect("equal", adjustable="datalim")
            fig.savefig(img_dir / f"{prefix}_{cl}_{theme}.png", transparent=True, bbox_inches="tight", pad_inches=0.02)
            plt.close(fig)


def run_view(view, store, gene_counts):
    print(f"\n=========== View: {view['title']} ===========")
    out = SITE / "views" / view["id"]
    data = out / "data"
    data.mkdir(parents=True, exist_ok=True)
    work = WORK / view["id"]
    work.mkdir(parents=True, exist_ok=True)
    members = view["members"]
    names = [m["name"] for m in members]

    # --- 1. DNA-embedding map --------------------------------------------------
    print("-> DNA-embedding map")
    ann = work / "annotations_gene_count.csv"
    gene_counts.to_csv(ann, index_label="Genome_ID")
    dna, dna_report = cluster_with_tool(view, "dna", [(m["name"], m["embeddings"]) for m in members],
                                        [], ann)
    ids = dna.index

    # --- 2. Functional map on eggNOG relative abundance --------------------------
    print("-> Functional map")
    ra, present, feats = store.relative_abundance(ids)
    inputs = []
    for m in members:
        rows = np.where((dna["Dataset"] == m["name"]).values)[0]
        path = work / f"relative_abundance_{m['id']}.csv"
        pd.DataFrame(ra[rows].toarray().astype(np.float32), index=ids[rows], columns=feats) \
            .to_csv(path, float_format="%.6g")
        inputs.append((m["name"], path))
    ann_f = work / "annotations_for_functional.csv"
    pd.DataFrame({"Gene_count": gene_counts.reindex(ids).values,
                  "DNA_embedding_cluster": dna["Cluster"].values}, index=ids).to_csv(ann_f, index_label="Genome_ID")
    fun, fun_report = cluster_with_tool(view, "functional", inputs, ["--metric", FUNCTIONAL_METRIC], ann_f)
    fun = fun.reindex(ids)

    # --- 3. Maps, each also coloured by the other clustering --------------------
    ann_d = work / "annotations_for_dna.csv"
    pd.DataFrame({"Gene_count": gene_counts.reindex(ids).values,
                  "Functional_cluster": fun["Cluster"].values}, index=ids).to_csv(ann_d, index_label="Genome_ID")
    make_maps(view, "dna", WORK / view["id"] / "dna" / "coordinates_and_clusters.csv", ann_d,
              f"{view['title']}: DNA-embedding map (DNABERT-S)", "Functional_cluster")
    make_maps(view, "functional", WORK / view["id"] / "functional" / "coordinates_and_clusters.csv", ann_f,
              f"{view['title']}: functional map (eggNOG)", "DNA_embedding_cluster")

    # --- 4. Tables ----------------------------------------------------------------
    genomes = pd.DataFrame({
        "Microbiome": dna["Dataset"], "Gene_count": gene_counts.reindex(ids).values,
        "DNA_UMAP_1": dna["UMAP_1"], "DNA_UMAP_2": dna["UMAP_2"], "DNA_cluster": dna["Cluster"],
        "Functional_UMAP_1": fun["UMAP_1"], "Functional_UMAP_2": fun["UMAP_2"], "Functional_cluster": fun["Cluster"],
    }, index=ids)
    genomes.to_csv(data / "genomes.csv", index_label="Genome_ID")
    (data / "clustering_report_dna.txt").write_text(dna_report, encoding="utf-8")
    (data / "clustering_report_functional.txt").write_text(fun_report, encoding="utf-8")
    cluster_summary(genomes["DNA_cluster"], genomes["Microbiome"], names).to_csv(data / "cluster_summary_dna.csv")
    cluster_summary(genomes["Functional_cluster"], genomes["Microbiome"], names) \
        .to_csv(data / "cluster_summary_functional.csv")

    top = cluster_top_cogs(ra, present, feats, store.desc, genomes["Functional_cluster"].values)
    top.to_csv(data / "functional_cluster_top_cog_terms.csv", index=False)

    both = (genomes["DNA_cluster"] >= 0) & (genomes["Functional_cluster"] >= 0)
    settings = {
        "n_genomes": int(len(genomes)), "n_ogs_filtered": int(len(feats)),
        "min_prevalence_tested": MIN_PREVALENCE_TESTED, "min_genomes_tested": MIN_GENOMES_TESTED,
        "functional_metric": FUNCTIONAL_METRIC, "min_cluster_prevalence": MIN_CLUSTER_PREVALENCE, "alpha": ALPHA,
        "min_samples_dna": MIN_SAMPLES["dna"], "min_samples_functional": MIN_SAMPLES["functional"],
        "dna_vs_functional_ari": float(adjusted_rand_score(genomes.loc[both, "DNA_cluster"],
                                                           genomes.loc[both, "Functional_cluster"])),
        "dna_vs_functional_share": float(both.mean()),
    }
    if len(names) == 2:
        a, b = names
        res = pair_test(ra, present, feats, store.desc, (genomes["Microbiome"] == b).values, a, b)
        res.to_csv(data / "comparison_mannwhitney.csv", index=False)
        settings["n_tested"] = int(len(res))
        print(f"-> {res['Significant_q<0.05'].sum()} of {len(res)} OGs significant at q < {ALPHA}")
    (data / "analysis_settings.json").write_text(json.dumps(settings, indent=2))

    mini_maps(fun, out / "img", "functional_cluster")
    print(f"-> DNA clusters: {genomes['DNA_cluster'].max() + 1}, functional clusters: "
          f"{genomes['Functional_cluster'].max() + 1}, agreement ARI {settings['dna_vs_functional_ari']:.2f}")


def main():
    cfg = load_config()
    all_views = views(cfg)
    wanted = sys.argv[1:]
    if wanted:
        unknown = set(wanted) - {v["id"] for v in all_views}
        if unknown:
            raise SystemExit(f"Unknown view(s) {sorted(unknown)}. Available: {[v['id'] for v in all_views]}")
        all_views = [v for v in all_views if v["id"] in wanted]

    WORK.mkdir(parents=True, exist_ok=True)
    store = EggnogStore(cfg["eggnog_features_dir"])
    for view in all_views:
        run_view(view, store, store.gene_counts())
    print("\nDone. Now run: python build_site.py")


if __name__ == "__main__":
    main()
