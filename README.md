# Human Microbiome Genome Atlas

Source for a static website that compares bacterial genomes from different human microbiomes by their DNABERT-S DNA embeddings and eggNOG gene-family content.

> **Status: early prototype.** Currently includes gut (UHGG) and vaginal (VMGC) genomes.

Visitors pick **one microbiome** to explore it on its own, or **two** to compare them. Every microbiome and every pair of microbiomes gets its own precomputed view. Each view shows the genomes on two maps, each clustered separately with UMAP + HDBSCAN:

| Map | Built from | Shows |
|---|---|---|
| **DNA-embedding map** | DNABERT-S genome embeddings (Euclidean distance) | A visual impression of sequence similarity; cluster sizes and, for pairs, their microbiome make-up |
| **Functional map** | Relative abundance of eggNOG gene families (Bray–Curtis dissimilarity) | Functional clusters with each cluster's **top 10 enriched COG families** |

Each map can also be coloured by the other map's clusters, and the page reports how well the two clusterings agree. Pair views add a **Mann-Whitney U test per eggNOG gene family** (relative abundance, Benjamini–Hochberg FDR) with a volcano plot and top-25 tables.

Outliers (genomes HDBSCAN did not assign to a cluster) are labelled "Outliers" on the site and have cluster `-1` in the CSV downloads.

## Layout

| Path | Purpose |
|---|---|
| `microbiomes.json` | **The list of microbiomes**: the only file to edit when adding data. |
| `atlas_config.py` | Reads the config; defines colours and the list of views. |
| `analysis.py` | Computes every view: runs the UMAP + HDBSCAN tool, relative abundance, statistics, mini-maps. Writes `docs/views/<view>/`. |
| `build_site.py` | Builds the HTML pages in `docs/` from those outputs. |
| `docs/` | The published website (GitHub Pages serves this folder). |
| `work/` | Intermediate files; not committed. |

## Adding a microbiome

1. **eggNOG annotations:** run eggNOG-mapper on the new genomes and add them to the shared feature store (`eggnog_features_dir` in `microbiomes.json`). With the current setup, add the new genome folder to `Random Forest\build_feature_matrix.py`, delete `features\file_list.csv`, and rerun it.
2. **Embeddings:** create a DNABERT-S embedding CSV for the new genomes (genome ID in the first column), using the same IDs as the eggNOG files.
3. **Config:** add an entry to `microbiomes.json`:
   ```json
   {
     "id": "oral",
     "name": "Oral",
     "catalogue": "Collection name and version",
     "citation": "Reference for the collection",
     "embeddings": "C:/path/to/oral_embeddings.csv",
     "color_slot": 3
   }
   ```
   `color_slot` (1–8) fixes the microbiome's colour across the whole site. Use the next unused slot.
4. **Rebuild:** existing views don't change, so compute only the new ones, then rebuild the pages:
   ```bash
   python analysis.py oral gut--oral vaginal--oral
   python build_site.py
   ```
   View ids are the microbiome id, or two ids joined by `--` in config order. `python build_site.py` lists any views that are missing.

## Rebuilding everything

```bash
python analysis.py      # a few minutes per view
python build_site.py    # seconds; rerun after editing page text or design
```

Requirements: the [interactive UMAP + HDBSCAN tool](../Interactive%20UMAP%20+%20HDBSCAN%20Tool/interactive-umap-hdbscan) and its requirements, plus `matplotlib`.

Settings are constants at the top of each script:

- **`analysis.py`:**
  - `MIN_PREVALENCE_TESTED`: a gene family must be present in at least this share of a view's genomes to be analysed (default 1%, and at least 5 genomes).
  - `MIN_CLUSTER_PREVALENCE`: a COG must be present in at least this share of a cluster's genomes to be ranked (default 10%).
- **`build_site.py`:**
  - `SITE_TITLE`: the site name.
  - `EFFECT_MIN`: the effect size needed to highlight a family (default 0.3).
  - `TOOL_REPO_URL`: link to the tool's GitHub repo in the footer.

## Publishing on GitHub Pages

1. Create a GitHub repository and upload this folder, excluding `work/`.
2. In the repository, go to **Settings → Pages**. Under **Build and deployment**, choose **Deploy from a branch**, select branch `main` and folder `/docs`, then click **Save**.
3. After a minute the site is live at `https://<username>.github.io/<repository>/`.

Each view adds about 2–3 MB: the maps load their plotting library from a CDN, so visitors need an internet connection. GitHub Pages sites can be up to 1 GB, so size won't limit how many microbiomes you add.
