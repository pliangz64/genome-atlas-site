"""Generate the static atlas website in docs/ from the per-view data written by analysis.py.

Usage:
    python build_site.py
"""
import html
import json
import re
from datetime import date

import numpy as np
import pandas as pd

from atlas_config import SITE, load_config, views

SITE_TITLE = "Human Microbiome Genome Atlas"
TOOL_REPO_URL = ""       # e.g. "https://github.com/<user>/interactive-umap-hdbscan"; hidden if empty
EFFECT_MIN = 0.3         # |rank-biserial| needed to count as a notable difference (with q < 0.05)
TOP_TABLE_ROWS = 25
TOP_COG_DISPLAY = 10     # COG families shown per functional cluster (up to 50 are saved in the CSV)
PLOTLY_CDN = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"

BASE_CSS = """
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --surface-2: #f3f2ee;
  --ink: #0b0b0b; --ink-2: #52514e; --muted: #6b6a66;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --neutral: #b9b8b1; --accent: #2a78d6;
  __LIGHT_MICROBIOMES__
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --surface-2: #222220;
    --ink: #ffffff; --ink-2: #c3c2b7; --muted: #9a9892;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --neutral: #55544f; --accent: #6da7ec;
    __DARK_MICROBIOMES__
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --surface-2: #222220;
  --ink: #ffffff; --ink-2: #c3c2b7; --muted: #9a9892;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --neutral: #55544f; --accent: #6da7ec;
  __DARK_MICROBIOMES__
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 16px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; }
a { color: var(--accent); }
.wrap { max-width: 1120px; margin: 0 auto; padding: 0 16px; }
header.site { background: var(--surface); border-bottom: 1px solid var(--border); }
header.site .wrap { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 24px; padding-top: 14px; padding-bottom: 14px; }
.brand { font-weight: 650; color: var(--ink); text-decoration: none; font-size: 17px; }
nav.main { display: flex; flex-wrap: wrap; gap: 4px; }
nav.main a { color: var(--ink-2); text-decoration: none; padding: 6px 10px; border-radius: 6px; font-size: 15px; }
nav.main a:hover { background: var(--surface-2); color: var(--ink); }
nav.main a[aria-current="page"] { background: var(--surface-2); color: var(--ink); font-weight: 600; }
main { padding: 28px 0 56px; }
h1 { font-size: 30px; line-height: 1.2; margin: 0 0 10px; letter-spacing: -0.01em; }
h2 { font-size: 22px; margin: 44px 0 10px; scroll-margin-top: 16px; }
h3 { font-size: 17px; margin: 0; }
.lede { font-size: 18px; color: var(--ink-2); max-width: 760px; margin: 0 0 20px; }
p { max-width: 780px; }
.muted { color: var(--muted); }
.small { font-size: 14px; }
.selector { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 14px 16px; margin: 0 0 28px; }
.selector form { display: flex; flex-wrap: wrap; align-items: center; gap: 10px 14px; }
.selector .label { font-weight: 600; margin-right: 4px; }
.chip { display: inline-flex; align-items: center; gap: 8px; padding: 6px 12px; border: 1px solid var(--border);
  border-radius: 999px; cursor: pointer; background: var(--surface); user-select: none; }
.chip input { accent-color: var(--ink); margin: 0; }
.chip:has(input:checked) { border-color: var(--ink-2); background: var(--surface-2); }
.selector button { font: inherit; font-weight: 600; padding: 7px 14px; border-radius: 8px; border: 0;
  background: var(--ink); color: var(--surface); cursor: pointer; }
.selector button:disabled { opacity: 0.4; cursor: not-allowed; }
.selector .hint { font-size: 14px; color: var(--muted); }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 18px 0; }
.tile { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 14px 16px; }
.tile .v { font-size: 28px; font-weight: 650; line-height: 1.15; }
.tile .k { color: var(--ink-2); font-size: 14px; margin-top: 2px; }
.swatch { display: inline-block; width: 10px; height: 10px; border-radius: 50%; margin-right: 6px; vertical-align: 0; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 12px; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px; }
a.card { color: inherit; text-decoration: none; display: block; }
a.card:hover { border-color: var(--axis); }
a.card h3 { color: var(--accent); margin-bottom: 4px; }
.subnav { display: flex; flex-wrap: wrap; gap: 6px; margin: 4px 0 8px; }
.subnav a { font-size: 15px; padding: 5px 10px; border-radius: 6px; background: var(--surface-2); color: var(--ink); text-decoration: none; }
.panel { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px; margin: 16px 0; }
.note { background: var(--surface-2); border-radius: 8px; padding: 12px 14px; font-size: 15px; color: var(--ink-2); max-width: 820px; }
.table-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 10px; background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: 14px; }
th, td { padding: 8px 10px; text-align: left; vertical-align: top; border-bottom: 1px solid var(--grid); }
th { color: var(--ink-2); font-weight: 600; background: var(--surface-2); white-space: nowrap; }
tr:last-child td { border-bottom: 0; }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.id { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; font-size: 13px; white-space: nowrap; }
.tabs { display: flex; flex-wrap: wrap; gap: 6px; margin: 8px 0 12px; }
.tabs button { font: inherit; font-size: 15px; padding: 7px 12px; border-radius: 8px; border: 1px solid var(--border);
  background: var(--surface); color: var(--ink-2); cursor: pointer; }
.tabs button[aria-selected="true"] { background: var(--ink); color: var(--surface); border-color: var(--ink); }
iframe.map { width: 100%; height: 80vh; min-height: 520px; max-height: 900px; border: 1px solid var(--border); border-radius: 10px; background: #fff; }
.chart { width: 100%; height: 520px; }
.legend { display: flex; flex-wrap: wrap; gap: 16px; font-size: 14px; color: var(--ink-2); margin: 4px 0 0; }
.cluster { display: grid; grid-template-columns: 200px minmax(0, 1fr); gap: 20px; align-items: start; }
.cluster img { width: 100%; height: auto; display: block; }
.bar { display: flex; gap: 2px; height: 10px; margin: 8px 0 4px; max-width: 420px; }
.bar span { display: block; height: 100%; }
.bar span:first-child { border-radius: 4px 0 0 4px; }
.bar span:last-child { border-radius: 0 4px 4px 0; }
.bar span:only-child { border-radius: 4px; }
footer { border-top: 1px solid var(--border); padding: 20px 0 32px; color: var(--muted); font-size: 14px; }
@media (max-width: 720px) {
  h1 { font-size: 25px; }
  .cluster { grid-template-columns: minmax(0, 1fr); }
  .cluster img { max-width: 220px; }
  .chart { height: 420px; }
}
"""

SELECTOR_JS = """
document.querySelectorAll('form.pick').forEach(form => {
  const order = JSON.parse(form.dataset.order);
  const root = form.dataset.root;
  const boxes = [...form.querySelectorAll('input[type=checkbox]')];
  const btn = form.querySelector('button');
  const hint = form.querySelector('.hint');
  let picked = boxes.filter(b => b.checked).map(b => b.value);
  function update() {
    btn.disabled = picked.length === 0;
    hint.textContent = picked.length === 0 ? 'Pick one microbiome to explore, or two to compare.'
      : picked.length === 1 ? 'Add a second microbiome to compare, or show this one alone.' : '';
  }
  boxes.forEach(b => b.addEventListener('change', () => {
    if (b.checked) {
      picked.push(b.value);
      if (picked.length > 2) {   // keep the two most recent choices
        const drop = picked.shift();
        boxes.find(x => x.value === drop).checked = false;
      }
    } else {
      picked = picked.filter(v => v !== b.value);
    }
    update();
  }));
  form.addEventListener('submit', e => {
    e.preventDefault();
    const id = [...picked].sort((a, b) => order.indexOf(a) - order.indexOf(b)).join('--');
    location.href = root + 'views/' + id + '/index.html';
  });
  update();
});
"""


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def esc(s):
    return html.escape("" if pd.isna(s) else str(s))


def pct(x):
    return f"{x:.0%}" if x >= 0.01 or x == 0 else "<1%"


def fmt_q(q):
    return "<1e-300" if q == 0 else (f"{q:.1e}" if q < 0.001 else f"{q:.3f}")


def desc_cell(s, limit=110):
    """Description table cell, shortened with the full text available on hover."""
    if pd.isna(s) or not str(s).strip():
        return "<td><span class='muted'>—</span></td>"
    s = str(s)
    if len(s) <= limit:
        return f"<td>{esc(s)}</td>"
    return f"<td title=\"{esc(s)}\">{esc(s[:limit].rsplit(' ', 1)[0])}…</td>"


def sw(m):
    return f'<span class="swatch" style="background:var(--m-{m["id"]})"></span>'


def tile(value, label, swatch=""):
    return f'<div class="tile"><div class="v">{value}</div><div class="k">{swatch}{label}</div></div>'


def page(path, title, body, root, current=None, extra_head=""):
    nav_items = [("index.html", "Home"), ("methods.html", "Methods"), ("downloads.html", "Downloads")]
    nav = "".join(f'<a href="{root}{href}"{" aria-current=\"page\"" if href == current else ""}>{label}</a>'
                  for href, label in nav_items)
    tool = f' · Built with the <a href="{TOOL_REPO_URL}">interactive UMAP + HDBSCAN tool</a>' if TOOL_REPO_URL else ""
    doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)} · {esc(SITE_TITLE)}</title>
<link rel="stylesheet" href="{root}assets/style.css">
{extra_head}
</head>
<body>
<header class="site"><div class="wrap">
  <a class="brand" href="{root}index.html">{esc(SITE_TITLE)}</a>
  <nav class="main" aria-label="Main">{nav}</nav>
</div></header>
<main><div class="wrap">
{body}
</div></main>
<footer><div class="wrap">Early prototype · generated {date.today():%d %B %Y}{tool}</div></footer>
<script src="{root}assets/selector.js"></script>
</body>
</html>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(doc, encoding="utf-8")


def selector(cfg, root, selected=()):
    chips = "".join(
        f'<label class="chip"><input type="checkbox" value="{m["id"]}"{" checked" if m["id"] in selected else ""}>'
        f'{sw(m)}{esc(m["name"])}</label>' for m in cfg["microbiomes"])
    order = json.dumps([m["id"] for m in cfg["microbiomes"]])
    return f"""<div class="selector">
<form class="pick" data-order='{order}' data-root="{root}">
  <span class="label">Microbiomes</span>{chips}
  <button type="submit">Show</button>
  <span class="hint"></span>
</form></div>"""


def stability_sentence(report):
    m = re.search(r"Stability: (.*)\n\s+(\d+) runs gave (\S+) clusters; .*?mean ([\d.]+), lowest ([\d.]+)", report)
    if not m:
        return ""
    verdict, runs, rng, mean, low = m.groups()
    return (f"{verdict} Across {runs} runs with different random seeds, the map gave {rng} clusters "
            f"(mean agreement {mean}, lowest {low}; Adjusted Rand Index).")


def load_view(view):
    d = SITE / "views" / view["id"] / "data"
    out = {
        "genomes": pd.read_csv(d / "genomes.csv", index_col=0),
        "top": pd.read_csv(d / "functional_cluster_top_cog_terms.csv"),
        "sum_dna": pd.read_csv(d / "cluster_summary_dna.csv", index_col=0),
        "sum_fun": pd.read_csv(d / "cluster_summary_functional.csv", index_col=0),
        "settings": json.loads((d / "analysis_settings.json").read_text()),
        "report_dna": (d / "clustering_report_dna.txt").read_text(encoding="utf-8"),
        "report_fun": (d / "clustering_report_functional.txt").read_text(encoding="utf-8"),
    }
    if len(view["members"]) == 2:
        res = pd.read_csv(d / "comparison_mannwhitney.csv")
        a, b = (m["name"] for m in view["members"])
        res["Effect"] = res[f"Rank_biserial_{b}_vs_{a}"]
        res["Notable"] = res["Significant_q<0.05"] & (res["Effect"].abs() >= EFFECT_MIN)
        out["res"] = res
    return out


# -----------------------------------------------------------------------------
# View page sections
# -----------------------------------------------------------------------------

def map_block(view, kind, other_label):
    """Tabbed interactive map for one clustering ('dna' or 'functional')."""
    tabs = [(f"maps/{kind}_cluster.html", "Colour by cluster")]
    if len(view["members"]) == 2:
        tabs.append((f"maps/{kind}_microbiome.html", "Colour by microbiome"))
    tabs.append((f"maps/{kind}_other.html", f"Colour by {other_label}"))
    buttons = "".join(f'<button role="tab" aria-selected="{str(i == 0).lower()}" data-src="{src}">{label}</button>'
                      for i, (src, label) in enumerate(tabs))
    overlay = " The smaller microbiome is drawn on top with a black outline." if len(view["members"]) == 2 else ""
    return f"""
<div class="tabs" role="tablist" aria-label="Colour the map by" data-frame="frame-{kind}">{buttons}</div>
<iframe class="map" id="frame-{kind}" src="{tabs[0][0]}" title="Interactive {kind} map" loading="lazy"></iframe>
<p class="small muted">Hover a point for its details. Click a legend entry to hide it; double-click to show only that
entry (double-click again to show all). Drag to zoom; double-click the plot to reset.{overlay} Grey points are outliers
that HDBSCAN did not assign to any cluster. <a href="{tabs[0][0]}" target="_blank" rel="noopener">Open full screen</a>.</p>"""


TABS_JS = """<script>
document.querySelectorAll('.tabs').forEach(group => group.querySelectorAll('button').forEach(b =>
  b.addEventListener('click', () => {
    group.querySelectorAll('button').forEach(x => x.setAttribute('aria-selected', x === b));
    document.getElementById(group.dataset.frame).src = b.dataset.src;
  })));
</script>"""


def cluster_label(cl):
    return "Outliers" if int(cl) == -1 else f"Cluster {cl}"


def composition_table(view, summary):
    """Compact table: genomes per cluster (and per microbiome for pairs)."""
    members = view["members"]
    pair = len(members) == 2
    head = "<th>Cluster</th><th class='num'>Genomes</th>" + (
        "".join(f"<th class='num'>{sw(m)}{esc(m['name'])}</th>" for m in members) if pair else "")
    order = [c for c in summary.index if c >= 0] + ([-1] if -1 in summary.index else [])
    rows = []
    for cl in order:
        row = summary.loc[cl]
        cells = "".join(f"<td class='num'>{row[m['name']]:,} ({row[m['name']] / row['Total']:.0%})</td>"
                        for m in members) if pair else ""
        rows.append(f"<tr><td>{cluster_label(cl)}</td><td class='num'>{row['Total']:,}</td>{cells}</tr>")
    return f"""<div class="table-wrap" style="max-width:640px"><table><thead><tr>{head}</tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>"""


def dna_section(view, vd):
    return f"""
<h2 id="dna">DNA-embedding map (DNABERT-S)</h2>
<p>Genomes are placed by DNABERT-S embeddings of their DNA sequence, so genomes with similar sequence composition sit
close together. Treat it as a visual impression of genetic similarity, not a measure of evolutionary distance.</p>
{map_block(view, 'dna', 'functional cluster')}
<p>{esc(stability_sentence(vd['report_dna']))}</p>
{composition_table(view, vd['sum_dna'])}"""


def functional_section(view, vd):
    members = view["members"]
    top, summary, settings = vd["top"], vd["sum_fun"], vd["settings"]
    cards = []
    for cl, row in summary.drop(index=-1, errors="ignore").iterrows():
        n = row["Total"]
        comp = ""
        if len(members) == 2:
            segs = "".join(f'<span style="width:{row[m["name"]] / n * 100:.2f}%;background:var(--m-{m["id"]})"></span>'
                           for m in members if row[m["name"]] > 0)
            text = " &nbsp; ".join(f'{sw(m)}{row[m["name"]]:,} {esc(m["name"].lower())} ({row[m["name"]] / n:.0%})'
                                   for m in members)
            comp = f'<div class="bar" aria-hidden="true">{segs}</div><div class="small muted">{text}</div>'
        t5 = top[(top["Cluster"] == cl) & (top["Rank"] <= TOP_COG_DISPLAY)]
        rows = "".join(
            f"<tr><td class='id'>{esc(r.COG)}</td>{desc_cell(r.Description)}"
            f"<td class='num'>{2 ** r.log2_fold_enrichment:,.0f}×</td>"
            f"<td class='num'>{pct(r.Prevalence_in_cluster)}</td><td class='num'>{pct(r.Prevalence_in_other_clusters)}</td></tr>"
            for r in t5.itertuples())
        table = (f"""<div class="table-wrap" style="margin-top:12px"><table>
      <thead><tr><th>COG</th><th>Description</th><th class="num">Enrichment</th><th class="num">In this<br>cluster</th><th class="num">In other<br>clusters</th></tr></thead>
      <tbody>{rows}</tbody></table></div>""" if rows else
                 '<p class="small muted">No COG family meets the prevalence threshold.</p>')
        cards.append(f"""
<section class="panel cluster" id="cluster-{cl}" aria-labelledby="h-cluster-{cl}">
  <picture>
    <source srcset="img/functional_cluster_{cl}_dark.png" media="(prefers-color-scheme: dark)">
    <img src="img/functional_cluster_{cl}_light.png" alt="Functional map with cluster {cl} highlighted" loading="lazy">
  </picture>
  <div>
    <h3 id="h-cluster-{cl}">Cluster {cl} · {n:,} genomes</h3>
    {comp}
    {table}
  </div>
</section>""")

    outliers = int(summary.loc[-1, "Total"]) if -1 in summary.index else 0
    single = len(cards) == 1
    return f"""
<h2 id="functional">Functional map (eggNOG gene families)</h2>
<p>Here genomes are placed by what their genes do. Each genome is described by the relative abundance of
{settings['n_ogs_filtered']:,} eggNOG gene families (each family's share of the genome's genes), and genomes with similar
profiles (Bray–Curtis dissimilarity) sit close together.</p>
{map_block(view, 'functional', 'DNA-embedding cluster')}
<p>{esc(stability_sentence(vd['report_fun']))} Agreement between these functional clusters and the DNA-embedding clusters:
Adjusted Rand Index {settings['dna_vs_functional_ari']:.2f} (1 = identical groupings, 0 = unrelated), computed on the
{settings['dna_vs_functional_share']:.0%} of genomes that are in a cluster on both maps.</p>
<h3 style="margin-top:28px">Top {TOP_COG_DISPLAY} COG gene families in each functional cluster</h3>
<p class="note"><strong>How they are chosen.</strong> Enrichment is the mean relative abundance of a COG family in the
cluster's genomes divided by its mean in all other clustered genomes of this view. Only families present in at least
{settings['min_cluster_prevalence']:.0%} of the cluster's genomes are ranked, so one-off genes can't top the list.
"In this cluster" is the share of the cluster's genomes carrying the family. The ranking is descriptive; no per-cluster
test is applied. {outliers:,} outlier genomes are excluded. The top 50 per cluster are in
<a href="data/functional_cluster_top_cog_terms.csv" download>functional_cluster_top_cog_terms.csv</a>.
{'With a single cluster there is nothing to compare against, so enrichment is not meaningful here.' if single else ''}</p>
{''.join(cards)}"""


def og_table(df, a, b):
    rows = "".join(
        f"<tr><td class='id'>{esc(r['OG'])}</td>{desc_cell(r['Description'])}"
        f"<td class='num'>{pct(r[f'Prevalence_{a}'])}</td><td class='num'>{pct(r[f'Prevalence_{b}'])}</td>"
        f"<td class='num'>{r[f'Mean_per_1000_genes_{a}']:.2f}</td><td class='num'>{r[f'Mean_per_1000_genes_{b}']:.2f}</td>"
        f"<td class='num'>{r['Effect']:+.2f}</td><td class='num'>{fmt_q(r['Q_value_BH'])}</td></tr>"
        for _, r in df.iterrows())
    if not rows:
        return '<p class="small muted">No gene family passes both thresholds.</p>'
    return f"""<div class="table-wrap"><table>
<thead><tr><th>Gene family</th><th>Description</th><th class="num">In {esc(a.lower())}<br>genomes</th><th class="num">In {esc(b.lower())}<br>genomes</th>
<th class="num">Per 1,000 genes<br>({esc(a.lower())})</th><th class="num">Per 1,000 genes<br>({esc(b.lower())})</th><th class="num">Effect<br>size</th><th class="num">q-value</th></tr></thead>
<tbody>{rows}</tbody></table></div>"""


def comparison_section(view, vd):
    ma, mb = view["members"]
    a, b = ma["name"], mb["name"]
    res, settings = vd["res"], vd["settings"]
    notable = res[res["Notable"]]
    up_a = notable[notable["Higher_in"] == a].sort_values("Effect").head(TOP_TABLE_ROWS)
    up_b = notable[notable["Higher_in"] == b].sort_values("Effect", ascending=False).head(TOP_TABLE_ROWS)

    pts = [[r.OG, (r.Description or "")[:90] if isinstance(r.Description, str) else "",
            round(r.Effect, 4), round(-np.log10(max(r.Q_value_BH, 1e-300)), 3),
            round(getattr(r, f"Prevalence_{a}"), 3), round(getattr(r, f"Prevalence_{b}"), 3),
            ("b" if r.Effect > 0 else "a") if r.Notable else "n"]
           for r in res.itertuples()]
    data_json = json.dumps(pts, separators=(",", ":"))

    return f"""
<h2 id="differences">What differs between {esc(a.lower())} and {esc(b.lower())} genomes?</h2>
<p>For each eggNOG gene family, a two-sided Mann-Whitney U test compares its relative abundance (its share of a genome's
genes) between the two microbiomes. P-values are corrected for all {settings['n_tested']:,} tests with the
Benjamini–Hochberg false discovery rate.</p>
<div class="tiles">
  {tile(f"{settings['n_tested']:,}", "Gene families tested")}
  {tile(f"{res['Significant_q<0.05'].sum():,}", "Significant (q &lt; 0.05)")}
  {tile(f"{(notable['Higher_in'] == a).sum():,}", f"Clearly higher in {esc(a.lower())} (effect ≥ {EFFECT_MIN})", sw(ma))}
  {tile(f"{(notable['Higher_in'] == b).sum():,}", f"Clearly higher in {esc(b.lower())} (effect ≥ {EFFECT_MIN})", sw(mb))}
</div>
<div class="panel">
  <h3>Every gene family tested</h3>
  <div id="volcano" class="chart" role="img" aria-label="Volcano plot of effect size against significance for all tested gene families"></div>
  <div class="legend">
    <span>{sw(ma)}Higher in {esc(a.lower())}</span><span>{sw(mb)}Higher in {esc(b.lower())}</span>
    <span><span class="swatch" style="background:var(--neutral)"></span>No clear difference</span>
  </div>
</div>
<p class="note"><strong>How to read this.</strong> The effect size (rank-biserial correlation, −1 to +1) says how
consistently one microbiome ranks higher: +1 means every {esc(b.lower())} genome has more of the family than every
{esc(a.lower())} genome, −1 the reverse, 0 no tendency. With many genomes even tiny differences reach significance, so
families are highlighted only with q &lt; 0.05 <em>and</em> an effect of at least {EFFECT_MIN}. "In {esc(a.lower())}
genomes" is the share of {esc(a.lower())} genomes carrying the family at all.</p>
<h3 style="margin-top:28px">{sw(ma)}Top {TOP_TABLE_ROWS} gene families higher in {esc(a.lower())} genomes</h3>
<p class="small muted">Sorted by effect size. All results: <a href="data/comparison_mannwhitney.csv" download>comparison_mannwhitney.csv</a>.</p>
{og_table(up_a, a, b)}
<h3 style="margin-top:28px">{sw(mb)}Top {TOP_TABLE_ROWS} gene families higher in {esc(b.lower())} genomes</h3>
<p class="small muted">Sorted by effect size.</p>
{og_table(up_b, a, b)}
<script>
const PTS = {data_json};
const EFFECT = {EFFECT_MIN};
function css(v) {{ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }}
function draw() {{
  const groups = {{n: ['No clear difference', '--neutral'], a: ['Higher in {esc(a.lower())}', '--m-{ma["id"]}'],
                  b: ['Higher in {esc(b.lower())}', '--m-{mb["id"]}']}};
  const traces = Object.entries(groups).map(([k, [name, col]]) => {{
    const p = PTS.filter(d => d[6] === k);
    return {{
      type: 'scattergl', mode: 'markers', name, showlegend: false,
      x: p.map(d => d[2]), y: p.map(d => d[3]),
      customdata: p.map(d => [d[0], d[1] || '—', (d[4] * 100).toFixed(1), (d[5] * 100).toFixed(1)]),
      hovertemplate: '<b>%{{customdata[0]}}</b><br>%{{customdata[1]}}<br>Effect size %{{x:+.2f}} · −log10 q %{{y:.1f}}' +
                     '<br>In %{{customdata[2]}}% of {esc(a.lower())}, %{{customdata[3]}}% of {esc(b.lower())} genomes<extra></extra>',
      marker: {{ size: k === 'n' ? 5 : 7, color: css(col), opacity: k === 'n' ? 0.6 : 0.85,
                line: {{ width: 1, color: css('--surface') }} }},
    }};
  }});
  const ink = css('--ink-2'), grid = css('--grid'), axis = css('--axis');
  const ymax = Math.max(...PTS.map(d => d[3])) * 1.05;
  const vline = x => ({{ type: 'line', x0: x, x1: x, y0: 0, y1: ymax, xref: 'x', yref: 'y', line: {{ color: axis, width: 1 }} }});
  Plotly.react('volcano', traces, {{
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: {{ family: 'system-ui, -apple-system, Segoe UI, sans-serif', color: ink, size: 13 }},
    margin: {{ l: 56, r: 12, t: 12, b: 48 }}, hovermode: 'closest',
    xaxis: {{ title: 'Effect size (← higher in {esc(a.lower())} · higher in {esc(b.lower())} →)', range: [-1.05, 1.05],
             gridcolor: grid, zeroline: true, zerolinecolor: axis }},
    yaxis: {{ title: '−log10 q-value', gridcolor: grid, zeroline: false, rangemode: 'tozero' }},
    shapes: [vline(-EFFECT), vline(EFFECT),
             {{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: -Math.log10(0.05), y1: -Math.log10(0.05), line: {{ color: axis, width: 1 }} }}],
  }}, {{ responsive: true, displaylogo: false }});
}}
draw();
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw);
</script>"""


def build_view(cfg, view, vd):
    members = view["members"]
    genomes = vd["genomes"]
    counts = "".join(tile(f"{(genomes['Microbiome'] == m['name']).sum():,}", f"{esc(m['name'])} genomes", sw(m))
                     for m in members)

    def n_clusters(summary):
        return int((summary.index >= 0).sum())

    pair = len(members) == 2
    subnav = ('<a href="#dna">DNA-embedding map</a><a href="#functional">Functional map and COG families</a>'
              + ('<a href="#differences">Differences</a>' if pair else ""))
    who = (f"{esc(members[0]['name'])} and {esc(members[1]['name'].lower())} genomes" if pair
           else f"{esc(members[0]['name'])} genomes ({esc(members[0]['catalogue'])})")
    lede = (f"{who} shown two ways: by their DNA sequence (DNABERT-S embeddings) and by their gene content (eggNOG gene "
            "families), each with its own clusters" + (", plus a test of which gene families differ between the two "
                                                        "microbiomes." if pair else "."))
    body = f"""
{selector(cfg, '../../', [m['id'] for m in members])}
<h1>{esc(view['title'])}</h1>
<p class="lede">{lede}</p>
<div class="tiles">
  {counts}
  {tile(n_clusters(vd['sum_dna']), "DNA-embedding clusters")}
  {tile(n_clusters(vd['sum_fun']), "Functional clusters")}
</div>
<nav class="subnav" aria-label="On this page">{subnav}</nav>
{dna_section(view, vd)}
{functional_section(view, vd)}
{comparison_section(view, vd) if pair else ''}
{TABS_JS}
"""
    head = f'<script src="{PLOTLY_CDN}"></script>' if pair else ""
    page(SITE / "views" / view["id"] / "index.html", view["title"], body, "../../", extra_head=head)


# -----------------------------------------------------------------------------
# Top-level pages
# -----------------------------------------------------------------------------

def build_home(cfg, all_views, data):
    ms = cfg["microbiomes"]
    cards = []
    for m in ms:
        g = data[m["id"]]["genomes"]
        cards.append(f'<a class="card" href="views/{m["id"]}/index.html"><h3>{sw(m)}{esc(m["name"])}</h3>'
                     f'<span class="small muted">{len(g):,} genomes · {esc(m["catalogue"])}</span></a>')
    pairs = [v for v in all_views if len(v["members"]) == 2]
    pair_cards = "".join(
        f'<a class="card" href="views/{v["id"]}/index.html"><h3>{esc(v["title"])}</h3>'
        f'<span class="small muted">{len(data[v["id"]]["genomes"]):,} genomes · DNA and functional maps, differences</span></a>'
        for v in pairs)
    body = f"""
<h1>{esc(SITE_TITLE)}</h1>
<p class="lede">Compare bacterial genomes from different human body sites. Each view shows the genomes on two maps: one
by DNA sequence (DNABERT-S embeddings) and one by gene content (relative abundance of eggNOG gene families), each with its
own clusters. Pick one microbiome to explore it on its own, or two to compare them.</p>
{selector(cfg, '')}
<h2>Microbiomes</h2>
<div class="cards">{''.join(cards)}</div>
<h2>Comparisons</h2>
<div class="cards">{pair_cards}</div>
<p class="small muted" style="margin-top:20px">{len(ms)} microbiomes · {len(all_views)} views. More microbiomes will be
added over time.</p>
"""
    page(SITE / "index.html", "Home", body, "", current="index.html")


def build_methods(cfg, any_settings):
    rows = "".join(f"<tr><td>{sw(m)}{esc(m['name'])}</td><td>{esc(m['catalogue'])}</td><td>{esc(m['citation'])}</td></tr>"
                   for m in cfg["microbiomes"])
    s = any_settings
    body = f"""
<h1>Methods</h1>
<h2>Genomes</h2>
<div class="table-wrap"><table><thead><tr><th>Microbiome</th><th>Genome collection</th><th>Reference</th></tr></thead>
<tbody>{rows}</tbody></table></div>
<p>Each genome represents one species-level group from its collection, not one sample.</p>
<h2>Views</h2>
<p>The site offers every microbiome on its own and every pair of microbiomes. Each view is computed separately, so a
single microbiome's clusters reflect its own structure, and map positions and cluster numbers are not comparable
between views.</p>
<h2>Gene families and relative abundance</h2>
<p>Genes were annotated with eggNOG-mapper and assigned to their root-level eggNOG orthologous group (OG). The relative
abundance of an OG in a genome is the number of its genes in that OG divided by the genome's total annotated genes,
which removes the effect of genome size. In each view, only OGs present in at least {s['min_prevalence_tested']:.0%} of
the view's genomes (and at least {s['min_genomes_tested']} genomes) are used.</p>
<h2>The two maps</h2>
<p>Each view has two maps, both made with UMAP (n_neighbors = 30, min_dist = 0.1, seed 42) and clustered on the 2-D map
with HDBSCAN (min_cluster_size = 5% of the view's genomes). Genomes in sparse regions are not assigned to a cluster and
are shown as outliers.</p>
<ul>
<li><strong>DNA-embedding map:</strong> genome embeddings from the DNA language model DNABERT-S, Euclidean distance,
HDBSCAN min_samples = {s['min_samples_dna']}. It reflects similarity of DNA sequence composition and is not a
phylogeny.</li>
<li><strong>Functional map:</strong> each genome's relative-abundance profile across the filtered OGs,
{s['functional_metric'].replace('braycurtis', 'Bray–Curtis')} dissimilarity, HDBSCAN min_samples =
{s['min_samples_functional']}.</li>
</ul>
<p>The min_samples values were chosen by sweeping the setting and checking cluster stability across random seeds. On the
DNA-embedding maps, values of 15 or more sometimes merged most genomes into one or two clusters; the functional maps were
most stable at around 10.</p>
<p>Stability was checked by repeating UMAP + HDBSCAN with five random seeds and comparing runs with the Adjusted Rand
Index (ARI) on genomes clustered in both runs. The verdict cut-offs (≥ 0.9 consistent, 0.75–0.9 mostly consistent) are
rule-of-thumb guides, not an established standard. Agreement between the DNA-embedding and functional clusters is also
reported as an ARI, on genomes clustered on both maps.</p>
<h2>Top COG families per functional cluster</h2>
<p>Among OGs that are COG groups, each functional cluster's families were ranked by log2 fold enrichment: mean relative
abundance in the cluster divided by the mean across all other clustered genomes in the view (pseudocount
10<sup>−6</sup>). Only COGs present in at least {s['min_cluster_prevalence']:.0%} of the cluster's genomes were ranked.
The top {TOP_COG_DISPLAY} are shown. This is descriptive; no per-cluster test is applied.</p>
<h2>Comparing two microbiomes</h2>
<p>For each OG, relative abundance in the two microbiomes was compared with a two-sided Mann-Whitney U test, and p-values
were corrected across all tested OGs with the Benjamini–Hochberg false discovery rate. Effect size is the rank-biserial
correlation. The site highlights OGs with q &lt; 0.05 and |effect| ≥ {EFFECT_MIN}.</p>
<h2>Limitations</h2>
<ul>
<li>Genomes share evolutionary history, so they are not fully independent observations. Results describe differences
between genome collections, not causal effects of the body site.</li>
<li>Relative abundance is compositional: an increase in one family's share lowers the others'.</li>
<li>Differences between microbiomes can partly reflect how each collection was built (isolates vs metagenome-assembled
genomes, completeness, annotation depth).</li>
<li>UMAP and HDBSCAN results depend on random seeds and settings; each view reports its own stability.</li>
</ul>
<h2>Software and method references</h2>
<ul class="small">
<li>Zhou Z, et al. DNABERT-S: Pioneering species differentiation with species-aware DNA embeddings. arXiv:2402.08777 (2024).</li>
<li>Cantalapiedra CP, et al. eggNOG-mapper v2. <em>Mol Biol Evol</em> 38:5825–5829 (2021).</li>
<li>Huerta-Cepas J, et al. eggNOG 5.0. <em>Nucleic Acids Res</em> 47:D309–D314 (2019).</li>
<li>McInnes L, Healy J, Melville J. UMAP. arXiv:1802.03426 (2018).</li>
<li>Campello RJGB, Moulavi D, Sander J. Density-based clustering based on hierarchical density estimates. PAKDD (2013).</li>
<li>Benjamini Y, Hochberg Y. Controlling the false discovery rate. <em>J R Stat Soc B</em> 57:289–300 (1995).</li>
</ul>
"""
    page(SITE / "methods.html", "Methods", body, "", current="methods.html")


def build_downloads(all_views):
    desc = {
        "genomes.csv": "Every genome in the view: ID, microbiome, gene count, and coordinates and cluster on both maps.",
        "cluster_summary_dna.csv": "Genomes per microbiome in each DNA-embedding cluster.",
        "cluster_summary_functional.csv": "Genomes per microbiome in each functional cluster.",
        "functional_cluster_top_cog_terms.csv": "Top 50 COG families per functional cluster by fold enrichment, with prevalence.",
        "comparison_mannwhitney.csv": "Mann-Whitney U results for every tested eggNOG gene family.",
        "clustering_report_dna.txt": "DNA-embedding map: settings and stability report.",
        "clustering_report_functional.txt": "Functional map: settings and stability report.",
    }
    sections = []
    for v in all_views:
        d = SITE / "views" / v["id"] / "data"
        rows = "".join(
            f"<tr><td><a href='views/{v['id']}/data/{f}' download>{f}</a></td><td>{text}</td>"
            f"<td class='num'>{(d / f).stat().st_size / 1024:,.0f} KB</td></tr>"
            for f, text in desc.items() if (d / f).exists())
        sections.append(f'<h2>{esc(v["title"])}</h2><div class="table-wrap"><table><thead><tr><th>File</th>'
                        f'<th>Contents</th><th class="num">Size</th></tr></thead><tbody>{rows}</tbody></table></div>')
    body = f"""
<h1>Downloads</h1>
<p class="lede">All data behind each view, as plain CSV and text files. Relative abundance is reported per 1,000 genes;
prevalence is the share of genomes carrying at least one gene in the family. In the CSV files, outliers have cluster
<code>-1</code> so the column stays numeric.</p>
{''.join(sections)}
"""
    page(SITE / "downloads.html", "Downloads", body, "", current="downloads.html")


def main():
    cfg = load_config()
    all_views = views(cfg)
    missing = [v["id"] for v in all_views if not (SITE / "views" / v["id"] / "data" / "genomes.csv").exists()]
    if missing:
        raise SystemExit(f"No analysis output for view(s) {missing}. Run: python analysis.py {' '.join(missing)}")

    light = " ".join(f"--m-{m['id']}: {m['light']};" for m in cfg["microbiomes"])
    dark = " ".join(f"--m-{m['id']}: {m['dark']};" for m in cfg["microbiomes"])
    css = BASE_CSS.replace("__LIGHT_MICROBIOMES__", light).replace("__DARK_MICROBIOMES__", dark)
    (SITE / "assets").mkdir(parents=True, exist_ok=True)
    (SITE / "assets" / "style.css").write_text(css.strip() + "\n", encoding="utf-8")
    (SITE / "assets" / "selector.js").write_text(SELECTOR_JS.strip() + "\n", encoding="utf-8")
    (SITE / ".nojekyll").write_text("")

    data = {v["id"]: load_view(v) for v in all_views}
    for v in all_views:
        build_view(cfg, v, data[v["id"]])
    build_home(cfg, all_views, data)
    build_methods(cfg, next(iter(data.values()))["settings"])
    build_downloads(all_views)
    print(f"Built {len(all_views)} views. Open docs/index.html in a browser.")


if __name__ == "__main__":
    main()
