"""Shared configuration for analysis.py and build_site.py: microbiomes, colours and views."""
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "docs"
WORK = ROOT / "work"

# Categorical palette (light, dark), in fixed slot order. A microbiome keeps its
# slot in every view, so its colour never changes when others are added.
PALETTE = {
    1: ("#2a78d6", "#3987e5"),  # blue
    2: ("#eb6834", "#d95926"),  # orange
    3: ("#1baf7a", "#199e70"),  # aqua
    4: ("#eda100", "#c98500"),  # yellow
    5: ("#e87ba4", "#d55181"),  # magenta
    6: ("#008300", "#008300"),  # green
    7: ("#4a3aa7", "#9085e9"),  # violet
    8: ("#e34948", "#e66767"),  # red
}


def load_config():
    cfg = json.loads((ROOT / "microbiomes.json").read_text(encoding="utf-8"))
    ids = [m["id"] for m in cfg["microbiomes"]]
    if len(set(ids)) != len(ids):
        raise SystemExit("microbiomes.json: microbiome ids must be unique")
    for m in cfg["microbiomes"]:
        m["light"], m["dark"] = PALETTE[m["color_slot"]]
    return cfg


def views(cfg):
    """Every single microbiome, then every pair, in config order.
    Each view: {'id', 'title', 'members': [microbiome dicts]}."""
    ms = cfg["microbiomes"]
    out = [{"id": m["id"], "title": m["name"], "members": [m]} for m in ms]
    out += [{"id": f"{a['id']}--{b['id']}", "title": f"{a['name']} vs {b['name']}", "members": [a, b]}
            for a, b in itertools.combinations(ms, 2)]
    return out
