"""Genus-level cell volume and growth-equivalent maintenance rate m.

Function form of the original micrm_shared/cell_size_m.py (same rules).

Cell diameter/length come from the Madin et al. (2020) condensed species trait
table (genus-level lookup; this is NOT inferred from genes).
Volume geometry: sphere V=pi d^3/6; rod V=pi d^2 (L-d)/4 + pi d^3/6.
Growth-equivalent maintenance (Lynch & Marinov 2015, 20 C):
    m = 0.01449 * V^-0.09   [1/h]
Genera without trait data receive a seeded random volume drawn uniformly from
the range of resolved genera ("sample"), or the maximum volume ("max").
"""

import math
import random
from collections import defaultdict

from common import read_csv

ALIASES = {
    "Escherichia-Shigella": ["Escherichia", "Shigella"],
    "Burkholderia-Caballeronia-Paraburkholderia": ["Burkholderia"],
    "Mediterraneibacter": ["Ruminococcus"],
    "Schaalia": ["Actinomyces"],
    "Segatella": ["Prevotella"],
    "Lancefieldella": ["Atopobium"],
    "Amylolactobacillus": ["Lactobacillus"],
    "Anaerobutyricum": ["Eubacterium"],
}


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def midpoint(lo, up):
    a, b = number(lo), number(up)
    if a is None and b is None:
        return None
    if a is None:
        return b
    if b is None:
        return a
    return 0.5 * (a + b)


def volume(d1, d2, shape):
    if d1 is None or d1 <= 0:
        return None
    sphere = math.pi * d1**3 / 6.0
    text = (shape or "").lower()
    if "cocc" in text or "sphere" in text or d2 is None or d2 <= d1:
        return sphere
    return math.pi * d1**2 * (d2 - d1) / 4.0 + sphere


def m_of(v):
    return 0.01449 * (v**-0.09)


def estimate(traits_path, genera, policy="sample", seed=42):
    species = defaultdict(list)
    for row in read_csv(traits_path):
        genus = (row.get("genus") or "").strip()
        if not genus:
            continue
        v = volume(midpoint(row.get("d1_lo"), row.get("d1_up")),
                   midpoint(row.get("d2_lo"), row.get("d2_up")), row.get("cell_shape"))
        if v is not None and v > 0:
            species[genus].append(v)
    stats = {g: sorted(vs)[len(vs) // 2] for g, vs in species.items()}
    genera = sorted(set(genera))
    resolved, source = {}, {}
    for genus in genera:
        if genus in stats:
            resolved[genus], source[genus] = stats[genus], "madin_genus"
            continue
        aliases = [a for a in ALIASES.get(genus, []) if a in stats]
        if aliases:
            resolved[genus] = sorted(stats[a] for a in aliases)[len(aliases) // 2]
            source[genus] = "madin_alias:" + ",".join(aliases)
    if not resolved:
        raise ValueError("No modeled genus has usable Madin cell-size data")
    rng = random.Random(seed)
    pool = sorted(resolved.values())
    rows = []
    for genus in genera:
        v = resolved.get(genus)
        if v is not None:
            m, src = m_of(v), source[genus]
        elif policy == "max":
            m, src = m_of(pool[-1]), "measured_volume_max"
        else:
            m, src = m_of(rng.uniform(pool[0], pool[-1])), "measured_volume_sample"
        rows.append(dict(genus=genus, volume_um3="" if v is None else v, m=m, source=src, resolved=int(v is not None)))
    return rows
