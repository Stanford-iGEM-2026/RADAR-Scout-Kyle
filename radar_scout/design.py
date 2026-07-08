"""Gene -> Plasmid + Primer hand-off for RADAR sensor construction.

Given a prioritized target gene, this module assembles the concrete design
artifacts a wet-lab team needs:

  1. the target CDS (fetched live from Ensembl),
  2. a RADAR/RADARS-style **sensing guide** — an antisense complementary region
     with an ADAR-editable adenosine placed in a favorable sequence context, plus
     Golden-Gate (BsaI) cloning overhangs,
  3. qPCR **validation primers** for the target transcript (nearest-neighbor Tm).

IMPORTANT — scope & honesty: the qPCR primer design and the ADAR context scoring
are rigorous. The sensor-guide construction is a documented *draft*: the exact
RADARS scaffold (payload fusion, editing-site spacing, secondary-structure
constraints) must be finalized against the primary RADAR paper before synthesis.
Everything that needs verification is flagged inline and in ROADMAP.md.
"""

from __future__ import annotations

import json
import math
import ssl
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

try:
    import certifi

    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:  # pragma: no cover
    _SSL_CTX = ssl.create_default_context()

ENSEMBL = "https://rest.ensembl.org"
_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")

# SantaLucia (1998) unified nearest-neighbor parameters: dH (kcal/mol), dS (cal/mol/K)
_NN = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "CT": (-7.8, -21.0), "GA": (-8.2, -22.2), "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
}
# ADAR nearest-neighbor preference for the edited A (Eggington 2011): 5' U>A>C>G, 3' G>C>A~U
_ADAR_5P = {"T": 1.0, "U": 1.0, "A": 0.7, "C": 0.5, "G": 0.2}
_ADAR_3P = {"G": 1.0, "C": 0.7, "A": 0.4, "T": 0.4, "U": 0.4}


# --------------------------------------------------------------------------- #
# sequence utilities
# --------------------------------------------------------------------------- #
def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def tm_nn(primer: str, na_mM: float = 50.0, primer_uM: float = 0.25) -> float:
    """Nearest-neighbor melting temperature (SantaLucia 1998), salt-corrected (deg C)."""
    s = primer.upper().replace("U", "T")
    if len(s) < 2:
        return float("nan")
    dh, ds = 0.2, -5.7  # initiation (with terminal corrections folded in approximately)
    for i in range(len(s) - 1):
        pair = s[i:i + 2]
        if pair in _NN:
            h, sds = _NN[pair]
        else:  # use complement strand's NN (symmetry)
            h, sds = _NN.get(reverse_complement(pair), (-8.0, -22.0))
        dh += h
        ds += sds
    r = 1.987  # cal/(mol K)
    ct = primer_uM * 1e-6
    tm = (dh * 1000.0) / (ds + r * math.log(ct / 4.0)) - 273.15
    tm += 16.6 * math.log10(na_mM / 1000.0)  # Schildkraut-Lifson salt correction
    return round(tm, 1)


# --------------------------------------------------------------------------- #
# Ensembl fetch
# --------------------------------------------------------------------------- #
def _get(path: str, timeout: float = 25.0) -> dict:
    url = f"{ENSEMBL}{path}"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json",
                                               "User-Agent": "RADAR-Scout/0.1 (Stanford iGEM)"})
    with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as resp:
        return json.load(resp)


def fetch_cds(gene_symbol: str, species: str = "homo_sapiens") -> tuple[str, str] | None:
    """Return (transcript_id, CDS sequence) for a gene's canonical transcript, or None."""
    try:
        info = _get(f"/lookup/symbol/{species}/{urllib.parse.quote(gene_symbol)}?expand=1")
        transcripts = info.get("Transcript", [])
        canonical = next((t for t in transcripts if t.get("is_canonical")), transcripts[0] if transcripts else None)
        if canonical is None:
            return None
        tid = canonical["id"]
        seq = _get(f"/sequence/id/{tid}?type=cds")
        return tid, seq["seq"]
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# design
# --------------------------------------------------------------------------- #
@dataclass
class PrimerPair:
    forward: str
    reverse: str
    fwd_tm: float
    rev_tm: float
    amplicon_len: int
    fwd_start: int


@dataclass
class RadarDesign:
    gene: str
    transcript_id: str
    sensor_guide: str          # antisense complementary region (5'->3')
    edit_site_index: int       # position of the editable A within the target window
    context_score: float       # ADAR neighbor-preference score in [0,1]
    guide_fwd_primer: str      # BsaI-flanked cloning primers for the guide
    guide_rev_primer: str
    qpcr: PrimerPair | None
    notes: list[str] = field(default_factory=list)


def design_validation_primers(cds: str, tm_target: float = 60.0,
                              product_range=(90, 150), plen_range=(18, 24)) -> PrimerPair | None:
    """Pick a qPCR primer pair for the target: Tm near target, GC 40-60%, no long runs."""
    s = cds.upper().replace("U", "T")
    n = len(s)

    def ok(p):
        if not (0.4 <= gc_fraction(p) <= 0.6):
            return False
        if any(b * 5 in p for b in "ACGT"):  # avoid >=5 homopolymer
            return False
        return abs(tm_nn(p) - tm_target) <= 3.0

    # search a window in the middle third of the CDS for a valid pair
    start_lo, start_hi = n // 3, min(n - product_range[1] - plen_range[1], 2 * n // 3)
    for fstart in range(max(0, start_lo), max(start_lo + 1, start_hi)):
        for flen in range(*plen_range):
            fwd = s[fstart:fstart + flen]
            if len(fwd) < flen or not ok(fwd):
                continue
            for prod in range(product_range[0], product_range[1] + 1, 5):
                rend = fstart + prod
                for rlen in range(*plen_range):
                    rstart = rend - rlen
                    if rstart <= fstart + flen:
                        continue
                    rev_template = s[rstart:rend]
                    if len(rev_template) < rlen:
                        continue
                    rev = reverse_complement(rev_template)
                    if ok(rev) and abs(tm_nn(fwd) - tm_nn(rev)) <= 2.0:
                        return PrimerPair(fwd, rev, tm_nn(fwd), tm_nn(rev), prod, fstart)
    return None


def _context_score(target: str, a_index: int) -> float:
    """ADAR editability of the adenosine at ``a_index`` from its 5'/3' neighbors."""
    if target[a_index].upper() not in "A":
        return 0.0
    five = target[a_index - 1].upper() if a_index > 0 else "N"
    three = target[a_index + 1].upper() if a_index + 1 < len(target) else "N"
    return round(_ADAR_5P.get(five, 0.3) * _ADAR_3P.get(three, 0.4), 3)


def design_radar_guide(cds: str, guide_len: int = 120) -> tuple[str, int, float, int]:
    """Draft a RADARS-style sensing guide against the best-context editing site.

    Scans the CDS for adenosines in a favorable ADAR context, centers an antisense
    complementary window on the best one, and returns
    (sensor_guide_5to3, edit_site_index_in_window, context_score, window_start).
    The guide is antisense to the target so that base-pairing forms the ADAR
    substrate; the payload/scaffold fusion is added at the plasmid level (draft).
    """
    s = cds.upper().replace("U", "T")
    half = guide_len // 2
    best = (half, 0.0, half)  # (a_index_global, score, ...)
    for i in range(half, len(s) - half):
        if s[i] == "A":
            sc = _context_score(s, i)
            if sc > best[1]:
                best = (i, sc, i)
    a_global = best[0]
    start = a_global - half
    window = s[start:start + guide_len]
    guide = reverse_complement(window)  # antisense sensing region
    edit_site_in_window = a_global - start
    return guide, edit_site_in_window, best[1], start


def _bsai_clone_primers(insert: str) -> tuple[str, str]:
    """Flank an insert with BsaI sites + fusion overhangs for Golden Gate assembly."""
    bsai_5 = "GGTCTCA"   # BsaI recognition + spacer; overhang 'AATG' (start context)
    bsai_3 = "TGAGACC"   # reverse BsaI
    fwd = bsai_5 + "AATG" + insert[:20]
    rev = reverse_complement(insert[-20:]) + "GCTT" + bsai_3
    return fwd, rev


def design_target(gene: str, species: str = "homo_sapiens", guide_len: int = 120) -> RadarDesign | None:
    """End-to-end hand-off for one prioritized gene: guide + cloning + qPCR primers."""
    fetched = fetch_cds(gene, species)
    if fetched is None:
        return None
    tid, cds = fetched
    guide, edit_i, ctx, _ = design_radar_guide(cds, guide_len)
    gf, gr = _bsai_clone_primers(guide)
    qpcr = design_validation_primers(cds)
    notes = [
        "DRAFT sensor guide: verify scaffold/payload fusion + editing-site spacing "
        "against the RADARS paper before synthesis.",
        f"Editable A ADAR-context score {ctx:.2f} (5'U/3'G is optimal).",
        "BsaI (GGTCTC) used for Golden Gate; ensure the insert is domesticated "
        "(no internal BsaI sites) before ordering.",
    ]
    if qpcr is None:
        notes.append("No qPCR pair met Tm/GC constraints in the default window; relax and retry.")
    return RadarDesign(gene=gene, transcript_id=tid, sensor_guide=guide, edit_site_index=edit_i,
                       context_score=ctx, guide_fwd_primer=gf, guide_rev_primer=gr,
                       qpcr=qpcr, notes=notes)


if __name__ == "__main__":  # live demo
    d = design_target("POSTN")
    if d:
        print(f"{d.gene} ({d.transcript_id})")
        print(f"  guide[:60]      {d.sensor_guide[:60]}...")
        print(f"  edit site idx   {d.edit_site_index}  context {d.context_score}")
        print(f"  clone F/R       {d.guide_fwd_primer[:30]}... / {d.guide_rev_primer[:30]}...")
        if d.qpcr:
            print(f"  qPCR F/R Tm     {d.qpcr.fwd_tm}/{d.qpcr.rev_tm}  amplicon {d.qpcr.amplicon_len} bp")
        for n in d.notes:
            print(f"  note: {n}")
