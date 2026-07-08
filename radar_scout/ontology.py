"""Free-text -> controlled-vocabulary harmonization via the EBI Ontology Lookup Service.

Maps a user's disease name to MONDO/DOID and a cell-type name to Cell Ontology (CL),
so that "breast cancer", "breast carcinoma", and "mammary cancer" all resolve to one
canonical term — the string the CELLxGENE Census actually filters on. This is the
front door of RADAR-Scout: it turns messy input into an automatable query.

Uses only the Python standard library (urllib) so it stays dependency-free and easy
to test. Network calls are best-effort; callers should handle ``None``.
"""

from __future__ import annotations

import json
import ssl
import urllib.parse
import urllib.request
from dataclasses import dataclass

OLS_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"

# The macOS python.org framework build ships without a CA bundle, so HTTPS
# verification fails against certifi-less defaults. Prefer certifi when present.
try:
    import certifi

    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:  # pragma: no cover
    _SSL_CTX = ssl.create_default_context()


@dataclass(frozen=True)
class OntologyHit:
    query: str
    id: str  # e.g. "MONDO:0006574" or "CL:0000057"
    label: str  # canonical label (what Census stores)
    ontology: str  # "mondo" | "doid" | "cl"
    iri: str

    def __str__(self) -> str:
        return f"{self.query!r} -> {self.label} [{self.id}]"


def _ols_search(query: str, ontology: str, rows: int = 8, timeout: float = 20.0) -> list[dict]:
    params = urllib.parse.urlencode(
        {"q": query, "ontology": ontology, "rows": rows, "type": "class"}
    )
    url = f"{OLS_SEARCH}?{params}"
    # EBI OLS rejects the default Python-urllib user-agent with 403; set a real one.
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": "RADAR-Scout/0.1 (Stanford iGEM)"}
    )
    with urllib.request.urlopen(req, timeout=timeout, context=_SSL_CTX) as resp:
        data = json.load(resp)
    return data.get("response", {}).get("docs", [])


def _best(query: str, docs: list[dict], ontology: str) -> OntologyHit | None:
    """Pick the best doc: exact label match first, then correct-prefix, then first."""
    if not docs:
        return None
    ql = query.strip().lower()
    prefix = ontology.upper() + ":"

    def make(d: dict) -> OntologyHit:
        return OntologyHit(query=query, id=d.get("obo_id", ""), label=d.get("label", ""),
                           ontology=ontology, iri=d.get("iri", ""))

    exact = [d for d in docs if d.get("label", "").lower() == ql]
    if exact:
        return make(exact[0])
    syn = [d for d in docs if any(s.lower() == ql for s in d.get("synonym", []) or [])]
    if syn:
        return make(syn[0])
    pref = [d for d in docs if str(d.get("obo_id", "")).startswith(prefix)]
    if pref:
        return make(pref[0])
    return make(docs[0])


def resolve_disease(name: str) -> OntologyHit | None:
    """Resolve a disease name to MONDO (falling back to DOID)."""
    for onto in ("mondo", "doid"):
        try:
            hit = _best(name, _ols_search(name, onto), onto)
        except Exception:
            hit = None
        if hit and hit.id:
            return hit
    return None


def resolve_cell_type(name: str) -> OntologyHit | None:
    """Resolve a cell-type name to Cell Ontology (CL)."""
    try:
        return _best(name, _ols_search(name, "cl"), "cl")
    except Exception:
        return None


def resolve(disease: str, cell_type: str) -> dict:
    """Normalize a (disease, cell_type) query. Returns canonical labels + IDs.

    ``disease_label`` / ``cell_type_label`` are the strings to feed into the
    Census ``obs_value_filter``.
    """
    d = resolve_disease(disease)
    c = resolve_cell_type(cell_type)
    return {
        "disease_query": disease, "disease_id": d.id if d else None,
        "disease_label": d.label if d else None,
        "cell_type_query": cell_type, "cell_type_id": c.id if c else None,
        "cell_type_label": c.label if c else None,
    }


if __name__ == "__main__":  # quick live demo
    for q in ["keloid", "breast cancer", "mammary cancer", "pulmonary fibrosis"]:
        print(resolve_disease(q))
    for q in ["fibroblast", "T cell", "CD3+ cell", "keratinocyte"]:
        print(resolve_cell_type(q))
