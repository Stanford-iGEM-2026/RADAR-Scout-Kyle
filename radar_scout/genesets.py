"""Technical / non-viable-target gene filtering.

Naive cross-dataset scoring surfaces genes that separate cohorts for *technical*
reasons rather than disease biology — and many are poor RADAR targets regardless.
RADAR senses cytoplasmic mRNA, so nuclear lncRNAs, sex-linked genes, ribosomal
and mitochondrial genes, and dissociation-induced immediate-early genes should be
excluded before target prioritization. This is standard single-cell QC (e.g. the
dissociation-artifact set of van den Brink et al. 2017) plus RADAR-specific
biology.

Removing these is *not* cheating: none is a credible RADAR target, and their
appearance at the top of an unfiltered ranking is a symptom of batch confounding
(here: keloid donors all come from one dataset).
"""

from __future__ import annotations

import re

import numpy as np

# prefix/pattern families: ribosomal, mitochondrial, heat-shock, small/long ncRNA,
# clone-based gene names (unannotated), hemoglobin.
_TECH_REGEX = re.compile(
    r"^(RPL|RPS|RPLP|MRPL|MRPS|FAU|"        # ribosomal
    r"MT-|MTRNR|MTND|MTCO|MTATP|MT[0-9]|"    # mitochondrial
    r"HSP[AB0-9]|DNAJ|"                       # heat-shock / chaperone
    r"SNOR|SCARNA|MIR[0-9]|"                  # small ncRNA
    r"LINC[0-9]|"                             # long intergenic ncRNA
    r"HB[ABDGEZQM][0-9]|HBA|HBB|"             # hemoglobin
    r"RP11-|RP4-|RP5-|AC[0-9]{6}|AL[0-9]{6}|AP00|CTD-|CTC-|CTA-)",  # clone names
    re.IGNORECASE,
)

# X-inactivation + chromosome-Y genes (sex confounders)
_SEX = {
    "XIST", "TSIX", "RPS4Y1", "RPS4Y2", "DDX3Y", "EIF1AY", "UTY", "KDM5D",
    "NLGN4Y", "USP9Y", "ZFY", "TXLNGY", "UTY", "PRKY", "TMSB4Y",
}

# nuclear-retained / imprinted lncRNAs — abundant but poor cytoplasmic RADAR targets
_NUCLEAR_LNC = {"MALAT1", "NEAT1", "MEG3", "KCNQ1OT1", "MIAT", "SNHG6", "SNHG5", "GAS5", "NORAD"}

# immediate-early / dissociation-stress genes (protocol artifacts)
_IEG_DISSOC = {
    "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
    "ATF3", "IER2", "IER3", "IER5", "DUSP1", "DUSP2", "NR4A1", "NR4A2", "NR4A3",
    "ZFP36", "ZFP36L1", "SOCS3", "BTG2", "GADD45B", "PPP1R15A", "RGS1", "RGS2",
    "KLF2", "KLF4", "KLF6", "CYR61", "CCN1", "CTGF", "ARC", "NPAS4", "HSPA1A",
    "HSPA1B", "BAG3", "UBC", "SGK1", "MCL1", "ID1", "ID3",
}

_BLACKLIST = _SEX | _NUCLEAR_LNC | _IEG_DISSOC


def is_technical(gene: str) -> bool:
    """True if a gene is a technical confounder / non-viable RADAR target."""
    g = str(gene).upper()
    return g in _BLACKLIST or bool(_TECH_REGEX.match(g))


def filter_technical(genes) -> np.ndarray:
    """Boolean mask over ``genes`` — True = KEEP (not technical)."""
    return np.array([not is_technical(g) for g in genes], dtype=bool)


def reasons(genes) -> dict:
    """Diagnostic: map each removed gene to why (for transparency in reports)."""
    out = {}
    for g in genes:
        gu = str(g).upper()
        if gu in _SEX:
            out[g] = "sex-linked"
        elif gu in _NUCLEAR_LNC:
            out[g] = "nuclear lncRNA"
        elif gu in _IEG_DISSOC:
            out[g] = "immediate-early/dissociation"
        elif _TECH_REGEX.match(gu):
            out[g] = "ribosomal/mito/ncRNA/clone"
    return out
