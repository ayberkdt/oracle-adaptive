"""Sayisal capraz denetim across the documents and the draft.

    Esikler, sayimlar ve maliyetler dosyalar arasinda ayni mi?

check_stale.py eskimis TERIM arar; bu betik eskimis SAYI arar. Ikisi ayri
kategori: bir esik dosyalar arasinda ayrisirsa terim taramasi bunu goremez.

Desenler bilerek DAR: gurultu ureten bir kapi yok sayilir, o yuzden asama
yaylari yalnizca asama/maliyet TABLO satirlarindan okunur, serbest metinden
degil. DECISIONS.md haric (tarihsel kayit eski sayilari korur).

Kullanim: python check_numbers.py   (cikis kodu 1 ise tutarsizlik var)
"""
import collections
import glob
import pathlib
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FILES = ([f for f in glob.glob("*.md") if f != "DECISIONS.md"]
         + glob.glob("paper/chapters/*.tex") + ["paper/main.tex"])

KAPPA = "κ"      # κ
GE = "≥"         # ≥
FHAT = "f̂"      # f̂
NDASH = "–"      # –
STAGE = r"\|\s*\*{0,2}M%d\*{0,2}\s*\|[^|]*\|\s*~?(\d+)\s*\|"

# Bazi esikler belirli dosyalarda BULUNMAK ZORUNDA: bir yerde dogru olmasi,
# baska yerde tamamen eksik olmasini mesru kilmaz.
REQUIRED = {
    "G5b esigi f-hat":     {"PREREG.md", "ROADMAP.md"},
    "sertifika esigi g_E": {"PREREG.md", "NOTATION.md"},
    "ablasyon sayisi":     {"PLAN.md"},
    "M3 yay (tablo)":      {"ROADMAP.md"},
    "M4 yay (tablo)":      {"ROADMAP.md"},
    "M5 yay (tablo)":      {"ROADMAP.md"},
    "M6 yay (tablo)":      {"ROADMAP.md"},
    "prob ek yuku (N=120/300/600)": {"NOTATION.md", "OUTCOMES.md"},
    # WP21 mimari dallanma esikleri (D116): tescil, yol haritasi ve
    # manuskript ayni sayiyi tasimak zorunda.
    "T6a esigi (geometri)": {"PREREG.md", "ROADMAP.md", "WP.md",
                             "paper/chapters/09_ablations.tex"},
    "T4 diz esigi":         {"PREREG.md", "paper/chapters/09_ablations.tex"},
}

PROBES = {
    # G4 kaldirildi (D106): bu sayi artik kapi degil, uyari esigi / sekil
    # referans cizgisi. Deger yine tek olmali, yoksa iki farkli cizgi cizilir.
    "kappa tani referansi": (KAPPA + r"[^\n]{0,25}?" + GE + r"\s*(0\.\d+)", "0.7"),
    "G5b esigi f-hat":     (FHAT + r"\s*" + GE + r"\s*(0\.\d+)", "0.15"),
    "sertifika esigi g_E": (r"g_E\s*<\s*(0\.\d+)", "0.10"),
    "ablasyon sayisi":     (r"(\d+)\s*ablasyon\b", "14"),
    "aday paneli":         (r"(\d+)\s*yörüngelik perilun-tabakalı panel", "26"),
    "ablasyon paneli":     (r"(\d+)\s*yörüngelik küçük panel", "8"),
    "M2 maliyeti":         (r"~(\d+)\s*propagasyon eşdeğeri", "270"),
    # D119: M1 artik entegrasyonsuz degil; sayi M2'ninkiyle karismasin diye
    # ayri bir ibareyle yaziliyor.
    "M1 STM maliyeti":     (r"~?(\d+)\s*(?:varyasyonel eşdeğeri|"
                            r"düz-propagasyon karşılığı)", "180"),
    "M3 yay (tablo)":      (STAGE % 3, "670"),
    "M4 yay (tablo)":      (STAGE % 4, "576"),
    "M5 yay (tablo)":      (STAGE % 5, "3060"),
    "M6 yay (tablo)":      (STAGE % 6, "1730"),
    # D141: tek bir yuzde degil, dereceye bagli uclu. Uc sayi da uc
    # dosyada ayni olmali; "%12-19" hala gecerli ama yalniz
    # dusuk-orta bant icin ve o ibareyle birlikte yaziliyor.
    "prob ek yuku (N=120/300/600)":
        (r"%(8\.5) ?/ ?%?(21) ?/ ?%?(42)", ("8.5", "21", "42")),
    # WP21 / D116
    # Anchor'lar dar: ROADMAP'teki alt esik (0.25) ayni satirda, onu YAKALAMAMALI.
    "T6a esigi (geometri)": (r"(?:≥ ?|more than )\$?(0\.60)", "0.60"),
    "T4 diz esigi":         (r"(?:reaching|ulaşan|≥) ?\$?(0\.90)", "0.90"),
}


def check_register_ranges() -> int:
    """Compare the README's stated D and Q ranges with the register.

    The document index quotes ranges like "D1-D121" and "Q1-Q15" in prose, and
    prose does not update itself.  This drifted silently through fourteen
    rounds (it still read "D1-D66, Q1-Q10" long after both had moved) because
    a terminology scan looks for stale *words* and a numeric scan looks for
    stale *thresholds*; a range that is merely out of date is neither.

    Returns
    -------
    int
        Number of inconsistencies found.
    """
    bad = 0
    decisions = pathlib.Path("DECISIONS.md").read_text(encoding="utf-8")
    readme = pathlib.Path("README.md").read_text(encoding="utf-8")

    for kind, pattern in (("D", r"\*\*(?:~~)?D(\d+)"), ("Q", r"\*\*(?:~~)?Q(\d+)")):
        found = [int(n) for n in re.findall(pattern, decisions)]
        if not found:
            bad += 1
            print(f"!! {kind} sicili: DECISIONS.md'de hic {kind} kaydi yok")
            continue
        top = max(found)
        claimed = re.findall(kind + r"1[^0-9][^0-9]?" + kind + r"(\d+)", readme)
        if not claimed:
            bad += 1
            print(f"!! README {kind} araligi: hic beyan yok "
                  f"(sicilde en buyuk {kind}{top})")
        elif any(int(c) != top for c in claimed):
            bad += 1
            print(f"!! README {kind} araligi: beyan {kind}1-{kind}"
                  f"{'/'.join(claimed)}, sicilde en buyuk {kind}{top}")
        else:
            print(f"OK {kind} sicili: {kind}1-{kind}{top}  (README ile uyumlu)")
    return bad


def check_derived_claims() -> int:
    """Compare the README's claims about the draft with the draft.

    Page count, placeholder count and note count are recomputed from the
    sources and from ``main.log`` rather than compared against a second
    hardcoded number.  A probe in :data:`PROBES` pins two literals to each
    other and goes stale when both move; this one cannot, because one side is
    always the artefact.

    Returns
    -------
    int
        Number of inconsistencies found.
    """
    log = pathlib.Path("paper/main.log")
    if not log.exists():
        print("-- taslak iddialari: paper/main.log yok, atlandi "
              "(once latexmk calistirin)")
        return 0

    sources = [*sorted(glob.glob("paper/chapters/*.tex")), "paper/main.tex"]
    body = "".join(pathlib.Path(f).read_text(encoding="utf-8") for f in sources)
    readme = pathlib.Path("README.md").read_text(encoding="utf-8")
    bs = re.escape(chr(92))

    pages = re.findall(r"\((\d+) pages", log.read_text(encoding="utf-8",
                                                       errors="replace"))
    actual = {
        "sayfa": int(pages[-1]) if pages else None,
        "ph": len(re.findall(bs + r"ph\{", body)),
        "dnote": len(re.findall(bs + r"dnote\{", body)),
    }
    claimed = {
        "sayfa": _first_int(readme, r"\*\*(\d+) sayfa\*\*"),
        "ph": _first_int(readme, r"`?" + bs + r"?ph`? sayısı (\d+)"),
        "dnote": _first_int(readme, r"`?" + bs + r"?dnote`? sayısı (\d+)"),
    }

    bad = 0
    for key, got in actual.items():
        want = claimed[key]
        if got is None or want is None:
            bad += 1
            print(f"!! taslak {key}: okunamadi (beyan={want}, olcum={got})")
        elif got != want:
            bad += 1
            print(f"!! taslak {key}: README {want} diyor, olcum {got}")
        else:
            print(f"OK taslak {key}: {got}  (kaynaktan doğrulandı)")
    return bad


def _first_int(text: str, pattern: str) -> int | None:
    """First integer captured by ``pattern``, or ``None`` if it does not match."""
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


def main() -> int:
    """Run every check and report the number of inconsistencies."""
    bad = check_register_ranges() + check_derived_claims()
    for name, (pat, want) in PROBES.items():
        hits = collections.defaultdict(list)
        for f in FILES:
            for m in re.finditer(pat, pathlib.Path(f).read_text(encoding="utf-8")):
                v = tuple(g for g in m.groups() if g)
                if v:
                    hits[v].append(f)
        exp = want if isinstance(want, tuple) else (want,)
        if not hits:
            bad += 1
            print(f"!! {name}: HIC ESLESME YOK "
                  f"(silinmis mi, yoksa desen mi bayat?)")
            continue
        seen = {f.replace("\\", "/") for fs in hits.values() for f in fs}
        missing = REQUIRED.get(name, set()) - seen
        if missing:
            bad += 1
            print(f"!! {name}: zorunlu dosyalarda yok: {sorted(missing)}")
        wrong = {v for v in hits if v != exp}
        if wrong:
            bad += 1
            print(f"!! {name}: beklenen {exp}")
            for v, fs in hits.items():
                mark = "<--" if v != exp else "   "
                print(f"   {mark} {v}  {sorted(set(fs))}")
        else:
            n = len({f for fs in hits.values() for f in fs})
            print(f"OK {name}: {exp[0] if len(exp) == 1 else exp}  ({n} dosya)")
    print()
    print("TUTARSIZLIK:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
