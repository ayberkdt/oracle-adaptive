"""Sayisal capraz denetim: esikler, sayimlar ve maliyetler dosyalar arasinda ayni mi?

check_stale.py eskimis TERIM arar; bu betik eskimis SAYI arar. Ikisi ayri
kategori: bir esik dosyalar arasinda ayrisirsa terim taramasi bunu goremez.

Desenler bilerek DAR: gurultu ureten bir kapi yok sayilir, o yuzden asama
yaylari yalnizca asama/maliyet TABLO satirlarindan okunur, serbest metinden
degil. DECISIONS.md haric (tarihsel kayit eski sayilari korur).

Kullanim: python check_numbers.py   (cikis kodu 1 ise tutarsizlik var)
"""
import io
import glob
import re
import sys
import collections

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
    "prob ek yuku":        (r"~?%(12)" + NDASH + r"(19)", ("12", "19")),
    # WP21 / D116
    # Anchor'lar dar: ROADMAP'teki alt esik (0.25) ayni satirda, onu YAKALAMAMALI.
    "T6a esigi (geometri)": (r"(?:≥ ?|more than )\$?(0\.60)", "0.60"),
    "T4 diz esigi":         (r"(?:reaching|ulaşan|≥) ?\$?(0\.90)", "0.90"),
}


def check_register_ranges() -> int:
    """Do the README's stated D and Q ranges match the actual register?

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
    decisions = io.open("DECISIONS.md", encoding="utf-8").read()
    readme = io.open("README.md", encoding="utf-8").read()

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


def main():
    bad = check_register_ranges()
    for name, (pat, want) in PROBES.items():
        hits = collections.defaultdict(list)
        for f in FILES:
            for m in re.finditer(pat, io.open(f, encoding="utf-8").read()):
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
