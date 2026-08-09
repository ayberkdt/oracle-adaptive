"""Eskimis-terim tarayicisi.

Her satir: (desen, neden eskimis, superseding karar, izin-verilen-baglam).
Bir eslesme, izin verilen baglam icinde gecerse atlanir (tarihsel kayit,
karsilastirma, ablasyon adi vb.).

Kullanim:  python check_stale.py
Cikis kodu 1 ise supheli eslesme var.
"""
import glob
import pathlib
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BS = chr(92)
HISTORICAL = {"DECISIONS.md"}   # karar gunlugu eski ifadeyi KORUR

FILES = (glob.glob("*.md")
         + glob.glob("paper/chapters/*.tex")
         + ["paper/main.tex", "paper/preamble.tex"])

STALE = [
    (r"piR/N|" + re.escape(BS) + r"pi R/N|πR/N",
     "yuzey yaricapi; uydu yaricapi olmali", "D63/D67",
     r"D3 |tarihsel"),
    (r"\bw_j\b",
     "kuadratur agirligi; omega_j olmali", "D76",
     r"omega"),
    (r"2kN|2k/N",
     "eski prob maliyeti", "D50",
     r"D50|D2 |es-konumlu|co-located|does not apply|yalniz|yalnız|C-lite|"
     r"Clite|increment|already|zaten|sanil|sanıl|wrong one|rather than"),
    (r"Kaula tail|Kaula kuyru|calibrated tail",
     "genlik artik olculen spektrumdan", "D53",
     r"D2 |abl-kaula|power law|guc yasas|güç yasas|tail criterion|tail rule"),
    (r"w\(t\)\s*=|‖Φ\(T,t\)B‖",
     "odunc alinmis skaler agirlik; K_i / diagonal kernel olmali", "D75/D110",
     r"borrowing|odunc|ödünç|alternative|keyfi|keyfî"),
    # D110: A-sens olcutu artik yerel cekirdek. Ham ayrik kosegen IZGARAYA BAGLI
    # (u_i ~ Dt_i => Dt_i^2). PLAN/WP/NOTATION bunu bir tur boyunca kacirdi.
    (r"A-sens[^\n]{0,60}(diag Q|köşegeni|kosegeni)|"
     r"u_i\(N\)ᵀ Q_ii u_i\(N\)|u_q\(N\)ᵀ Q_qq u_q\(N\)|"
     r"mathbf u_i\^\{" + re.escape(BS) + r"top\}" + re.escape(BS) +
     r"mathbf Q_\{ii\}",
     "A-sens olcutu yerel cekirdek K_i; ham ayrik diag Q izgaraya bagli",
     "D110",
     r"literally|harfiyen|kullanilmaz|kullanılmaz|would make|degil|değil|"
     r"yerel|local|K_i|" + re.escape(BS) + r"mathbf K|Delta t_i.{0,20}Delta"),
    # D106: kappa tani, kapi degil. "gate/kapi" kelimesi kappa'ya baglanmamali.
    # Beyaz listeye D106/D74 KOYMA: celiskinin yasadigi yer tam da D106'yi
    # anlatan paragrafti; karar numarasini muaf tutmak kapiyi kor eder.
    # Yalnizca ACIK olumsuzlama muaftir. Satir sonu ayirici degil (metin sarili).
    (r"(kapı|kapi|gate)[^.]{0,40}κ|"
     r"κ[^.]{0,40}(kapısı|kapisi|kapıdır|kapidir|üzerindedir|gate is|gated)",
     "kappa kapi degil, tani (G4 kaldirildi)", "D106",
     r"kapı değil|kapi degil|not a gate|kapı olmaktan|kapı olması değil|"
     r"değil, tanı|degil, tan|yeni bir kapı|tek başına|tek basina|"
     r"uyarıdır|durdurma değil"),
    (r"zaman-indeksli|time-indexed",
     "denetleyici plani faz-indeksli", "D52",
     r"abl-timeindex|D12 |D52|A-sens|A-sign|kiyas|kıyas|vs|yerine|replaced|"
     r"benchmark|degil|değil|yanlis|yanlış|olurdu|would"),
    (r"τ_dec|tau_\{" + re.escape(BS) + r"mathrm\{dec\}\}",
     "tau_corr olarak yeniden adlandirildi", "D63",
     r"Delta t"),
    (r"float32",
     "float64 varsayilan", "D61",
     r"D61|D44|parity|paritesi|opsiyon|gerek yok|gerek kalmad|gerek birak|getirmez"),
    (r"~240|240 varyasyonel|240 propagasyon",
     "M2 maliyeti ~270", "D79",
     r"D46|D79"),
    (r"~430|430 yay",
     "M3 maliyeti ~670", "D92",
     r"D45|D47|D92|erken|Erken"),
    (r"10 ablasyon",
     "ablasyon sayisi 14", "D79",
     r"D46|D79"),
    (r"B2,C \+ B\+ = B_tot = B2,F",
     "eslestirme: capa + aday kalibrasyonu", "D91",
     r"D23|D91"),
    # D142: butce bir TAVAN. Esitlik yazmak hem A-sign'i kotulestirir hem FW
    # sertifikasini gecersiz kilar; J derecede monoton degil.
    (r"B2 ?\+ ?B\+ ?= ?B_tot|B_2 ?\+ ?B_\+ ?" + re.escape(BS) + r";=" +
     re.escape(BS) + r"; ?B_\{" + re.escape(BS) + r"mathrm\{tot\}\}",
     "butce bir tavan; = degil <= olmali", "D142",
     r"D107|D142|tavan|ceiling|boyutsal|dimension"),
    (r"M7'de değil|M7.de degil",
     "asama numaralari kaydi", "D42",
     r"D8 "),
    # D75/D110: olu skaler agirlik `w`. Onceki desen yalniz "w(t) =" ariyordu;
    # ciplak $w$ ve `w d(t,N)` uc yerde hayatta kalmisti.
    (r"\$w\$|\bw\\,d\(|,w,\\lambda|\\mathbf c,w\b",
     "olu skaler agirlik w; K_i olmali", "D75/D110",
     r"borrowing|odunc|ödünç|keyfi|keyfî|would be|olurdu"),
    # D106 (ingilizce): "declared gate", "the gate is kappa" vb.
    (r"κ[^.]{0,30}(declared gate|the gate)|"
     r"(declared|hard|stopping) gate[^.]{0,30}(κ|\\kappa)|"
     r"\\kappa[^.]{0,30}(declared gate|as the gate)|Gate G4|G4 lives|gate G4",
     "kappa kapi degil, tani; G4 kaldirildi", "D106",
     r"not a gate|no stopping rule|diagnostic|warning and not|leading indicator"),
    # D110: ciplak Q_ii de olcut/agirlik olarak gecmemeli.
    (r"\\mathbf Q_\{ii\}",
     "tasinan agirlik K_i (yerel cekirdek), ham Q_ii degil", "D110",
     r"literally|raw diagonal|would make|extra " + re.escape(BS) + r"?Delta"),
    (r"cannot be tabulated|copy of the field|cannot work",
     "sert dil; cozunurluk-maliyet ifadesi olmali", "D30/D62",
     r"not an impossibility|degil|değil"),
]


def main() -> int:
    """Scan every non-historical document for superseded wording."""
    issues = 0
    for f in FILES:
        if f in HISTORICAL:
            continue
        t = pathlib.Path(f).read_text(encoding="utf-8")
        for pat, why, dec, ok in STALE:
            for m in re.finditer(pat, t):
                ctx = re.sub(r"\s+", " ", t[max(0, m.start() - 200):m.start() + 200])
                if re.search(ok, ctx, re.I):
                    continue
                issues += 1
                line = t[:m.start()].count("\n") + 1
                print(f"[{dec}] {f}:{line}  ({why})")
                print(f"    ...{ctx}...")
    print()
    print("SUPHELI ESLESME:", issues)
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
