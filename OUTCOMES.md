# OUTCOMES — hangi sonuç ne kadar iyi, ve neyi yanlış okuruz

## 0. Birincil ölçüler

Semboller [NOTATION.md](NOTATION.md)'de. Kısaca:

- **ρ(komparatör, aday) = E_komparatör / E_aday**, eşleştirme
  ortak tavan `B2 + B+ ≤ B_tot` altında (D142); `B_tot` = `F-op(β)`'nın gerçekleşen işi ve
  **adaylar ona kalibre edilir**, komparatör kaymaz (D91). Yörünge başına oranların
  **medyanı**; `ρ > 1` aday lehine.
- **`ρ̂` ve `ρ` ayrıdır:** `ρ̂` öngörülen hatadan (`Ê = √J`, M1–M2, propagasyon
  yok), `ρ` propagate edilmiş hatadan (M3+). Her tabloda hangisi olduğu yazılır.
- **f = (E_F-op − E_C-plan) / (E_F-op − E_A-sign)**, yalnızca
  `E_A-sign < E_F-op` ve fark çözünürlük zarfını geçtiğinde tanımlı; aksi
  hâlde `N/A` (NOTATION §3) — uçurulabilir
  denetleyicinin, tavanın açtığı boşluğun ne kadarını aldığı.
- **Çözünen sayım** — `M_res > 1`; kararsızlar hiçbir tarafa sayılmaz.
- **Tekrar** — kaç bağımsız popülasyonda aynı yön.

### Komparatör ayrımı — bant sınırlarına gömülü

| Komparatör | Uçurulabilir | İddiadaki rolü |
|---|---|---|
| `F-op` | ✔ | **Yenmek şart.** Yenemezse yöntem yok. |
| `F-env` | ✘ (post-hoc alt zarf) | **Bonus.** Yenmek bandı yükseltir; yenememek başarısızlık değil. |
| `R-int` | ✔ | Mevcut en iyi uçurulabilir üye. Yenmek şart (aksi hâlde önceki makalenin sonucundan ileri gidilmiyor). |

Bu ayrım kritiktir: `F-env` her yörüngenin hatasını gördükten sonra en iyi
sabit dereceyi seçiyor. Önceki makalede doygun dereceden medyan **2.75 kat**
daha iyiydi. Uçurulabilir bir yöntemin onu yenmesi beklenmez.

---

## 1. Sonuç bandı tablosu

### 🟢 ÇOK İYİ — "yöntem çalışıyor"

- `C-plan`, eşleştirme denklemi altında **`F-op`'u (H4) ve `R-int`'i (H4b)**
  yener; çözünen karşılaştırmaların çoğunluğunda, **≥3 bağımsız popülasyonda**.
- β ∈ [0.5, 1.5] boyunca, **β = 0.5 dahil** (H5 doğrulanmış).
- Medyan ρ ≥ 2, yakalama oranı f ≥ 0.33.
- Zorlanmış varyasyonel enstrüman propagate sıralamayı öngörüyor (mekanizma
  doğrulanmış).
- Geniş-eliptik popülasyonda da kazanıyor (rejime hapsolmamış).
- Ek yük `B+` bütçeye dahilken kazanç ayakta.

**Yükseltici (zorunlu değil):** `F-env`'i de yenerse sonuç bandın üst ucudur —
"sabit aileye karşı post-hoc en iyi seçim bile geride kalıyor". Bu, JGCD'de
manşet cümlesi olur.

**Ne demek:** "Yörünge-farkında derece tahsisi, sabit dereceyi eşit
gerçekleşen işte yeniyor." Önceki makalenin negatif sonucunu yapıcıya çeviren
tam makale.

### 🟩 İYİ — "sınırlı ama gerçek bir kazanç"

Yukarıdakilerin çoğu var, biri–ikisi eksik. Tipik biçimler:

- β ≈ 1'de kazanıyor, **β = 0.5'te kaybediyor** (H5 reddedilmiş, H4 sağlam).
- `F-op`'u yeniyor, `R-int` ile berabere.
- 2 popülasyonda kazanıyor, 3.'de kararsız.
- **G7 geçti ama M5 kaynak nedeniyle koşulmadı** → bu bandın alt ucu: kazanç
  β ∈ [0.75, 1.25] boyunca gerçek, ama tek popülasyonda gösterilmiş (D64).
- ρ 1.2–2, f 0.1–0.33.

**Ne demek:** yayınlanabilir tam makale, iddia rejim-sınırlı yazılır:
"kritik derecenin bütçesi civarında, düşük-orta radyal açıklıklı alçak ay
yörüngelerinde". Önceki makalenin dilini devral — nerede *geçerli* olduğunu
söyle.

### 🟡 ORTA — "ödül var, uçurulabilir kısmı yok"

- `A-sign` büyük bir kazanç gösteriyor (masada ödül var) ama hiçbir
  uçurulabilir denetleyici çözünür bir kısmını yakalayamıyor (f < 0.1).
- Ya da `κ` düşük çıkıyor ve **G5b** bunu doğruluyor: sadeleşme taşınamıyor,
  L3 yalnızca L1 seviyesine iniyor. (`κ` tek başına kapı değil — D106.)
- Ya da G1 düşüyor: `ρ̂(R-int, A-sign) < 1.5` → ödül zaten küçük.
- Ya da **G5b** düşüyor: hiçbir adayın öngörülen yakalama oranı `f̂ ≥ 0.15`
  değil → M3'e propagasyon harcanmaz.
- Ya da **G7** düşüyor: `C-plan` β = 0.75 ve β = 1.25'ten yalnız birinde (veya
  hiçbirinde) çözünen çoğunluğu alıyor → kazanç β=1'e özgü, M5 koşulmaz, sonuç
  bu bantta kalır (D64).
- Ya da uçurulabilirlik ikamelerinden biri (Kaula genliği, pilot yay, pilot
  STM) kazancı yiyor.

**Ne demek:** yayınlanabilir, biçim değişir → *yörünge-seviyesi tahsis kıyası
ve uçurulabilir kuralların ona neden ulaşamadığı*. Kısa makale / teknik not.
Önceki makaleyi **güçlendirir**: "kuralı düzeltmek yetmiyor, bilgi kısıtı
esas". `πr/N` dekorelasyon argümanı bu bandın teorik omurgası olur.

M1 kapıları bu bandı **iki gecelik maliyetle** tespit etmek için var.

### 🟠 KÖTÜ — "ödül yok"

- `A-sign`, tam duyarlılık ve işaret bilgisiyle bile, eşit `B2`'de `F-op`'u
  çözünür biçimde yenemiyor.
- Ya da `A-sign`'ın kazancı **propagasyonda hayatta kalmıyor** (H1b reddi):
  doğrusallaştırılmış modelde sömürülen sadeleşme gerçek yörüngede yok oluyor.
  WP17 (`A-sign+`) bunu ayırır.

**Ne demek:** **güçlü bir negatif iddia** ve soruyu kapatır: "tam duyarlılık
bilgisine sahip bir tahsis bile bu bütçelerde sabit dereceyi yenemiyor."
Teknik not olarak yayınlanabilir; yöntem çıkmaz. Önceki makalenin sonuç
bölümünü doğrudan tamamlar.

### ⚫ SONUÇ DEĞİL — enstrüman arızası

Bunlar bant değil, **koşunun geçersiz olduğunun işareti**. "Kötü sonuç" diye
raporlanamaz; tasarım değiştirilip tekrarlanır.

- Her şey sayısal çözünürlük tabanının altına düşüyor, hiçbir karşılaştırma
  karara bağlanmıyor → **G2**; bütçe/yay uzunluğu seçimi gözden geçirilir.
- Çizelge `Δt_acc` ile kararsız → **G3**; yakınsak ızgarada tekrarlanır.
- Varyasyonel tahmin ile propagate sonuç işaret olarak uyuşmuyor → **G6**;
  uygulama hatası aranır.

### 🔴 SAHTE-İYİ — iyi görünen artefaktlar

Sonuç değil, hata. Çoğu önceki kampanyada bilfiil yaşandı. Kontrolleri
ön-tescil edilmeden hiçbir kazanç iddiası yazılmaz.

| Artefakt | Belirtisi | Zorunlu kontrol |
|---|---|---|
| **Referans-komşuluğu** | ρ ≈ 10³–10⁵; hata mikrometre ölçeğinde | Referans dereceyi 300 → 600 çıkarıp tekrar (önceki makalede B020: 1.57×10⁵ → 18.9) |
| **Derece tavanı teması** | Politika referans derecesine ulaşıyor, defekt tanım gereği sıfır | Tavan denetimi + sansür; tavana değen yörüngeler medyandan çıkar |
| **Sızıntı** | Çağrı başına eşit, gerçekte %29–%120 fazla | Birincil bütçe `B2`; dört basamak ayrı |
| **Ek yükün gizlenmesi** | Pilot yay / prob / IFBDA maliyeti dipnotta | `B+` bütçeye **yazılır**, her tabloda görünür |
| **Bayat yön** ★★ | Prob sınırdan önce alınmış, karar sonraki aralık için — yön çoktan dekorele | İleri prob (D26); `abl-lookahead` gerekliliği ölçer |
| **Zayıf komparatör** | Yalnızca doygun sabit dereceye karşı kazanmak | `R-int` de komparatör; `F-env` ayrıca raporlanır |
| **Öngörü sızması** | "Uçurulabilir" denetleyici referans yayı/alanı okuyor | Her `C-*` kaydı kullandığı bilgi kaynağını yazar; kod yolu denetimi |
| **Sadeleşmeye aşırı uyum** | Tek referans yayda muhteşem, faz kaydırınca yok | Faz MC; farklı başlangıç epoku |
| **Ayrıklaştırma artefaktı** ★ | Kazanç `Δt_acc` ile kayboluyor; 120 s'de işaretli integral aliasing | **WP4** — M1'de, geç aşamada değil |
| **Tutarlılık yanılsaması** ★ | Pilot yay planı gerçek yayda geçersiz ama tabloda öyle görünmüyor | **WP6** — `T_coh` ölçülür; `v̂` uçuş anında probe edilir |
| **Prob maliyetinin eksik sayılması** ★★★★ | Ek yük `2k/N` ≈ %5 sanılıyor; gerçekte ileri prob yeni bir konumda tam sentez gerektiriyor, gerçek yük yay integrali `[∫dt/τ_corr]/N_RHS` ve **dereceyle doğrusal** (D141): %8.5 / %21 / %42 @ N=120/300/600 | **D50** — maliyet ölçülür, modellenmez; `κ_eff(n_probe)` eğrisi ve `C-lite` uç noktası WP5'te |
| **Bayat plan hizası** ★★★★ | Pilot yayın ~100 s zaman kayması perilun geçişi mertebesinde; zaman-indeksli plan `c` ve `K_i`'yi tam da baskın oldukları yerde yanlış epoklara koyuyor | **D52** — plan faz-indeksli; WP6 kararlılığı faz-indeksli ölçüyor |
| **Eleyici döngüselliği** ★★★★ | `A-sign` `J`'yi minimize ediyor, eleyici `J`'yi hesaplıyor; önceki 100/100 kalibrasyon optimize edilmemiş politikalar içindi | **D54** — G6 taşıyıcı kapı hâline geldi; başarısızlığı D41 kaçışını tetikler |
| **Güç yasası kuyruğu** ★★★★ | `γ` uydurulmuş Kaula yasasından; ama ay spektrumu tek güç yasası değil (`p_spec` 2.13 vs `p_fit` 1.76) | **D53** — `γ` gerçek `P_n` spektrumundan (1D tablo, uçurulabilir) |
| **Optimizatörün tavana yaslanması** ★★★★ | `A-sign`, bütçe izin verdiğinde dereceyi referansa itip defekti tanım gereği sıfırlıyor; sabit politikada bu teşvik yok | **D57** — çizelge başına tavan teması kesri, yüksek temaslılar orandan çıkar |
| **Tek bütçe noktası** ★★★ | β=1'de kazanıyor ama komşu bütçelerde kaybediyor; önceki makalede iç üyenin bütün sonucu β=0.5'te ters dönmüştü | **G7 / WP20a** — popülasyon genişletmeden önce bütçe probu |
| **Hata düzeyi karışması** ★★★ | `Ê` (öngörülen) ile `E` (propagate) aynı sembolle yazılıp M1 sonucu M3 sonucu sanılıyor | NOTATION §3'te ayrı semboller; her hipotez satırında düzey yazılı |
| **Anahtarlama serbestliği** ★★ | `A-sign` kazanıyor ama sebebi daha çok bilgi değil, dereceyi daha sık değiştirebilmesi | Karar uzayı `K_dec` boyutlu (D16); `C-plan` ile aynı karar ızgarasında; `Δt_dec` taraması |
| **Yerel minimum şansı** ★ | Tek başlangıçtan çıkan çizelge tesadüfen iyi | Çok başlangıçlı koordinat inişi; yayılım raporlanır |
| **Sertifikasız "oracle"** ★ | Üst sınır çözümüne tavan denmesi | **WP7**; medyan `g_E < 0.10` sağlanmazsa "oracle" kelimesi kullanılmaz (D29) |
| **Yuvarlama şansı** | Hata derecede monoton değil; komşu dereceler mertebe farkı | Derece merdiveni taraması |
| **Bin/ızgara bağımlılığı** | Sonuç çizelgenin biçimine (faz-indeksli vs irtifa-binli) veya karar ızgarasına ait | WP12 |
| **Çoklu karşılaştırma** | 6 aday × 5 bütçe × 5 popülasyon → biri hep kazanır | Eleme kuralı önceden; M3'e **tam iki** aday |

★ = Tur 2 · ★★ = Tur 3 · ★★★ = Tur 5 (iç denetim) · ★★★★ = Tur 6 (astrodinamik denetimi).

---

## 2. Risk sicili

| # | Risk | Etki | Azaltım |
|---|---|---|---|
| R1 | Koordinat inişi kötü yerel minimuma takılır | `A-sign` tavan olmaktan çıkar | Çok başlangıçlı; yayılım raporlanır; WP7 alt sınırı bağımsız kontrol |
| R2 | Sertifika boşluğu büyük kalır | "oracle" dili kullanılamaz | Dil kuralı peşinen bağlı (PREREG OA-02); *erişilebilir tahsis* dili |
| R3 | Doğrusallaştırma düşük perilunde (31 km) geçersiz | Genlikler yanlış | Genlik yorumlanmaz, yalnız işaret; gradyan derecesi kontrolü |
| R4 | Kaula genliği `d̂` gerçek `d`'yi kötü temsil eder | Uçurulabilir varyant çöker | `abl-kaula` ayrı ölçülür; `p_fit` yeniden kalibre edilebilir |
| R5 | **Prob yön isabeti düşük (`κ < 0.7`)** | Sadeleşme taşınamaz → 🟡 | WP5; yedek seçenekler; null model kontrolü |
| R6 | **Pilot yay tutarlılık ufku kısa** | "Tablo yükle ve uygula" mimarisi geçersiz | WP6; `v̂` uçuş anında probe; gerekirse WP11 alıcı-ufuk |
| R7 | `B+` kazancı yer | Net kazanç sıfır | Ek yük M5'te ölçülür ve bütçeye dahil edilir |
| R8 | Değişken derece entegratörü fazladan adıma zorlar | Gerçekleşen iş şişer | `C-fb` bütçeyi durum olarak taşır; adım sayısı raporlanır |
| R9 | **`Δt_acc` yakınsamaması, tabloları büyütür** | ~87 MB/yörünge, 128 yörünge | M1'de sıkıştırma/akış çözülür, M4'te değil |
| R10 | 60 günlük yayda çizelge geçersizleşir | Yapıcı sonuç uzun yayda çöker | WP11 + WP18 |
| R11 | Aday sayısı çokluğu (forking paths) | Sahte kazanan | Ön-tescil + sabit eleme + post-hoc etiketi |
| R12 | Eşzamanlı koşular birbirini aç bırakır | Duvar saati patlar | İşçi sayısı kısılır, BelowNormal öncelik |
| R13 | Figür üretimi eski veriyle üzerine yazar | Bayat şekil | Yalnızca `make_figures_oa*.py` |
| R15 | **Prob ek yükü kazancı yer** — dereceyle doğrusal, `N=600`'de **%42**'ye çıkıyor (D141), artı pilot yay %8 | Net kazanç sıfır veya negatif | `κ_eff(n_probe)` eğrisi; `C-lite` varyantı; ek yük `B_tot`'tan düşülür (NOTATION §2) |
| R16 | **Prob perilunda en pahalı, tam da en gerekli olduğu yerde** (`τ_corr` orada en kısa) | Yöntem kendi tatlı noktasında en pahalı | Ölçülür; gerekirse prob yalnızca perilun penceresinde açılır, apolunda kapalı — ayrı varyant olarak sınanır |
| R17 | Eleyici döngüselliği (D54) | M2'nin elemesi yanıltıcı olabilir | G6 taşıyıcı; kaçış D41 |
| R14 | `A-sign` bir *referans-yay* optimumu çıkar | Tavan iddiası zayıflar | WP17 sabit nokta; sonuç öyleyse öyle adlandırılır |

---

## 3. Erken uyarı göstergeleri

M1 biter bitmez bakılacak beş sayı:

| # | Gösterge | Kötü değer | Sonuç |
|---|---|---|---|
| 1 | `ρ̂(R-int, A-sign)` medyanı | < 1.5 | 🟡 |
| 2 | `G_sign / G_sens` (log-iyileşme adımları) | ≈ 0 | Sadeleşme önemsiz; makale zayıflar ama biter |
| 3 | `A-sign` kazancı / sayısal zarf | < 3 | ⚫ enstrüman arızası |
| 4 | `κ` (prob yön isabeti) | < 0.7 | 🟡 |
| 5 | `T_coh` / yay uzunluğu | ≪ 1 | WP11 zorunlu; mimari değişir |
