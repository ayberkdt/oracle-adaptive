# ROADMAP — aşamalar, kapılar, maliyet

İş paketlerinin tanımı [WP.md](WP.md)'de; buradaki tablo onları zamana ve
karar kapılarına diziyor.

**Tasarım ilkesi:** ucuz önce, pahalı sonra, her aşamanın sonunda yazılı bir
durdurma kriteri. Önceki makalenin en pahalı dersi, propagasyonların
enstrümanın ölçemediği bütçelerde harcanmasıydı (β=3'te 57 karşılaştırmanın
56'sı kararsız). İkinci dersi: β=1 tek bir noktadır ve orada kazanmak
yöntemin kazandığı anlamına gelmez.

Bu yüzden **popülasyon genişletmeden önce bütçe probu var** (M4/G7). Tek bir
bütçede kazanan ama komşusunda kaybeden bir yöntem için üç popülasyon
propagate etmek, önceki kampanyanın yaptığı hatanın aynısı olur.

---

## Aşama haritası

| Aşama | WP | Yay | Duvar saati | Çıkış kapısı |
|---|---|---|---|---|
| **M0** | WP0 | 0 | saatler | kabul kontrolü (+ geri düşüş) |
| **M1** | WP1–WP6 | ~10 pilot + **~180 varyasyonel eşdeğeri** (D119) + **3.1 M sentezlik tablo** (D137) | 2–3 gece | **G1, G2, G3** (+`κ` tanısı) |
| **M2** | WP7–WP13 | ~270 varyasyonel çözüm | 3–4 gece | **G5, G5b** |
| **M2b** | WP21 (T1–T7) | **0** | 1–2 gece | mimari dallanma (kapı değil) |
| **M3** | WP14 + WP17a | ~670 | 3–4 gece | **G6** |
| **M4** | WP20a (bütçe probu) | ~576 | 2–3 gece | **G7** |
| **M5** | WP15 | ~3060 | 6–9 gece | — |
| **M6** | WP20b (tam ızgara) | ~1730 | 5–7 gece | — |
| **M7** | WP16 + WP17b | ~150 (seri) | 2 gece + boşta makine | — |
| **M8** | WP18 | ~40 uzun yay | 2–3 gece | — |
| **M9** | WP19 | 0 | — | konsey GEÇER |

Kapılarda durursa maliyet kesildiği yerde kalır: M1'de ~3 gece, M2'de ~7 gece,
M4'te ~13 gece. Planın esas amacı bu.

---

## M0/M1 öncesi açık borçlar

Tam liste, koşu başlamadan önce. Bir madde burada yoksa kapanmış demektir;
kapanınca buradan silinir ve DECISIONS'a geçer.

### Koşuyu bloke edenler

| # | Borç | Neden bloke | Nerede |
|---|---|---|---|
| **B1** | **Kümülatif-dereceye-göre çekirdek giriş noktası** yok | M1'in en büyük kalemini 60 kat düşürüyor (`N=600`: 8.4 sa → 8 dk) **ve** §7'nin prob maliyet modelini açıyor. Doğrulama şartı: bantların toplamı = sabit-derece sonucu, makine hassasiyetinde | D120/D138 |
| **B2** | **Q13 kapalı değil:** kıyasın STM'i hangi gradyan derecesinde? | Düşük-perilun popülasyonu var ve önceki kampanya derece-120 gradyanın 31 km'de yetersiz olduğunu **ölçtü**. Ana kıyas için `N_G = N_ref`, sonra 8–14 yörüngelik tabakalı panelde `120→300→600` yakınsama kontrolü | Q13 |
| **B3** | **`lunaris` commit SHA'sı pinlenmedi**; ayrıca özel (`_` önekli) `_compute_sh_acceleration_dual_numba` API'sine bağımlıyız | Yeniden üretilebilirlik; özel bir API sessizce değişirse kampanya iki kod yolundan gelir | D135 |
| **B4** | **WP0 kabul kontrolü hiç koşulmadı** | Bütün kampanyanın "tek kod yolu" kuralı buna dayanıyor | WP0 |
| **B5** ★ | **Bütçe `≤` mi `=` mi — sertifika geçerliliği.** Problem `Σ W_q N_q² ≤ B` diye tanımlı ama Algoritma 1 `W−B` üzerinde bisekt ediyordu ve FW uygulaması her yinelemede `|W(θ)−B|`'yi makine hassasiyetinde sıfıra zorluyordu. **`J` derecede monoton değil** — makalenin bütün fikri sadeleşme; daha yüksek derece bir çiftin bir terimini küçültüp `J`'yi **artırabilir**. Dolayısıyla `W* < B` gayet mümkün. LP'nin gerçek optimumu boşluk bırakıyorsa ve onu zorla yüzeye itersen **LMO'yu çözmemiş olursun** ve `L_FW ≤ J*` garantisi düşer. **Düzeltildi:** KKT (`λ≥0`, `B−W≥0`, `λ(B−W)=0`), `λ=0` önce sınanıyor; LMO `≤B` üzerinde çözülüyor ve LP'nin kendi optimallik sertifikası doğrulanıyor; eşleştirme "aynı gerçekleşen-iş **tavanı**" oldu ve her politikanın harcadığı iş hatasının yanında raporlanıyor. | **KAPANDI (D142)** |
| **B6** | `build_decision_grid`'in eşleşme kontrolü `K_dec × M` yayın çarpımı kuruyordu: `N=300`'de 2.6×10⁸ karşılaştırma ve ~4 GB float geçici. Birim testte görünmez, gerçek yayda öldürücü. **Düzeltildi:** ikili arama, `O(K log M)` ve `O(K)` RAM. | **KAPANDI (D143)** |
| **B7** | T3 ham `A_i` üzerinde iz/özdeğer/özyön ölçüyordu; `A_i` durum üzerinde **karışık birimli** bir kuadratik form, yani verdikt m/s ↔ km/s değişiminde değişebilir. **Düzeltildi:** ölçüm `Ã_i = Sᵀ A_i S` üzerinde (`tda.stm.nondimensionalise_form`, **eşlenik** dönüşüm — operatör benzerliği değil). | **KAPANDI (D144)** |
| **B8** | T3a/T3b kanonik dallanma kuralı PREREG'de vardı, ROADMAP ve WP'de yoktu; `p=3` hiç tanımlı değildi. **Düzeltildi:** kural üç dosyada aynı — `p ≤ 2` **ve** medyan asal açı `< 15°` ⇒ `C-rank`; aksi hâlde tam durum; `p=3` açıkça tam duruma gidiyor ve `C-rank` gelecek-iş satırına yazılıyor. | **KAPANDI (D145)** |

### Yazılmamış kod

| Modül | İçerik | Bağımlı olduğu |
|---|---|---|
| `archive.py` | salt-okunur arşiv erişimi + WP0 kabul kontrolü | — |
| `tables.py` | `Δa(t,N)` memmap tablosu, şema, provenance | archive, field, grids |
| `kernel.py` | `A_i`, `K_i`, `c_i`, `J` — önek/sonek, Neumaier | stm, grids, tables |
| `allocate/` | `separable`, `descent`, `frankwolfe` (LP-LMO, D132), `rounding` | kernel |
| `policies.py` | `F-op`, `F-env`, `R-rad`, `R-int` | grids |
| `probe.py`, `controller/` | bant probu, IFBDA, `C-plan`/`C-lite`/`C-fb` | spectrum, kernel |
| `analysis/` | WP21'in T1–T7'si | allocate, kernel |
| `io/` | manifest, digest zinciri, ön-tescil hash'i | — |

### Ölçülmemiş sayılar (şu an tahmin)

| Nicelik | Şu anki değer | Nasıl kapanır |
|---|---|---|
| `M` (hücre sayısı) | ~51 100 @ `N=300`, **vekil yörüngede** | gerçek arşiv yaylarında WP1 |
| Tablo boyutu | ~74 MB/yörünge | `M` ölçülünce |
| **Prob ek yükü** | **%12–19 yazılı ama dereceyle doğrusal** (D141) | WP5/WP16, derece bandı başına |
| `inverse_residual`, `κ₂(Φ)` | hiç ölçülmedi | WP1, yay başına |
| `float32` vs `float64` verdict paritesi | koşulmadı | WP1 |
| Çekirdek zamanı bu makinede | arşivin makinesinden alındı | WP16 |

### Manuskript borçları

- 21 float'ın 13'ü metinden referanslanmıyor (WP19 teslim listesinde).
- §2'de iki placeholder: kuadratik-sırt çantası / tamsayı-QP literatürü ve
  çok-fidelite paragrafı.
- **Q14** yazılmadı: ileri prob, en ucuz rakibi olan *bir önceki devirde aynı
  fazda ölçüm*e karşı gerekçelendirilmemiş (yer izi kayması ≈33 km vs
  `πr/N` = 48/19/10 km).
- **Q15** yazılmadı: hızlı alan değerlendirme literatürüne karşı konum;
  §7.1'in çözünürlük argümanı alanın kendisi için de kurulmalı.
- 296 `\ph` yer tutucusu, 11 `\dnote`.

### Depo

- **`LICENSE` yok.** Public depo; karar bekliyor.
- `paper/main.log` izlenmediği için CI'da taslak-iddia denetimi atlanıyor
  (kasıtlı; lokalde `latexmk` sonrası koşuyor).

---

## M0 — İskelet ve kabul

`../codebase` arşivi salt-okunur içe alınır; bir arşiv değeri aynı kod
yolundan **birebir** yeniden üretilir.

**Geri düşüş (D39).** Birebir eşleşme olmazsa proje durmaz, ama körlemesine de
devam etmez. Sıra: (1) farkın kaynağını teşhis et (ortam, BLAS, çekirdek
sürümü, derleyici); (2) fark sınırlıysa **ilan edilmiş bir tolerans** ile
devam et ve toleransı manifeste yaz; (3) **hiçbir karşılaştırma iki kod yolunu
karıştırmaz** — bir kampanyanın bütün sayıları tek yoldan gelir. Teşhis
edilemeyen fark → dur.

---

## M1 — Kazanç var mı, ve ölçülebilir mi

Projenin değerini bu aşama söylüyor. **Hiç *politika* propagasyonu yok**,
dolayısıyla buradaki bütün hatalar öngörülmüş hatadır: `Ê = √J`, `E` değil
(NOTATION §3).

**Ama entegrasyonsuz değil (D119).** `Q` tamamen `Φ`'ye bağlı ve arşiv `Φ`'yi
saklamıyor; WP1 her yörünge için 42 durumlu varyasyonel sistemi koşmak
zorunda. Merkezi-fark gradyanı çağrı başına 6 sentez eklediği için bir
varyasyonel yay ≈ **7× düz propagasyon** (eşit gradyan derecesinde).
26 yörüngelik panelde bu **~180 düz-propagasyon karşılığı** — M1'in bir
maliyeti; **daha büyüğü tablo kurulumu** (D137): `M × |𝒩|` sentez, `N=300`'de
panel başına ~1.3 saat, `N=600`'de **~8.4 saat** (8 çekirdek). Kümülatif
çekirdek giriş noktası (D120) bunu ~60 kat düşürür.

**Manifeste iki sayı yazılır:** varyasyonel eşdeğer ve tablo sentez sayısı.
Gradyan derecesi referansın altına çekilirse birincisi düşer — bu tam olarak
Q13'ün ölçtüğü ödünleşim; ikincisi ise `N_ref` ve `|𝒩|` ile belirlenir ve
yalnızca D120 ile düşer.

| Kapı | Ölçüt | Geçilemezse |
|---|---|---|
| **G1** | `ρ̂(R-int, A-sign)` tasarım medyanı ≥ 1.5 | 🟡 orta senaryosu: makale bir tahsis-kıyası + negatif sonuç notuna döner |
| **G2** | `A-sign` öngörülen kazancı / sayısal zarf ≥ 3. **Zarf, aynı yörüngenin arşivlenmiş R14 kaydından alınır** — M1'de yeni propagasyon olmadığı için kendi zarfını üretemez | enstrüman ölçemeyecek; bütçe veya yay uzunluğu gözden geçirilir. **Sonuç değil, enstrüman arızası** |
| **G3** | `A-sign` çizelgesi ızgara inceliğiyle kararlı — uyarlanabilir ızgarada ölçüt `n_s = τ_corr/Δt_i`, düzgün ızgarada `Δt_acc` (WP4, D68) | tüm kampanya yakınsak ızgarada koşulur; 120 s yalnızca karar ızgarası kalır |
| ~~**G4**~~ | **kapı değil, tanı (D106).** `κ` medyanı ve korunan kesir raporlanır; `κ < 0.7` bir **uyarıdır**, durdurma değil | Yöntem G5b'ye kadar yaşar; gerçek uçurulabilirlik kapısı `f̂ ≥ 0.15` |

**M1'in ürettiği dört sayı** (üçü kapı, biri tanı):

1. `ρ̂(R-int, A-sign)` medyanı — **G1**
2. `G_sign / G_sens` — sadeleşme, ağırlığın üstüne ne katıyor
3. `κ` ve korunan kesir — **tanı**, kapı değil (D106)
4. `T_coh` — pilot planın ömrü; kısa çıkarsa WP11 zorunlu

---

## M2 — Denetleyici tasarımı ve eleme (propagasyonsuz)

WP7 (sertifika), WP8 (pilot STM), WP9 (IFBDA), WP10 (`C-fb`), WP11 (koşullu),
WP12 (biçim kontrolleri), WP13 (varyasyonel eleme).

**Maliyet uyarısı — bu aşama ucuz değil.** Her varyasyonel çözüm bir
propagasyona mal oluyor. 6 aday × 26 yörüngelik panel = 156 çözüm, artı **14**
ablasyon × 8 yörüngelik küçük panel = 112 çözüm → **~270 propagasyon
eşdeğeri**. Ablasyonların küçük panelde koşulması bilinçlidir ve panel
büyüklüğü her tabloda yazılır.

| Kapı | Ölçüt | Geçilemezse |
|---|---|---|
| **G5** | Eleme kuralı **koşudan önce** yazılır; **tam iki** uçurulabilir aday M3'e geçer. Post-hoc aday seçimi yasak. | — (usul kuralı) |
| **G5b** | En iyi adayın **öngörülen** yakalama oranı `f̂ ≥ 0.15`, **`B+` dahil** (prob ek yükü **dereceyle doğrusal**, D141: vekil yayda `N=120`→%8.5, `N=300`→%21, `N=600`→%42) | 🟡 bandı: masada ödül var ama uçurulabilir kısmı yok. M3'e ~670 yay harcanmaz; makale kıyas + bilgi-sınırı biçiminde yazılır |

`f̂ = 0.15` eşiği, H3'ün hedefi olan `f ≥ 0.33`'ün altındadır — bilerek. Amaç
umut vaat eden bir adayı erken elememek, yalnızca hiç yakalayamayan bir aday
için propagasyon harcamamak.

---

## M2b — Mimari dallanma (WP21, 0 yay)

**Statü ayrımı (D116).** `A-sign` ile `C-plan` aynı güven düzeyinde değil:

| Katman | Güven |
|---|---|
| Yörünge hatasını işaretli, STM-ağırlıklı fonksiyonel olarak kurmak | çok yüksek |
| Problemin ayrıştırılamaz olduğu | çok yüksek |
| `Q` / çekirdek formülasyonu | yüksek |
| Koordinat inişinin **en iyi çözücü** olduğu | orta → **T1/T2 ölçer** |
| `A-sign`'ın gerçek doğrusal-olmayan optimuma yakınlığı | bilinmiyor → sertifika + `A-sign+` |
| Bant probunun **en iyi bilgi edinme** yöntemi olduğu | düşük/orta → **T6/T7 ölçer** |
| IFBDA'nın **en iyi uçurulabilir denetleyici** olduğu | orta → **T3/T4/T5 ölçer** |
| `C-plan`'ın küresel olarak en iyi yöntem olduğu | **hayır** |

`C-plan`'ın görevi ilk **yorumlanabilir, denetlenebilir, düşük-karmaşıklıklı**
denemedir. M1/M2 bunun yerine başka bir aile gerektiğini söylerse mimariyi
oraya çevirmek bilimsel olarak doğrudur ve **bir başarısızlık değildir**.

**Dallanma tablosu — koşudan önce yazıldı.** Gözlemler WP21'in yedi testinden
ve M1/M2 kapılarından gelir; hiçbiri ek yay istemez.

| Gözlem | Eşik | Ne söyler | Mimari |
|---|---|---|---|
| G1 düşer | `ρ̂(R-int,A-sign) < 1.5` | masada ödül yok | **dur** — makale tahsis-kıyası + negatif sonuç |
| **T3a+T3b** birlikte | `p ≤ 2` **ve** medyan asal açı `< 15°` | sadeleşme durumu `p` skaler **ve** altuzay yay boyunca tutarlı | **`C-rank`** — `p` skaler, çevrimiçi güncellenir; IFBDA yerine |
| **T3a** geçer, **T3b** düşer | `p ≤ 2` ama açı `≥ 15°` | rank düşük ama baskın yön dönüyor; 1-B bilgi her epokta başka yönde | **tam durum** — `C-plan`/IFBDA |
| **T3a** | `p = 3` | ara bölge; indirgeme var ama zayıf | **tam durum**, `C-rank` gelecek-iş satırına yazılır |
| **T3a** | `p ≥ 4` | tam 6×6 eşleşme gerekli | `C-plan`/IFBDA doğru tasarım |
| **T3** ölçüm tabanı | `Ã_i = Sᵀ A_i S` | `A_i` karışık birimli bir kuadratik form; ham iz/özdeğer/özyön birim sistemine bağlı | boyutsuzlaştırılmış formda ölçülür, ham `A_i`'de değil |
| **T4** ufuk dizi | `H* ≤ 1 devir` | eşleşme kısa menzilli | kısa-ufuklu MPC ucuz; ama `C-plan` da yeter |
| **T4** | `H*` yayın yarısını aşar | eşleşme yay-boyu | çevrimdışı plan zorunlu; MPC pahalı |
| **T5** `c` duyarlılığı | profil uzaklığı > G3 eşiği | dondurulmuş `c` yetmez | çevrimiçi eşdurum: adjoint/MPC veya `C-rank` |
| **T6a** düzgün-durum uyumu | `≥ 0.60` | kazanç doku değil **geometri** | prob gereksiz; ucuz durum-geri-besleme politikası |
| **T6a / T6b** | `T6a ≤ 0.25`, T6b yükseliyor | kazanç gerçekten dokuda | prob veya vektör-vekil; **T7 Pareto'su** seçer |
| **T1** boşluk ayrışması | iniş terimi baskın | çözücü zayıf | `S-round`/MIQP; kanonik çözücü değişebilir |
| **T1** | gevşetme terimi baskın | çözücü iyi | çözücüye yatırım yok; `g_E` muhafazakârlık olarak raporlanır |
| `T_coh` (WP6) | `≪ T` | plan bayatlıyor | alıcı-ufuk gerçekten gerekli (WP11 tetiklenir) |
| `A-sign+` ≠ `A-sign` | ciddi fark | doğrusallaştırma bağlayıcı | yeniden-planlama mimarisi |

**Usul kuralı.** Bu tablo M3'ten **önce** okunur ve seçilen aile yazılı olarak
ilan edilir. Dallanma post-hoc yapılamaz; M3'ün ~670 yayı ilan edilen aileye
harcanır. Seçilmeyen dallar "koşulmadı" olarak raporlanır, sessizce düşmez.

**Ne yapılmaz.** MPC, öğrenilmiş politika, MIQP ve vektör-vekil bu aşamada
**uygulanmaz** — yalnızca gerekip gerekmedikleri ölçülür. Uygulama kararı
tabloya bakılarak verilir ve bedeli o zaman bütçelenir.

---

## M3 — İlk propagasyon: onay/ret

Tasarım A, β = 1, 64 yörünge. **Propagate edilen tam liste — beş politika:**

| Politika | Neden |
|---|---|
| Aday 1, Aday 2 | M2'den geçen iki uçurulabilir denetleyici |
| **`A-sign`** | H1b'nin ve `f`'in paydası; propagate edilmeden ikisi de ölçülemez |
| `F-op` | iddianın şartı olan komparatör |
| **`R-int`** | mevcut en iyi uçurulabilir üye; H4b'nin komparatörü |

`F-env` propagate edilmez: sabit aile taramasından alt zarf olarak kurulur.

Buna ek olarak **WP17a**: `A-sign+` sabit-nokta testinin 16 yörüngelik indirgenmiş
hâli. Bu M5'ten M3'e alındı (D40), çünkü `A-sign`'ın bir *yörünge seviyesi*
kıyas mı yoksa bir *referans-yay* optimumu mu olduğunu belirliyor — ve bu,
makalenin merkezî nesnesinin ne olduğu sorusudur. Popülasyon genişletmeden
önce bilinmesi gerekir.

**Bütçe ve eşleştirme (D91).** `B_tot` = `F-op(β)`'nın o yaydaki
**gerçekleşen** işi; β=1'de arşivlenmiş kritik-derece koşusu, yani çapa hazır.
**Komparatör değil, her aday çapaya kalibre edilir** — çünkü `f`'in payı ve
paydası aynı `E_F-op`'u kullanmak zorunda. `A-sign` ve `R-int`'in ek yükü yok,
tam `B_tot` harcar; `C-*` için `B2 + B+ = B_tot`.

**Yay sayısı — kalibrasyon dahil.** Gerçekleşen iş ancak propagasyondan sonra
bilindiği için her aday çapaya **2–3 propagasyonda** yakınsıyor (%2 bandı).
Dört kalibre edilecek aday (Aday 1, Aday 2, `A-sign`, `R-int`) × 64 yörünge ×
~2.5 ≈ **640**, artı `F-op` çapası (arşivden, 0), artı WP17a (~32)
≈ **~670**. (Erken planlarda "~430" yazılmıştı; kalibrasyon iterasyonu
sayılmamıştı — D91.)

| Kapı | Ölçüt | Geçilemezse |
|---|---|---|
| **G6** | Varyasyonel tahmin ile propagate sonuç **işaret olarak** uyuşuyor | Uygulama hatası aranır → M2'ye dön. **Kaçış (D41):** iki teşhis turundan sonra hâlâ uyuşmuyorsa, eleyici bu denetleyici sınıfı için *kalibre değil* ilan edilir; M2'nin eleme sonucu iptal edilir ve M3 tek aday ile, yalnızca propagasyona dayanarak yürür. Maliyet sonucu açıkça yazılır. |

---

## M4 — Bütçe probu (yeni aşama, D42)

**Neden burada:** β=1 tek bir noktadır. Önceki makalede iç üyenin bütün
yapıcı sonucu β=0.5'te ters dönmüştü. Üç popülasyona 1800 yay harcamadan
önce, kazancın komşu bütçelerde de durup durmadığı bilinmelidir.

Tasarım A'nın **perilun-tabakalı 16 yörüngelik** alt kümesi (aynı sıralama
kuralı), β ∈ {0.50, 0.75, 1.25, 1.50}.

**Yay sayısı (D98 formülü):** 4 kalibre aday × 16 × 4 × ~2 (β=1'den sıcak
başlatma) = 512, artı **çapa** — `F-op(β)` β≠1'de arşivde yok, 16 × 4 = 64
propagasyon → **~576**.

| Kapı | Ölçüt | Geçilemezse |
|---|---|---|
| **G7** | `C-plan`, **hem** β = 0.75 **hem** β = 1.25'te çözünen çoğunluğu alıyor | Kazanç β=1'e özgüdür → M5 **iptal**, sonuç **🟡** bandı (D64). Geçtiği hâlde M5 kaynak nedeniyle koşulmazsa **🟩'nin alt ucu**: gerçek ama tek popülasyonda gösterilmiş |

---

## M5 — Popülasyon genişletme

WP15. En büyük compute kalemi ve **yalnızca G7 geçilirse** koşulur.
**Üç kalibre aday (D108):** `C-plan`, `R-int`, `F-op` çapası. `A-sign` burada
propagate edilmez; sonucu, `f` ve H3'ün yalnızca tasarım A/B/C'de
raporlanması.

| Popülasyon | Neden zorunlu |
|---|---|
| Tasarım B | tekrar — yön tek karıştırmaya ait değil |
| Tasarım C | bağımsız üçüncü kapsam tasarımı |
| Geniş-eliptik (Kaguya benzeri) | önceki makalede radyal kuralın **kazandığı** rejim |
| Beş geometri stratumu | geometri bağımlılığı |
| Düşük perilun 31–50 km | doğrusallaştırmanın en zayıf olduğu yer; yalnız işaret okunur |

Bir yön ≥3 bağımsız popülasyonda tekrar etmiyorsa iddia rejim-sınırlı yazılır.

Eşzamanlı koşularda işçi sayısı kısılır ve süreç BelowNormal önceliğe alınır.

---

## M6 — Tam bütçe ızgarası

WP20b. Tasarım A'da 64 yörüngenin tamamı, β ∈ {0.50, 0.75, 1.25, 1.50}
(β=1 M3'ten). Tasarım B yalnızca çaprazlamaya en yakın tek değerde — önceki
makalenin sahne kuralı. **β = 3 propagate edilmez**, enstrümanı tüketiyor.

Kritik hücre **β = 0.5** (H5): önceki makalede her uçurulabilir alternatifin
kaybettiği yer.

---

## M7 — Maliyet muhasebesi

WP16 (maliyet merdiveni + `B+` + ölçülen zaman) ve WP17b (`A-sign+` tam
panel). `B+` dipnota atılmaz.

**Güncel `B+` modeli (D50/D58) — prob baskın kalem:**

| Kalem | Tahmin | Kaynak |
|---|---|---|
| Pilot yay | %8 | önceki makalede ölçülmüş |
| **Çevrimiçi ileri prob** | **~%12–19** (yay integrali, D67) | `[∫dt/τ_corr]/N_RHS`; WP5 ölçer |
| Anahtarlama (5040 sınır) | ~%2 | önceki kampanyadan: ~0.45 çağrı/anahtarlama |
| Konum öngörücüsü (Kepler) | ~0 | analitik |
| IFBDA planlaması | çevrimdışı, küçük | — |

`C-lite` varyantında prob ~%2'ye iner, karşılığında yön yalnız aralık
başında bilinir. Maliyet–isabet eğrisi WP5'ten gelir ve **M2'nin eleme
ölçütüne girer**: aday sıralaması `B1 + B+` altında yapılıyor, dolayısıyla
%21'lik bir prob ek yükü G5b eşiğini doğrudan zorlaştırıyor.

---

## M8 — Uzun yay ve dayanıklılık

WP18. 60 günlük panel önceki makalenin yapıcı sonucunun çöktüğü yer;
alıcı-ufuklu varyant (WP11) varsa asıl sınavı burasıdır.

---

## M9 — Manuskript

`konsey` → `literatur` → `gonderim`. Hedef JGCD.

---

## Ön-tescil planı

| Tescil | Kapsam | Ne zaman |
|---|---|---|
| `oa01_preregistration.json` | H1a, H2, H6, H7; G1–G4 eşikleri; `Ê`/`E` ayrımı; ızgara, prob ve pencere taramaları; sansür ve çözünürlük kuralı | M1 öncesi |
| `oa02_preregistration.json` | aday listesi, eleme kuralı, G5b eşiği, M3'e geçen aday sayısı, sertifika metriği ve dil kuralı | M2 öncesi |
| `oa03_preregistration.json` | propagasyon matrisi (beş politika), komparatör katmanları, bütçe eşleştirme denklemi, G6 kaçışı, G7 eşiği, popülasyon listesi, ızgara genişletme kuralı | M3 öncesi |

---

## Hipotezler

Hangi hata düzeyinde olduğu **her satırda** yazılıdır: `Ê` öngörülen
(propagasyonsuz), `E` propagate edilmiş.

| # | Hipotez | Düzey | Bütçe | Aşama |
|---|---|---|---|---|
| **H1a** | `A-sign`, eşit `B1`'de `F-op`'u yener | `Ê` | `B1` | M1 |
| **H1b** | Aynı kazanç eşit `B2`'de propagasyonda hayatta kalır | `E` | `B2` | M3 |
| **H2** | `G_sign > G_sens` (NOTATION §3) | `Ê` | `B1` | M1 |
| **H3** | `f ≥ 0.33`, `f`'in tanımlı olduğu yörüngelerde | `E` | `B2`+`B+` | M3 |
| **H4** | `C-plan`, `F-op`'u β ∈ [0.75, 1.25]'te yener | `E` | `B2`+`B+` | M3–M6 |
| **H4b** | `C-plan`, **`R-int`**'i aynı aralıkta yener | `E` | `B2`+`B+` | M3–M6 |
| **H5** | Aynısı β = 0.5'te *(en riskli; ayrı, reddi H4'ü götürmez)* | `E` | `B2`+`B+` | M6 |
| **H6** | `abl-phi` kazancın yarısından fazlasını kaybeder | `Ê` | `B1` | M2 |
| **H7** | `C-plan`, yalnız-genlik + rastgele-işaret null modelini yener | `Ê` | `B1` | M2 |

**H4b neden ayrı bir hipotez (D43):** OUTCOMES'un 🟢 bandı `R-int`'in de
yenilmesini şart koşuyordu ama bunu test eden tescilli bir hipotez yoktu.
Bir bant sınırının, ön-tescilde karşılığı olmayan bir iddiaya dayanması
kabul edilemez.

`F-env`'i yenmek **hiçbir hipotezin şartı değildir**; ayrı ve bonus olarak
raporlanır.

Reddedilen hipotezler yeniden yazılmaz, **reddedilmiş olarak raporlanır**.

---

## Kaba maliyet

**Maliyet formülü (D98) — her aşamada aynı uygulanır.** D91'den sonra bir aday
tek propagasyonla bitmiyor; çapaya kalibre olması gerekiyor:

```
aşama yayı = (kalibre aday) × (yörünge) × (β noktası) × n_cal
           + (çapa: F-op(β) propagasyonları; β=1 arşivden, 0)
n_cal = 2.5 ilk kalibrasyon  |  2.0 komşu β'ya sıcak başlatma
```

Bu formül **yalnızca M3'e** uygulanmıştı; M4, M5 ve M6 "politika başına tek
propagasyon" aritmetiğinde kalmış ve β≠1'deki çapa hiç sayılmamıştı. Aşağısı
formülün tutarlı uygulanmış hâli.

| Aşama | Yay / çözüm | Duvar saati | Risk |
|---|---|---|---|
| M0 | 0 | saatler | yok |
| M1 | ~10 pilot | 2–3 gece | düşük |
| M2 | ~270 varyasyonel | 3–4 gece | orta |
| M3 | ~670 | 3–4 gece | orta |
| M4 | ~576 | 2–3 gece | orta |
| M5 | ~3060 | 6–9 gece | **yüksek** |
| M6 | ~1730 | 5–7 gece | yüksek |
| M7 | ~150 seri | 2 gece (+boşta makine) | orta |
| M8 | ~40 uzun yay (her biri ~8.5×) | 2–3 gece | yüksek |

**Toplam ~6230 yay** (M3 672 + M4 576 + M5 3060 + M6 1730 + M7 150 + M8 40),
kabaca **dört–altı haftalık gece koşusu**, kapılar geçilirse. G7'de durursa M5
ve M6 düşer ve toplam ~1400'e iner — kapının koruduğu miktar budur.

**M5 kapsamı kapatıldı (D108).** Popülasyon aşamasında **üç** aday kalibre
edilir: `C-plan`, `R-int` ve `F-op` çapası. `A-sign` popülasyon düzeyinde
propagate **edilmez** — RQ1 ve yakalama oranı `f` tasarım A, B ve C'de
kuruluyor; geniş-eliptik, stratumlar ve düşük-perilun yalnızca **verdikti**
(`C-plan` vs `F-op` vs `R-int`) test ediyor. ~920 yay tasarruf.
**Bedeli açıkça yazılır:** `f` o popülasyonlarda raporlanamaz, ve H3 yalnızca
A/B/C üzerinde sınanır.

**Depolama (D44, D55/D61 ile revize).** Düzgün 10 s'lik ızgarada tablo yörünge
başına ~87 MB olurdu; **uyarlanabilir ızgara bunu ~9 MB'ye indiriyor** (D55),
dolayısıyla **float32'ye gerek yok — `float64` saklanır** (D61). Kalan çözüm
**bellek-eşlemeli sıralı erişim**. Koordinat inişinin ihtiyaç duyduğu önek/sonek toplamları
`M × 6` boyutunda (~3 MB) ve bellekte tutulur; `u_i(N)` tablosu her süpürmede
diskten sırayla okunur. Tabloyu tamamen belleğe almak gerekmiyor. Bu, M1'de
çözülür — M5'te değil.