# PREREG — ön-tescil taslakları

Bu dosya **taslaktır**. Her tescil ilgili aşamadan önce JSON'a yazılıp
hash'lenir ve o andan sonra değiştirilmez. Sonradan eklenen her şey
**post-hoc** olarak etiketlenir ve ayrı raporlanır.

Önceki makalenin dersi: skorlama kuralının ne zaman yazıldığı, kuralın
kendisinden daha çok sorgulanıyor (R37 skorlama düzeltmesi post-hoc kalmak
zorunda kaldı). Bu yüzden skorlama kuralı burada, koşudan önce, açıkça yazılı.

---

## OA-01 — Tahsis kıyası (M1 öncesi)

### Hipotezler
- **H1a.** `A-sign`, eşit `B1`'de (kısıt `Σ_q W_q N_q² ≤ B = B1·T`, D60)
  `F-op`'tan küçük **öngörülen** hata `Ê = √J` verir, her tasarımın
  yörüngelerinin çoğunluğunda. **Bu bir `Ê` hipotezidir**: M1'de propagasyon
  yoktur ve `E` ölçülemez (NOTATION §3).
- **H2.** `G_sign > G_sens`, log-iyileşme adımlarıyla
  `G_sens = log(Ê_A-force / Ê_A-sens)`, `G_sign = log(Ê_A-sens / Ê_A-sign)` —
  **öngörülen** hata, M1'de propagasyon yok.
  (İşaretli sadeleşme, duyarlılık ağırlığından daha çok kazandırır.) Ham hata
  farkı yörünge ölçeğine bağlı olduğundan kullanılmaz.
- **H6.** `abl-phi` (`(T−t)` ağırlığı, `Φ` yok) `A-sign` kazancının yarısından
  fazlasını kaybeder.

### Sabitlenen tanımlar
- Amaç: NOTATION §3'teki `E²`, konum seçicisi `H_r = [I₃ 0]` açık.
- **Hata düzeyi:** M1 ve M2'nin bütün nicelikleri `Ê = √J` (öngörülen);
  `E` (propagate) yalnızca M3 ve sonrasında. Oranlar `ρ̂` / `ρ` ile ayrılır ve
  her tabloda hangisi olduğu yazılır.
- **Karar uzayı:** `(N_1,…,N_Kdec)`, karar ızgarasında (NOTATION §4).
  Biriktirme ızgarasında `M` boyutlu bir optimizasyon **yapılmaz**; yapılırsa
  `A-sign` uçurulabilir adaydan yalnızca bilgi değil, anahtarlama serbestliği
  bakımından da ayrışır ve karşılaştırma geçersiz olur.
- Bütçe: `β = B1/N_crit²` (yalnız `B1` düzeyinde), ızgara
  β ∈ {0.5, 0.75, 1.0, 1.25, 1.5, 3.0}. Gerçekleşen normalizasyon `β₂` ayrı
  raporlanır ve `β` yerine kullanılmaz.
- Politikalar ve komparatörler: NOTATION §1.
- Derece kümesi `𝒩`: `rev14_oracle.degree_grid` ile aynı + karşılaştırılan
  politikaların kendi dereceleri eklenir.
- Çözünürlük: `M_res > 1`; kararsızlar hiçbir tarafa sayılmaz.
- Sansür: politika referans derecesine ulaşırsa o hücre sansürlenir,
  komparatör referans derecede kıstırılmaz.
- İstatistik: **p-değeri yok.** Yön, bağımsız tasarımlarda tekrar ile kurulur.
- Oran istatistiği: yörünge başına oranların **medyanı**.

### Taramalar (koşudan önce ilan)
- Izgara yakınsaması iki parametrede taranır (D68):
  **düzgün** ızgarada `Δt_acc ∈ {240, 120, 60, 30, 10}` s, ve **uyarlanabilir**
  ızgarada korelasyon-zamanı başına örnek sayısı
  `n_s = τ_corr(t)/Δt_i ∈ {0.5, 1, 2, 4, 8}`. Kampanya uyarlanabilir ızgarayla
  koşacağı için yakınsama ölçütü `n_s`'tir; mutlak `Δt_acc` yalnızca önceki
  makaleyle karşılaştırılabilirlik için taranır. Panel: 8–14 yörünge,
  perilun-yayılı.
- `Δt_dec ∈ {60, 120, 300}` s.
- Prob derinliği `k ∈ {1, 2, 3, 4}`.
- Aday penceresi yarı-genişliği `δ ∈ {2, 4, 8}` — **aday-ızgara indisi**
  cinsinden (D32), derece cinsinden değil. İndüklenen derece açıklığı
  `Δ_span` her satırda ayrıca raporlanır.
- Prob **ileri** alınır (D26): sınırdan önceki hiçbir yön bilgisi yeniden
  kullanılmaz. `abl-lookahead` bunun gerekliliğini ölçer.
- **Prob maliyeti bir eğri olarak raporlanır (D50/D51):** `κ_eff(n_probe)` ve
  karşılığındaki ölçülen maliyet. İki uç koşudan önce ilan edilmiştir:
  `C-plan` (tam kapsayan ileri prob) ve `C-lite` (eş-konumlu tek prob).
  Maliyet modeli `1/(τ_corr·r_RHS)`; ölçüm buna uyuyor mu ayrıca raporlanır.
- **Genlik tamamlaması `γ` gerçek `P_n` spektrumundan (D53)**; güç yasası
  varyantı `abl-kaula` olarak ayrıca koşulur.
- **Plan faz-indekslidir (D52)**; zaman-indeksli varyant `abl-timeindex`
  olarak koşulur ve pilot yay zaman kaymasının bedelini ölçer.
- **Biriktirme ızgarası uyarlanabilir (D55)**; düzgün ızgara varyantı WP4'te
  ayrıca koşulur ve `B1` zaman-ağırlıklı tanımlanır.
- `λ` monotonluk kontrolü (D33): yoğun log ızgarada `W(λ)` süpürülür; monoton
  değilse bisekt yerine ızgara üzerinde en küçük `J`'li olurlu nokta.
- Prob konumu öngörücüsü: iki-cisim propagasyonu; `abl-predictor` ile
  düşük-dereceli mikro-propagasyona karşı ölçülür (D31).
- Koordinat inişi başlangıçları: `F-op`, `R-rad`, `R-int`, `A-sens` + 8
  rastgele tohum (tohumlar aritmetik kuralla türetilir, tek çekiliş).

### Karar kapıları (eşikler önceden)
- **G1:** `ρ̂(R-int, A-sign)` tasarım medyanı ≥ 1.5 (öngörülen hata; M1'de
  propagasyon yok).
- **G2:** `A-sign` öngörülen kazancı / sayısal zarf ≥ 3. **Zarf, aynı
  yörüngenin arşivlenmiş R14 kaydından alınır** — M1'de yeni propagasyon
  olmadığı için kendi zarfı üretilemez. Bu, zarfın dışsal ve önceden var olan
  bir nicelik olması bakımından önceki makalenin skorlama düzeltmesiyle aynı
  gerekçeye dayanır.
- **G3:** ızgara yakınsaması. Kampanya uyarlanabilir ızgarayla koştuğu için
  ölçüt **`n_s = τ_corr/Δt_i`** (D68); düzgün ızgaradaki `Δt_acc` taraması
  yalnızca karşılaştırılabilirlik için. Çizelge profil uzaklığı eşiği tescilde
  sayı olarak yazılacak.
- ~~**G4**~~ **kapı olmaktan çıkarıldı, tanı oldu (D106).** `κ` medyanı ve
  korunan kesir raporlanır; `κ < 0.7` **uyarıdır, durdurma değil**. Gerekçe:
  `κ` korunan kesrin kendisi değil, dolayısıyla ona sert bir kapı bağlamak
  erken yanlış-negatif üretir — `κ = 0.6` olup hatanın tamamı `z`'ye dik
  olabilir (kesir iyi), ya da `κ = 0.9` olup küçük açısal hata tam `z`
  doğrultusunda projeksiyonu bozabilir. Uçurulabilirlik kararı zaten **G5b**
  (`f̂ ≥ 0.15`) ve o doğrudan "ödülün ne kadarı yakalanıyor"u ölçüyor.
  **`κ` yine de raporlanır** (D74): sadeleşme teriminin korunan kesri `⟨v̂,z⟩/⟨Δa,z⟩`
  (`z = (ΦB)ᵀc`) ve `ΦB` açıları korumaz. Kesir, `A-sens` çizelgesinden
  türetilen `c` ile aynı koşuda hesaplanıp `κ`'nın yanında raporlanır.
  **İkisi de tanısaldır**; devam/durdurma kararı G5b'deki `f̂` üzerinden verilir.
  `κ`'yı önde gelen gösterge olarak tutmamızın nedeni kapı olması değil,
  propagasyonsuz ve çizelgeden bağımsız ölçülebilmesidir.

### Panel seçim kuralı
Perilun irtifasına göre sıralanıp eşit aralıklı rütbelerden seçim —
önceki makalenin ölçülen-zaman ve varyasyonel panelleriyle **aynı kural**.
Sonuç üzerinde tabakalama yok.

### Raporlama taahhüdü
Reddedilen hipotez yeniden yazılmaz. Ölçülemez hücreler (G2) sayılır ve
raporlanır, sessizce düşürülmez. Koordinat inişinin başlangıçlar arası
yayılımı her zaman raporlanır.

---

## OA-02 — Denetleyici eleme (M2 öncesi)

### Sabitlenen liste
Elenecek uçurulabilir adaylar, koşudan önce tam liste:
`C-plan(J=1)`, `C-plan(J=2)`, `C-plan(J=3)`, **`C-lite`**, `C-fb`, `C-tgo`,
[koşullu] `C-rh` (WP11 tetiklenirse).

`C-lite` (D51) listeye eklendi: NOTATION ve PLAN'de ilan edilmişti ama aday
listesinde yoktu. Prob maliyet–isabet eğrisinin ucuz ucu olduğu için elemeye
girmesi zorunlu — aksi hâlde eğrinin bir ucu hiç propagasyon şansı bulamaz.

### Eleme kuralı
- Skorlayıcı: zorlanmış varyasyonel çözüm. **Paneller koşudan önce ilan
  edilmiştir:** adaylar 26 yörüngelik perilun-tabakalı panelde, **14 ablasyon
  8 yörüngelik küçük panelde** (D79). Toplam ~270 propagasyon eşdeğeri.
  Ablasyon paneli küçüktür ve her tabloda büyüklüğü yazılır; ablasyon satırları
  birbiriyle karşılaştırılabilir, merdiven tablosuyla **değil**.
- Sıralama ölçütü: `B2` **tahmini** değil — M2'de `B2` yok — `B1 + B+` altında
  öngörülen `E`'nin yörünge-medyanı.
- **Tam iki aday** M3'e geçer. Beraberlik hâlinde daha düşük `B+` olan geçer;
  o da eşitse daha düşük `J`.
- Panel sonrası aday **eklenemez**.

### Sertifika metriği ve dil kuralı
Metrik dondurulmuştur (D29): `g_J = (J_desc − L_FW)/J_desc`,
`g_E = 1 − √(L_FW/J_desc)`; `L_FW = max{0, en iyi FW sınırı}`.
**Eşik hata uzayındadır: medyan `g_E < 0.10`** (≈ `g_J 0.19`). `L_FW = 0`
kalan yörüngeler "boş sertifika" olarak ayrı sayılır ve eşiğe girmez.
Bu eşik sağlanırsa `A-sign` metinde "oracle" olarak adlandırılabilir. Aksi hâlde "linearized trajectory-aware allocation
benchmark" ve *erişilebilir tahsis, sınır değil* dili kullanılır. Bu karar
sonuca değil, boşluk sayısına bağlıdır ve burada peşinen bağlanmıştır.

### Karar kapısı G5b
En iyi adayın **öngörülen** yakalama oranı `f̂ ≥ 0.15`. Altında M3'e propagasyon
harcanmaz; proje 🟡 bandına geçer ve makale kıyas + bilgi-sınırı biçiminde
yazılır. Eşik, H3'ün hedefi olan `f ≥ 0.33`'ün bilerek altındadır: amaç umut
vaat eden bir adayı erken elememek, yalnızca hiç yakalayamayan bir aday için
~670 yay harcamamaktır.

### Hipotezler
- **H7.** `C-plan`, yalnız-genlik + rastgele-işaret null modelini yener (`Ê`).

### Mimari dallanma testleri — WP21 / T1–T7 (D116)

**Statü beyanı, koşudan önce.** `A-sign` yörünge-hatası fonksiyonelinden
türetildiği için kıyas olarak savunulur. **`C-plan`'ın en iyi uçurulabilir
çözüm olduğu iddia edilmez**; ilk yorumlanabilir, denetlenebilir,
düşük-karmaşıklıklı denemedir. Aşağıdaki eşikler M3 öncesinde okunur ve
seçilen denetleyici ailesi **yazılı olarak** ilan edilir; dallanma post-hoc
yapılamaz.

| Test | Ölçülen | Önceden ilan edilen eşik |
|---|---|---|
| **T1** | `J_desc − J*` ve `J* − L_FW` (kaba kuvvet/MIQP küresel optimumla) | iniş terimi toplamın **yarısından** fazlaysa çözücü değiştirilir |
| **T2** | `S-round` çizelgesinin `J`'si | çok-başlangıçlı inişten iyiyse **kanonik çözücü değişir**, iniş kontrole düşer |
| **T3a** | `A_i` izinin %95'i için gereken özyön sayısı `p` | `p ≤ 2` **gerekli koşul**, tek başına yeterli değil |
| **T3b** | baskın altuzayın zaman tutarlılığı: ardışık epoklar arasında asal açı `∠(U_p(t_i), U_p(t_{i+1}))` | medyan asal açı, bir korelasyon zamanı boyunca birikmiş hâliyle **< 15°** ise altuzay tutarlı; `C-rank1` dalı **ancak T3a ve T3b birlikte** geçerse açılır |
| **T4** | `ρ̂(H)/ρ̂(T)`, diz noktası `H*` | `≥ 0.90`'a ulaşan en küçük `H` = `H*`; `H* ≤ 1 devir` → kısa ufuk yeter |
| **T5** | pilot–referans `c` pertürbasyonu altında çizelge profil uzaklığı | G3 eşiğinin üstündeyse çevrimiçi eşdurum gerekir |
| **T6a** | düzgün-durum politikasının koruduğu kazanç kesri (A'da uydur, **B'de ölç**) | `≥ 0.60` → prob gereksiz; `≤ 0.25` → kazanç dokuda |
| **T6b/T7** | vekil-model parametre sayısı vs korunan kesir | Pareto eğrisi olarak raporlanır, tek eşik yok |

**T6'nın model sınıfı da tescillidir.** Aksi hâlde sonuçlar görüldükten sonra
ağaç derinliği, düğüm sayısı, öznitelik normalizasyonu ve etkileşim terimleri
oynanabilir; o zaman test yanlışlama testi olmaktan çıkar. Koşudan önce
sabitlenen liste:

- **Model sınıfı iki tanedir ve başkası denenmez:** (i) öznitelik başına
  **kübik spline**, düğüm sayısı `{4, 8}`, toplamsal, etkileşim yok;
  (ii) **karar ağacı**, `max_depth ∈ {3, 5, 7}`, başka hiperparametre yok.
- **Öznitelikler T6a'da tam olarak** `(h, r, |v|, faz, kalan süre, harcanmış
  bütçe)`; normalizasyon her öznitelik için tasarım A'nın çeyrekler açıklığına
  bölme, sabit.
- **Seçim yalnız tasarım A'nın çapraz doğrulamasında** yapılır (5 kat,
  yörünge düzeyinde bölünmüş — aynı yörüngenin epokları iki kata dağılmaz).
- **Verdikt yalnız tasarım B'de** okunur ve **bir kez** okunur.

**Yanlışlama taahhüdü.** T6a `≥ 0.60` çıkarsa bu makalenin **bilgi iddiasının
reddi**dir: kazanç doku değil geometri kaynaklıdır, bant probu gereksizdir ve
bu sonuç yeniden yorumlanmadan, reddedildiği hâliyle raporlanır.

**Koşulmayan dallar.** MPC, öğrenilmiş politika, MIQP çözücü ve vektör-vekil
denetleyici bu kampanyada **uygulanmaz**; yalnızca gerekip gerekmedikleri
ölçülür. Seçilmeyen dallar "ölçüldü, seçilmedi" olarak raporlanır.

---

## OA-03 — Propagasyon kampanyası (M3 öncesi)

### Propagasyon matrisi
Koşudan önce tam olarak yazılır; asimetri (hangi tasarım hangi β'da) burada
ilan edilir ve tablolarda "yok" değil "koşulmadı" olarak görünür.

- M3: tasarım A, β = 1, 64 yörünge, **beş politika**: Aday 1, Aday 2,
  `A-sign`, `F-op`, `R-int`. `F-env` propagate edilmez; sabit aile
  taramasından alt zarf olarak kurulur.
- **M4 (bütçe probu):** tasarım A'nın 16 yörüngelik perilun-tabakalı alt
  kümesi, β ∈ {0.50, 0.75, 1.25, 1.50}, beş politika. G7 buradan çıkar.
- **M5 (popülasyon genişletme, G7'ye koşullu):** tasarım B ve C, geniş-eliptik
  popülasyon, beş stratum, düşük perilun. **Üç kalibre aday (D108):**
  `C-plan`, `R-int`, `F-op` çapası. `A-sign` burada propagate edilmez, bu
  yüzden **`f` ve H3 yalnızca tasarım A/B/C'de raporlanır**; diğer
  popülasyonlar H4/H4b verdiktini taşır.
- **M6 (tam bütçe ızgarası):** tasarım A'nın 64 yörüngesi,
  β ∈ {0.50, 0.75, 1.25, 1.50}; tasarım B yalnızca çaprazlamaya en yakın tek
  değerde. **β = 3 propagate edilmez.**

*(Aşama numaraları D42 ile kaydı: eski M4 = popülasyon genişletme artık M5;
M4 yeni bütçe probudur.)*

`A-sign` listede olmak zorundadır: H1b onun propagate edilmiş hatasını test
eder ve `f`'in paydası odur. Propagate edilmezse iki tescilli nicelik de
ölçülemez.

**Arşiv yeniden kullanımı.** Önceki kampanyanın `F-op` / `R-int` yayları
yalnızca gerçekleşen işleri yeni `B2` sözleşmesi altında eşleştirilebiliyorsa
yeniden kullanılır; eşleşemiyorsa yeniden propagate edilir. Her yeniden
kullanılan yay kayıtta böyle işaretlenir.

### Bütçe eşleştirme sözleşmesi (NOTATION §2, D91)
```
her aday:  B2 + B+ = B_tot ,   B_tot := F-op(β)'nın GERÇEKLEŞEN işi
```
- **Çapa önceden ilan edilir:** β=1'de arşivlenmiş kritik-derece koşusu.
- **Komparatör sabit, adaylar kalibre edilir.** Tersi (`f`'in payı ve paydası
  için iki farklı `F-op`) `f`'i tanımsız yapardı.
- **Tolerans:** `|B2+B+ − B_tot|/B_tot ≤ %2`; aday başına 2–3 propagasyon.
  Ulaşılan eşleşme her verdiktin yanında raporlanır.
- Sabit komparatöre hayali `B+` eklenmez; denetleyici kendi ek yükü kadar daha
  az gravite bütçesiyle çalışır. `A-sign` ve `R-int`'in ek yükü yoktur.

1. Birincil: yukarıdaki denklem.
2. İkincil: `B1` (nominal çağrı başına), yalnızca önceki makaleyle
   karşılaştırılabilirlik için.
3. Doğrulama: `B3` (ölçülen seri çekirdek zamanı), 14 yörüngelik panel,
   0.90–1.10 bandı, 3 tekrar, boşta makine.

### Komparatör katmanları
- `F-op` — iddianın şartı.
- `F-env` — bonus; yenilememesi başarısızlık sayılmaz. Post-hoc alt zarf
  olduğu her tabloda yazılır.

### Hipotezler
- **H1b.** H1a'nın kazancı eşit `B2`'de propagasyonda hayatta kalır
  (propagate edilmiş `E_A-sign` ile ölçülür).
- **H3.** `f ≥ 0.33`, **tasarım A/B/C'de ve `f`'in tanımlı olduğu yörüngelerde**
  (D108: diğer popülasyonlarda `A-sign` propagate edilmiyor): `f` yalnızca
  `E_A-sign < E_F-op` ve bu fark çözünürlük zarfını geçtiğinde tanımlıdır
  (NOTATION §3). Tanımsız yörüngeler `N/A` olarak sayılır ve raporlanır;
  sıfıra veya büyük bir sayıya çevrilmez.
- **H4.** `C-plan`, `F-op`'u β ∈ [0.75, 1.25]'te yener (çözünen çoğunluk).
- **H4b.** `C-plan`, **`R-int`**'i aynı aralıkta yener. *(Ayrı hipotez, çünkü
  OUTCOMES'un 🟢 bandı bunu şart koşuyor ve bir bant sınırının ön-tescilde
  karşılığı olmayan bir iddiaya dayanması kabul edilemez.)*
- **H5.** Aynısı β = 0.5'te. *(Ayrı hipotez: reddi H4'ü götürmez.)*

### Karar kapıları G6 ve G7
- **G6 (işaret uyumu) ve kaçışı:** varyasyonel tahmin ile propagate sonuç
  işaret olarak uyuşmazsa uygulama hatası aranır. **İki teşhis turundan sonra**
  hâlâ uyuşmuyorsa, eleyici bu denetleyici sınıfı için *kalibre değil* ilan
  edilir; M2'nin eleme sonucu iptal edilir ve kampanya tek adayla, yalnızca
  propagasyona dayanarak yürür. Maliyet sonucu açıkça yazılır. Sonsuz döngü
  yasak.
- **G7 (bütçe probu) ve bant kuralı (D64).** M4'te, tasarım A'nın 16
  yörüngelik perilun-tabakalı alt kümesinde:
  - **GEÇER** = `C-plan` **hem** β = 0.75 **hem** β = 1.25'te çözünen
    çoğunluğu alıyor. M5 koşulur.
  - **KALIR** = ikisinden yalnız biri, veya hiçbiri. Kazanç β=1'e özgüdür,
    M5 **koşulmaz**, sonuç **🟡** bandında raporlanır.
  - G7 geçtiği hâlde M5 kaynak nedeniyle koşulmazsa, sonuç **🟩 bandının alt
    ucu**: gerçek ama tek popülasyonda gösterilmiş.

  Bu üç durum burada tek kuralla bağlanmıştır; ROADMAP ve OUTCOMES bunu
  yansıtır ve aralarında yorum farkı bırakılmaz.

### Genişletme kuralı
Ek bütçe noktası ancak tasarım medyanlarının çaprazlaması iki komşu ızgara
değeri arasına düşerse propagate edilir; o durumda **tam olarak** o iki değer,
önce tasarım A'da. Başka hiçbir gerekçeyle bütçe eklenmez. Bu kural dışında
eklenen her nokta **post-hoc** etiketlenir ve hipotez sayımlarına girmez.

### Zorunlu kontroller (sonuç ne olursa olsun koşulur)
- Referans derecesi 300 → 600 (referans-komşuluğu artefaktı).
- Derece tavanı denetimi + sansür sayımı.
- Faz-kaydırma MC (sadeleşmeye aşırı uyum).
- Çizelge biçimi (WP12): **faz-indeksli** plan (D52) vs irtifa-binli;
  zaman-indeksli varyant ayrıca `abl-timeindex` olarak koşulur.
- **Bilgi-kaynağı denetimi (D62 ile keskinleştirildi):** her `C-*` kaydı
  kullandığı bilgiyi yazar. **Yasak olan** referans yay ve alanın yay boyunca
  **noktasal değerlendirmeleri** (`Δa_ref(t)`); **serbest olan** küresel /
  çevrimdışı model üst-verisi (`P_n` derece varyansları). "Referans alan
  görünürse geçersiz" biçimindeki eski ifade yanlıştı — bir propagatör
  katsayıları tanım gereği taşır.

---

## Ortak taahhütler

1. Hiçbir toplu sonuca bakılmadan hash'lenir.
2. Kabul kontrolü: yeni koşu, arşivlenmiş bir değeri birebir üretmeden kabul
   edilmez.
3. Manifest zinciri: her kampanya öncekini digest'e çivileyerek taşır.
4. Post-hoc her şey etiketlenir ve iki okuma (tescilli / post-hoc) birlikte
   arşivlenir.
5. Duvar saatiyle sınırlanan koşularda, saat dolduğunda ne kalırsa raporlanır;
   panel seçim kuralı iç içe (nested) olduğundan yarım kalan koşu yanlı bir
   önek değil, kapsayıcı bir panel bırakır.
