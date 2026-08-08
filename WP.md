# WP — iş paketi ve koşu sicili

Yapılacak **her şey** burada. Bir iş bu listede yoksa yapılmaz; yapılması
gerekiyorsa önce buraya eklenir. Amaç, süreç içinde neyin neden yapıldığını
kaybetmemek.

Sütunlar: **Soru** (ne cevaplıyor) · **Girdi** · **Çıktı** · **Maliyet** ·
**Bağımlılık** · **Kapı** (geçilemezse ne olur).

---

## Özet

| WP | Ad | Propagasyon | Aşama | Kritiklik |
|---|---|---|---|---|
| WP0 | İskelet + kabul kontrolü | 0 | M0 | zorunlu |
| WP1 | `Δa(t,N)` tablolaması + `Φ` içe alma | 0 | M1 | zorunlu |
| WP2 | `A-sens` (L1) çözücü | 0 | M1 | zorunlu |
| WP3 | `A-sign` (L2) koordinat inişi | 0 | M1 | **çekirdek** |
| WP4 | Izgara yakınsaması (`Δt_acc`, `Δt_dec`) | 0 | M1 | **çekirdek** |
| WP5 | Bant probu — yön isabeti `κ` | 0 | M1 | **çekirdek** |
| WP6 | Pilot yay tutarlılık ufku `T_coh` | ~10 pilot yay | M1 | **çekirdek** |
| WP7 | Frank–Wolfe alt sınırı / sertifika | 0 | M2 | yüksek |
| WP8 | Düşük-dereceli pilot STM | 0 | M2 | yüksek |
| WP9 | IFBDA denetleyicisi `C-plan` | 0 | M2 | **çekirdek** |
| WP10 | Bütçe geri beslemesi `C-fb` | 0 | M2 | orta |
| WP11 | Alıcı-ufuklu asgari varyant | 0 | M2 | **koşullu** (WP6) |
| WP12 | Çizelge biçimi kontrolleri | 0 | M2 | orta |
| WP13 | Varyasyonel eleme | ~270 çözüm | M2 | zorunlu |
| WP14 | İlk propagasyon kampanyası | ~645 yay | M3 | zorunlu |
| WP17a | `A-sign+` indirgenmiş sabit nokta | ~32 yay | **M3** | **çekirdek** |
| WP20a | Bütçe probu (16 yörünge × 4 β) | ~576 yay | **M4** | **çekirdek** |
| WP15 | Popülasyon genişletme | ~3060 yay | M5 | koşullu (G7) |
| WP20b | Tam bütçe ızgarası | ~1730 yay | M6 | zorunlu |
| WP16 | Maliyet merdiveni + `B+` + ölçülen zaman | ~150 seri yay | M7 | zorunlu |
| WP17b | `A-sign+` tam panel | ~64 yay | M7 | orta |
| WP18 | Uzun yay + dayanıklılık kontrolleri | ~40 uzun yay | M8 | orta |
| WP21 | Mimari seçim testleri T1–T7 | **0** | M2b | **çekirdek** |
| WP19 | Manuskript + denetim skilleri | 0 | M9 | zorunlu |

---

## Ayrıntılar

### WP0 — İskelet ve kabul kontrolü
- **Soru:** kod yolumuz arşivlenmiş bir sayıyı birebir üretiyor mu?
- **Girdi:** `../codebase/metrics/*` (salt-okunur), `environment-py312.yml`
- **Çıktı:** `metrics/oa00_admissibility.json`
- **Kapı ve geri düşüş (D39):** birebir eşleşme yoksa proje durmaz ama
  körlemesine devam etmez. (1) farkın kaynağını teşhis et (ortam, BLAS,
  çekirdek sürümü); (2) fark sınırlıysa **ilan edilmiş bir tolerans** ile devam
  et, toleransı manifeste yaz; (3) **hiçbir karşılaştırma iki kod yolunu
  karıştırmaz** — bir kampanyanın bütün sayıları tek yoldan gelir. Teşhis
  edilemeyen fark → dur.

### WP1 — Tablolama
- **Soru:** `Δa(t,N)` **vektör** olarak (yalnızca norm değil) ve `Φ(t_0,t_i)`
  hangi ızgarada, hangi maliyetle elde edilir?
- **Not:** önceki `rev14_oracle.py` yalnızca `d = ‖Δa‖²` saklıyordu. Buradaki
  her şey vektöre bağlı → tablo formatı yeniden tanımlanmalı. Bellek maliyeti
  `M × |𝒩| × 3` float64; `M` aşağıdaki inceltme hesabından çıkar.
- **Uyarlanabilir ızgara (D55, D129 ile yeniden kalibre).** `τ_corr = πr/(Nv)`
  tek paylaşılan inceltme derecesinde **5.6 kat** değişiyor (`N=300`: perilun
  9.5 s, apolun 53.6 s). Izgara `τ_corr`'e göre inceltilir; düzgün-en-ince
  ızgaraya karşı kazanç **2.5 kat**, eskiden sanılan 10 kat değil.
  Sonuç: `B1` **zaman-ağırlıklı** ortalama olmak zorunda.
- **İnceltme, karar aralıklarının İÇİNDE yapılır (D130).** Böylece her karar
  sınırı bir biriktirme kenarıdır; hiçbir hücre derece anahtarlamasını
  kesmez, `W_q` tam aralık uzunluğudur ve boş karar aralığı yapısal olarak
  imkânsızdır.
- **CPU maliyeti — M1'in baskın kalemi ve şimdiye kadar bütçelenmemişti (D137).**
  Tablo `M × |𝒩|` sentez istiyor: `N=300`, `M≈51 100`, `|𝒩|≈60` ile yörünge
  başına **3.07 milyon sentez**. Arşivin ölçtüğü çekirdek süresiyle
  (`e3_kernel_timing`: 454 µs @ N=300, 1508 µs @ N=600) 26 yörüngelik panel,
  8 işlemci paralel:

  | `N_ref` | sentez/yörünge | panel (8 çekirdek) | kümülatif çekirdekle |
  |---|---|---|---|
  | 120 | 1.23 M | ~0.1 sa | <0.01 sa |
  | 300 | 3.07 M | **~1.3 sa** | ~0.02 sa |
  | 600 | 6.14 M | **~8.4 sa** | ~0.14 sa |

  Yani M1'in "2–3 gece"sine sığıyor, ama STM entegrasyonundan (~180
  düz-propagasyon karşılığı) **daha büyük** ve planda hiç yoktu.
- **Bu, D120'yi maliyet iddiasından çıkarıp fizibilite kalemine taşıyor.**
  Kümülatif-dereceye-göre çekirdek giriş noktası `M × |𝒩|`'i `M`'e indiriyor:
  tek Legendre geçişi bütün adayların `a_{≤N}`'ini kısmi toplamla veriyor.
  **~60 kat**, yani `N=600`'de 8.4 saat → 8 dakika. Aynı giriş noktası prob
  maliyet modelini de açıyor. Öncelik: **WP1'den önce**.

- **Telafili toplama (D56).** Önek/sonek toplamlarında Kahan/Neumaier.
  `S_j` altmış bin işaretli terimin toplamı ve sadeleşme yöntemin ta kendisi;
  telafisiz toplamak gereksiz risk.
- **Depolama çözümü, WP1'in teslimi (D44, D60, D129 ile revize):** tablo
  `M × |𝒩| × 3 × 8` bayt; `N=300` ve `|𝒩|≈60` ile **~74 MB/yörünge**,
  `N=600`'de ~147 MB, 26 yörüngelik panelde disk ~2 GB. Tablo **`float64`
  saklanır** — küçük olduğu için değil, **bellek-eşlemeli ve sıralı** okunduğu
  için; RAM'e sığması gerekmiyor. **Gerçek `M` WP1'in ilk teslimi olarak
  ölçülür**, tahmin edilmez. Erişim yine bellek-eşlemeli ve
  sıralı; önek/sonek toplamları (`M × 6`, ~0.3 MB) bellekte.
- **Zorunlu test (D60):** `float32` vs `float64` **verdict parity** — aynı
  yörüngelerde çizelge ve verdict değişiyor mu. Makalenin merkezindeki sonuç
  işaretli sadeleşme; depolama hassasiyetine teslim edilmez. Kahan/Neumaier
  toplama, float32'ye yazarken kaybedilen bitleri **geri getirmez** — ikisi
  farklı sorunlar. Parity testi geçerse float32 yalnızca bir opsiyon olarak
  kayda geçer, varsayılan değil.
- **Çıktı:** `metrics/oa01_tables/<orbit>.npy` (memmap) + şema JSON'u

### WP2 — `A-sens`
- **Soru:** yalnızca duyarlılık ağırlığı ne kadar kazandırıyor? (RQ2'nin paydası)
- **Yöntem (D75, D110):** `argmin_N[ Σ_{i∈I_q} Δt_i · Δa_i(N)ᵀ K_i Δa_i(N) + λ W_q N² ]`,
  `λ` bütçeye bisekte. Ağırlık **yerel duyarlılık çekirdeği**
  `K_i = (1/T)·Bᵀ Φ(t₀,t_i)ᵀ A_i Φ(t₀,t_i) B`, yani amacın kendi çekirdeğinin
  köşegen üzerindeki değeri `K(t_i,t_i)` — ödünç alınmış bir skaler değil.
  **Ayrık `u_iᵀ Q_ii u_i` biçimi kullanılmaz (D110):** `u_i ∝ Δt_i` olduğundan
  o ifade `Δt_i²` taşır ve ızgara inceldikçe sistematik olarak küçülür;
  uyarlanabilir ızgarada `Δt_i` yay boyunca değiştiği için bu ayrıca yay içi
  bir eğilim de üretir. Eşdeğer yazım: `Σ_i u_iᵀ Q_ii u_i / Δt_i`.
  `K_i` aynı `A_i` sonek dizisinden okunur, ek maliyet yok.
  `A-force` aynı biçim ama `Φ`'siz (ağırlık yerine `‖Δa‖²`).
- **Onaylama testi:** `n_s` iki katına çıkarıldığında `A-sens` çizelgesi
  değişmemeli (profil uzaklığı G3 eşiğinin altında). Eski `diag Q` biçimi bu
  testi tanım gereği geçemez — regresyon nöbetçisi olarak da kullanılabilir.
- **Çıktı:** çizelgeler + öngörülen `E`

### WP3 — `A-sign` ★
- **Soru:** RQ1 — masada ne kadar kazanç var?
- **Yöntem:** PLAN §1.2 önek/sonek yapısı üzerinde blok koordinat inişi,
  çok başlangıçlı (`F-op`, `R-rad`, `R-int`, `A-sens`, +8 rastgele tohum)
- **Karar uzayı: `K_dec` boyutlu, `M` boyutlu DEĞİL** (PLAN §1.2b). Blok-içi
  terim yürüyen önek toplamıyla `O(m)`, süpürme yine `O(M·|𝒩|)`. Kodda bunun
  bir onaylama testi olsun: `Δt_dec` iki katına çıkarıldığında karar
  değişkeni sayısı yarıya inmeli.
- **λ monotonluk kontrolü (D33).** Bisektten önce `λ` yoğun log ızgarada
  süpürülür; en-iyi-başlangıç çözümünün nominal işi `W(λ)` azalan mı diye
  bakılır. Monotonluk **global** minimizörün özelliği; biz yerel çözüm
  alıyoruz, basin değişebilir. Monoton değilse bisekt bırakılır, ızgarada en
  küçük `J`'li olurlu nokta seçilir, bu yörüngeler sayılıp raporlanır.
- **Raporlanacak:** başlangıçlar arası yayılım (yerel minimum riski)
- **Kapı G1:** `ρ(R-int, A-sign)` tasarım medyanı **< 1.5** → 🟡 senaryosu
  (bkz. OUTCOMES); denetleyici geliştirmeye büyük compute harcanmaz.
- **Kapı G2:** öngörülen kazanç yörüngenin sayısal zarfının altındaysa o hücre
  ölçülemez olarak işaretlenir.

### WP4 — Izgara yakınsaması ★
- **Soru:** `A-sign`'ın kazancının ne kadarı ayrıklaştırma artefaktı?
- **Yöntem, iki parametre (D68):** kampanya **uyarlanabilir** ızgarayla
  koşacağı için asıl yakınsama ölçütü **korelasyon-zamanı başına örnek sayısı**
  `n_s = τ_corr(t)/Δt_i ∈ {0.5, 1, 2, 4, 8}`; mutlak `Δt_acc ∈ {240,120,60,30,10}` s
  yalnızca **düzgün** ızgarada ve önceki makaleyle karşılaştırılabilirlik için
  taranır. 8–14 yörüngelik perilun-yayılı panelde. Üç şey ölçülür:
  1. amaç değerinin yakınsaması,
  2. **çizelgenin kendisinin** yakınsaması (`N*_{240}` vs `N*_{30}` arasındaki
     derece-profili uzaklığı),
  3. verdict kararlılığı.
  Ayrıca `Δt_dec ∈ {60, 120, 300}` s ayrı taranır.
- **Öngörü (yazılıp sonra test edilecek):** genlik 120 s'de yakınsak (önceki
  makale ölçtü), **işaretli terim değil**. `πr/N` dekorelasyonu bunu
  gerektiriyor.
- **Kapı:** çizelge `Δt_acc` ile anlamlı biçimde kayıyorsa, tüm kampanya
  yakınsak ızgarada koşulur ve 120 s yalnızca karar ızgarası olarak kalır.
- **Bu WP geç aşamada değil, M1'dedir.** Sonradan bulunursa her şey
  tekrarlanır. (Aşama numaraları D42 ile kaydı; bu satır numaraya bağlı
  yazılmıyor.)

### WP5 — Bant probu ★
- **Soru:** atlanan kuvvetin **yönü** ucuza tahmin edilebilir mi? (RQ3'ün kilidi)
- **Yöntem:** `v̂ = −γ(N,h) Σ_{n=N+1}^{N+k} a_n(t)`, `k ∈ {1,2,3,4}`.
  **`γ` alanın gerçek derece-varyans spektrumundan `P_n` gelir (D53), uydurulmuş
  güç yasasından değil** — ay spektrumu tek bir güç yasası değil (`p_spec=2.13`
  vs `p_fit=1.76`). `P_n` ~1800 sayılık 1D tablo; uçurulabilirliği zedelemez.
  Ablasyon `abl-kaula` artık iki varyantı karşılaştırıyor: spektrum vs güç yasası;
  `κ = cos∠(v̂, Δa_ref)` irtifa × `N` × konum üzerinde
- **Ayrıca:** prob maliyetinin ölçümü (teorik vs gerçek çekirdek zamanı), ve
  prob kadansı `Δt_probe` için gereken tazeleme sıklığı
- **Prob maliyeti — asıl kalem (D50/D51). `2k/N` yanlıştı.** Paylaşılan yığın
  ve pencere yalnızca *bir noktadaki adaylar arası* maliyeti düşürüyor;
  *noktalar arası* maliyeti düşürmüyor. İleri prob yeni bir konumda alındığı
  için `P_{n,m}` rekürsiyonu sıfırdan `N+k`'ya koşmak zorunda → **her prob
  noktası ≈ bir tam sentez.** Maliyet tabanı
  `n_probe/n_RHS = 1/(τ_corr · r_RHS)`, **`Δt_dec`'ten bağımsız**. Ama bütçeye
  giren sayı **yay integralidir** (D67): `[∫dt/τ_corr]/N_RHS`. Yerel `τ_corr`
  ile küresel `r_RHS` çarpılmaz.
  **Bu integral, ızgara hücre sayısıyla aynıdır (D141):** `n_probe = M/n_s`.
  Dolayısıyla dereceyle **doğrusal** ve tek bir yüzde değil — vekil eksantrik
  yayda `N_RHS≈121k` ile `N=120`→**%8.5**, `N=300`→**%21**, `N=600`→**%42**.
  Eski "%12–19" hangi dereceye ait olduğu yazılmadan taşınıyordu. **WP5/WP16
  bunu derece bandı başına ölçer ve tek sayı olarak raporlamaz.**
- **Bu yüzden WP5'in birincil çıktısı bir eğri, bir yüzde değil:**
  `κ_eff(n_probe)` — aralık boyunca ortalanmış etkin yön isabeti, prob nokta
  sayısının fonksiyonu — ve karşısında ölçülen maliyet. İki uç ilan edilmiştir:
  `C-plan` (tam kapsayan ileri prob, ~%21+) ve `C-lite` (aralığın ilk RHS
  çağrısı yükseltilmiş derecede, eş-konumlu, ~%2, yön yalnız aralık başında).
  `2k/N` yalnızca `C-lite`'a uygulanır — orada prob gerçekten eş-konumlu.
- **Ölçülecek:** (a) paylaşılan yığından kısmi toplamla adayların çıkması
  gerçekten ek sentez gerektirmiyor mu; (b) `δ ∈ {2,4,8}` indis
  yarı-genişliği ve indüklenen `Δ_span`; (c) `κ_eff` vs `n_probe` eğrisi;
  (d) ölçülen çekirdek zamanının teorik modele uyumu.
- **İleri prob geçerliliği (D26) ve konum öngörücüsü (D31).** Prob, sınırdan
  önce değil, önündeki aralık boyunca alınır; konumları **iki-cisim
  propagasyonu** verir (alan değerlendirmesi yok). Üç kontrol:
  (a) iki aday derecenin bir aralıkta ayrışması `½‖Δa‖Δt_dec²` gerçekten
  `πr/N`'in kaç mertebe altında — döngüselliği kıran sayı budur;
  (b) Kepler öngörücüsünün hatası `½‖a_pert‖Δt_dec²` toleransın içinde mi;
  (c) `abl-predictor` — `N_plan`'da düşük-dereceli mikro-propagasyon `κ`'yı ve
  nihai hatayı değiştiriyor mu. Öngörücü maliyeti `B+`'a yazılır.
- **Ablasyon:** `abl-lookahead` — ileri prob yerine yalnız sınırda tek prob.
  İleri probun gerçekten gerekli olduğunu gösterir (veya göstermez).
- **`κ` yanında korunan kesir de raporlanır (D74).** `κ` ivme uzayında bir yön
  isabetidir; sadeleşme teriminin korunan kesri ise `⟨v̂,z⟩/⟨Δa,z⟩` ile
  `z = (ΦB)ᵀc` — `ΦB` açıları korumadığı için ikisi aynı sayı **değil**.
  `κ` önde gelen gösterge (propagasyonsuz, çizelgeden bağımsız); korunan kesir
  `A-sens` çizelgesinden türetilen `c` ile aynı koşuda hesaplanır. Yeni kapı
  değil, aynı ölçümün ikinci çıktısı.
- **Tanı, kapı değil (D106):** `κ` medyanı **< 0.7** bir **uyarıdır** —
  sadeleşmenin uçurulabilir biçimde taşınmasının zor olacağına işaret eder ve
  yedek seçeneklerin (PLAN §3.1) denenmesini tetikler. **Tek başına hiçbir
  sonuç bandı belirlemez ve WP5'i durdurmaz**; akış her hâlükârda G5b'ye devam
  eder. Sonuç bandını `f̂` verir: `f̂ < 0.15` ise uçurulabilir yakalama
  başarısızdır, L3 yalnızca L1 seviyesine iner ve proje 🟡 bandına kayar.
  Düşük `κ` ancak G5b de başarısızlığı teyit ederse o bandın *gerekçesi* olur.
- **Negatif kontrol:** yalnız-genlik + rastgele işaret null modeli. Kazanç
  bunun üstünde değilse işaret bilgisi taşınmıyor demektir.

### WP6 — Tutarlılık ufku ★
- **Soru:** pilot yay, sadeleşme planı için ne kadar süre geçerli?
- **Yöntem:** N=40 gevşek-toleranslı pilot yay ile referans yay arasındaki
  konum farkını zamanın fonksiyonu olarak ölç; `πr/N_plan(t)` ile karşılaştır →
  `T_coh`. Ayrıca `c_i`'nin pilot ve referans yay arasındaki kararlılığını
  ölç (bunun `v̂`'den çok daha dayanıklı olması beklenir).
- **Kararlılık FAZ-İNDEKSLİ ölçülür, zaman-indeksli değil (D52).** N=40 pilot
  yay 7 günde ~100 km along-track sapıyor = **~60–100 s zaman kayması**;
  perilun geçişi birkaç dakika. Zaman-indeksli karşılaştırma, gerçekte
  dayanıklı olan bir niceliği dayanıksız gösterir. İkisi de raporlanır ve
  fark, faz indekslemesinin gerekçesi olarak kullanılır.
- **Neden kritik:** `T_coh` yaydan çok kısaysa "tabloyu yükle ve uygula"
  mimarisi geçersizdir; WP11 zorunlu hâle gelir.
- **İkinci ölçüm (D26):** kısa ufuk. Geminin kendi durum tahmininin bir karar
  aralığı (~120 s) boyunca hatası, `πr/N` ile karşılaştırılır. İleri probun
  geçerliliği buna dayanıyor; "besbelli küçük" diye geçilmez, ölçülür.
- **Çıktı:** `T_coh` dağılımı + `c_i` kararlılık metriği + kısa-ufuk tahmin
  hatası / `πr/N` oranı

### WP7 — Sertifika
- **Soru:** `A-sign` gerçekten tavan mı, yoksa iyi bir çözüm mü?
- **Yöntem:** dışbükey zarf gevşetmesi + Frank–Wolfe. LMO, `A-sens/L1` ile
  **aynı yapıdadır** (karar aralığı başına ayrışır, tek çarpan) ama aynı
  problem değildir: amacı `⟨∇_q J, u_q(N)⟩`, ve `∇J = 2Qu` mevcut gevşek
  çözüme bağlıdır. `O(M·|𝒩|)`/iterasyon.
- **Kritik uygulama koşulu:** LMO **tam** çözülmelidir. `λ` bisektini en yakın
  ayrık bütçede kesip bırakmak sertifikayı geçersiz kılar. Gevşetilmiş
  problemde kesirli karışım serbest olduğundan bütçe tam karşılanabilir:
  kritik çarpanda iki komşu derece arasında kesirli ağırlık kurulur. Her
  iterasyonda LMO'nun bütçeyi makine hassasiyetinde karşıladığı doğrulanır ve
  kayda yazılır.
- **Metrik (D29):** `g_J = (J_desc − L_FW)/J_desc`,
  `g_E = 1 − √(L_FW/J_desc)`. **Eşik hata uzayında: medyan `g_E < 0.10`.**
  İkisi de raporlanır.
- **Dejenere hâller:** `L_FW = max{0, en iyi sınır}` (FW sınırı erken
  iterasyonda negatif olabilir, `J ≥ 0`). İterasyon bütçesi bitince
  `L_FW = 0` ise sertifika **boş**: boşluk raporlanmaz, ayrı sütunda sayılır.
- **Dil kuralı:** eşik sağlanırsa "oracle" denebilir; değilse NOTATION §6'daki
  uzun ad + "erişilebilir tahsis" dili

### WP8 — Düşük-dereceli pilot STM
- **Soru:** `Φ`, referans alan olmadan yeterince iyi elde edilebilir mi?
- **Yöntem:** pilot yay üzerinde `Φ̇ = A₄₀(t)Φ` entegrasyonu; referans `Φ` ile
  karşılaştırma (özdeğer/operatör normu ve `w(t)` üzerindeki etki)
- **Not:** Kepler/HCW **ana yöntem değil**, yalnızca `abl-stm` ablasyonu.
  Ay mascon dinamiğinde düşük-dereceli gerçek STM daha güvenilir ve zaten
  pilot yay koşuluyor.

### WP9 — IFBDA `C-plan` ★
- **Soru:** RQ3 — kazancın ne kadarı uçurulabilir biçimde yakalanır?
- **Yöntem:** PLAN §3.2'nin ileri–geri iterasyonu; `J ∈ {1,2,3}`
- **Mimari (WP6'nın dayattığı):** `c`, `w`, `λ`, `N_plan` uçuş öncesi pilot
  yaydan; `v̂` uçuş anında ileri bant probuyla
- **Plan FAZ ile indekslenir (D52):** `(devir indeksi, devir içi faz)`,
  mutlak zamanla değil. Pilot yayın ~100 s zaman kayması perilun geçişi
  mertebesinde; zaman-indeksli plan `c` ve `K_i`'yi tam da baskın oldukları
  yerde yanlış hizalar. Gemi kendi fazını biliyor. Kalan ufuk devir
  numarasıyla korunuyor, yani yarıçapın tek değerli fonksiyonuna dönmüyoruz.
- **Pilot yay üretimi ayrı bir kalemdir.** WP6 yalnızca ~10 tanesini
  tutarlılık ölçümü için koşuyor. `C-plan`, propagate edilen **her** yörünge
  için bir pilot yay ister: M3'te 64, M5'te ~370. Her biri gerçek yayın
  medyan %8'i → M3'te ~5, M5'te ~30 yay eşdeğeri. Küçük ama sıfır değil;
  maliyet tablosunda ve `B+`'da ayrı satır.
- **Çıktı:** **faz-indeksli** plan (devir, faz) → `c`, `K_i`, `λ`, `N_plan`
  + `B+` muhasebesi

### WP10 — `C-fb`
- **Soru:** bütçe sızıntısı tasarımla kapanır mı?
- **Yöntem:** harcanan `Σ N²`'yi RHS çağrılarında say, `λ`'yı çevrimiçi ayarla
- **Ölçüt:** `B2 / B1` oranı 1'e ne kadar yakın (önceki makalede 1.29, β=0.5'te 2.20)

### WP11 — Alıcı-ufuklu asgari varyant *(koşullu)*
- **Tetikleyici:** WP6'da `T_coh` ≪ yay uzunluğu
- **Kapsam:** yalnızca periyodik yeniden planlama; tam MPC değil
- **Aksi hâlde:** future work'te kalır

### WP12 — Çizelge biçimi kontrolleri
- **faz-indeksli** plan (D52) vs irtifa-binli (önceki makalenin konvansiyonu)
- zaman-indeksli varyant ayrıca `abl-timeindex` olarak koşulur
- `Δt_dec` = 60/120/300 s
- **Soru:** kazanç çizelgenin biçimine mi ait, yoksa binli hâlde de duruyor mu

### WP13 — Varyasyonel eleme
- **Soru:** hangi iki uçurulabilir aday propagasyona değer?
- **Gerekçe:** önceki makalede varyasyonel enstrüman çözünen 100
  karşılaştırmanın 100'ünde propagate işareti doğru verdi → kalibre bir eleyici
- **Yöntem:** 6 aday × 26 yörüngelik perilun-tabakalı panel = 156 çözüm;
  **14** ablasyon × 8 yörüngelik küçük panel = 112 çözüm → **~270 propagasyon
  eşdeğeri**. Ablasyonların küçük panelde koşulması bilinçlidir; panel
  büyüklüğü her tabloda yazılır. (Ablasyon sayısı D71 ile 10'dan 14'e çıktı —
  `abl-predictor`, `abl-timeindex`, `abl-window`, `abl-lookahead` eklendi;
  maliyet buna göre güncellendi — D79.)
- **Uygulama: tek çözüm, çok kanal (D59).** `rev42` zaten böyle çalışıyor —
  bütün politikalar aynı referans yayı ve aynı `G`'yi paylaşıyor, yalnızca
  `Δa` farklı. Dolayısıyla 270 **ayrı** çözüm olarak uygulanmamalı: yörünge
  başına tek artırılmış entegrasyon, içinde tüm aday kanalları. Adım kontrolü
  paylaşılıyor, marjinal maliyet yalnızca kanal başına `Δa` değerlendirmesi.
- **DÖNGÜSELLİK UYARISI (D54).** `A-sign` `J`'yi minimize ediyor; varyasyonel
  eleyici de `J`'yi hesaplıyor. Önceki makalenin 100/100 kalibrasyonu `J`'ye
  göre **optimize edilmemiş** politikalar içindi. `J`'ye göre optimize edilmiş
  bir adayı `J` ile skorlamak, `J`'nin `E`'yi öngördüğünü test etmez —
  varsayar. **Bu yüzden G6 artık taşıyıcı bir kapı**, teyit değil, ve
  başarısızlığı D41'in kaçışını tetikler.
- **Maliyet uyarısı:** bu aşama "ucuz" değil. Her varyasyonel çözüm bir
  propagasyona mal oluyor (önceki makalenin ölçümü). Erken planlarda "~26
  çözüm" yazılmıştı; doğrusu ~270.
- **Kapı G5:** eleme kuralı ve geçen aday sayısı (tam iki) **koşudan önce**
  yazılır.
- **Kapı G5b:** en iyi adayın **öngörülen** yakalama oranı `f̂ ≥ 0.15`. Altında
  M3'e ~670 yay harcanmaz; proje 🟡 bandına geçer. Eşik, H3'ün hedefi olan
  0.33'ün bilerek altındadır — umut vaat eden adayı erken elememek için.

### WP14 — İlk propagasyon
- Tasarım A, β = 1, 64 yörünge, **beş politika propagate edilir:**
  Aday 1, Aday 2, **`A-sign`**, `F-op`, **`R-int`**
- **Eşleştirme (D91): komparatör değil, adaylar çapaya kalibre edilir.**
  `B_tot` = `F-op(β)`'nın gerçekleşen işi (β=1'de arşivden). Gerekçe: `f`'in
  payı ve paydası aynı `E_F-op`'u kullanmak zorunda; komparatör her adaya göre
  kayarsa iki farklı `F-op` örneği çıkar ve `f` tanımsızlaşır.
- **Yay sayısı — kalibrasyon iterasyonu dahil.** Gerçekleşen iş ancak
  propagasyondan sonra bilinir (nominal kalibrasyon medyan %29 ıskalıyor),
  dolayısıyla her aday **2–3 propagasyonda** %2 bandına yakınsar.
  4 kalibre aday × 64 × ~2.5 ≈ **640**, artı pilot yaylar (~5 eşdeğer),
  artı `F-op` çapası (arşivden 0) ≈ **~645**; WP17a ile M3 toplamı ~670.
  (Erken planlarda "~200", sonra "~430" yazılmıştı; ikisi de kalibrasyon
  iterasyonunu saymıyordu.)
- `A-sign` zorunlu: H1b onun propagate hatasını test ediyor, `f`'in paydası o
- `R-int` zorunlu: 🟢 bandının şartlarından biri onu yenmek
- `F-env` propagate edilmez; sabit aile taramasından alt zarf olarak kurulur
- **Arşiv yeniden kullanımı:** eski `F-op`/`R-int` yayları yalnızca gerçekleşen
  işleri yeni `B2` sözleşmesi altında eşleştirilebiliyorsa kullanılır; aksi
  hâlde yeniden propagate edilir. Her yeniden kullanılan yay işaretlenir
- Birincil bütçe **`B2` (gerçekleşen)**; `B1` ikincil raporlanır
- **Kapı:** varyasyonel tahmin ile propagate işaret uyuşmuyorsa mekanizma
  değil uygulama hatası aranır → WP13'e dön

### WP15 — Popülasyon genişletme *(koşullu — G7)*
- **Tetikleyici:** M4/WP20a'nın G7 kapısı geçilmeli. Geçilmezse bu WP
  **koşulmaz** ve makale tek popülasyonlu, rejim-sınırlı yazılır.
- **Kalibre aday: 3 (D108).** `C-plan`, `R-int`, artı `F-op` çapası.
  `A-sign` popülasyon düzeyinde propagate **edilmez**.
- **Yay sayısı (D98 formülü):** 3 × ~368 × ~2.5 = 2760, artı çapa (`F-op`
  tasarım A dışında arşivde yok) ~304 → **~3060**. Pilot yaylar
  (~368 × %8) ≈ 30 yay eşdeğeri, ayrıca.
- **Bedeli açıkça raporlanır:** `f` (yakalama oranı) bu popülasyonlarda
  hesaplanamaz — paydası `E_A-sign`. H3 yalnızca tasarım A/B/C'de sınanır;
  diğer popülasyonlar verdikt (H4, H4b) testi.

| Popülasyon | Neden zorunlu |
|---|---|
| Tasarım A, B | tekrar — yön tek karıştırmaya ait değil |
| Tasarım C | bağımsız üçüncü kapsam tasarımı |
| Geniş-eliptik (Kaguya benzeri) | önceki makalede radyal kuralın **kazandığı** rejim |
| Beş geometri stratumu | geometri bağımlılığı |
| Düşük perilun 31–50 km | doğrusallaştırmanın en zayıf olduğu yer |

### WP16 — Maliyet merdiveni
- `B1 / B2 / B3` dört basamak (önceki makalenin `tab:cost-ladder` biçimi)
- **`B+` kalemleri, güncel model (D50/D58):**
  | Kalem | Tahmin | Kaynak |
  |---|---|---|
  | Pilot yay | %8 | önceki makalede ölçülmüş |
  | **Çevrimiçi ileri prob** | **~%12–19** (yay integrali, D67) | **en büyük kalem**; `[∫dt/τ_corr]/N_RHS`, WP5 ölçer |
  | Anahtarlama (5040 sınır) | ~%2 | önceki kampanyadan kalibre: ~0.45 çağrı/anahtarlama |
  | Konum öngörücüsü (Kepler) | ~0 | analitik, alan değerlendirmesi yok |
  | IFBDA planlaması | çevrimdışı, küçük | `J` süpürmeleri |
  `C-lite` varyantında prob kalemi ~%2'ye iner, karşılığında yön yalnızca
  aralık başında bilinir. Maliyet–isabet eğrisi WP5'ten gelir.
- **`A-sign` tavan teması denetimi (D57):** optimizatör, bütçenin izin verdiği
  yerde dereceyi referansa itip defekti tanım gereği sıfırlamaya *teşviklidir*;
  sabit politikanın böyle bir teşviki yok. Her `A-sign` çizelgesi için tavan
  teması kesri raporlanır, yüksek temaslı yörüngeler orandan çıkarılır.
- **Eşleştirme (D91):** `B2 + B+ = B_tot`; çapa `F-op(β)`'nın gerçekleşen işi,
  adaylar ona %2 bandında kalibre edilir. Sabit komparatöre
  hayali ek yük eklenmez
- 14 yörüngelik ölçülen seri çekirdek zamanı paneli, 3 tekrar, 0.90–1.10 bandı

### WP17a — `A-sign+` indirgenmiş (M3'te) ★
- **Soru:** referans-yay doğrusallaştırması bağlayıcı bir sınırlama mı?
- **Neden M3'te (D40):** bu, makalenin **merkezî nesnesinin ne olduğunu**
  belirliyor — `A-sign` bir *yörünge seviyesi* kıyas mı, yoksa bir
  *referans-yay* optimumu mu? Popülasyon genişletmeden önce bilinmeli;
  sonradan öğrenilirse M5'in 1800 yayı yanlış adlandırılmış bir nesne
  üzerinde harcanmış olur.
- **Yöntem:** perilun-tabakalı 16 yörünge; çizelgeyi propagate et → yeni yayda
  yeniden tablola → yeniden çöz → tekrar propagate (2 iterasyon) ≈ 32 yay
- **Sonuç yorumu:** çizelge veya hata anlamlı biçimde kayıyorsa `A-sign` bir
  *referans-yay optimumu*dur, makale böyle adlandırır, ve `f`'in paydası
  `A-sign+` ile de raporlanır

### WP17b — `A-sign+` tam panel (M7'de)
- WP17a'nın 32 yörüngeye genişletilmiş hâli; yalnızca WP17a bir kayma
  gösterirse zorunlu, aksi hâlde doğrulama

### WP20a — Bütçe probu (M4) ★
- **Soru:** kazanç β=1'e mi özgü?
- **Yöntem:** tasarım A'nın perilun-tabakalı 16 yörüngelik alt kümesi
  (panel seçim kuralı aynı), β ∈ {0.50, 0.75, 1.25, 1.50}, beş politika
  → **~576 yay** (4 kalibre aday × 16 × 4 × ~2 = 512, artı çapa 16 × 4 = 64;
  D98 formülü)
- **Kapı G7:** `C-plan` β = 0.75 ve β = 1.25'te de çözünen çoğunluğu alıyor mu.
  Almıyorsa **M5 iptal**; makale tek popülasyonlu ve rejim-sınırlı yazılır
- **Neden popülasyon genişletmeden önce:** önceki makalede iç üyenin bütün
  yapıcı sonucu β=0.5'te ters dönmüştü. Tek bütçede kazanan bir yöntem için
  üç popülasyona 1800 yay harcamak aynı hatanın tekrarı olur

### WP20b — Tam bütçe ızgarası (M6)
- Tasarım A'nın 64 yörüngesi, β ∈ {0.50, 0.75, 1.25, 1.50}; tasarım B
  yalnızca çaprazlamaya en yakın tek değerde (önceki makalenin sahne kuralı)
- **β = 3 propagate edilmez** — enstrümanı tüketiyor, yalnızca `Ê` düzeyinde

### WP18 — Uzun yay ve dayanıklılık
- 60 günlük panel (8 yörünge, mevcut referanslar)
- Faz-kaydırma MC — sadeleşme sömürüsüne aşırı uyum testi
- Referans derecesi kontrolü 300 → 600 — referans-komşuluğu artefaktı
- Bin kontrolü

### WP21 — Mimari seçim testleri (T1–T7) ★

**Neden var (D116).** `A-sign` ile `C-plan` aynı statüde değil. `A-sign`
bir sezgiselden değil, doğrudan yörünge-hatası fonksiyonelinden çıkıyor —
ona güveniyoruz. `C-plan` ise **ilk yorumlanabilir, denetlenebilir,
düşük-karmaşıklıklı deneme**; en iyi uçurulabilir çözüm olduğu iddia
edilmiyor. Bu WP, M1/M2'yi "`C-plan`'ın ön koşulu" olmaktan çıkarıp
**hangi denetleyici ailesinin mantıklı olduğunu seçen keşif aşaması** yapar.

**Ortak özellik: yedi testin hiçbiri propagasyon istemiyor.** Hepsi WP1'in
`Δa`/`Φ` tabloları, WP3'ün `A_i` sonek dizisi ve WP7'nin FW yinelemesi
üzerinde koşuyor. Maliyet duvar saati; yay bütçesi **0**.

---

**T1 — Boşluk ayrıştırması (küçük örneklerde küresel optimum).**
`g_E` tek başına iki farklı şeyi karıştırıyor: inişin optimumdan uzaklığı
ve gevşetmenin gevşekliği. Şu an bunları ayıramıyoruz, dolayısıyla büyük bir
`g_E` "çözücü zayıf" mı "sertifika muhafazakâr" mı bilmiyoruz.
Kısa segmentlerde (`K_dec = 8`, `|𝒩| = 6` → 1.7 M kombinasyon) **kaba kuvvetle**
gerçek `J*` bulunur; `K_dec ≈ 30`, `|𝒩| ≈ 10` için dışbükey MIQP
(çözücü yoksa kaba kuvvet tabanı yeter, bağımlılık **isteğe bağlı**).
- **Çıktı:** `J_desc − J*` (iniş terimi) ve `J* − L_FW` (gevşetme terimi).
- **Tescilli okuma:** iniş terimi baskınsa daha güçlü çözücü değer; gevşetme
  terimi baskınsa çözücüye yatırım yapılmaz ve `g_E` muhafazakârlık olarak
  raporlanır.

**T2 — `S-round`: gevşetme + yuvarlama, ikinci çözücü.**
FW zaten gevşek optimum `θ*`'ı üretiyor; şu an yalnızca alt sınır için
kullanılıyor. Aynı yinelemeden üç adımla bir **çizelge** de çıkar:
aralık başına `argmax θ`, ardından karma-tamsayılı optimal kontrolden gelen
**toplaya-yuvarla (sum-up rounding)**, ardından bir tur koordinat cilası.
- **Maliyet:** WP7'nin yan ürünü, ek koşu yok.
- **Tescilli okuma:** `S-round` çok-başlangıçlı inişi geçerse **kanonik
  çözücü değişir** ve iniş bir kontrol olarak raporlanır. Geçmezse iniş
  lehine bir kanıt satırı olur.

**T3 — Eşleşmenin etkin rankı ve altuzay tutarlılığı → olası `C-rank` ★.**

**Ölçüm boyutsuzlaştırılmış formda yapılır.** `A_i` durum üzerinde bir
kuadratik form ve `(r,v)` blokları farklı birim taşıyor; ham izi, özdeğerleri
ve özyönleri birim sistemine bağlı — hızı m/s yerine km/s yazınca verdikt
değişebilir. Doğru dönüşüm **eşlenik**tir (benzerlik değil):
`Ã_i = Sᵀ A_i S` ile `S = diag(L,L,L,L/T,L/T,L/T)`
(`tda.stm.nondimensionalise_form`). Aynı sorunu STM koşul sayısında zaten
yakalamıştık; burası da aynı sınıf.

**İki alt test, `C-rank` dalı ikisi birden geçerse açılır:**
- **T3a** — `Ã_i` izinin %95'i için gereken özyön sayısı `p`. Gerekli koşul.
- **T3b** — baskın altuzayın **zaman tutarlılığı**: ardışık epoklar arasında
  asal açı `∠(U_p(t_i), U_p(t_{i+1}))`, bir korelasyon zamanı boyunca
  birikmiş. `p ≤ 2` olması **tek bir** 1-B altuzay olduğu anlamına gelmez;
  baskın özyön yay boyunca dönebilir ve o zaman 6-B yapı gerçekten
  küçülmemiş olur.

`M_j`'nin rankı tam olarak 3 (`H_r` üç bileşen seçiyor), ama `A_i`'nin
**özdeğer dağılımı** ölçülmedi. STM'in baskın yönü along-track seküler
büyüme olduğundan `A_i`'nin izinin çoğunun 1–2 özyönde toplanması **kuvvetle
muhtemel**. Ölçüm bedava: `A_i` zaten var.
- **Yöntem:** yay boyunca `A_i`'nin özdeğerleri; izin %95'i için gereken `p`;
  `A_i`'yi ilk `p` özçiftine kesip (a) `J` hatası, (b) yeniden çözülen
  çizelgenin profil uzaklığı.
- **Neden önemli:** `p ≤ 2` çıkarsa sadeleşme durumu **`p` skalerdir**:
  `s⁽ᵐ⁾ = e_mᵀS`, eşleşme terimi `Σ_m σ_m (e_mᵀu_i) s⁽ᵐ⁾`. O zaman denetleyici
  6×6 makine ve dondurulmuş `c` tablosu taşımaz; `p` skaleri **çevrimiçi,
  kendi uçtuğu yörüngede** günceller. Bu hem IFBDA'dan ucuz hem pilot
  kaymasına dayanıklı — yani danışmanın istediği "çevrimiçi eşdurum" fikrinin
  ucuz hâli. `p ≥ 4` çıkarsa tam `Q` gerekli ve `C-plan` doğru tasarımdır.
- **Aynı ölçüm DP/Viterbi sorusunu da kapatıyor:** DP'nin durumu `(S, bütçe)`
  olduğu için 6 boyutta patlıyor; `p ≤ 2` ise `p+1` boyutta yaşayabilir.
  Gelecek-iş satırı artık argüman değil, sayı olur.

**T4 — Ufuk yeterliliği eğrisi ★.**
Şu an eşleşmenin **ne kadar ileriye uzandığını bilmiyoruz** — bu hem MPC
sorusunun cevabı hem de makalenin kendi başına eksik olan bir fizik sonucu.
`A_i` yerine ufukla kesilmiş `A_i^H = Σ_{i ≤ j ≤ i+H/Δt} ω_j M_j` konur,
tahsis her sınırda yalnız `[t_q, t_q+H]` üzerinden çözülür, `N_q` uygulanır,
ilerlenir — yani **tam bilgili alıcı-ufuk `A-sign`**.
- **Süpürme:** `H ∈ {⅛, ¼, ½, 1, 2, 4}` devir ve tam yay.
- **Çıktı:** `ρ̂(H)/ρ̂(T)` eğrisi ve diz noktası `H*`.
- **Tescilli okuma:** `H*` ≤ 1 devir ise kısa ufuklu MPC ucuz ve yeterli;
  `H*` yayın yarısını aşıyorsa çevrimdışı plan zorunlu ve MPC pahalı.

**T5 — `c` duyarlılığı (dondurulmuş plan yeter mi).**
WP6 pilot–referans farkını zaten ölçüyor. Bu fark `c`'ye bir pertürbasyon
olarak uygulanır ve çizelge yeniden çözülür.
- **Tescilli okuma:** profil uzaklığı G3 eşiğinin altındaysa dondurulmuş `c`
  meşrudur (`C-plan` tamam); üstündeyse çevrimiçi eşdurum güncellemesi
  (adjoint/MPC veya T3'ün `C-rank`'i) gerekli.
- **Maliyet:** WP6'nın yan ürünü.

**T6 — Politika ifade edilebilirliği `A-fit` ★ (makalenin kendi tezinin
yanlışlama testi).**
`A-sign` çizelgelerine bir durum-geri-besleme politikası uydurulur ve
kazancın ne kadarını koruduğu ölçülür. **Tasarım A'da uydurulur, tasarım
B'de değerlendirilir** — kıyastan bilgi sızmasını kapatan tek şey budur.
- **T6a — yalnız düzgün durum:** `(h, r, |v|, faz, kalan süre, harcanmış bütçe)`.
  Küre üzerinde konum **yok**.
- **T6b — + konum:** `(φ, λ)` artan çözünürlükte (T7'nin vekil-modeli).
- **Tescilli okuma:** T6a kazancın **≥ 0.60**'ını koruyorsa kazanç
  doku değil **geometri** demektir — bant probu gereksizdir ve ucuz bir
  durum-geri-besleme politikası doğru mimaridir. T6a **%25'in altında**
  kalıp T6b yükseliyorsa kazanç gerçekten dokudadır; makalenin bilgi iddiası
  **ölçülmüş** olur, savunulmuş değil.
- **Bu bir "NN eğittik" iddiası değildir**, bir kontroldür; en basit
  ifade eden (ağaç/spline) seçilir ve karmaşıklığı raporlanır.

**T7 — Vekil-model maliyet–koruma Pareto'su.**
Mevcut `abl-probe`'un genişletilmiş hâli: tek bir vekil değil, artan parametre
sayısında bir aile. `πr/N` argümanını **iddia olmaktan çıkarıp ölçülmüş bir
eğri** yapar. Prob ek yükü %12–19 çıkarsa ilk bakılacak alternatif budur.

---

**Nerede duruyor.** T1/T2 WP7'nin, T5 WP6'nın yan ürünü; T7 `abl-probe`'un
genişlemesi. **Net yeni iş T3, T4, T6.** Üçü de M1/M2 tablolarında koşuyor,
üçü de 0 yay. Kapı değiller: çıktıları **mimari dallanma tablosuna**
(ROADMAP) girer.

### WP19 — Manuskript
- **Gönderim kapıları — ikisi de çıkış kodu 0 vermeli:**
  `python check_stale.py` (eskimiş **terim**, D94) ve
  `python check_numbers.py` (eskimiş **sayı/eşik**, D101).
  `DECISIONS.md` ikisinde de bilerek hariç.
- **`check_stale.py` (D94).** Çıkış kodu 0 olmalı;
  1 ise eskimiş terim var ve gönderim durur. `DECISIONS.md` bilerek hariç.
- **Gönderim öncesi tutarlılık kontrolleri (D73, D79):**
  ablasyon sayısı × panel = M2 maliyeti satırı (üç dosyada birden);
  aday listesi PREREG = manuskript; kapı listesi ROADMAP = README = manuskript.
  Kapsam her büyüdüğünde bu üç eşitlik yeniden doğrulanır.
- **Gönderim öncesi sıfırlanacak sayaçlar:** `\ph{}` sayısı 0,
  `\dnote{}` sayısı 0, **metinden referanslanmayan float sayısı 0** (şu an
  21'in 13'ü referanssız — referans cümleleri `\ph{one paragraph}`
  bölgelerinin içine gelecek), `UNVERIFIED` bib kaydı 0.
- `konsey` → `literatur` → `gonderim` sırasıyla
- Figürler **yalnızca** `make_figures_oa*.py` ile; ortak betik çağrılmaz
- Commit mesajlarına asistan imzası eklenmez
