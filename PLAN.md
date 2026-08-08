# PLAN — Yörünge-farkında derece tahsisi

> Semboller, politika kodları ve maliyet tanımları [NOTATION.md](NOTATION.md)'de
> dondurulmuştur. Bu dosya bilimsel içeriği anlatır; koşu sırası
> [ROADMAP.md](ROADMAP.md), iş paketleri [WP.md](WP.md).

## 0. Üç araştırma sorusu

İlk makale bu üçüne odaklanır, fazlasına değil.

- **RQ1.** Yörünge-farkında bilgi altında, sabit/radyal tahsise göre ne kadar
  **erişilebilir kazanç** var?
- **RQ2.** Bu kazancın ne kadarı **duyarlılık ağırlığından**, ne kadarı
  **işaretli sadeleşmeden** geliyor?
- **RQ3.** Referans yay olmadan ve alanın **yay boyunca noktasal
  değerlendirmeleri** olmadan (yalnız `P_n` gibi küresel/çevrimdışı üst-veriyle)
  bunun anlamlı bir kısmı yakalanabilir mi?

§7'deki genişletme fikirleri (MPC, dağıtımsal tahsis, alan transferi, GRAIL
kovaryansı) **future work**'tür ve roadmap'in zorunlu hattında değildir.

---

## 1. Problem

### 1.1 Amaç ve neden ayrıştırılamaz

`codebase` makalesi amacın tam biçimini yazdı
([08_discussion.tex:189](../codebase/chapters/08_discussion.tex:189)):

```
E² = (1/T) ∫₀^T ‖ H_r ∫₀^t Φ(t,τ) B Δa(τ) dτ ‖² dt ,     H_r = [I₃ 0]
```

İç integral 6 boyutlu bir **durum** pertürbasyonudur; metrik bir **konum**
RMS'i olduğundan norm alınmadan önce konum bloğu seçilmelidir. `H_r` her
ifadede açık taşınır.

Norm, **işaretli integralden sonra** alınıyor. Sonuçları:

| Amaç biçimi | Ayrıştırılabilir | Yörünge sıralamasını verir |
|---|---|---|
| `Σ_i ‖Δa_i‖²` (kuvvet defekti) | ✔ | ✘ — önceki makalenin ana bulgusu |
| `Σ_i w_i ‖Δa_i‖²` (Φ-ağırlıklı büyüklük) | ✔ | ✘ — `D_I` basamağı hâlâ ters sıralıyor (45/64, 44/64) |
| `‖Σ_i v_i‖²` (işaretli) | **✘** | hedef |

Önceki makale ikinci satırın yetmeyeceğini peşinen söylemişti: *"replacing
d(t,N) by a Φ-weighted version would correct the weighting while still summing
magnitudes, bounding the gap from one side without closing it."* Yeni
makalenin teknik çekirdeği üçüncü satırdır.

### 1.2 Ayrıklaştırma: çözülebilir yapı

Biriktirme ızgarası `t_1..t_M` (`Δt_acc`, bkz. NOTATION §4). Konum seçicisi
`H_r = [I₃ 0]` açıkça tanımlanır — boyut tartışması bırakmamak için:

Izgara **düzgün varsayılmaz** (D55/D60); adım `Δt_i`, dış integral bir
kuadratür, ağırlıklar `ω_j ≥ 0` toplamı `T`:

```
u_i(N) = Φ(t_0, t_i) B Δa(t_i, N) Δt_i
S_j    = Σ_{i≤j} u_i
δr_j   = H_r Φ(t_j, t_0) S_j
M_j    = Φ(t_j,t_0)ᵀ H_rᵀ H_r Φ(t_j,t_0)
J(u)   = (1/T) Σ_j ω_j Sⱼᵀ M_j Sⱼ = uᵀ Q u
```

`Q_ik = (1/T) A_{max(i,k)}` ile `A_i = Σ_{j≥i} ω_j M_j`; bloklar **6×6**
(`u_i` bir durum pertürbasyonu, seçim `M_j`'nin içinde). **Ağırlıklar yalnızca
sonek dizisinin içine giriyor, yapı bozulmuyor** — makalenin matematik katkısı
uyarlanabilir ızgarada aynen ayakta. Düzgün ızgarada `ω_j = Δt`, `T = MΔt` ve
eski `1/M` formuna indirgeniyor. **`ω_j ≥ 0` şarttır** — `Q ⪰ 0` ve dolayısıyla
dışbükey gevşetmenin tamamı buna bağlı; dikdörtgen/yamuk sağlar, yüksek
dereceli Newton–Cotes sağlamaz (D76).

**İç kuadratür kuralı serbest değil — önek yapısı onu belirliyor (D124).**
`S_j = Σ_{i≤j} u_i`'nin **düz** önek toplamı olması, `i` epokunun iç
ağırlığının `j`'den **bağımsız** olmasını gerektirir:

- **Yamuk kullanılamaz.** Uç katsayısı `(t_j − t_{j−1})/2`, yani `j`'ye bağlı;
  `Q_ik = A_max(i,k)` yerine "o artı `j`-bağımlı düzeltme" kalırdı.
- **Dikdörtgen yapıyı korur ama birinci mertebedir.** `n_s = 2`'de (Nyquist
  sınırı) bu gerçek bir doğruluk kaybı, yuvarlama ayrıntısı değil.
- **Uygulanan: orta nokta + kenar şeması.** Biriktirme düğümleri hücre
  **orta noktalarında**, dış integral hücre **kenarlarında** örnekleniyor.
  İç ağırlık tam hücre genişliği → `j`'den bağımsız → yapı **aynen** korunur;
  kural orta nokta → **ikinci mertebe**; ve `S_j` tam olarak `j+1`. kenara
  kadarki iç integralin **orta-nokta ayrıklaştırmasının tam önek toplamıdır**.
  (Sürekli integrale hâlâ bir yaklaşımdır — ikinci mertebeden; "tam" olan,
  ayrıklaştırmanın önek toplamı olmasıdır.) Dış integralde önek toplamı
  olmadığı için orada yamuk serbest.

Yani yukarıdaki formüllerde `t_i` **hücre orta noktası**, `t_j` **hücre
kenarı**dır; ikisi çakışmaz. Uygulama: `tda.grids.AccumulationGrid`.

```
c_i := Σ_k Q_ik u_k = (1/T) [ A_i · Σ_{k≤i} u_k  +  Σ_{k>i} A_k u_k ]
```

Terminal amaç `J_T = ‖H_r Φ(t_M,t_0) Σ_i u_i‖²`, yani
`M_j = δ_{jM} Φ(t_M,t_0)ᵀ H_rᵀ H_r Φ(t_M,t_0)`. Taşıma matrisi düşürülmez:
`u_i`'ler `t_0`'a taşınmıştır.

### 1.2b Karar değişkeni karar ızgarasındadır ★

Amaç `M` biriktirme örneğinde toplanıyor, ama **derece `M` kez seçilemez.**
`g(i) = q`, `I_q = {i : g(i) = q}` ile

```
N_i = N_{g(i)} ,   W_q = Σ_{i∈I_q} Δt_i ,   min over (N_1,…,N_{K_dec})
kısıt:  Σ_q W_q N_q² ≤ B = B₁·T
```

**Ağırlık `W_q`, örnek sayısı `|I_q|` değil (D60).** Uyarlanabilir ızgarada
sayım kullanmak, perilunu yalnızca daha yoğun örneklendiği için bütçede fazla
ağırlıklandırır. Koordinat güncellemesinde ceza `λ W_q N²`; FW gevşetmesinde
de `Σ_q W_q Σ_N θ_qN N²`.

Bu, uygulanabilirlik kısıtı değil, **karşılaştırmanın geçerlilik koşuludur**:
`M` boyutlu uzayda çözülen bir `A-sign`, uçurulabilir denetleyiciden yalnızca
daha çok *bilgi* kullandığı için değil, daha çok *anahtarlama serbestliği*
olduğu için de kazanır — ve o zaman ölçtüğümüz şey artık yörünge-farkındalık
değil, ızgara inceliğidir.

Gruplamanın maliyeti yok. Karar değişkeni `q` için diğer bloklarla etkileşim

```
c_i^{(−q)} = Σ_{k∉I_q} Q_ik u_k          (önek/sonek toplamlarından blok-içi kısım çıkarılarak)
```

ve blok-içi terim, `I_q = {i_1 < … < i_m}` için

```
Σ_{i,k∈I_q} u_iᵀ A_max(i,k) u_k
  = Σ_p u_{i_p}ᵀ A_{i_p} u_{i_p} + 2 Σ_p ( Σ_{p'<p} u_{i_{p'}} )ᵀ A_{i_p} u_{i_p}
```

yani blok içinde tek bir yürüyen önek toplamıyla `O(m)`. Toplam süpürme
maliyeti yine **O(M·|𝒩|)**.

**Sonuç:** ayrıştırılamayan problem, arşivlenmiş `Φ` ve tablolanmış
`Δa(t,N)` üzerinde, **tek bir propagasyon harcamadan** optimize edilebilir.

### 1.3 İki ızgara

`Δa` yaklaşık `πr/N` ölçeğinde dekorele olur — **uydu yarıçapında**, ve bu
yörünge boyunca sabit değil: `τ_corr = πr/(Nv)`, tek paylaşılan inceltme
derecesinde perilun/apolun oranı **5.6 kat** (D129; `N=300`'de 9.5 s ve 53.6 s).
Eski "~50 kat" iki ucu farklı derecede okumaktan geliyordu. Karar değişkeni `N(t)` ise yörünge zaman ölçeğinde
değişir. Bunlar aynı ızgara olamaz:

- **Biriktirme ızgarası** `τ_corr`'a göre **uyarlanabilir** inceltilir; düzgün
  ve yeterince ince bir ızgara apolunda saf israf, yeterince kaba olanı
  perilunda aliasing. Yakınsama ölçütü `n_s = τ_corr/Δt_i` (D68).
- **Karar ızgarası** `Δt_dec` kaba kalabilir (uygulanabilirlik için kalmalı).

Önceki makale kuvvet defektinin **büyüklüğünün** 120 s'de yakınsak olduğunu
ölçtü (10 s'ye inince %1.2–1.8 değişim). **İşaretli terim için aynı şey
gösterilmiş değildir.** WP4 bunu ölçer ve `A-sign`'ın bir kısmının
ayrıklaştırma artefaktı olup olmadığına karar verir.

---

## 2. Yöntem merdiveni

| # | Kod | Ne ekler | Uçurulabilir | Sorusu |
|---|---|---|---|---|
| L0a | `F-op` | — | ✔ | taban |
| L0b | `R-rad` | yarıçap | ✔ | önceki makalenin kaybeden ucu |
| L0c | `R-int` | yarıçap | ✔ | mevcut en iyi uçurulabilir üye |
| L0d | `A-force` | referans alan | ✘ | kuvvet seviyesi tavan (O26) |
| **L1** | `A-sens` | yerel duyarlılık çekirdeği `K_i = K(t_i,t_i)` | ✘ | **RQ2 payda** |
| **L2** | `A-sign` | tüm `Q` (köşegen dışı = işaretli sadeleşme) | ✘ | **RQ1** |
| L2+ | `A-sign+` | kendi yörüngesinde sabit nokta | ✘ | doğrusallaştırma bağlayıcı mı |
| **L3** | `C-plan`, `C-fb`, `C-lite` | ileri prob + `P_n` spektrumu + pilot STM | ✔ | **RQ3** |
| L4 | ablasyonlar | — | — | hangi bileşen kazandırıyor |

> **Hata düzeyi (D37).** Bu bölümün bütün hataları **öngörülmüş** hatadır,
> `Ê = √J`; propagasyon yok. `E` (propagate edilmiş) yalnızca M3 ve sonrasında
> geçerlidir ve oranlar `ρ̂` / `ρ` ile ayrılır. Bkz. NOTATION §3.

### 2.1 L1 — `A-sens`

**Ağırlık ödünç alınmaz, türetilir (D75).** Ama **köşegen kısıtlaması ayrık
`Q`'da değil, çekirdekte alınır (D110).** Sürekli hâlde `J` bir çekirdeğe karşı
çift integraldir, `J = ∫∫ Δa(τ)ᵀ K(τ,σ) Δa(σ) dτ dσ`, ve `τ=σ` kümesinin ölçüsü
sıfırdır — "köşegen dışını at" harfiyen alınabilecek bir tanım değildir. Ayrık
`u_iᵀ Q_ii u_i` ifadesi `u_i ∝ Δt_i` yüzünden `Δt_i²` taşır: ızgara inceldikçe
sistematik olarak küçülür, uyarlanabilir ızgarada ise `Δt_i` yay boyunca
değiştiğinden yay içi sahte bir eğilim üretir. İyi tanımlı olan, **çekirdeğin
köşegen üzerindeki değeridir**:

```
K(t,t) = (1/T) ∫_t^T Bᵀ Φ(τ,t)ᵀ H_rᵀ H_r Φ(τ,t) B dτ
K_i    = (1/T) Bᵀ Φ(t₀,t_i)ᵀ A_i Φ(t₀,t_i) B          (cocycle özdeşliğiyle)

A-sens:  argmin_N [ Σ_{i∈I_q} Δt_i · Δa(t_i,N)ᵀ K_i Δa(t_i,N) + λ W_q N² ]
A-force: argmin_N [ Σ_{i∈I_q} ‖Δa(t_i,N)‖² Δt_i               + λ W_q N² ]
```

`A-sens` böylece `O(Δt)` mertebesinde ve ızgara-değişmez. Eşdeğer yazım:
`Σ_i u_iᵀ Q_ii u_i / Δt_i`. `K_i` aynı `A_i` sonek dizisinden okunduğu için
ek maliyet yok; `Q` türetimi, önek/sonek yapısı ve FW sertifikası etkilenmez.

`w(t) = ‖Φ(T,t)B‖²` gibi bir skaler ödünç almak keyfî bir ithal olurdu ve
ayrıca **terminal** bir proxy'yi **yay-RMS** amacına takmak demekti. Gerek yok.

Merdiven böylece tek bir yörünge çekirdeğinin ayrıştırması oluyor:
`A-force` çekirdeği yok sayar, `A-sens` yerel köşegenini `K_i` tutar, `A-sign`
epoklar arası tam eşleşmeyi tutar. **RQ2 tek bir soruya iniyor: kazancın ne
kadarı köşegen dışında yaşıyor?**

### 2.2 L2 — `A-sign` ★

```
min_{(N_1,…,N_Kdec) ∈ 𝒩^Kdec}  uᵀ Q u      s.t.  Σ_q W_q N_q² ≤ B = B₁·T
```

Karar uzayı `K_dec` boyutlu (§1.2b). `B1` çağrı başına ortalama, kısıt ise
işin zaman integrali: `B = B₁·T`, ve aralık ağırlığı `W_q` (D60).

**Üst sınır (erişilebilir çözüm): blok koordinat inişi.** L1'den başla, §1.2b'nin
`c_i^{(−q)}` ve blok-içi terimiyle karar aralığı karar aralığı güncelle, `λ`
dış döngüde bisekte. Çok başlangıçlı: `F-op`, `R-rad`, `R-int`, `A-sens` ve
rastgele tohumlar; yayılım raporlanır.

### 2.3 Sertifika — "oracle" kelimesini hak etmek

Sertifika yokken "trajectory-level oracle" demek hakeme kapı bırakır. Ama
sertifika **hesaplanabilir**:

`Q` PSD ve kısıt karar aralığı başına ayrık bir kümeye ait. **Dışbükey zarf
gevşetmesi** — karar aralığı `q` için `θ_{qN} ≥ 0`, `Σ_N θ_{qN} = 1`, katkı ve
maliyet aynı konveks kombinasyonla — bir dışbükey QP verir ve optimumu
**geçerli bir alt sınırdır** (her tamsayı çizelge, gösterge `θ` ile gevşetmede
olurlu ve aynı amacı verir).

Çözüm yolu: **Frank–Wolfe**. Lineer alt problem (LMO)

```
s ∈ argmin_θ ⟨ ∇J(u), u(θ) ⟩   s.t.  bütçe
```

`L1/A-sens` ile **aynı yapıdadır** — karar aralığı başına ayrışır ve tek bir
çarpanla çözülür — ama **aynı problem değildir**:
LMO'nunki `⟨∇_q J, u_q(N)⟩` (L1'inki ise yerel çekirdek biçimi
`Σ_{i∈I_q} Δt_i Δa_iᵀ K_i Δa_i`, D75/D110), ve
`∇J = 2Qu` mevcut gevşek çözüme bağlıdır.
Metinde "aynı yapıda", "L1'in kendisi" değil. Gradyan önek/sonek yapısıyla
`O(M)`; iterasyon başına `O(M·|𝒩|)`.

**LMO tam çözülmelidir.** `λ` bisektini en yakın ayrık bütçede kesip bırakmak
sertifikayı geçersiz kılar. Gevşetilmiş problemde kesirli karışım serbest
olduğundan bütçe **tam** karşılanabilir: kritik çarpanda iki komşu derece
arasında kesirli ağırlıkla eşitlik kurulur. Uygulamada bu açıkça çözülür ve
LMO'nun bütçeyi makine hassasiyetinde karşıladığı her iterasyonda doğrulanır.

Alt sınır FW iterasyonunun **kendisinden** değil, her adımda bedavaya çıkan
**FW dualite boşluğundan** okunur:

```
J(u) + ⟨∇J(u), s − u⟩  ≤  J*_gevşek  ≤  J*_tam  ≤  J(u_koordinat)
```

→ Koordinat inişi **üst sınır**, FW dualite boşluğu **alt sınır**, ikisi
arasındaki fark **sertifikalı boşluk**.

**Sertifikalı boşluk — tanım (NOTATION §6):**

```
g_J = (J_desc − L_FW)/J_desc ,   g_E = 1 − √(L_FW/J_desc) = 1 − √(1−g_J)
L_FW = max{0, iterasyonlar boyunca en iyi FW sınırı}
```

`J = E²` olduğundan ikisi aynı sayı değil. **Eşik hata uzayında:** medyan
`g_E < 0.10` (≈ `g_J 0.19`). `L_FW = 0` kalırsa sertifika boştur, o yörünge
ayrı sütunda sayılır.

**Dil kuralı:** medyan `g_E` 0.10'un altına inerse `A-sign` "oracle" olarak
adlandırılabilir. Aksi hâlde NOTATION §6'daki uzun ad kullanılır ve — önceki
makalenin O26'da yaptığı gibi — *erişilebilir tahsis, sınır değil* denir;
bildirilen boşluklar muhafazakârdır.

### 2.4 L2+ — doğrusallaştırma bağlayıcı mı

`A-sign`, arşivlenmiş referans yay üzerindeki `Φ` ve `Δa`'yı kullanıyor.
Politika değişince yörünge değişir → `Φ` değişir → yerel defekt değişir.
Sabit-nokta testi: `A-sign` çizelgesini propagate et, yeni yay üzerinde
`Φ` ve `Δa`'yı yeniden tablola, yeniden çöz, tekrar propagate (2 iterasyon).

Çizelge ve hata anlamlı biçimde kayıyorsa, `A-sign` bir **referans-yay
optimumu**dur ve makale bunu böyle adlandırır.

---

## 3. L3 — uçurulabilir denetleyici

`A-sign` kıyası, bir propagatörde bulunmayan üç bilgiyi kullanıyor. Üçünü de
elemek gerekiyor:

| `A-sign`'ın kullandığı | İkame | Durum |
|---|---|---|
| Referans yay (irtifa + faz geçmişi) | ucuz pilot yay, N=40 gevşek tolerans | ✔ önceki makalede fiyatlanmış: medyan **%8** gravite işi |
| Sayısal `Φ` | **pilot yay üzerinde düşük-dereceli varyasyonel denklem** `Φ̇ = A₄₀(t)Φ` | plan; Kepler/HCW ablasyona indi (§3.4) |
| Referans alanın **noktasal** değerlendirmeleri `Δa_ref(t)` | `P_n` spektrumu (**genlik**, D53) + ileri bant probu (**yön**) | §3.1 — **asıl açık problem** |

### 3.1 Açık problem: atlanan kuvvetin **yönü** ★

L2'nin L1'den güçlü olmasının tek sebebi işaret/yön bilgisi. Bir kuyruk
ölçütü yalnızca genlik verir (`d̂ ≈ ‖Δa‖²`). Sadeleşme terimi `2⟨v̂(t,N), c_i⟩` için
**yön** şart. Bu, makalenin gerçek mühendislik icadı olabilir.

**Önerilen çözüm — bant probu.** Atlanan kuyruk
`Δa(t,N) = −Σ_{n>N} a_n(t)` ve `(R/r)^n` sönümü yüzünden ilk birkaç atlanan
bant baskın. Öyleyse

```
v̂(t,N) = −γ(N,r) · Σ_{n=N+1}^{N+k} a_n(t)
```

`k ∈ {1,2,3,4}`; `γ` **alanın ölçülen derece-varyans spektrumundan** gelir
(D53), uydurulmuş bir güç yasasından değil:

```
γ(N,r) = [ Σ_{n>N} σ_a²(n;r) / Σ_{n=N+1}^{N+k} σ_a²(n;r) ]^{1/2}
```

Yön probdan, genlik spektrumdan. `P_n` ~1800 sayılık 1D tablo; taşımak
uçurulabilirliği zedelemez.

**Maliyet — `2k/N` yalnızca eş-konumlu prob için geçerli (D50).**

Zaten değerlendirme yapılan bir noktada bantları `N+k`'ya uzatmak
`(N+k)² − N² ≈ 2kN` işlem, yani `2k/N` ≈ %5 (N=120, k=3): ALF rekürsiyonu
artımlı, prob onu miras alıyor.

**Ama ileri prob eş-konumlu değil.** Yeni bir konumda `P_{n,m}` rekürsiyonu
sıfırdan koşmak zorunda; tek bir bandı `O(N)` işlemle almanın standart bir
yolu yok. Yani **her ileri prob noktası ≈ bir tam sentez.**

Aralığı kapsayan prob için `n_probe ≈ Δt_dec/τ_corr`, aralıktaki çağrı
`n_RHS ≈ Δt_dec·r_RHS`, dolayısıyla

```
prob ek yükü ≈ n_probe/n_RHS = 1/(τ_corr · r_RHS)      — Δt_dec'ten BAĞIMSIZ
```

Ama bu **yerel** bir orandır. Bütçeye giren sayı yay integralidir (D67):

```
ek yük = [ ∫₀^T dt/τ_corr(t) ] / N_RHS,toplam = [ ∫ N(t)v(t)/(πr(t)) dt ] / N_RHS
```

`τ_corr` en kısa olduğu yerde uydu en az zaman geçiriyor ve orada `r_RHS` de
yüksek — iki etki zıt yönde. Arşiv yaylarından: dairesele yakın **~%19**,
perilun 50/apolun 1000 **~%14**, perilun 30/apolun 2500 **~%12**. Eksantrik
yörüngeler uzun apolun kolları neredeyse hiç prob istemediği için **daha
ucuz**. Yani düşük-orta derece bandında **~%12–19**; taslağın önceki sürümündeki "%21–53" yerel `τ_corr`
ile küresel `r_RHS`'yi çarpan bir kategori hatasıydı. Prob yine de `B+`'ın
baskın kalemi.

İki özellik maliyeti düşürüyor ama **bir prob noktasının içinde**, noktalar
arasında değil:

1. **Paylaşılan yığın.** Bantlar bir kez `N_max+k`'ya kadar hesaplanır; her
   adayın `v̂`'si kısmi toplamla çıkar, ek sentez yok.
2. **Pencere.** `𝒩_q = {N_{j−δ},…,N_{j+δ}}`, `δ` aday-ızgara indis
   yarı-genişliği (D32); indüklenen derece açıklığı `Δ_span` yığın derinliğini
   belirler.

**Sonuç bir eğri, tek sayı değil (D51).** WP5 `κ_eff(n_probe)` üretiyor:
aralık boyunca ortalanmış etkin yön isabeti, prob nokta sayısına karşı, ve
karşısında ölçülen maliyet. İki uç ilan edilmiş:

| Varyant | Prob | Maliyet | Yön |
|---|---|---|---|
| `C-plan` | `n_probe` ayrı ileri nokta | ~%12–19 (yay integrali) | aralığı kapsar |
| `C-lite` | aralığın ilk RHS çağrısı yükseltilmiş derecede | ~%2 | yalnız aralık başı |

`C-lite` eş-konumlu olduğu için `2k/N` **ona** gerçekten uygulanır.

**Anahtarlama ucuz (D58).** Önceki kampanyada radyal politika ~16 000
anahtarlamada yalnızca %7 fazla çağrı harcadı → ~0.45 çağrı/anahtarlama.
`Δt_dec = 120 s`'de 5040 sınır → ~%2. Pahalı olan prob.

**Genlik `γ` gerçek spektrumdan (D53).** Uydurulmuş güç yasasından değil:
`γ(N,r) = [Σ_{n>N} σ_a²(n;r) / Σ_{n=N+1}^{N+k} σ_a²(n;r)]^{1/2}`, `P_n`'den.
Ay spektrumu tek bir güç yasası değil (`p_spec=2.13` vs `p_fit=1.76`), ve
`P_n` yalnızca ~1800 sayılık 1D tablo. Makalenin bilgi ifadesi de bu:
**taşınması ucuz olan spektrum (1D), pahalı olan doku (2D alan).**

**Neden ucuz bir vekil-model kolay bir çıkış yolu değil — dikkatli ifade.**
Kuyruk yönü `πr/N` açısal ölçeğinde değişiyor; N=300'de ~18 km. Bu ölçeği
çözmeyen bir `(r,φ,λ)` vekil-modeli yönü **temsil edemez**; çözen bir model ise
atlanan katsayılarla karşılaştırılabilir sayıda serbestlik derecesi taşır ve
onun depolama/değerlendirme maliyetini üstlenir.

Bu bir **çözünürlük ve maliyet** ifadesidir, herhangi bir fonksiyon sınıfı
hakkında imkânsızlık iddiası değil. Metinde "hiçbir vekil bunu yapamaz" gibi
kategorik cümle kullanılmaz — hakem "neden bir sinir ağı öğrenemesin?" diye
semantik tartışma açar. `abl-probe` vekili gerçekten kurar ve neyi koruduğunu,
neye mal olduğunu ölçer.

**Ölçülecek (WP5, propagasyonsuz):** `κ = cos∠(v̂_probe, Δa_gerçek)`,
`k`, `N` ve irtifanın fonksiyonu olarak. Karar eşiği: `κ` medyanı ~0.7'nin
altındaysa sadeleşme terimi uçurulabilir biçimde taşınamaz ve L3
tasarımı yalnızca L1 seviyesine iner (proje 🟡 bandına kayar).

**Yedek seçenekler** (WP5'te ablasyon olarak, `κ` yetersizse):
bant-artım vektörü ile ekstrapolasyon; `(r,φ,λ,N) → Δa` düşük-ranklı vekil
(alanı taşımak demektir, "referans alansız" iddiasını zayıflatır — bu yüzden
yedek); yalnız-genlik + rastgele-işaret null modeli (kazancın işaretten
geldiğinin negatif kontrolü).

### 3.2 Açık problem: `c(t)` tek bir geri STM geçişinden çıkmaz

Danışman haklı: `c_i = Σ_k Q_ik u_k` diğer epoklarda seçilen derecelere
bağlı, yani `c_i = c_i(N_1,…,N_{i-1},N_{i+1},…)`. Pilot yaydan tek bir geri
adjoint geçişi bunu vermez.

**Çözüm — İteratif İleri–Geri Derece Tahsisçisi (IFBDA).** Uçuş öncesi,
pilot yay üzerinde:

```
N⁽⁰⁾ ← L1 çizelgesi (duyarlılık ağırlıklı, işaretsiz)
for j = 1..J:                                  # J = 2..3
    ileri:  v̂(t, N⁽ʲ⁻¹⁾) biriktir → S_i
    geri:   A_i sonek toplamları → c_i
    güncelle: N⁽ʲ⁾(t) = argmin_N [ ‖v̂(t,N)‖²_Q + 2⟨v̂(t,N), c_i⟩ + λN² ]
    λ ← bütçeye bisekte
çıktı: faz-indeksli plan (devir, faz) → c, w, λ, N_plan
```

Yakınsama, iterasyon sayısı duyarlılığı ve `J`'nin maliyeti `B+`'a yazılır.

### 3.2b İleri prob — nedensellik düzeltmesi ★

Bir önceki sürümde "aralık içindeki problar bir sonraki sınırda kullanılır"
yazıyordu. **Bu yanlıştı** ve kendi `πr/N` argümanımızla çatışıyordu: 120 s
önce ölçülmüş bir yön çoktan dekorele olmuş olur.

Doğrusu: her karar sınırında **önündeki aralık boyunca ileri probe edilir.**

```
t_q sınırında:
  1. [t_q, t_{q+1}) için kendi kısa-ufuk durum tahminini üret
  2. n_probe = ⌈Δt_dec/Δt_probe⌉ noktada paylaşılan bant yığınını değerlendir
  3. her aday N için Σ_{i∈I_q} u_i(N) topla, N_q seç, aralık boyunca sabit tut
```

Meşruiyeti **iki farklı koherans ölçeğinden** geliyor:

| Ufuk | Durum tahmini hatası | `πr/N` karşısında |
|---|---|---|
| Bir karar aralığı (~120 s) | geminin kendi kısa-yay tahmini | **çok altında → geçerli** |
| Bütün yay (pilot vs uçulan) | `T_coh` sonrası yüzlerce km | **çok üstünde → geçersiz** |

Uzun ufuktaki başarısızlık plan/prob ayrımını dayatıyor (§3.3); kısa ufuktaki
başarı ileri probu mümkün kılıyor. Sınırdan önceki hiçbir yön bilgisi yeniden
kullanılmaz. Derecenin aralık boyunca sabit tutulması ayrıca `C-plan`'ı
`A-sign` ile aynı anahtarlama serbestliğinde tutar.

**Açık kalan risk:** eğer sadeleşme teriminin işe yaraması için `Δt_dec`'in
`Δt_probe`'a doğru küçülmesi gerekiyorsa, bu bir bulgudur ve provoke ettiği
entegratör maliyetiyle birlikte raporlanır — `Δt_dec` sonuca bakılarak
seçilmez. Kararı `abl-grid` taraması verir.

### 3.3 Açık problem: pilot yay tutarlılık ufku ★

Bu ne planda ne danışmanda vardı ve L3'ün mimarisini belirliyor.

Kuyruk yönü `πr/N` ölçeğinde dekorele oluyorsa (N=120 → ~45 km) ve N=40
gevşek-toleranslı pilot yay yedi günde referanstan bunun çok üstünde
sapıyorsa, **pilot yayda probe edilmiş yön gerçek yayda geçersizdir**.

Ama iki nicelik farklı davranır:

| Nicelik | Doğası | Pilot yay hatasına dayanıklı mı |
|---|---|---|
| `v̂(t,N)` yerel atlanan vektör | yerel alan dokusu, `πr/N`'de dekorele | **✘** |
| `c_i` sadeleşme yönü | *biriken* yer değiştirmenin fonksiyonu; önceki makale defektin koherent olduğunu ölçtü (büyüme üsteli 1.9–2.6, sıfır-ortalama nullü 1.5) | **✔ muhtemelen** |
| `w(t)`, `λ` | düzgün, yörünge zaman ölçeğinde | ✔ |

**Bunun dayattığı mimari:**

- `c`, `w`, `λ`, `N_plan` → **uçuş öncesi**, pilot yaydan (IFBDA ile),
  **faz-indeksli** (D52: pilot yayın ~100 s zaman kayması perilun geçişi
  mertebesinde, zaman indeksi kayar).
- `v̂` → **uçuş anında ileri**, geminin kendi kısa-ufuk tahmini üzerinde
  (maliyet §3.1, ~%12–19 tam kapsama / ~%2 `C-lite`).

Yani `C-plan` bir "tablo yükle ve uygula" denetleyicisi değil: planlanmış bir
sadeleşme hedefi + çevrimiçi yerel prob. Bu hem uçurulabilir hem de öngörü
sızması içermiyor.

**Ölçülecek (WP6):** `T_coh` — pilot yayın konum hatasının `πr/N`'i aştığı
süre; ve `c_i`'nin pilot yay ile referans yay arasındaki kararlılığı. `T_coh`
yaydan çok kısaysa, periyodik yeniden planlama (alıcı-ufuklu) *isteğe bağlı
değil zorunlu* hâle gelir ve minimum kapsamda plana alınır (WP11).

### 3.4 `C-fb` — gerçekleşen iş geri beslemesi

Önceki makalenin en zayıf noktası sızıntıydı: çağrı başına eşitlenen bütçe
entegratörle temas edince %29 (β=0.5'te %120) aşıyordu. `C-fb` harcanan
`Σ N²`'yi RHS çağrılarında sayar ve `λ`'yı çevrimiçi ayarlar → sızıntı
**tasarımla** kapanır. Önceki makalenin en büyük ölçüm zaafı, yeni yöntemin
satış argümanı olur.

### 3.5 Ablasyonlar (L4)

Hepsi ucuz; varyasyonel tahminle skorlanır, propagasyon yok.

**Kanonik liste — 14 ablasyon.** Manuskript §9 tablosu bununla birebir aynı
olmak zorunda; sayı M2 maliyetini belirliyor (14 × 8 = 112 çözüm, D79) ve
gönderim öncesi tutarlılık kontrolüne tabi (WP19).

| # | Kod | Ne kapatılır | Ne öğrenilir |
|---|---|---|---|
| 1 | `abl-sign` | işaretli terim atılır; `C-plan` yerel çekirdek `K_i` seviyesine iner | **RQ2'nin uçurulabilir taraftaki cevabı** |
| 2 | `abl-phi` | `Φ` yerine sadece `(T−t)` ufuk ağırlığı → `C-tgo` | duyarlılık gerçekten gerekli mi |
| 3 | `abl-probe` | prob yerine `(r,φ,λ,N)` vekil-modeli | §3.1'in çözünürlük argümanını doğrudan ölçer |
| 4 | `abl-null` | prob yönü, genlik sabit tutulup rastgeleleştirilir | **negatif kontrol — kazancı yok etmeli** |
| 5 | `abl-kaula` | spektrum `P_n` yerine uydurulmuş güç yasası | genlik ikamesinin bedeli (D53) |
| 6 | `abl-pilot` | pilot yay yerine referans yay | öngörü bedeli |
| 7 | `abl-stm` | düşük-dereceli pilot STM yerine Kepler | kapalı-form yeter mi *(ana yöntem değil)* |
| 8 | `abl-J` | IFBDA iterasyonu `J = 1` | kaç geçiş gerekli |
| 9 | `abl-k` | prob derinliği `k = 1` | prob derinliğinin bedeli/kazancı |
| 10 | `abl-window` | aday penceresi `δ` daraltılıp genişletilir | plan/prob ayrımının bedeli |
| 11 | `abl-grid` | karar ızgarası `Δt_dec` taranır | sadeleşmenin istediği kadans vs entegratör maliyeti |
| 12 | `abl-lookahead` | ileri prob yerine yalnız sınırda tek prob | ileri probun gerekliliği (D26) |
| 13 | `abl-predictor` | iki-cisim öngörücü yerine düşük-dereceli mikro-propagasyon | öngörücü seçiminin bedeli (D31) |
| 14 | `abl-timeindex` | faz-indeksli plan yerine zaman-indeksli | pilot yayın zaman kaymasının bedeli (D52) |

Adlandırma küçük harfle sabit (`abl-stm`, `abl-phi`) — manuskriptte `\code{}`
içinde geçtiği için büyük/küçük harf farkı iki ayrı ablasyon gibi görünüyordu.

---

## 4. Neden bu makale, ne söylüyor

Anlatı hattı önceki makaleden kesintisiz devam ediyor:

1. Önceki makale: kuvvet defekti yörünge hatasını sıralamıyor; sıralamayı
   `Φ` eşleşmesi belirliyor.
2. Bu makale, adım 1: o hâlde tahsisi doğrudan `Φ`-ağırlıklı işaretli amaç
   üzerinde yap → `A-sign`, ve ne kadar kazanç *masada* olduğunu ölç (RQ1).
3. Adım 2: kazancın kaynağını ayır — ağırlık mı, işaret mi (RQ2).
4. Adım 3: referans alan/yay olmadan ne kadarı yakalanır (RQ3). Cevap bant
   probu + IFBDA + pilot STM; `πr/N` dekorelasyonu bunun neden gerekli
   olduğunu söylüyor.
5. Kalibrasyon devamlılığı: `d̂`'nin dayandığı `p_fit = 1.76` önceki makalenin
   kendi katkısı.

---

## 5. Makale iskeleti

1. Giriş — force–trajectory boşluğu, bırakılan yön, RQ1–3.
2. Ayrıştırılamayan tahsis problemi — §1, neden Φ-ağırlıklı gevşetme yetmez.
3. Algoritma — `Q`'nun önek/sonek yapısı, koordinat inişi, Frank–Wolfe alt
   sınırı, sertifikalı boşluk.
4. Erişilebilir kazanç — `A-sens` ve `A-sign`; RQ1 ve RQ2.
5. Uçurulabilir denetleyici — bant probu, dekorelasyon argümanı, IFBDA,
   tutarlılık ufku, bütçe geri beslemesi.
6. Kampanya — popülasyonlar, bütçe ızgarası, gerçekleşen-iş eşleştirmesi.
7. Ablasyonlar — RQ3'ün ayrıştırılması.
8. Sınırlar — geometri, yay uzunluğu, ızgara, çözünürlük tabanı, sertifika
   boşluğu.

Hedef dergi: **JGCD** (mevcut kampanya takımı zaten JGCD'ye kalibre).

---

## 6. Şekil ve tablo listesi (ön taslak)

| # | İçerik |
|---|---|
| F1 | Merdiven: `F-op → R-int → A-force → A-sens → A-sign` hata/bütçe |
| F2 | Sertifikalı boşluk: FW alt sınırı ile koordinat inişi üst sınırı |
| F3 | Prob isabeti `κ` — `k`, `N`, irtifa |
| F4 | Dekorelasyon ölçeği `πr/N` ve prob kadansı |
| F5 | `C-plan` vs `F-op` yörünge yörünge, eşit `B2`'de (önceki makalenin β=1 paneli biçiminde) |
| F6 | Yakalama oranı `f` — bütçe ızgarası boyunca |
| F7 | Ablasyon şelalesi — hangi bileşen ne kadar kaybettiriyor |
| F8 | Varyasyonel parite: öngörülen vs propagate edilen oran |
| T1 | Politika tanımları ve bilgi kullanımı (NOTATION §1'in makale hâli) |
| T2 | Maliyet merdiveni `B1/B2/B3/B+` |
| T3 | `A-sign` sertifika boşluğu, tasarım tasarım |
| T4 | Ana verdict tablosu: çözünen sayımlar, `ρ`, popülasyon popülasyon |
| T5 | Ablasyon tablosu |
| T6 | Izgara ve bin kontrolleri |

---

## 7. Future work (ilk makalede **yok**)

Bilinçli olarak kapsam dışı; kaybolmamak için burada kayıtlı:

1. Alıcı-ufuklu / MPC tahsis — WP11 yalnızca `T_coh` gerektirirse minimum
   kapsamda girer.
2. Anahtarlama cezası ve Viterbi/trellis formülasyonu (bit tahsisi bağı).
3. Sağlam / dağıtımsal tahsis: faz MC topluluğu üzerinde ortalama amaç.
4. Alan transferi (Mars, Dünya) — R16 altyapısı mevcut.
5. GRAIL katsayı kovaryansı ile riske göre tahsis.
6. Sürekli derece + kesirli maliyet modeli.
7. Amaç seçimi çeşitlemesi: terminal durum, kovaryans izi, sadece in-track.
8. Operasyonel duvar-saati gösterimi (LRO / Kaguya benzeri).
