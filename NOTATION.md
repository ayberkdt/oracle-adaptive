# NOTATION — dondurulmuş adlar, semboller, tanımlar

Bu dosya kampanya başlamadan **dondurulur**. Önceki makalede en çok zaman
kaybettiren şey ölçüm konvansiyonlarının karışmasıydı (üç ayrı maliyet
tanımı, iki ayrı üstel, iki ayrı skorlama düzeyi). Burada her ad bir kez
tanımlanır; kod, tablo, metin ve ön-tescil aynı adı kullanır.

---

## 1. Politikalar (tahsis üreticileri)

| Kod | Ad (metinde) | Bilgi kullanımı | Uçurulabilir |
|---|---|---|---|
| `F-op` | operational fixed comparator | bütçe | ✔ |
| `F-env` | fixed-family lower envelope | bütçe + **her yörüngenin sonucu** (post-hoc) | ✘ |
| `R-rad` | budget-calibrated radial endpoint | yarıçap | ✔ |
| `R-int` | interior member (k=0.5) | yarıçap | ✔ |
| `A-force` | force-level allocation benchmark (önceki O26) | referans yay + referans alan | ✘ |
| `A-sens` | sensitivity-weighted allocation benchmark (L1) | + yerel duyarlılık çekirdeği `K_i = K(t_i,t_i)` (türetilmiş, D75/D110) | ✘ |
| `A-sign` | **linearized trajectory-aware allocation benchmark** (L2) | + `Q`'nun köşegen dışı blokları | ✘ |
| `A-sign+` | aynısı, politikanın kendi yörüngesinde sabit nokta (L2+) | + propagasyon | ✘ |
| `C-plan` | planned trajectory-aware schedule (L3b) | pilot yay + ileri prob + spektrum `P_n` | ✔ |
| `C-lite` | ucuz varyant: eş-konumlu tek prob (D51) | aynısı, `n_probe = 1` | ✔ |
| `C-fb` | budget-feedback variant (L3c) | + gerçekleşen iş durumu | ✔ |
| `C-tgo` | time-to-go only ablation (L4-abl-Φ) | kalan süre | ✔ |

> **`F-env` özel durumu.** Post-hoc alt zarftır; her yörüngenin hatasını
> gördükten sonra en iyi sabit dereceyi seçer. Uçurulabilir bir yöntem
> **değildir**. `C-plan`'ın onu yenmesi *bonus*, yenememesi *başarısızlık
> değil*. Bu ayrım OUTCOMES'ta bant sınırlarına gömülüdür.

---

## 2. Maliyet tanımları

Üç düzey artı bir ek yük kalemi. Hiçbiri diğerinin yerine geçmez ve her
sayının yanında hangisi olduğu yazılır.

| Kod | Tanım | Ne zaman ölçülebilir |
|---|---|---|
| `B1` | nominal çağrı başına iş, `⟨N²⟩` referans çıktı ızgarasında | M1 (propagasyonsuz) |
| `B2` | gerçekleşen toplam kuadratik iş, `Σ_{RHS} N_k²` | ancak propagasyondan sonra |
| `B3` | ölçülen seri gravite-çekirdek zamanı | boşta makinede, 3 tekrar |
| `B+` | **denetleyici ek yükü**, `B2` ile **aynı birimde** (aşağı) | analitik + ölçülen |

### `B+`'ın birimi — dönüşüm kuralı (D107)

`B2 = Σ_{RHS} N²` kuadratik-iş birimindedir. `B+` ise iki farklı türde iş
içeriyor ve **toplanabilmesi için ikisi de aynı birime çevrilmelidir**:

| Kalem | Tür | Dönüşüm |
|---|---|---|
| Pilot yay | gravite değerlendirmesi | **doğrudan** `Σ N²` |
| Çevrimiçi ileri prob | gravite değerlendirmesi | **doğrudan** `Σ N²` (bant yığını dahil) |
| Konum öngörücüsü (Kepler) | analitik, alan yok | ölçülen zamandan çevrilir |
| IFBDA planlama süpürmeleri | CPU, gravite değil | ölçülen zamandan çevrilir |
| Çevrimiçi `argmin`, defter tutma | CPU, gravite değil | ölçülen zamandan çevrilir |

Gravite dışı kalemler için **eşdeğer kuadratik iş**:

```
B+,eş = t_+ / ( t_çekirdek / B2 )
```

`t_çekirdek` = o yayın ölçülen seri gravite-çekirdek zamanı, `B2` = aynı yayın
kuadratik işi; oran bir "birim iş başına saniye" katsayısı verir ve `t_+`
(gravite dışı CPU zamanı) bununla bölünerek `N²` birimine çevrilir. Katsayı
mimariye bağlı olduğundan `B3` ile aynı koşulda ölçülür (boşta makine, üç
tekrar) ve her tabloda hangi makinede ölçüldüğü yazılır.

**Neden şart:** aksi hâlde `B2 + B+` boyutsal olarak saniye ile `ΣN²`'yi
topluyor olurdu. Hakem bunu ilk sorar.

### Normalizasyon — iki ayrı β, karıştırılmaz

`β` **yalnızca `B1` düzeyinde** tanımlıdır; tahsis hedefidir, gerçekleşen iş
değil:

```
β  = B1 / N_crit²                       (nominal tahsis hedefi)
β₂ = B2 / (K_ref · N_crit²)             (gerçekleşen, çağrı sayısına normalize)
```

`K_ref` = **`F-op` komparatörünün** aynı yörüngede, aynı `B_tot`'ta yaptığı
RHS çağrı sayısı — konvansiyon budur, başka bir politikanınki değil. `B2`'yi
doğrudan `N_crit²`'ye bölmek çağrı sayısını da orana sokar ve iki farklı şeyi
tek sayıda toplar — yapılmaz.

### Eşleştirme sözleşmesi — çapa, kim ayarlanır, tolerans (D91)

Karşılaştırmalar **toplam bütçe** `B_tot` üzerinde eşleştirilir:

```
her aday:           B2 + B+  =  B_tot        (B+ yalnız C-* için, diğerlerinde 0)
sabit komparatör:   B2,F     =  B_tot
```

Denetleyici, kendi ek yükü yüzünden gravite değerlendirmesine **daha az** bütçe
bırakır. Sabit komparatöre hayali bir `B+` **eklenmez**.

**Çapa (anchor).** `B_tot` her yörünge ve her `β` için **önceden ilan edilmiş
bir sayıdır**: `F-op(β)`'nın o yayda **gerçekleşen** toplam kuadratik işi.
`β = 1`'de bu, kritik derecenin arşivlenmiş yayıdır — yani çapa zaten var ve
yeniden hesaplanmıyor.

**Kim ayarlanır: aday, komparatör değil.** Önceki makale komparatörü adaya
göre kaydırıyordu; burada bu **yapılamaz**, çünkü `f`'in payı ve paydası aynı
`E_F-op`'u kullanmak zorunda:

```
f = (E_F-op − E_C-plan) / (E_F-op − E_A-sign)
```

Komparatör her adaya göre kayarsa payda ve payda **iki farklı** `F-op` örneği
olur ve `f` tanımsızlaşır. Bu yüzden `F-op` sabit çapa, **her aday ona
kalibre edilir**.

**Nasıl.** Bir çizelgenin gerçekleşen işi ancak propagasyondan sonra bilinir
(nominal `B1` kalibrasyonu `B2`'yi vermez — önceki makalenin ölçtüğü medyan
%29'luk sızıntı tam olarak budur). Dolayısıyla her aday için:

```
propagate → B2 (+B+) ölç → λ'yı ayarla → tekrarla,  |B2+B+ − B_tot|/B_tot ≤ %2 olana dek
```

Tipik olarak aday başına **2–3 propagasyon**. Bu, M3 ve sonrasının yay
sayısını doğrudan belirliyor (ROADMAP) ve gerçekleşen eşleştirmenin bedelidir;
"tek propagasyon başına bir aday" varsayımı yanlış olurdu.

**Kural:** `C-plan`, `C-lite` ve `C-fb` için raporlanan her bütçe `B+`'ı
içerir; ek yük dipnota atılmaz. `A-sign` ve `R-int` uçurulabilir olmadıkları
veya ek yükleri olmadığı için tam `B_tot`'u gravite değerlendirmesine harcar.

---

## 3. Amaç ve metrikler

`H_r = [I₃ 0]` konum seçicisi **açık yazılır**; iç integral 6 boyutlu bir durum
pertürbasyonudur ve konum bloğu seçilmeden norm alınamaz (birim karışır):

```
E² = (1/T) ∫₀^T ‖ H_r ∫₀^t Φ(t,τ) B Δa(τ) dτ ‖² dt      (birincil amaç)
```

Ayrık hâlde `M_j = Φ(t_j,t_0)ᵀ H_rᵀ H_r Φ(t_j,t_0)`, dolayısıyla `M_j`, `A_i`
ve `Q_ik` blokları **6×6**'dır (`u_i` bir durum pertürbasyonu; seçim `M_j`'nin
içinde). Terminal özel hâli:

```
J_T = ‖ H_r Φ(t_M,t_0) Σ_i u_i ‖²
```

— taşıma matrisi düşürülmez, çünkü `u_i`'ler `t_0`'a taşınmıştır.

### Amaç çekirdeği `K(τ,σ)` ve yerel köşegeni `K_i` (D110)

`J` sürekli hâlde bir çekirdeğe karşı **çift integraldir**:

```
J = ∫∫ Δa(τ)ᵀ K(τ,σ) Δa(σ) dτ dσ
```

`τ = σ` kümesinin ölçüsü sıfır olduğundan "köşegen dışını at" ayrık `Q`'da
harfiyen alınabilecek bir işlem **değildir**. İyi tanımlı olan, **çekirdeğin
köşegen üzerindeki değeridir**:

```
K(t,t) = (1/T) ∫_t^T Bᵀ Φ(τ,t)ᵀ H_rᵀ H_r Φ(τ,t) B dτ
K_i    = K(t_i,t_i) = (1/T) Bᵀ Φ(t₀,t_i)ᵀ A_i Φ(t₀,t_i) B
```

İkinci satır cocycle özdeşliği `Φ(τ,t_i) = Φ(τ,t₀)Φ(t₀,t_i)` ile birinciden
çıkar; yani `K_i` **`A_i` sonek dizisinden okunur**, ek maliyeti yoktur.
Blok boyutu `A_i` gibi 6×6 değil, `B`'nin sütun sayısı kadar (ivme uzayı, 3×3).

**Neden ham `Q_ii` değil.** `u_i ∝ Δt_i` olduğundan `u_iᵀ Q_ii u_i` ifadesi
`Δt_i²` taşır: ızgara inceldikçe sistematik olarak sıfıra gider ve
uyarlanabilir ızgarada `Δt_i` yay boyunca değiştiğinden yay içi sahte bir
ağırlık eğilimi üretir. `A-sens` bu yüzden `Σ_i Δt_i Δa_iᵀ K_i Δa_i` biçiminde
kurulur — eşdeğer olarak `Σ_i u_iᵀ Q_ii u_i / Δt_i`, yani fazladan bir `Δt`
çıkarılmış hâli. `Q` türetimi, önek/sonek yapısı ve FW sertifikası
etkilenmez; bu yalnızca **köşegen kısıtlamasının** tanımıdır.

`K` ile `S_j` (önek toplamı) ayrı sembollerdir; çekirdek için `S` kullanılmaz.

### Bütçe ile kısıtın ilişkisi

`B1 = ⟨N²⟩` çağrı başına **ortalama**; optimizasyon kısıtı ise toplam:

```
ana biçim (her zaman):   Σ_q W_q N_q² ≤ B = B1 · T ,   W_q = Σ_{i∈I_q} Δt_i
düzgün ızgara özel hâli: Σ_i N_{g(i)}² ≤ M · B1          (W_q = |I_q|·Δt)
```

`B1` bir **ortalama**, kısıt ise bir toplam; bağlantı bir kez kurulur.
Uyarlanabilir ızgarada toplam **zaman-ağırlıklı** olmak zorundadır, yoksa
perilun (örnek yoğunluğu yüksek olduğu için) sistematik olarak fazla
ağırlıklandırılır. Toplamın karar aralıkları yerine biriktirme örnekleri
üzerinden yazılması ayrıca, aralıklar eşit uzunlukta olmadığında da doğru
ağırlıklandırır.

| Sembol | Tanım |
|---|---|
| `E` | **propagate edilmiş** yedi günlük konum RMS hatası, referansa karşı |
| `Ê` | **birinci mertebe öngörülen** hata, `Ê = √J` — propagasyon yok |
| `ρ(X, Y)`, `ρ̂(X, Y)` | `E_X/E_Y` (propagate) ve `Ê_X/Ê_Y` (öngörülen) — **X komparatör, Y aday**; `> 1` Y lehine. Yörünge başına oranların **medyanı** (medyanların oranı değil) |
| `f` | yakalama oranı (aşağıda; geçerlilik koşuluna tabi) |
| `G_sens`, `G_sign` | log-iyileşme adımları (aşağıda) |
| `M_res` | çözünürlük marjı; `M_res > 1` ise karşılaştırma **çözünür** |
| `κ` | prob yön isabeti: `cos∠(v̂_probe, Δa_ref)` |
| `T_coh` | tutarlılık ufku: pilot–referans konum farkının `πr/N_plan(t)`'yi aştığı ilk an. Eşik zamanla değişir, çünkü `N_plan` değişir |

### `Ê` ile `E` asla karıştırılmaz (D37)

M1 ve M2'de **hiç propagasyon yok**, dolayısıyla oradaki her hata öngörülmüş
hatadır:

```
Ê = √J    (birinci mertebe, referans yay üzerinde, M1–M2)
E         (propagate edilmiş yedi günlük RMS, M3 ve sonrası)
```

- `Ê` kullanan nicelikler: G1, G2, H1a, H2, H6, H7, `G_sens`, `G_sign`,
  merdiven tablosu, sertifika (`J` zaten `Ê²`).
- `E` kullanan nicelikler: H1b, H3, H4, H4b, H5, `f`, kampanya verdiktleri.
- `ρ` her ikisiyle de yazılabilir ama **hangisi olduğu her tabloda belirtilir**:
  `ρ̂(X,Y)` öngörülen, `ρ(X,Y)` propagate edilmiş.

İkisinin karışması, önceki makalenin ölçüm-konvansiyonu tuzaklarının tam
muadili olurdu: bir sayı iki farklı şeyi gösterir ve fark yalnızca hangi
aşamada üretildiğine bakılarak anlaşılır.

**Oran yönü — her yerde açık yazılır.** `A-sign / R-int` gibi kısaltmalar
kullanılmaz; bir okur bunu hata oranı sanıp yönü ters çevirebilir. Doğru
biçim `ρ(R-int, A-sign)`: paydada aday, payda komparatör, `> 1` aday lehine.

**Yakalama oranı ve geçerlilik koşulu.**

```
f = (E_F-op − E_C-plan) / (E_F-op − E_A-sign)
```

`f` **yalnızca** şu iki koşul birlikte sağlandığında tanımlıdır:
1. `E_A-sign < E_F-op` (payda pozitif), **ve**
2. bu fark o yörüngenin çözünürlük zarfını geçiyor (`M_res > 1`).

Aksi hâlde `f = N/A` raporlanır, sıfır veya çok büyük bir sayı değil. Payda
sıfıra yaklaştıkça `f` patlar; koşul tam olarak bunu engellemek için var.

**`f`'in üç ayrıntısı, koşudan önce bağlanmış (D38):**

1. **Bütçe asimetrisi kasıtlıdır.** `A-sign` uçurulabilir değil, ek yükü yok:
   tam `B_tot` alır. `C-plan` ise `B_tot − B+` alır. Yani `f`, uçurulabilirliğin
   bedelini *içerir*; bu bir kusur değil, ölçmek istediğimiz şeyin ta kendisi.
2. **`f > 1` mümkündür** ve hata değildir: `A-sign` referans yay üzerinde
   birinci mertebe optimumdur, propagate edilince gerçek yayda daha kötü
   olabilir. Çıkarsa raporlanır ve WP17 ile birlikte okunur.
3. **Payı ve paydası aynı düzeyden gelir:** ikisi de `E` (propagate),
   `Ê` değil. `f` M3'ten önce hesaplanamaz.

**İyileşme adımları — ham fark değil, log-oran.** Ham hata farkı yörünge
ölçeğine bağlı olduğundan medyan alınamaz. Bunun yerine:

```
G_sens = log( Ê_A-force / Ê_A-sens )     (ağırlıklandırmanın kattığı)
G_sign = log( Ê_A-sens  / Ê_A-sign )     (işaretli eşleşmenin üstüne kattığı)
```

İkisi de pozitifse iyileşme var. H2 buna göre `G_sign > G_sens` biçiminde
yazılır.

---

## 4. Izgaralar — ikisi ayrıdır

Önceki makalede tek ızgara vardı; burada ayrılması zorunlu.

| Izgara | Sembol | Ne için | Varsayılan |
|---|---|---|---|
| **Biriktirme ızgarası** | `Δt_acc`, `M` örnek | `Δa`, `Φ`, `S_j` integrali | **uyarlanabilir** (aşağı) |
| **Karar ızgarası** | `Δt_dec`, `K_dec` aralık | `N`'nin parçalı sabit olduğu aralık — **karar değişkeni burada** | 120 s (kontrol edilir) |
| **Prob ızgarası** | `Δt_probe` | kuyruk yön probunun tazelenme sıklığı | `≈ τ_corr` (aşağı) |

**Biriktirme ızgarasının iki düğüm kümesi var (D124).** Hücre **kenarları**
`e_0 = 0 < … < e_M = T` ve hücre **orta noktaları** `m_i = (e_i+e_{i+1})/2`.
`Δa`, `Φ(t_0,·)` ve `u_i` **orta noktalarda**; amaç fonksiyonu ve `M_j`
**kenarlarda** örneklenir. İç ağırlık `Δt_i = e_{i+1} − e_i` (orta nokta
kuralı, ikinci mertebe, `j`'den bağımsız → önek yapısı korunur); dış ağırlık
`ω_j` kenarlarda yamuk. Bir propagasyonun ziyaret etmesi gereken epok kümesi
ikisinin birleşimidir (`2M+1` nokta). Gerekçe: PLAN §1.2.

Dekorelasyon zaman ölçeği **uydu yarıçapıyla** yazılır (yüzey yarıçapıyla
değil) ve yörünge boyunca sabit değildir:

```
τ_corr = π r / (N v)
```

`τ_corr` ızgara inceltmesinde **en yüksek aday dereceyle** hesaplanır: ızgara
bütün adaylarca paylaşılıyor ve en ince doku en yüksek truncation'a ait; daha
düşük bir dereceyle inceltmek kampanyanın kurulduğu karşılaştırmaları
örtüşmeye (aliasing) sokar. İnceltme **eşdağıtım** ile (D125).

Referans yörünge: perilun 50 km (`r=1788` km), apolun 2500 km (`r=4238` km),
`a=3013` km, `e=0.407`, periyot 4.12 sa. **Izgara tek bir paylaşılan derecede
inceltildiği için tablonun her satırı aynı `N`'i kullanır** — bu düzeltilmiş
tablodur (D129).

| İnceltme derecesi | `τ_corr` perilun | `τ_corr` apolun | oran | 7 günde hücre (`n_s=2`) |
|---|---|---|---|---|
| `N = 120` | 23.8 s | 133.9 s | 5.6× | ~20 500 |
| `N = 300` | 9.5 s | 53.6 s | 5.6× | ~51 100 |
| `N = 600` | 4.8 s | 26.8 s | 5.6× | ~102 300 |

**Beş buçuk kat değişiyor, elli değil.** Oran `(r_a/v_a)/(r_p/v_p)`'dir ve
dereceden bağımsızdır; eski tablodaki "~50 kat" perilunda `N=300`, apolunda
`N=20` kullanmaktan geliyordu — paylaşılan tek ızgarada böyle bir seçenek yok.
Yine de 120 s'lik düzgün bir ızgara perilunda işaretli integrali örtüşmeye
sokar; asıl mesele budur ve o değişmedi.

### Biriktirme ızgarası uyarlanabilirdir (D55)

Düzgün 10 s'lik ızgara 7 günlük yayda 60 481 örnek üretir ve çoğu hiçbir şeyin
olmadığı apolundadır. Izgara `τ_corr`'e göre inceltilir (perilunda yoğun,
apolunda seyrek). İki sonucu var:

1. Depolama **2.5 kat** düşer — düzgün-en-ince ızgaraya karşı, tek paylaşılan
   derecede (D129). Mutlak boyut inceltme derecesiyle **doğrusal** büyür:
   `M × |𝒩| × 3 × 8` bayt, `|𝒩| ≈ 60` ile **`N=300`'de ~74 MB/yörünge**,
   `N=600`'de ~147 MB. Eski "~87 → ~9 MB" hesabı hem 50 katlık sahte orandan
   hem de düşük bir dereceden geliyordu.
   **Sonuç `float64` kararını (D61) değiştirmiyor** ama gerekçesini
   değiştiriyor: tablo artık "küçük olduğu için" değil, **bellek-eşlemeli ve
   sıralı okunduğu için** `float64` kalıyor; RAM'e sığması gerekmiyor. 26
   yörüngelik panelde disk ~2 GB. `float32` vs `float64` verdict parity testi
   yine zorunlu ve artık daha da gerekli.
2. `B1` artık **zaman-ağırlıklı** ortalamadır:
   ```
   B1 = Σ_i N_{g(i)}² Δt_i / T ,     B = Σ_i N_{g(i)}² Δt_i / Δt_ref
   ```
   Düzgün ızgarada eski tanıma indirgenir. Ağırlıksız ortalama, uyarlanabilir
   ızgarada perilunu sistematik olarak fazla ağırlıklandırır.

Yakınsama testi (WP4) hem düzgün hem uyarlanabilir ızgarada koşulur; kampanya
yakınsak olanla yürür.

### Eşleme ve karar değişkeni — dondurulmuş

`g(i) = q` biriktirme örneği `i`'nin ait olduğu karar aralığını versin,
`I_q = {i : g(i) = q}`. O hâlde

```
N_i = N_{g(i)} ,     karar değişkeni:  (N_1, …, N_{K_dec}) ,   K_dec ≪ M
```

**Optimizasyon her zaman `K_dec` boyutlu uzayda yapılır, `M` boyutlu uzayda
değil.** Aksi hâlde `A-sign`, uçurulabilir yöntemden yalnızca daha çok bilgi
kullandığı için değil, dereceyi 10–30 saniyede bir değiştirebildiği için de
avantaj kazanır ve karşılaştırma geçersiz olur.

### Prob–karar ilişkisi — İLERİ prob (düzeltildi, D26)

**Önceki kural yanlıştı:** "aralık içindeki problar bir sonraki sınırda
kullanılır" demek, 120 s önce ölçülmüş bir yönle karar vermek demekti —
`πr/N` argümanının tam olarak imkânsız dediği şey.

Doğru kural: her karar sınırında **önündeki aralık boyunca ileri probe edilir.**

```
t_q sınırında:
  1. önümüzdeki [t_q, t_{q+1}) aralığı için kendi kısa-ufuk durum tahminini üret
  2. n_probe = ⌈Δt_dec / Δt_probe⌉ noktada bant yığınını değerlendir
  3. aday derecelerin u-katkısını topla, N_q'yu seç, aralık boyunca sabit tut
```

Bu meşrudur çünkü **iki farklı koherans ölçeği** var:

| Ufuk | Durum tahmini hatası | `πr/N` karşısında |
|---|---|---|
| Bir karar aralığı (~120 s) | iki-cisim öngörücüsü, metre mertebesi (D31) | **km ölçeğinin çok altında → ileri prob geçerli** |
| Bütün yay (pilot vs uçulan, 7 gün) | yüzlerce km, `T_coh` sonrası | **çok üstünde → önceden hesaplanmış `v̂` geçersiz** |

Uzun ufuktaki başarısızlık plan/prob ayrımını dayatıyor; kısa ufuktaki başarı
ileri probu mümkün kılıyor. Sınırdan önceki hiçbir yön bilgisi yeniden
kullanılmaz.

### Prob konumlarını ne öngörüyor? (D31)

`t_q`'da `N_q` henüz seçilmedi, dolayısıyla entegratör o aralığı geçmedi. Prob
noktaları nereden geliyor?

**Döngü yok, çünkü prob konumu dereceye neredeyse hiç bağlı değil.** İki aday
derece bir aralıkta en fazla `½‖Δa‖ Δt_dec²` kadar ayrışır — `‖Δa‖ ~ 1e-6 m/s²`
ve 120 s için **milimetre** mertebesi, `πr/N` ~ 9–45 km'ye karşı. Yani
öngörücünün `N_q`'yu bilmesi gerekmiyor.

**Öngörücü:** mevcut durumun **iki-cisim (Kepler) propagasyonu** — hiç alan
değerlendirmesi yok, analitik. Hatası ihmal ettiği perturbasyondan geliyor:
`½‖a_pert‖ Δt_dec²` ≈ **metre** mertebesi, yine km toleransının çok içinde.
Maliyeti `B+`'a yazılır.

**Ablasyon `abl-predictor`:** Kepler yerine `N_plan`'da düşük-dereceli
mikro-propagasyon; `κ` ve nihai hata farkı ölçülür. Seçim varsayılmaz.

### Aday penceresi ve paylaşılan bant yığını

Çevrimiçi `argmin` birden çok aday dereceyi karşılaştırıyor; her aday için ayrı
prob **gerekmez**:

- **Paylaşılan yığın:** bantlar `a_n` bir kez `N_max+k`'ya kadar hesaplanır;
  her aday `N`'nin `v̂`'si kısmi toplamla çıkar, ek alan değerlendirmesi yok.
- **Pencere:** `𝒩 = {N_1 < N_2 < …}` sıralı aday ızgarası, `j(q)` planın bu
  ızgaradaki indisi olsun. Aday kümesi

  ```
  𝒩_q = { N_{j(q)−δ}, …, N_{j(q)+δ} }
  ```

  **`δ` bir aday-ızgara İNDİS yarı-genişliğidir, derece değil (D32).** Böylece
  ızgara inceltilse/kabalaştırılsa da denetleyici karmaşıklığı (`2δ+1` aday)
  sabit kalır. Maliyeti belirleyen ise indis sayısı değil, indüklenen **derece
  açıklığı** `Δ_span = N_{j+δ} − N_{j−δ}`:

  ```
  yığın: N_{j−δ}+1 … N_{j+δ}+k
  ```

### Genlik tamamlaması `γ` — spektrumdan, güç yasasından değil (D53) ★

`γ(N,h)` = tam atlanan kuyruğun probe edilen bantlara oranı. Bunu uydurulmuş
bir Kaula güç yasasından almak **sistematik hata taşır**: önceki makale ay
spektrumunun tek bir güç yasası olmadığını ölçtü (`p_spec = 2.13` ile
`p_fit = 1.76` farkının sebebi tam olarak bu).

Bunun yerine alanın **gerçek derece-varyans spektrumu** `P_n` kullanılır:

```
γ(N,h) = √( Σ_{n>N}   σ_a²(n;r) )  /  √( Σ_{n=N+1}^{N+k} σ_a²(n;r) )
```

`σ_a(n;r)` önceki makalenin Eq. (degree-rms) ifadesi, `P_n`'den. `P_n`
yalnızca ~1800 sayılık **1D tablo** — yüklemek uçurulabilirliği hiç zedelemez.

Bu, makalenin bilgi argümanını da keskinleştiriyor:
**taşınması ucuz olan spektrum (1D); taşınması pahalı olan doku (2D alan).**

  Her iki sayı da (indis `δ` ve derece `Δ_span`) maliyet tablosunda yan yana
  raporlanır.

### Prob maliyeti — `2k/N` yanlıştı (D50) ★

Paylaşılan yığın ve pencere, **bir noktadaki adaylar arası** maliyeti düşürür.
Ama **noktalar arası** maliyeti düşürmez ve asıl kalem orada:

> Bir ileri prob noktası **yeni bir konumdur**. Orada `P_{n,m}` rekürsiyonu
> sıfırdan `N+k`'ya kadar koşmak zorundadır — tek bir bandı `O(N)` işlemle
> almanın standart bir yolu yoktur. Dolayısıyla **her prob noktası ≈ bir tam
> sentez**, `2k/N`'lik bir artım değil.

`2k/N` yalnızca prob, zaten değerlendirme yapılan bir noktada olsaydı
geçerliydi — ileri prob tanımı gereği öyle değil (D26 ile D50 birbirini
kısıtlıyor).

**Maliyet — yay integrali, nokta değeri değil (D67 ile düzeltildi).**
Aralık başına `n_probe ≈ Δt_dec/τ_corr`, `n_RHS ≈ Δt_dec·r_RHS`, dolayısıyla
*yerel* ek yük `1/(τ_corr·r_RHS)` ve `Δt_dec`'ten bağımsız. Ama bir yay için
raporlanacak sayı bu değil, **yay integrali**:

```
ek yük = [ ∫₀^T dt/τ_corr(t) ] / N_RHS,toplam
       = [ ∫₀^T N(t)v(t)/(πr(t)) dt ] / N_RHS,toplam
```

`τ_corr` en kısa olduğu yerde (perilun, yüksek `N`) uydu **en az zaman
geçiriyor**, ve orada `r_RHS` de yüksek. İki etki zıt yönde; yerel bir `τ_corr`
ile küresel bir `r_RHS`'yi çarpmak kategori hatasıdır.

Arşivlenmiş referans yaylardan hesaplanan yay-integrali değerler:

| Geometri | derece profili | ek yük |
|---|---|---|
| Dairesele yakın 50–100 km | N ≈ 120 | **~%19** |
| Perilun 50 / apolun 1000 km | N 200→40 | **~%14** |
| Perilun 30 / apolun 2500 km | N 300→20 | **~%12** |

Yani düşük-orta derece bandında **~%12–19**, önceki taslakta yazdığım "%21–53" değil. **Dereceye bağlı (D141).** Bu aralık, ölçüldüğü politika derecesi yazılmadan taşınıyordu. `n_probe = M/n_s` olduğu için ek yük `N` ile doğrusal: vekil yayda `N=120`→%8.5, `N=300`→%21, `N=600`→%42. Aşağıdaki sayı düşük-orta derece bandına aittir ve **tek sayı olarak raporlanmayacaktır**; WP5/WP16 derece bandı başına ölçer. `r_RHS(t)`'nin
perilunda yükselmesi integrale henüz katılmadı, dolayısıyla bu değerler hâlâ
üst taraftan sapıyor. Kesin sayı WP5/WP16'da ölçülür; formül yukarıdaki.

**Nitel sonuç değişmiyor:** prob hâlâ `B+`'ın baskın kalemi (pilot yay %8,
anahtarlama ~%2) ve hâlâ ölçülmesi gereken bir ödünleşim.

### Bunun sonucu bir eğri, tek bir sayı değil (D51)

Sorun kaçınılmaz olduğuna göre **ölçülür**. WP5 `κ_eff(n_probe)` üretir:
aralık boyunca ortalanmış etkin yön isabeti, prob noktası sayısının
fonksiyonu. Makale bir **maliyet–isabet ödünleşim eğrisi** raporlar.

İki uç ilan edilmiştir:

| Varyant | Prob | Maliyet | Yön bilgisi |
|---|---|---|---|
| `C-plan` | `n_probe` ayrı ileri nokta | yay integrali, **dereceye bağlı** (D141): %8.5 / %21 / %42 @ N=120/300/600 | aralığı kapsar |
| `C-lite` | aralığın **ilk RHS çağrısı** yükseltilmiş derecede (eş-konumlu) | ~%2 | yalnız aralık başı, `τ_corr` kadar geçerli |

`C-lite` eş-konumlu olduğu için `2k/N` ona **gerçekten** uygulanır.

### Anahtarlama ucuz, prob pahalı (D58)

Anahtarlama maliyeti önceki kampanyanın **ölçülmüş** %7 fazla RHS çağrısına
dayanıyor. Anahtarlama başına maliyet (~0.5 çağrı) ise **kestirim üstüne
kestirim**: %7'yi, radyal politikanın 10 km binlerle yaptığı anahtarlama
sayısına (~1–2×10⁴, tarafımdan tahmin) bölerek çıkıyor. Ölçülen sayı sağdaki
oran değil, %7'nin kendisi.

Bu kaydıyla: `Δt_dec = 120 s` ile 5040 anahtarlama → **~%2 mertebesi**.
Anahtarlama sayısı arşivden doğrudan sayılabilir ve WP16'da öyle yapılır.
Sonuç yine de sağlam: pahalı olan prob, anahtarlama değil.

`δ` ilan edilmiş bir parametredir, taranır, ve denetleyiciyi `A-sign`'a göre
kısıtlar — bu kısıt yakalama oranında görünür, gizlenmez.

---

## 5. Çizelge biçimi

| Biçim | Tanım | Kim kullanıyor |
|---|---|---|
| irtifa-binli | `N = N(h)`, 10 km binler | `R-rad`, `R-int` (önceki makalenin konvansiyonu) |
| zaman-indeksli | `N = N(t)`, karar ızgarasında | `A-sens`, `A-sign` (kıyaslar; referans yayda çözülüyor) |
| **faz-indeksli** | `N = N(devir, faz)` | **`C-plan`, `C-fb`, `C-lite`** |

Bu bilinçli bir sapmadır: yörünge-farkında tahsisin bütün noktası, derecenin
yarıçapın tek değerli fonksiyonu **olmaması**. Ama derece **zamanın** da
fonksiyonu olamaz — bkz. aşağı.

### Denetleyici planı neden faz-indeksli (D52) ★

N=40 gevşek-toleranslı pilot yay 7 günde referanstan ~100 km along-track
sapıyor; LLO'da 1.6 km/s ile bu **~60–100 s zaman kayması** demek. Perilun
geçişi birkaç dakika sürüyor — **aynı mertebe.**

`c` ve `K_i` perilun civarında baskın olduğundan, zaman-indeksli bir plan
onları yanlış epoklara hizalar: denetleyici, perilun geçişi için hesaplanmış
sadeleşme hedefini perilundan önce veya sonra uygular.

Faz indeksi bunu tamamen kaldırır: plan "**şu devirde, yörüngenin şu
fazında**, sadeleşmeyi şu yöne çevir" der. Gemi kendi fazını biliyor. Devirler
arası 100 s = ~0.014 devir, ihmal edilebilir.

Yarıçapın tek değerli fonksiyonuna geri dönmüyoruz: indeks **faz artı devir
numarası**, yani kalan ufuk hâlâ içeride. WP6 `c`'nin kararlılığını
**faz-indeksli** biçimde ölçer, zaman-indeksli değil.

---

## 6. Terminoloji

- Referans alan / referans yörünge → **"reference"**, "truth" değil.
  (Model-göreli bir ölçüm yapıyoruz, gerçek Ay alanına mesafe ölçmüyoruz.)
- **Bilgi sınırı (D62).** "Referans alan yok" demek yanlış olurdu: bir
  propagatör katsayıları tanım gereği taşır. Doğru ayrım **küresel/çevrimdışı
  model üst-verisi** ile **yörüngeye özgü değerlendirmeler** arasında:

  | | Denetleyiciye açık | Yasak |
  |---|---|---|
  | Alan | derece varyansları `P_n` (1D tablo) | `Δa_ref(t)` yay boyunca |
  | Yörünge | kendi hesapladığı ucuz pilot yay | referans yay |
  | Duyarlılık | düşük-dereceli `Φ` | tam `Φ` |

  Metinde "without the reference trajectory or pointwise reference-field
  evaluations" biçimi kullanılır, "reference field unavailable" değil.
- `A-sign` için **"oracle"** kelimesi, ancak sertifikalı boşluk eşiği
  sağlanırsa kullanılır (aşağıda). Aksi hâlde **"linearized trajectory-aware
  allocation benchmark"**.

### Sertifikalı boşluk — tanım dondurulmuş

```
g_J = (J_desc − L_FW) / J_desc            (amaç uzayında)
g_E = 1 − √(L_FW / J_desc) = 1 − √(1−g_J)  (hata uzayında)
```

`J = E²` olduğu için ikisi aynı sayı değildir. **Eşik hata uzayında
tanımlıdır:** medyan `g_E < 0.10` (kabaca `g_J ≈ 0.19`). İkisi de raporlanır.

`L_FW = max{0, iterasyonlar boyunca en iyi FW sınırı}` — FW dualite sınırı
geçerlidir ama erken iterasyonlarda negatif çıkabilir ve `J ≥ 0` olduğundan
sıfırla kırpılır. İterasyon bütçesi bittiğinde `L_FW = 0` ise sertifika
**boştur**: o yörünge için boşluk raporlanmaz, ayrı bir sütunda sayılır,
sessizce düşürülmez.

### Sayısal: telafili toplama (D56)

`S_j` altmış bin işaretli terimin toplamı ve **sadeleşme yöntemin ta kendisi**,
yani koşul sayısı yüksek. Önek/sonek toplamlarında **Kahan/Neumaier telafili
toplama** kullanılır. Çift duyarlıkta hata muhtemelen zaten kabul edilebilir
seviyede kalır, ama telafisiz toplamak, sonucu tam da sadeleşmeye dayanan bir
yöntemde gereksiz risktir.

### Diğer terim kuralları

- "optimal" kelimesi hiçbir politikaya yörünge seviyesinde uygulanmaz.
- Sonuçlar model-göreli ve test edilen tasarımlarla sınırlıdır.
- **`κ ≥ 0.7` eşiği — gerekçe düzeltildi (D74).** Önceki taslakta "`κ` doğrudan
  sadeleşme teriminin korunan kesridir" yazıyordu. **Bu doğru değil.** Terim
  `2⟨u_i(N), c_i⟩` ve `u_i = Φ B Δa Δt`, dolayısıyla

  ```
  ⟨Φ B v̂, c⟩ = ⟨v̂, z⟩ ,   z := (Φ B)ᵀ c
  ```

  Korunan kesir `⟨v̂,z⟩/⟨Δa,z⟩`'dir, yani `v̂` ile `Δa` arasındaki açıya değil,
  **ikisinin `z` ile yaptığı açılara** bağlıdır. `Φ B` bir lineer dönüşüm ve
  açıları korumaz. Uç örnek: `Δa ⊥ z` ise gerçek katkı sıfırdır ve `κ` ne
  olursa olsun prob sahte bir katkı üretebilir.

  Doğru okuma: `κ = 1` korunan kesri 1 yapar (genlik de doğruysa), ve düşük `κ`
  korunmayı bozar — yani `κ` **gerekli-benzeri bir gösterge**, kesrin kendisi
  değil. Eşik ilan edilmiş ve keyfîdir; `κ`'yı önde gelen gösterge olarak
  kullanıyoruz çünkü propagasyonsuz ve çizelgeden bağımsız ölçülebiliyor.
  Korunan kesrin kendisi (`⟨v̂,z⟩/⟨Δa,z⟩`, `A-sens` çizelgesinden türetilen
  `c` ile) WP5'te `κ`'nın **yanında** raporlanır — yeni bir kapı değil, aynı
  ölçümün ikinci çıktısı.

---

## 7. Dosya adlandırma

```
oa<NN>_<konu>.py              betik            (oa01_tabulate.py)
metrics/oa<NN>_*.json         ham çıktı
metrics/oa<NN>_*_table.tex    manuskript tablosu
metrics/oa<NN>_preregistration.json
figures/fig_oa_<konu>.pdf
make_figures_oa<NN>.py        kampanyaya özel — ortak betik ASLA çağrılmaz
```

Kampanya kodları `OA1..OAn`, iş paketleri `WP0..WP20` (WP17 ve WP20 a/b ayrık) (bkz. [WP.md](WP.md)).
