# DECISIONS — karar günlüğü

Verilmiş kararlar burada; bir daha tartışılmaz. Karar değişirse eski satır
silinmez, **üstü çizilir ve gerekçesiyle yenisi eklenir.** Amaç aynı tartışmayı
üç hafta sonra baştan yapmamak.

Durum: `KABUL` (uygulanacak) · `AÇIK` (karar bekliyor) · `DEĞİŞTİ`

---

## Tur 1 — danışman geri bildirimi (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D1** | L2'nin adı `A-sign` = **"linearized trajectory-aware allocation benchmark"**. "Oracle" kelimesi ancak WP7 sertifikası eşiği sağlarsa kullanılır. *(Eşik D29 ile hata uzayında sayısallaştırıldı: medyan `g_E < 0.10`.)* | KABUL | Sertifikasız "trajectory-level oracle" hakeme açık kapı. Dil kuralı PREREG OA-02'de peşinen bağlandı, sonuca göre seçilemez. |
| **D2** | `C-plan`, atlanan kuvvetin **yönünü** bant probuyla tahmin eder: `v̂ = −γ(N,h)·Σ_{n=N+1}^{N+k} a_n`. Yön probdan, genlik Kaula kuyruğundan. | KABUL | Kaula yalnız genlik veriyor; `2⟨v̂,c⟩` terimi yön istiyor. Bant maliyeti `≈2kN` vs `N²` → N=120,k=3'te ~%5. **→ D50 ile revize: bu yalnızca eş-konumlu prob için geçerli.** |
| **D3** | Kuyruk yönü için düzgün yüzey vekil-modeli (lookup/surrogate) **ana yöntem değil**, yalnızca yedek. | KABUL | Kuyruk `πR/N` ölçeğinde dekorele oluyor (N=300 → 18 km). Düzgün vekil bunu taşıyamaz. Bu, yöntemin gerekçesi ve makalede söylenmeye değer. |
| **D4** | `c(t)` tek bir geri STM geçişinden gelmiyor; `C-plan` bir **ileri–geri iterasyon** (IFBDA), `J = 2..3`. | KABUL | `c_i` diğer epokların derecelerine bağlı. Tek geçiş yanlış. |
| **D5** | Mimari ikiye ayrılır: `c(t)`, `w(t)`, `λ` **uçuş öncesi** pilot yaydan; `v̂(t,N)` **uçuş anında** gerçek durumda probe edilir. | KABUL | WP6'nın tutarlılık ufku problemi. `c` biriken ve koherent (dayanıklı), `v̂` yerel doku (dayanıksız). |
| **D6** | H1 ikiye ayrıldı: **H1a** allocation-space (`B1`, M1'de), **H1b** propagated (`B2`, M3'te). | KABUL | M1'de propagasyon yok, `B2` bilinemez. Önceki makalenin sızıntı sonucuna da doğrudan bağlanıyor. |
| **D7** | Komparatörler adlandırma ile ayrıldı: `F-op` (uçurulabilir, **yenmek şart**) ve `F-env` (post-hoc alt zarf, **bonus**). | KABUL | `F-env` uçurulabilir değil; ona kaybetmek başarısızlık sayılamaz. OUTCOMES bantları buna göre yumuşatıldı. |
| **D8** | Izgara yakınsaması (`Δt_acc ∈ {240..10}` s) **M1'de**, M7'de değil. | KABUL | 120 s'de işaretli integral aliasing riski var; sonradan bulunursa her şey tekrarlanır. |
| **D9** | **İki ızgara** ayrıldı: biriktirme (`Δt_acc`, ince) ve karar (`Δt_dec`, kaba). | KABUL | `Δa` ~10 s ölçeğinde değişiyor, `N(t)` yörünge ölçeğinde. Aynı ızgara olamaz. |
| **D10** | `Φ` için ana yöntem **pilot yay üzerinde düşük-dereceli varyasyonel denklem** `Φ̇ = A₄₀(t)Φ`. Kepler/HCW yalnızca `abl-STM` ablasyonu *(→ D71 ile `abl-stm` olarak küçük harfe sabitlendi)*. | KABUL | Ay mascon dinamiğinde HCW zayıf; pilot yay zaten koşuluyor, ek maliyet küçük. |
| **D11** | İlk makale **L0→L1→L2→L3** ve RQ1–RQ3 ile sınırlı. MPC, dağıtımsal tahsis, alan transferi, GRAIL kovaryansı → future work. | KABUL | Kapsam kontrolü. PLAN §7'de kayıtlı, kaybolmuyor. |
| **D12** | `A-sign` çizelgeleri **zaman-indekslidir**, irtifa-binli değil. | KABUL | Yörünge-farkında tahsisin noktası, derecenin yarıçapın tek değerli fonksiyonu olmaması. Kendi kontrolü var (WP12). |
| **D13** | `B+` (pilot yay + prob + IFBDA + çevrimiçi argmin) bütçeye **dahildir**, dipnota atılmaz. | KABUL | Aksi hâlde hakem haklı olarak reddeder. |
| **D14** | `A-sign+` (WP17) sabit-nokta testi plana alındı. | KABUL | `A-sign` referans yay üzerinde çözülüyor; politika değişince `Φ` ve `Δa` değişiyor. Bağlayıcı mı, ölçülür. |
| **D15** | Terminoloji: "reference", "truth" değil. "Optimal" hiçbir politikaya yörünge seviyesinde uygulanmaz. | KABUL | Önceki makalenin konvansiyonu; model-göreli ölçüm yapıyoruz. |

---

## Tur 2 — danışman geri bildirimi (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D16** | `A-sign` **karar ızgarasında** optimize edilir: `N_i = N_{g(i)}`, karar uzayı `K_dec` boyutlu. Biriktirme ızgarasında `M` boyutlu optimizasyon yasak. | KABUL | Aksi hâlde `A-sign`, uçurulabilir adaydan yalnızca daha çok *bilgi* kullandığı için değil, daha çok *anahtarlama serbestliği* olduğu için de kazanır; ölçtüğümüz şey ızgara inceliği olur. Gruplamanın maliyeti yok: blok-içi terim yürüyen önek toplamıyla `O(m)`, süpürme yine `O(M·\|𝒩\|)`. |
| **D17** | M3 propagasyon matrisine **`A-sign` ve `R-int` eklendi** (beş politika). `F-env` propagate edilmez, alt zarf olarak kurulur. | KABUL | H1b propagate edilmiş `E_A-sign` gerektiriyor ve `f`'in paydası odur; ikisi de listede yoktu. `R-int` OUTCOMES'ta 🟢 bandının şartı. |
| **D18** | Arşiv propagasyonları yalnızca gerçekleşen işleri yeni `B2` sözleşmesi altında eşleştirilebiliyorsa yeniden kullanılır; aksi hâlde yeniden koşulur ve her yeniden kullanım işaretlenir. | KABUL | Eski kampanya farklı bir eşit-iş konvansiyonunda koşuldu; sessiz yeniden kullanım farkı yönteme yazar. |
| **D19** | FW dili: LMO, `A-sens/L1` ile **aynı yapıda**, "L1'in kendisi" değil. Ayrıca LMO **tam** çözülür (kritik çarpanda kesirli karışımla), her iterasyonda doğrulanır. | KABUL | Amaç fonksiyonları farklı (`⟨∇_q J, u_q(N)⟩` vs `w_q d_q(N)`) ve `∇J` mevcut gevşek çözüme bağlı. Bisekti en yakın ayrık bütçede kesmek sertifikayı geçersiz kılar. |
| **D20** | Konum seçicisi `H_r = [I₃ 0]` açık yazılır; `M_j = Φᵀ H_rᵀ H_r Φ`. | KABUL | `\|_pos` gösterimi boyutsal olarak tartışmaya açık; hakem karşısında temiz olmalı. |
| **D21** | Oran her yerde `ρ(komparatör, aday)` biçiminde yazılır; `A-sign / R-int` gibi kısaltmalar kullanılmaz. | KABUL | Kısaltma hata oranı sanılıp yön ters okunabiliyor. |
| **D22** | H2 log-oranla yazılır: `G_sens = log(E_A-force/E_A-sens)`, `G_sign = log(E_A-sens/E_A-sign)`, hipotez `G_sign > G_sens`. | KABUL | Ham hata farkı yörünge ölçeğine bağlı, medyanı anlamsız. |
| **D23** | `f` yalnızca `E_A-sign < E_F-op` **ve** fark çözünürlük zarfını geçtiğinde tanımlıdır; aksi hâlde `N/A`. `β` yalnızca `B1` düzeyinde; gerçekleşen için ayrı `β₂`. Eşleştirme tek denklemle dondu: `B2,C + B+ = B_tot = B2,F`. | KABUL | Payda sıfıra giderse `f` patlar, negatifse anlamsızlaşır. `B2`'yi doğrudan `N_crit²`'ye bölmek çağrı sayısını orana sokar. Sabit komparatöre hayali `B+` eklemek yeni yöntemi kayırır. |

---

## Tur 3 — danışman geri bildirimi (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D24** | Sürekli amaç `H_r` ile yazılır: `E² = (1/T)∫‖H_r ∫Φ B Δa dτ‖² dt`. | KABUL | İç integral 6 boyutlu bir **durum** pertürbasyonu; konum RMS'i için blok seçilmeden norm alınamaz, alınırsa birim karışır. Ayrık kısımda zaten doğruydu, sürekli hâli senkron değildi. |
| **D25** | `M_j`, `A_i`, `Q_ik` blokları **6×6**. Terminal özel hâli `J_T = ‖H_r Φ(t_M,t_0) Σ u_i‖²`, `M_M = Φ(t_M,t_0)ᵀH_rᵀH_rΦ(t_M,t_0)`. | KABUL | `u_i` bir durum pertürbasyonu ve seçim `M_j`'nin içinde; "3×3 blok" yanlıştı. Taşıma matrisi düşürülemez, `u_i`'ler `t_0`'a taşınmış. |
| **D26** | **İleri prob.** Her karar sınırında, önündeki aralık boyunca, geminin kendi kısa-ufuk durum tahmini üzerinde `n_probe` noktada probe edilir. Sınırdan önceki yön bilgisi yeniden kullanılmaz. | KABUL | Eski kural (aralık içinde topla → sonraki sınırda kullan) `πR/N` argümanıyla çatışıyordu: 120 s önceki yön dekorele. İki koherans ölçeği farklı — bir karar aralığında durum tahmini `πR/N`'in çok altında, yedi günlük pilot yayda çok üstünde. |
| **D27** | Aday kümesi bir **pencere**: `𝒩_q = {N_plan(t_q) ± δ}`, ve tüm adaylar **tek paylaşılan bant yığınının** kısmi toplamlarından çıkar. Artımlı maliyet `2N(δ+k)`, `2kN` değil. **→ D50 ile revize: bu bir prob noktasının içi; noktalar arası maliyet ayrı ve baskın.** | KABUL | Çok adaylı çevrimiçi prob maliyeti sayılmamıştı; her aday için ayrı prob olsaydı ek yük %5 değil çok daha fazla olurdu. Pencere ayrıca plan/prob ayrımını somutlaştırıyor: plan kaba tahsisi, prob pencere içi düzeltmeyi yapıyor. `δ` taranır ve kısıtı yakalama oranında görünür. |
| **D28** | `B = M·B1` bağlantısı bir kez kurulur; kısıt toplam üzerinden yazılır (`Σ_i N_{g(i)}² ≤ B`), `B1` çağrı başına ortalama kalır. **→ D60 ile genelleştirildi: ana biçim `Σ_q W_q N_q² ≤ B1·T`; bu satır onun düzgün-ızgara özel hâli.** `β` yalnız `B1` düzeyinde, `β₂` ayrı — nomenklatüre de girdi. | KABUL | `B1 = ⟨N²⟩` ile `Σ N²` aynı sayı değildi. Toplamın biriktirme örnekleri üzerinden yazılması, karar aralıkları eşit uzunlukta olmadığında da doğru ağırlıklandırır. |
| **D29** | Sertifika metriği dondu: `g_J = (J_desc − L_FW)/J_desc`, `g_E = 1 − √(L_FW/J_desc)`. **Eşik hata uzayında: medyan `g_E < 0.10`.** `L_FW = max{0, en iyi sınır}`; `L_FW = 0` kalırsa sertifika boş, ayrı sütunda sayılır. | KABUL | `J = E²` olduğundan %20 J-boşluğu ile %20 E-boşluğu aynı şey değil; eşiğin hangi uzayda olduğu söylenmemişti. FW sınırı erken iterasyonda negatif çıkabilir ve kırpılmalı. |
| **D30** | Dil frenlendi: "no smooth surrogate can supply it" / "a copy of the field" gibi kategorik ifadeler kaldırıldı. Yerine **çözünürlük ve maliyet** ifadesi: bu ölçeği çözmeyen bir vekil yönü temsil edemez, çözen bir vekil karşılaştırılabilir serbestlik derecesi taşır. | KABUL | Kategorik imkânsızlık iddiası, hakemi "neden bir sinir ağı öğrenemesin?" semantik tartışmasına davet ediyor. `abl-probe` vekili zaten kuruyor ve ölçüyor; iddia ölçüme bırakılır. |

---

## Tur 4 — danışman geri bildirimi (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D31** | Prob konumlarını **mevcut durumun iki-cisim (Kepler) propagasyonu** öngörür. Hiç alan değerlendirmesi yok, analitik; maliyeti `B+`'a yazılır. `abl-predictor` ile `N_plan`'da düşük-dereceli mikro-propagasyona karşı ölçülür. | KABUL | `t_q`'da `N_q` seçilmemiş, entegratör aralığı geçmemiş — "propagatör zaten üretiyor" demek döngüseldi. Döngüyü kıran şey: iki aday derece bir aralıkta en fazla `½‖Δa‖Δt²` ≈ **mm** ayrışıyor, `πR/N` ise **km**. Öngörücünün `N_q`'yu bilmesi gerekmiyor. Kepler hatası da `½‖a_pert‖Δt²` ≈ metre, toleransın çok içinde. |
| **D32** | `δ` bir **aday-ızgara indis yarı-genişliği**, derece değil: `𝒩_q = {N_{j−δ},…,N_{j+δ}}`. Maliyeti belirleyen indüklenen derece açıklığı `Δ_span = N_{j+δ} − N_{j−δ}`; ikisi de tabloda yan yana. | KABUL | `𝒩` seyrek bir ızgara olabilir (20,30,40,…); "δ=4" dört derece mi dört indis mi belirsizdi. İndis olarak sabitlemek, ızgara değişse de denetleyici karmaşıklığını (`2δ+1` aday) sabit tutar. |
| **D33** | `λ` bisektinden önce **monotonluk kontrolü**: yoğun log ızgarada `W(λ)` süpürülür ve en-iyi-başlangıç çözümünün işi azalan mı diye bakılır. Monoton değilse bisekt bırakılır, ızgara üzerinde en küçük `J`'li olurlu nokta seçilir; bu yörüngeler sayılıp raporlanır. | KABUL | `W(λ)` monotonluğu **global** minimizör için geçerli; biz çok-başlangıçlı yerel çözüm alıyoruz ve basin değişebilir. Kör bisekt sessizce yanlış nokta verebilir. Kontrol çok ucuz. |
| **D34** | Abstract'taki süpürme ifadesi: "`K_dec` karar aralığı üzerinde, `M` biriktirme örneği kullanılarak, `O(M·|𝒩|)`". | KABUL | D16 ile tam tutarlılık; "M epok üzerinde süpürme" eski formülasyonun kalıntısıydı. |
| **D35** | Kalan sert dil temizlendi: §7.1 başlığı "cannot be tabulated" → "tabulating it is not cheap"; §9 "copy of the field" → ölçülen serbestlik-derecesi karşılaştırması; intro amber-branch `dnote`'unda "impossibility" → "measured cost statement". | KABUL | D30'un manuskript dallarına ve başlıklara da uygulanması. |
| **D36** | README'deki amaç fonksiyonu `H_r` ile senkronlandı. | KABUL | "Bir sembol bir yerde tanımlanır" kuralı; D24'ten sonra README bayat kalmıştı. |

---

## Tur 5 — iç denetim (2026-08-08)

Danışman geri bildirimi değil; sistemin baştan sona kendi kendine denetimi.
On üç sorun bulundu, dördü ciddi.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D37** | **`Ê` (öngörülen, `√J`) ile `E` (propagate) ayrı semboller.** M1–M2'nin bütün nicelikleri `Ê`; `E` yalnızca M3+. Oranlar `ρ̂` / `ρ`. Her hipotez satırında düzey yazılı. | KABUL | **Ciddi.** M1'de hiç propagasyon yok, ama G1/H1a/H2 "hata" diyordu ve NOTATION `E`'yi propagate edilmiş RMS olarak tanımlıyordu. Aynı sembol iki farklı şeyi gösteriyordu — önceki makalenin ölçüm-konvansiyonu tuzaklarının birebir muadili. |
| **D38** | `f`'in üç ayrıntısı bağlandı: (1) `A-sign` uçurulabilir olmadığı için ek yüksüz tam `B_tot` alır, `C-plan` `B_tot − B+` alır — asimetri kasıtlı ve ölçülen şeyin kendisi; (2) `f > 1` mümkündür ve hata değildir; (3) pay ve payda ikisi de `E`, `f` M3'ten önce hesaplanamaz. | KABUL | Bütçe asimetrisi hiçbir yerde yazılı değildi; hakem "aynı bütçede mi karşılaştırıyorsunuz?" diye sorar. `f > 1` de mümkün: `A-sign` referans yayda birinci mertebe optimum, propagate edilince gerçek yayda daha kötü olabilir. |
| **D39** | WP0'a **geri düşüş** eklendi: birebir eşleşme olmazsa teşhis → ilan edilmiş tolerans → tek kod yolu kuralı; teşhis edilemezse dur. | KABUL | Eski kapı "eşleşme yoksa hiçbir sayı üretilmez" diyordu, yani projeyi öldüren ve alternatifi olmayan bir kapıydı. Bir plan ne yapılacağını söylemeli. |
| **D40** | `A-sign+` sabit-nokta testinin indirgenmiş hâli (16 yörünge) **M5'ten M3'e** alındı (WP17a). | KABUL | Bu test makalenin **merkezî nesnesinin ne olduğunu** belirliyor: `A-sign` bir yörünge-seviyesi kıyas mı, referans-yay optimumu mu? Sonradan öğrenilirse M5'in ~1800 yayı yanlış adlandırılmış bir nesne üzerinde harcanmış olur. |
| **D41** | G6'ya **kaçış** eklendi: iki teşhis turundan sonra hâlâ işaret uyuşmuyorsa eleyici bu sınıf için kalibre değil ilan edilir, M2'nin elemesi iptal, kampanya tek adayla propagasyona dayanır. | KABUL | Eski metin "→ M2'ye dön" diyordu; sonsuz döngü riski vardı. |
| **D42** | **Yeni aşama M4: bütçe probu.** Tasarım A'nın 16 yörüngelik tabakalı alt kümesi, β ∈ {0.50, 0.75, 1.25, 1.50}, beş politika (~320 yay). **Kapı G7:** komşu bütçelerde de kazanıyor mu. Geçilmezse M5 popülasyon genişletmesi iptal. Aşamalar M4–M8 → M5–M9 kaydı. | KABUL | **Ciddi.** Eski sıra β=1'de kazanır kazanmaz üç popülasyona ~1800 yay harcıyordu. Önceki makalede iç üyenin bütün yapıcı sonucu β=0.5'te ters dönmüştü — aynı hatanın tekrarı olurdu. Prob 320 yaya mal oluyor ve 1500 yayı riske atmaktan koruyor. |
| **D43** | **H4b** eklendi: `C-plan`, `R-int`'i de yener. | KABUL | OUTCOMES'un 🟢 bandı `R-int`'in yenilmesini şart koşuyordu ama bunu test eden tescilli hipotez yoktu. Bir bant sınırı, ön-tescilde karşılığı olmayan bir iddiaya dayanamaz. |
| **D44** | Depolama çözümü WP1'in **teslimi**: float32 + bellek-eşlemeli sıralı erişim. Önek/sonek toplamları `M × 6` (~3 MB) bellekte, `u_i(N)` tablosu her süpürmede diskten sırayla. | KABUL | "Sıkıştırma veya akış gerekir" yazıyordu ama hiçbir WP sahiplenmiyordu ve çözüm belirsizdi. Koordinat inişinin erişim deseni sıralı olduğu için memmap yeterli; tabloyu belleğe almak gerekmiyor. |
| **D45** | **G5b** eklendi: M2'nin gerçek durdurma kriteri, en iyi adayın öngörülen yakalama oranı `f̂ ≥ 0.15`. | KABUL | G5 bir usul kuralıydı (kural önceden yazılır), pas/kal kriteri değildi. Bütün adaylar kötüyse M3'e 430 yay boşuna harcanacaktı. Eşik H3'ün 0.33'ünün bilerek altında. |
| **D46** | M2 maliyeti düzeltildi: "~26 varyasyonel çözüm" → **~240 propagasyon eşdeğeri** (6 aday × 26 + 10 ablasyon × 8). Ablasyonlar küçük panelde koşulur, panel büyüklüğü her tabloda yazılır. | KABUL | Altı kat hata. Her varyasyonel çözüm bir propagasyona mal oluyor; aday ve ablasyon sayısıyla çarpılmamıştı. |
| **D47** | M3 yay sayısı düzeltildi: "~200" → **~430** (5 politika × 64 + eşit-iş komparatörleri + pilot yaylar + WP17a). | KABUL | Beş politikaya geçilince (D17) sayı güncellenmemişti. |
| **D48** | Pilot yay üretimi ayrı maliyet kalemi: `C-plan` propagate edilen **her** yörünge için bir pilot yay ister (M3'te 64, M5'te ~370), her biri %8 → ~5 ve ~30 yay eşdeğeri. | KABUL | WP6 yalnızca 10 pilot yayı tutarlılık ölçümü için sayıyordu; denetleyicinin kendi ihtiyacı hiçbir yerde yoktu. |
| **D49** | Küçük düzeltmeler: `T_coh` eşiği zamanla değişir (`πR/N_plan(t)`, `N` sabit değil); `κ = cos∠(v̂, Δa_ref)` ("gerçek" değil, terminoloji kuralı); iki-koherans tablosunda "metre altı" → "metre mertebesi" (D31 ile tutarlı); `K_ref` = **`F-op`'un** çağrı sayısı (hangi politika olduğu belirsizdi); `κ ≥ 0.7` eşiğine gerekçe yazıldı (kosinüs, sadeleşme teriminin korunan kesridir); NOTATION §6'da kopmuş madde işaretleri onarıldı. | KABUL | — |

---

## Tur 6 — hesaplamalı astrodinamik denetimi (2026-08-08)

Fizik ve sayısal taraf, kendi iddialarım yeniden hesaplanarak. **Dört ciddi
hata**, biri maliyet modelini kökten değiştiriyor.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D50** | **Prob maliyeti yeniden modellendi. `2k/N` YANLIŞ.** İleri prob yeni bir konumda alınıyor; orada ALF rekürsiyonu sıfırdan `N+k`'ya kadar koşmak zorunda. `2k/N` yalnızca prob, zaten değerlendirme yapılan bir noktada olsaydı geçerliydi — ileri prob tanımı gereği öyle değil. Doğru model: **her prob noktası ≈ bir tam sentez.** Aralığı tam kapsayan prob için maliyet tabanı `n_probe/n_RHS = 1/(τ_corr · RHS_hızı)` ve **`Δt_dec`'den bağımsızdır**: N=120'de **~%21**, N=300'de **~%53**. **→ D67 ile düzeltildi: bunlar yerel değerler; yay integrali ~%12–19.** | KABUL | **En ciddi bulgu.** Tek bir bant için bile `P_{n,m}` rekürsiyonu tüm alt derecelerden geçmek zorunda; tek bir dereceyi `O(N)` işlemle almanın standart bir yolu yok. `B+` %5 değil, ~%20–50. |
| **D51** | Bunun sonucu bir **tasarım eğrisi**, tek bir sayı değil. WP5 artık `κ_eff(n_probe)` ölçüyor: aralık boyunca ortalanmış etkin yön isabeti, prob noktası sayısının fonksiyonu. Makale bir **maliyet–isabet ödünleşim eğrisi** raporluyor, bir yüzde değil. Ayrıca ucuz varyant `C-lite`: aralığın ilk RHS çağrısı yükseltilmiş derecede yapılır (eş-konumlu, ~%2), yön yalnızca aralık başında bilinir. | KABUL | Sorun kaçınılmaz olduğuna göre ölçülmeli. İki uç: tam kapsama (~%12–19, D67) veya tek eş-konumlu prob (~%2, kısmi yön). Aradaki eğri makalenin mühendislik katkısı. |
| **D52** | **Plan `zaman` ile değil `faz` ile indekslenir.** `c(t)`, `w(t)`, `λ`, `N_plan` → `(devir indeksi, devir içi faz)`. | KABUL | **Ciddi.** N=40 gevşek pilot yay 7 günde ~100 km along-track sapıyor = **~60–100 s zaman kayması**. Perilun geçişi birkaç dakika sürüyor — aynı mertebe. Zaman-indeksli bir plan, `c` ve `w`'yi perilun civarında yanlış epoklara hizalar; oysa bunlar tam da perilunda baskın. Faz indeksi bu duyarlılığı tamamen kaldırır (devirler arası 100 s = 0.014 devir, ihmal edilebilir). Yarıçapın tek değerli fonksiyonu olmaya da geri dönmüyoruz: indeks faz **artı** devir numarası. |
| **D53** | **`γ` genlik tamamlaması, uydurulmuş güç yasasından değil, alanın gerçek derece-varyans spektrumundan** (`P_n`, 1D tablo) alınır. | KABUL | Önceki makale ay spektrumunun tek bir güç yasası **olmadığını** ölçtü (`p_spec=2.13` vs `p_fit=1.76` farkının sebebi tam olarak bu). Güç yasasıyla kuyruk ekstrapolasyonu sistematik hata taşır. `P_n` yalnızca ~1800 sayı — yüklemek uçurulabilirliği hiç zedelemez. Ayrıca makalenin bilgi argümanını **keskinleştiriyor**: taşınması ucuz olan **spektrum (1D)**, taşınması pahalı olan **doku (2D alan)**. |
| **D54** | **Döngüsellik uyarısı: eleyici ile amaç aynı fonksiyonel.** `A-sign` `J`'yi minimize ediyor, varyasyonel eleyici de `J`'yi hesaplıyor. Önceki makalenin 100/100 kalibrasyonu, `J`'ye göre **optimize edilmemiş** politikalar içindi. Bu yüzden **G6 artık taşıyıcı bir kapı**, teyit değil. | KABUL | `J`'ye göre optimize edilmiş bir politikayı `J` ile skorlamak, `J`'nin `E`'yi öngördüğü varsayımını test etmez. G6'nın rolü yükseldi ve OUTCOMES'a sahte-iyi satırı olarak girdi. |
| **D55** | `Δt_acc` **uyarlanabilir**: `τ_corr = πr/(N v)` yörünge boyunca ~50 kat değişiyor (perilunda 11–28 s, apolunda ~1500 s). Perilun-inceltilmiş ızgara kullanılır. `B1` o zaman **zaman-ağırlıklı** ortalamadır: `B1 = Σ N_i² Δt_i / T`. | KABUL | Düzgün 10 s'lik ızgara 7 günlük yayda 60481 örnek üretiyor ve çoğu hiçbir şeyin olmadığı apolunda. Uyarlanabilir ızgara depolamayı ~10 kat düşürüyor (D44'ün sorununu büyük ölçüde çözüyor) ve fizik olarak doğru olan da bu. |
| **D56** | Önek toplamlarında **telafili toplama** (Kahan/Neumaier) kullanılır. | KABUL | `S_j` 60 bin işaretli terimin toplamı ve **sadeleşme yöntemin ta kendisi** — koşul sayısı yüksek. Çift duyarlıkta bile telafisiz toplama gereksiz risk; önceki makale zaten `kahan1965pracniques` atıflı. |
| **D57** | `A-sign`'ın **referans dereceye doygunlaşma eğilimi** ayrıca denetlenir: optimizatör, bütçenin izin verdiği yerde dereceyi referansa itip defekti tanım gereği sıfırlamaya *teşvikli*; sabit bir politikanın böyle bir teşviki yok. Tavan teması kesri her `A-sign` çizelgesi için raporlanır, yüksek temaslı yörüngeler orandan çıkarılır. | KABUL | Önceki makalenin tavan denetimi sayım yapıyordu; burada optimizatörün eğilimi daha güçlü ve ayrı bir kontrol gerektiriyor. |
| **D58** | Anahtarlama maliyeti **ölçülmüş veriden** kalibre edildi, tahminden değil: önceki kampanyada radyal politika 10 km binlerle ~16 000 anahtarlama yapıp yalnızca **%7 fazla RHS çağrısı** harcadı → anahtarlama başına ~0.45 çağrı. `Δt_dec = 120 s` ile 5040 anahtarlama → **~%2**. Yani anahtarlama ucuz; pahalı olan prob. | KABUL | Kaba tahminim (anahtarlama başına ~5–12 çağrı) gerçek ölçümün 10–25 katıydı. Doğru sayı önceki kampanyanın kendi verisinden çıkıyor. |
| **D59** | Küçük düzeltmeler: aday ayrışması `½‖Δa‖Δt²` **mm değil ~cm** (2–7 cm, `‖Δa‖`'ya göre) — sonuç değişmiyor, km ölçeğinin hâlâ 5–6 mertebe altında; dekorelasyon `πr/N` (uydu yarıçapı), `πR/N` değil; varyasyonel eleme **tek çözümde çok kanallı** yapılır (rev42 zaten öyle: tüm politikalar aynı referans ve `G`'yi paylaşır, yalnızca `Δa` farklı) — 240 ayrı çözüm olarak uygulanmamalı. | KABUL | — |

---

## Tur 7 — son düzeltmeler, sonra DONDURULDU (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D60** | **Türetim baştan ağırlıklı yazıldı.** `u_i = Φ B Δa_i Δt_i`; `J = (1/T) Σ_j w_j S_jᵀ M_j S_j`; `A_i = Σ_{j≥i} w_j M_j`; `Q_ik = A_max(i,k)/T`; ceza `λ W_q N²` (`W_q = Σ_{i∈I_q} Δt_i`); kısıt `Σ_q W_q N_q² ≤ B = B₁T`; FW gevşetmesinde de `W_q`. | KABUL | **Majör tutarsızlık.** D55 uyarlanabilir ızgarayı getirdi ama ana türetim düzgün ızgara varsayımında kalmıştı; `|I_q|` örnek sayısı ağırlık yerine kullanılıyordu — bu, perilunu yalnızca daha yoğun örneklendiği için bütçede fazla ağırlıklandırırdı. Önek/sonek yapısı bozulmuyor: ağırlıklar yalnızca sonek dizisinin içine giriyor. Düzgün ızgarada eski forma indirgeniyor. |
| **D61** | **`float64` saklanır; `float32` yalnızca zorunlu parity testiyle opsiyon.** | KABUL | Uyarlanabilir ızgara depolamayı ~87 → ~9 MB/yörünge'ye indirdiği için float32'nin gerekçesi kalmadı. Kahan toplama, float32'ye **yazarken** kaybedilen biti geri getirmez — ikisi ayrı sorun. Makalenin merkezindeki işaretli sadeleşme sonucunu depolama hassasiyetine teslim etmek gereksiz risk. |
| **D62** | **RQ3 dili düzeltildi.** "reference field unavailable" değil; **"without the reference trajectory or pointwise reference-field evaluations"**. Ayrım: `P_n` = küresel/çevrimdışı model üst-verisi (açık); `Δa_ref(t)` = yörüngeye özgü (yasak). | KABUL | Denetleyici artık `P_n` spektrumunu kullanıyor (D53) ve bir propagatör katsayıları zaten taşır. Eski ifade literal olarak yanlıştı. Yeni ayrım hem doğru hem de makalenin bilgi iddiasını netleştiriyor. |
| **D63** | **`τ_dec` → `τ_corr`.** | KABUL | `Δt_dec` (denetleyici karar kadansı) ile `τ_dec` (alan dokusu koherans süresi) kesin karışırdı. |
| **D64** | **G7 bant kuralı tek yerde donduruldu.** GEÇER = hem β=0.75 hem β=1.25'te çözünen çoğunluk → M5 koşulur. KALIR = yalnız biri veya hiçbiri → M5 iptal, **🟡**. GEÇTİ ama M5 kaynak nedeniyle koşulmadı → **🟩'nin alt ucu**. | KABUL | ROADMAP G7 düşmesini 🟩, OUTCOMES 🟡 sayıyordu. Ön-tescilde tek karar olmalı; üçü de artık aynı kuralı yansıtıyor. |
| **D65** | Senkronizasyon: PREREG'de G1 → `ρ̂`, H2 → `Ê`; abstract'ta "calibrated Kaula tail" → ölçülen spektrum, "sweep over M epochs" → `K_dec` karar aralığı; `πR/N` → `πr/N` (uydu yarıçapı) her yerde. | KABUL | — |
| **D66** | **KAPSAM DONDURULDU.** Bu turdan sonra yeni kontrol, yeni ablasyon, yeni kapı, yeni WP **eklenmez**. Bulunan her yeni fikir doğrudan PLAN §7 future work'e yazılır. | KABUL | Sistem şu an 7 kapı, M0–M9, ~270 varyasyonel çözüm, ~3700 propagate yay *(→ D98 ile ~7150; kalibrasyon iterasyonu sayılmamıştı)*, 14 ablasyon, iki denetleyici varyantı içeriyor. Her hakem ihtimalini gönderim öncesinde çözmeye çalışmak makaleyi bir teze çevirir. Bundan sonrası ölçüm. |

---

## Tur 8 — dondurma sonrası kusur taraması (2026-08-08)

Yeni kapsam yok (D66 geçerli). Yalnızca mevcut sistemdeki **hatalar**.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D67** | **Prob ek yükü yay integrali olarak yeniden hesaplandı: `~%12–19`, `%21–53` değil.** Formül `[∫dt/τ_corr(t)] / N_RHS,toplam`. | KABUL | **Kendi hatam.** *Yerel* `τ_corr`'u (perilun, yüksek `N`) *küresel* `r_RHS` ile çarpmışım — kategori hatası. `τ_corr`'un en kısa olduğu yerde uydu en az zaman geçiriyor ve orada `r_RHS` de yüksek; iki etki zıt yönde. Arşiv yaylarından hesap: dairesele yakın %19, perilun 50/apolun 1000 %14, perilun 30/apolun 2500 %12 — eksantrik yörüngeler uzun apolun kolları neredeyse hiç prob istemediği için **daha ucuz**. Nitel sonuç değişmiyor: prob hâlâ `B+`'ın baskın kalemi. |
| **D68** | **Uyarlanabilir ızgaranın yakınsama parametresi tanımlandı:** `n_s = τ_corr(t)/Δt_i ∈ {0.5,1,2,4,8}`. Mutlak `Δt_acc` taraması yalnızca düzgün ızgarada ve önceki makaleyle karşılaştırılabilirlik için. | KABUL | D55 ızgarayı uyarlanabilir yaptı ama WP4 hâlâ **mutlak** `Δt_acc` tarıyordu — uyarlanabilir ızgarada mutlak adımın tek bir değeri yok, dolayısıyla G3'ün ölçütü tanımsız kalmıştı. |
| **D69** | Anahtarlama başına maliyet (~0.5 çağrı) **kestirim üstüne kestirim** olarak etiketlendi; ölçülen sayı %7'nin kendisi. Anahtarlama sayısı arşivdeki derece geçmişlerinden **sayılır**, tahmin edilmez (WP16). | KABUL | %7'yi kendi tahmin ettiğim anahtarlama sayısına bölüp "ölçülmüş 0.45" gibi sunmuştum. |
| **D70** | `C-lite` ön-tescilli aday listesine eklendi (OA-02). | KABUL | D51'de ilan edilmiş, NOTATION/PLAN/manuskriptte var, ama eleme listesinde yoktu — maliyet–isabet eğrisinin ucuz ucu hiç propagasyon şansı bulamayacaktı. |
| **D71** | Ablasyon adları hizalandı: `abl-bin` → `abl-timeindex`; `abl-predictor`, `abl-window`, `abl-lookahead` manuskript tablosuna eklendi. | KABUL | Dokümanlarda ilan edilen dört ablasyon manuskript tablosunda yoktu; `abl-bin`/`abl-timeindex` aynı şeyin iki adıydı. |
| **D72** | Nomenklatüre eklendi: `Δt_i`, `w_j`, `W_q`, `T`, `Δ_span`, `n_probe`, `r_RHS`, `N_RHS`, `P_n`, `L_FW`, `κ_eff`, `τ_corr`. §3'teki "`O(M)` süpürme" ifadesi düzeltildi: çapraz terimler `O(M)`, süpürme `O(M·\|𝒩\|)` — ikisi ayrı büyüklük. | KABUL | `W_q` matematikte 10, `w_j` 6 kez geçiyordu ama tanım listesinde yoktu. Ağırlıklı türetime geçişin artığı. |

| **D73** | Bütçe kısıtı tek ana biçime indirgendi: `Σ_q W_q N_q² ≤ B = B1·T`; `Σ_i N² ≤ M·B1` bunun düzgün-ızgara özel hâli olarak etiketlendi. | KABUL | NOTATION ve PREREG iki biçimi eşdeğermiş gibi yan yana listeliyordu; uyarlanabilir ızgarada yalnız biri doğru. |

---

## Tur 9 — denklem denklem doğrulama (2026-08-08)

Her denklemi kendim türeterek. Eq. (5)–(11) **doğru çıktı**; dört başka kusur.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D74** | **`κ` eşiğinin gerekçesi yanlıştı, düzeltildi.** "`κ` doğrudan sadeleşme teriminin korunan kesridir" demiştim. Değil: terim `2⟨u_i,c_i⟩` ve `u_i = ΦBΔaΔt`, dolayısıyla korunan kesir `⟨v̂,z⟩/⟨Δa,z⟩` (`z = (ΦB)ᵀc`) — `v̂` ile `Δa` arasındaki açıya değil, **ikisinin `z` ile yaptığı açılara** bağlı. `ΦB` açıları korumaz. `Δa ⊥ z` ise gerçek katkı sıfır ve prob sahte katkı üretebilir. `κ` artık **önde gelen gösterge** olarak sunuluyor; korunan kesir WP5'te yanında raporlanıyor (yeni kapı değil, aynı ölçümün ikinci çıktısı). | KABUL | Matematiksel olarak yanlış bir gerekçeyi bir kapıya dayamak, hakemin ilk yakalayacağı şeylerden. |
| **D75** | **`A-sens` türetilmiş tanıma geçti: `diag Q`.** Eskiden `w(t) = ‖Φ(T,t)B‖²` ödünç alınmış skalerdi — hem keyfî bir ithal hem de **terminal** bir proxy'yi **yay-RMS** amacına takmak. Köşegen dışı bloklar silinince geriye tam olarak `Q_ii = A_i/T` kalıyor. | KABUL | **Rigor kazancı.** Merdiven artık tek bir matrisin ayrıştırması: `A-force` = `Q` yok, `A-sens` = `diag Q`, `A-sign` = tüm `Q`. RQ2 tek soruya iniyor: *kazancın ne kadarı `Q`'nun köşegeni dışında yaşıyor?* Ek hesap maliyeti sıfır — `A_i` zaten hesaplanıyor. |
| **D76** | Kuadratür ağırlığı `w_j → ω_j` yeniden adlandırıldı ve **`ω_j ≥ 0`** şartı yazıldı. | KABUL | `w` hem kuadratür ağırlığı hem duyarlılık ağırlığı olarak kullanılıyordu — aynı bölümde, iki anlamda. Ayrıca `Q ⪰ 0` (dışbükey gevşetmenin dayandığı şey) `ω_j ≥ 0` gerektiriyor: dikdörtgen/yamuk sağlar, yüksek dereceli Newton–Cotes **sağlamaz**. Bu şart hiçbir yerde yazılı değildi. |
| **D77** | Algoritma 1 listesi düzeltildi: çok-başlangıç döngüsü **bisektin içine** açıkça yazıldı (eski liste "bir başlangıçtan başlat" diyordu, tek başlangıç kodlanırdı) ve dönüş değeri eklendi. Sıcak başlatma **kasten** kullanılmıyor: bisekt yoluna bağımlılık yaratırdı — D33'ün monotonluk kontrolünün yakalamaya çalıştığı şey tam olarak bu. Ayrıca Eq. (10)'a eksik `1/T` çarpanı kondu. | KABUL | Metin çok-başlangıç diyordu, sözde-kod tek başlangıç. Bir uygulayıcı sözde-kodu kodlar. |
| **D78** | Sessizce başarısız olmuş üç yama yakalandı ve elle düzeltildi: merdiven tablosu başlığı hâlâ `E`/`ρ` kullanıyordu (`Ê`/`ρ̂` olmalı), `G_sens`/`G_sign` denklemi `E` kullanıyordu, G1 `dnote`'u `A-sign/R-int` kısaltmasını kullanıyordu. | KABUL | Betikle yapılan `replace` çağrıları `assert` olmadan sessizce geçmişti. **Bundan sonra LaTeX'e yapılan her betik yaması `assert` ile korunuyor.** |

---

## Tur 10 — okunmamış bölümlerin taranması (2026-08-08)

Daha önce uçtan uca okumadığım §2, §5, §8, §10 ve doküman maliyet sayıları.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D79** | **M2 maliyeti yine yanlıştı: ~240 → ~270.** Ablasyon sayısı D71/D74 ile 10'dan **14**'e çıkmış ama maliyet güncellenmemişti (6 aday × 26 + 14 × 8 = 268). | KABUL | Sayı üç dosyada birden bayattı. Kapsam büyüdükçe maliyet satırının otomatik güncellenmemesi tekrar eden bir kusur — WP19 gönderim sayaçlarına "ablasyon sayısı × panel = M2 maliyeti" tutarlılığı eklendi. |
| **D80** | **§8, var olmayan bir maddeye atıf yapıyordu:** "Section~\ref{sec:setup}'un kaçış maddesi". §5'te böyle bir madde yoktu — G6 kaçışı yalnızca ROADMAP/PREREG'deydi. §5'e **"Ön-tescil ve durdurma kuralları"** alt bölümü yazıldı; bütçe kapısı (G7) ve eleme kaçışı (G6) artık manuskriptte. | KABUL | LaTeX referansı geçerliydi (etiket vardı) ama işaret ettiği içerik yoktu — mekanik tarama bunu yakalayamaz. |
| **D81** | **`B₁`'in iki farklı tanımı vardı:** §5 "referans çıktı ızgarasında `⟨N²⟩`", §3.5 "zaman-ağırlıklı `Σ N²Δt/T`". Tek tanıma indirgendi: **zaman ortalaması**; düzgün çıktı ızgarasında düz ortalamaya indirgenir (önceki makalenin raporladığı biçim). | KABUL | Uyarlanabilir ızgaraya geçince ikisi ayrıştı ve §5 güncellenmemişti. |
| **D82** | §10'daki dört kusur: "vekil **çalışamaz**" (D30 ile yumuşatılmıştı, burada kalmış); "**üç** tasarım seçimi" (gerçekte **beş**); "**çağrı başına** prob" (aralık başına, ve D50'den sonra maliyet modeli değişti); Limits listesinde **prob ek yükü** ve **eleyici döngüselliği** yoktu. Dördü de düzeltildi. | KABUL | §10 son üç turun hiçbirinde okunmamıştı. |
| **D83** | §2'de **ayrıştırılamayan tahsis literatürüne konumlanma yoktu.** Ayrıştırılabilir hâl (Everett, Shoham) atıflı, ama makalenin çekirdeği tek kısıtlı bir tamsayı kuadratik program ve bu bir literatür sınıfı. Konumlandırma paragrafı eklendi, atıf `\ph{}` ile literatür turuna bırakıldı (uydurma kaynak yazılmadı). | KABUL | Hakem "bu bilinen bir problem sınıfı değil mi?" diye sorar; sormadan önce söylemek gerek. |
| **D84** | §10 ve §11'e türetilmiş merdiven dili taşındı: "üç kriter tek bir matrisin okumaları; ağırlık zaten köşegende, eksik olan hiçbir zaman ağırlık değildi". | KABUL | D75 §3/§6'da uygulanmıştı ama tartışma ve sonuç eski çerçeveyi anlatıyordu. |

---

## Tur 11 — kalan bölümler ve çapraz atıf denetimi (2026-08-08)

§9, §11, `main.tex` ve PREREG uçtan uca okundu; ayrıca her `D<n>` çapraz
atıfının işaret ettiği karara gerçekten karşılık gelip gelmediği denetlendi.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D85** | **Dört yanlış çapraz atıf düzeltildi: G7 bant kuralı için `D62` yazılmış, doğrusu `D64`.** (`D62` = RQ3 dili.) PREREG'de bir, OUTCOMES'ta iki, ROADMAP'te bir. | KABUL | Karar günlüğünün bütün değeri atıfların doğru olmasında. Otomatik denetim eklendi: her `D<n>` atıfı, işaret ettiği kararın konusuyla birlikte listeleniyor. |
| **D86** | **Abstract'ta iki bayat yer.** "sweep over $M$ epochs" hâlâ duruyordu (D34'ün yaması daha önce bir `assert` hatasına takılıp yazılmamış); "texture that cannot be tabulated at all" ise D30/D62'nin temizlediği sert dilin kalıntısıydı. İkisi de düzeltildi; ayrıca merdiven `diag Q` / köşegen-dışı diliyle yeniden yazıldı. | KABUL | Abstract makalenin en çok okunan kısmı ve üç turdur bayattı. |
| **D87** | **§9 "hepsi aynı panelde" diyordu** — oysa adaylar 26, ablasyonlar 8 yörüngelik panelde (D79). Panel büyüklükleri ve ablasyon satırlarının **merdiven tablosuyla karşılaştırılamaz** olduğu açıkça yazıldı. Ayrıca `abl-kaula`'nın tarifi tersti: artık birincil `P_n` spektrumu, ablasyon güç yasası (D53). | KABUL | "Aynı panel" ifadesi hem yanlıştı hem de M2 maliyet modeliyle çelişiyordu. |
| **D88** | **§11 "ne klasik çarpan yöntemi ne de duyarlılık-ağırlıklı onarımı *uygulanır*" diyordu.** Yanlış: ikisi de gayet uygulanır (ayrıştırılabilirler), yalnızca **yanlış fonksiyoneli** optimize ederler. D75'in dili sonuca da taşındı. Sertifikalı boşluk için `\ph{}` eklendi. | KABUL | D75 §3/§6'da uygulanmıştı, §10 ve §11 eski çerçeveyi anlatmaya devam ediyordu. |
| **D89** | **PREREG'de altı tutarsızlık:** G3 hâlâ `Δt_acc` diyordu (D68'den sonra ölçüt `n_s`); aşama numaraları eskiydi (M4 = popülasyon genişletme yazıyordu, oysa D42 ile M4 = bütçe probu, M5 = genişletme); bilgi-kaynağı denetimi "referans alan görünürse geçersiz" diyordu (D62 ile `P_n` serbest); çizelge biçimi kontrolü zaman-indeksli diyordu (D52 ile faz-indeksli); G4'te `κ`'nın gösterge olduğu yazılı değildi (D74); ablasyon paneli büyüklüğü ilan edilmemişti (D79). | KABUL | PREREG dondurulacak dosya; içindeki her tutarsızlık koşu sırasında karar boşluğu demek. |
| **D90** | Aşama numarasına bağlı yazılmış cümleler ("M7'de değil") numaradan bağımsız hâle getirildi ("geç aşamada değil"). | KABUL | D42'nin yeniden numaralandırması bu cümleleri sessizce yanlış yapmıştı; bir daha olmaması için numara referansı kaldırıldı. |

---

## Tur 12 — bütüncül okuma (2026-08-08)

Sistem bir bütün olarak okundu: tutarsızlık değil, **eksik tanım** arandı.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D91** | **Eşleştirme denklemi hangi tarafın ayarlandığını söylemiyordu.** Üç parça eklendi: (1) **çapa** — `B_tot` her yörünge ve `β` için `F-op(β)`'nın *gerçekleşen* işi, β=1'de arşivden; (2) **kim ayarlanır** — komparatör değil, **her aday çapaya kalibre edilir**; (3) **nasıl ve ne pahasına** — gerçekleşen iş ancak propagasyondan sonra bilindiğinden aday başına **2–3 propagasyon**, %2 bandına. | KABUL | **En büyük eksik.** Önceki makale komparatörü adaya göre kaydırıyordu; burada bu **yapılamaz**, çünkü `f = (E_F-op − E_C-plan)/(E_F-op − E_A-sign)` payda ve paydada **aynı** `E_F-op`'u gerektiriyor. Komparatör kayarsa iki farklı `F-op` örneği çıkar ve `f` tanımsızlaşır. Denklem üç turdur "simetrik görünen ama eksik" hâldeydi ve kimse hangi tarafın hareket ettiğini soramamıştı. |
| **D92** | **M3 yay sayısı düzeltildi: ~430 → ~670.** Kalibrasyon iterasyonu hiç sayılmamıştı. 4 kalibre aday × 64 × ~2.5 ≈ 640, artı pilot yaylar ve WP17a. | KABUL | D91'in doğrudan sonucu. Bu, "~200 → ~430 → ~670" olarak üçüncü düzeltme; her seferinde sebep aynı: kapsam büyüyünce maliyet satırı elle güncellenmiyor. WP19'un tutarlılık kontrolü artık bunu da kapsıyor. |
| **D93** | PLAN'de dört bayat yer: `πR/N` → `πr/N`; §1.3 hâlâ tek/düzgün ızgara anlatıyordu (D55 uyarlanabilir yaptı); merdiven tablosunda `C-lite` yok ve "Kaula" yazıyordu (D53 ile spektrum); FW bölümünde L1'in amacı hâlâ `w_q d_q(N)` (D75 ile `u_qᵀQ_qq u_q`). | KABUL | PLAN bilimsel içeriğin kanonik dosyası; içindeki bayat tanım koşuya yanlış girer. |

---

## Tur 13 — otomatik eskimis-terim taramasi (2026-08-08)

Elle okumak yerine, eskiyen her terim icin bir desen + izin-verilen-baglam
yazilip butun dosyalara uygulandi (`check_stale.py`). 75 supheli eslesmenin
cogu DECISIONS'in tarihsel kayitlariydi; **on bir gercek eskimis yer** cikti.

| # | Karar | Durum | Gerekce |
|---|---|---|---|
| **D94** | **`check_stale.py` eklendi ve gonderim kapisi oldu.** On dort eskimis terim deseni (`πR/N`, `w_j`, `2k/N`, Kaula-genlik, odunc `w(t)`, zaman-indeksli plan, `τ_dec`, float32, ~240, ~430, 10-ablasyon, eski eslestirme, asama-numarasina bagli cumleler, sert dil) + izin-verilen-baglam. `DECISIONS.md` bilerek haric — karar gunlugu eski ifadeyi korur. | KABUL | On iki turdur elle bulduklarimin cogu ayni birkac kategoriye giriyordu. Kategori basina bir desen yazmak, her turda yeniden okumaktan hem hizli hem guvenilir. Cikis kodu 1 ise gonderim durur. |
| **D95** | Tarayicinin buldugu **on bir gercek eskimis yer** duzeltildi: PLAN'de `w_j`→`ω_j` (4 yer) ve `ω_j ≥ 0` sarti; PLAN'de genlik hala "Kaula kuyrugu" (D53 ile spektrum) ve `γ(N,h)`→`γ(N,r)`; PLAN/WP/NOTATION/OUTCOMES'ta 18 adet `πR/N`→`πr/N`; ROADMAP'te depolama hala float32 (D61 ile float64); WP9 ciktisi "zaman-indeksli cizelge" (D52 ile faz-indeksli); WP12 kontrolu zaman-indeksli (faz-indeksli olmali). | KABUL | Hepsi bir kararla degistirilmis ama tek dosyada guncellenmis terimlerdi. |
| **D96** | **Iki yanlis karar-atfi daha:** NOTATION'da "Bilgi siniri (D61)" — dogrusu **D62** (D61 float64); WP'de "ablasyon sayisi D71/D74 ile 14'e cikti" — **D74** ablasyon eklememisti (`κ` gerekcesiydi), dogrusu yalniz **D71**. | KABUL | D85'te dort yanlis atif bulunmustu; bu ikisi ayni kategorinin kalintisi. Atif konusu esleme kontrolu artik denetimin parcasi. |
| **D97** | Plan artik `w` skaleri tasimiyor: denetleyicinin tasidigi agirlik `Q_ii` (D75). NOTATION/WP/OUTCOMES'taki `w(t)` kalintilari temizlendi. | KABUL | D75 §3/§6'da uygulanmisti; uc dokumanda skaler `w` kalmisti. |

---

## Tur 14 — sayısal çapraz denetim (2026-08-08)

Eskimiş *terim* taraması temizdi, bu yüzden eskimiş *sayı* taraması yazıldı:
her kritik eşik ve sayım dosyalar arasında karşılaştırıldı.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D98** | **D91'in kalibrasyon maliyeti yalnızca M3'e uygulanmıştı.** M4, M5 ve M6 hâlâ "politika başına tek propagasyon" aritmetiğindeydi, ve **β≠1'de çapanın kendisinin propagate edilmesi gerektiği hiç sayılmamıştı** (`F-op(β)` yalnızca β=1'de arşivde). Maliyet formülü yazıldı ve üç aşamaya tutarlı uygulandı: M4 320→**576**, M5 1200–1800→**3980**, M6 960→**1730**. Kampanya toplamı ~3700 → **~7150 yay**. | KABUL | **Ciddi.** Kaynak planı iki kat yanlıştı. Formül artık dosyada yazılı olduğu için bir sonraki kapsam değişikliğinde elle yeniden hesaplanması gerekmiyor. |
| **D99** | M5 kampanyanın **%56'sı**. Kalibre aday sayısını 4→3 indirmek (`A-sign`'ı popülasyon düzeyinde bırakmak; RQ1 zaten tasarım A ve B'de kuruluyor) ~920 yay tasarruf ediyor. **Bu bir kapsam kararı olarak işaretlendi, plan tek başına almadı** — koşu başlamadan karara bağlanmalı. | AÇIK | Freeze (D66) kapsam *eklemeyi* yasaklıyor; kapsam *azaltmak* kullanıcının kararı, planın değil. |
| **D100** | Tarayıcının bağlam penceresi ±90'dan ±200'e genişletilmişti; bunun gerçek kusur saklayıp saklamadığı denetlendi. **Saklamamış** — atlanan dört eşleşmenin dördü de açıkça "eski değer şuydu, yanlıştı" diyen açıklayıcı bağlam. Kapı dürüst. | KABUL | Kendi kapısını gevşeterek temizleyen bir denetim, denetim değildir; bu yüzden ayrıca doğrulandı. |
| **D101** | Sayısal çapraz denetim eklendi: eşikler (`κ≥0.7`, `f̂≥0.15`, `f≥0.33`, `g_E<0.10`, %2 bandı), sayımlar (14 ablasyon, 26/8 panel, 9 hipotez), maliyetler. Beş uyarı çıktı, **beşi de yanlış pozitif** (regex artefaktı ve açıkça geri çekilmiş `%21–53`). | KABUL | Terim taraması sayı tutarsızlığını yakalayamıyordu; ikisi ayrı kategori. |

---

## Tur 15 — henüz hiç bakılmamış dört kategori (2026-08-08)

İki kapı da temiz olduğu için **kapıların kapsamadığı** yerlere bakıldı:
kaynakça kullanımı, makro kullanımı, şekil adlandırması ve **birinci turda
yazılan açık soruların hâlâ geçerli olup olmadığı**.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D102** | **Açık soruların kendisi eskimişti.** `Q3` hâlâ mutlak `Δt_acc` soruyordu — D68 ile ölçüt `n_s = τ_corr/Δt_i` oldu. `Q6` (vektör tablosu depolaması) D44/D55/D61 ile **zaten cevaplanmıştı** ama açık listesinde duruyordu. Q3 güncellendi, Q6 üstü çizilip kapatıldı. | KABUL | Açık sorular listesi, koşu sırasında "neye karar verilmedi" diye bakılacak yer. İçinde cevaplanmış soru olması, o listeyi güvenilmez yapar. On dört turdur hiç denetlenmemişti. |
| **D103** | **`Q11` eklendi:** M5'te kalibre aday 4 mü 3 mü (D99). Karar tablosunda AÇIK işaretliydi ama açık sorular listesinde yoktu. | KABUL | İki yerde tutulan bir listenin ikisinin de aynı şeyi söylemesi gerekiyor. |
| **D104** | **Üç kaynakça kaydı atıfsızdı** — `joe2008sobol`, `montenbruck2000satellite`, `vallado2013fundamentals`. Bu bir temizlik meselesi değil, **eksik atıf** sinyaliydi: tasarımlar karıştırılmış Sobol ama §5 bunu atıflandırmıyordu, ve pertürbasyon modeli konvansiyonlarının kaynağı yoktu. Üçü de yerlerine kondu; Sobol atfı ayrıca "noktalar bağımlı, bu yüzden p-değeri yok" cümlesini kaynağına bağlıyor. | KABUL | Atıfsız kayıt bibtex çıktısını bozmaz, o yüzden derleme temiz görünüyordu — ama içerik boşluğuydu. |
| **D105** | Ölü makro `\Rtwo` kaldırıldı (tanımlıydı, hiç kullanılmıyordu). Şekil adlandırması denetlendi: 6 şekil ortamı, 6 dosya adı, hepsi `fig_oa_*` konvansiyonunda (NOTATION §7). | KABUL | — |

---

## Tur 16 — dört tanım düzeltmesi, sonra GERÇEKTEN freeze (2026-08-08)

Danışmanın son turu: yeni deney yok, mevcut tanımların matematiksel düzeltmesi.

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D106** | **`κ` kapı olmaktan çıktı, tanı oldu. G4 kaldırıldı.** `κ` medyanı ve korunan kesir raporlanır; `κ < 0.7` **uyarıdır, durdurma değil**. | KABUL | D74 `κ`'nın korunan kesir **olmadığını** kanıtlamıştı, ama ön-tescil hâlâ ona sert bir kapı bağlıyordu — mantıksal tutarsızlık ve **erken yanlış-negatif riski**: `κ=0.6` olup hatanın tamamı `z`'ye dik olabilir (kesir iyi), ya da `κ=0.9` olup küçük açısal hata tam `z` doğrultusunda projeksiyonu bozabilir. Uçurulabilirlik kararı zaten **G5b** (`f̂ ≥ 0.15`) ve o doğrudan doğru şeyi ölçüyor. Kapı sayısı 7→6. |
| **D107** | **`B+`'ın birimi tanımlandı.** Pilot yay ve prob doğrudan `Σ N²`; gravite dışı CPU (öngörücü, IFBDA süpürmeleri, `argmin`, defter tutma) **ölçülen orandan** çevrilir: `B+,eş = t_+ / (t_çekirdek/B2)`. | KABUL | `B2 + B+ = B_tot` denklemi saniyeyle `ΣN²`'yi topluyordu. Boyutsal olarak tutarsızdı ve hakemin ilk soracağı şeydi. |
| **D108** | **D99 kapatıldı: M5'te üç kalibre aday** — `C-plan`, `R-int`, `F-op` çapası. `A-sign` popülasyon düzeyinde propagate edilmez. M5 3980→**3060**, toplam 7150→**6230**. **Bedeli açıkça yazıldı:** `f` ve H3 yalnızca tasarım A/B/C'de raporlanır; diğer popülasyonlar H4/H4b verdiktini taşır. | KABUL | Kullanıcı kararı. RQ1 ve yakalama oranı zaten A/B/C'de kuruluyor; geniş-eliptik, stratumlar ve düşük-perilun **verdikt** testi, kıyas testi değil. |
| **D109** | `check_numbers.py` sertleştirildi: **hiç eşleşme yoksa artık hata** (eskiden sessizce geçiyordu — bir eşik bütün dosyalardan silinse kapı 0 dönerdi), ve kritik eşikler için **zorunlu dosya** listesi eklendi (bir yerde doğru olması, başka yerde eksik olmasını meşru kılmaz). | KABUL | Gönderim kapısının sessizce geçmesi, kapı olmamasından kötüdür. |

---

## Tur 17 — senkronizasyon yaması (2026-08-08)

| # | Karar | Durum | Gerekçe |
|---|---|---|---|
| **D110** | **`A-sens`'in kanonik tanımı PLAN/WP/NOTATION'da da `K_i`'ye geçirildi.** Üç dosya hâlâ `argmin Σ u_iᵀ Q_ii u_i + λW_qN²` (ızgaraya bağlı) yazıyordu; manuskript D106 turunda `Σ_i Δt_i Δa_iᵀ K_i Δa_i`'ye geçmişti. `A-sens` bilgi tanımı "yerel duyarlılık çekirdeği `K(t,t)`" oldu; `c` ile birlikte taşınan ağırlık her yerde `K_i`. Makaledeki slogan da "ignore the trajectory kernel / retain its local diagonal / retain the full cross-epoch coupling" olarak hassaslaştırıldı. | KABUL | **Terminoloji değil, algoritma sorunu:** WP'yi takip ederek kod yazan biri makaledekinden farklı bir ölçüt uygular. `check_stale.py`'nin bunu yakalayamamış olması D111'in gerekçesi. |
| **D111** | **QA betikleri iki yeni kanona göre güncellendi.** `check_stale.py`: (a) `w(t)` satırının hedefi `diag Q` → `K_i / diagonal kernel`; (b) yeni desen — tarihsel olmayan dosyalarda `A-sens = diag Q` veya ham `u_iᵀ Q_ii u_i` ölçüt olarak kullanılamaz; (c) yeni desen — `κ`'yı kapı/gate diye anan cümle. `check_numbers.py`: `"G4 esigi kappa"` → `"kappa tani referansi"` (G4 artık yok, sayı yalnızca uyarı eşiği/şekil referans çizgisi). | KABUL | D110'daki üç dosyalık kayma betikten sessizce geçti. Kapının yakalayamadığı bir kanon, kanon değildir. |
| **D124** ★ | **Önek yapısı iç kuadratür kuralını belirliyor: orta nokta + kenar şeması.** `S_j = Σ_{i≤j} u_i`'nin *düz* önek toplamı olması, `i` epokunun iç ağırlığının `j`'den bağımsız olmasını gerektiriyor. Trapez bunu ihlal ediyor (uç katsayısı `(t_j−t_{j−1})/2`) — yani `Q_ik = A_max(i,k)` yapısı bozulurdu. Dikdörtgen yapıyı koruyor ama **birinci mertebe**, ve `n_s = 2`'de (Nyquist) bu gerçek bir kayıp. Çözüm bedava: **biriktirme düğümleri hücre orta noktalarında, dış integral hücre kenarlarında.** İç ağırlık tam hücre genişliği → `j`'den bağımsız → yapı **aynen** korunuyor; kural orta nokta → **ikinci mertebe**; ve `S_j` tam olarak `j+1`. kenara kadarki integral, yaklaşığı değil. Dış integralde önek yok, trapez serbest. İki ağırlık kümesi de negatif değil → `Q ⪰ 0` korunuyor. Manuskript §3.3'e paragraf eklendi; `tda.grids` bunu uyguluyor. | KABUL | Makale `S_j`'yi düz önek toplamı olarak yazıyordu ama hangi kuralı zorunlu kıldığını söylemiyordu; uygulayan kişi trapez seçip yapıyı sessizce bozabilirdi. |
| **D125** | Izgara inceltmesi **eşdağıtım** (equidistribution, de Boor 1973) ile: yoğunluk `ρ = n_s/τ_corr` kırpılıp integre ediliyor, kenarlar `S(t)`'nin tam sayıları geçtiği yere konuyor. **Genişlikler değil yoğunluk kırpılıyor** — sonradan kırpma kenarları kaydırır ve komşuları sınır dışına atar. Tek artık: düğüm sayısı tam sayı olmak zorunda, `⌈S(T)⌉/S(T)` yeniden ölçeklemesi her hücreyi en fazla `1/M` oranında sıkıyor; sınır bu toleransla geçerli ve öyle belgelendi. Ayrıca `τ_corr` **en yüksek** aday dereceyle hesaplanmalı (ızgara bütün adaylarca paylaşılıyor, en ince doku en yüksek truncation'da). | KABUL | Alternatif (son kenarı yapıştırmak) tek bir hücreyi tam bir adım sınır dışına bırakırdı — daha kötü takas. |
| **D126** | Depo: `github.com/ayberkdt/oracle-adaptive` (public, `main`). Derleme çıktıları (`.aux`, `.log`, `.ruff_cache`, `__pycache__`) hariç; `paper/main.pdf` **bilerek izleniyor** (taslak inceleme nesnesi, okur TeX kurmak zorunda kalmasın). `.gitattributes` ile satır sonları LF'e normalleniyor. Commit'lerde **yalnızca** `Ayberk <141725692+ayberkdt@users.noreply.github.com>`; asistan imzası/etiketi yok, geçmiş bu kurala karşı tarandı. | KABUL | Ön-tescilin herkese açık ve zaman damgalı olması ön-tescilin **amacı**; koşulardan önce yayınlanması artı. |
| **D122** | **Kod denetimi turu: iki gerçek hata, üç sessiz zayıflık.** (1) `test_symplectic_defect_detects_a_broken_matrix` **geçmeyecek bir test**ti: `Φ = I + 1e-3·e₀e₃ᵀ` aslında simplektik (serbest sürüklenme kesmesi). Genel ölçüt türetildi — `E = e_a e_bᵀ` ancak `J e_a ∝ e_b` iken Hamiltoniyen, ve rank-1 için ikinci mertebe `EᵀJE` daima sıfır — üç pertürbasyon buna göre seçildi, kesmenin simplektikliği ayrı bir test olarak kaydedildi. (2) Beş `np.allclose` **varsayılan `rtol=1e-5`** ile çağrılmıştı; 1.8e6 m yarıçapta bu 18 metrelik tolerans demek, yani test bozuk bir propagatörde de geçerdi. `rtol=0` veya `array_equal` yapıldı. (3) Arşiv yükleyicileri `rev10`'da değil `rev3_common`'da; (4) `Ω_ay` yuvarlanmış literal yerine arşivin ifadesiyle yazıldı; (5) `estimate_syntheses_per_rhs` sayım veriyordu, maliyet değil — `estimate_rhs_work_factor` eklendi (`1+6(n_g/N)²`), Q13'ün ödünleşimini sayıya çeviriyor. | KABUL | "Geçen test" ile "doğru test" aynı şey değil; ikisi de statik okumayla yakalandı, koşuyla değil. |
| **D123** | **README dizini on dört tur bayat kalmıştı** (`D1–D66`, `Q1–Q10`, `WP0–WP20`, `G1–G7`, `M0–M9`). İki kapı da bunu göremiyordu: `check_stale` bayat *sözcük*, `check_numbers` bayat *eşik* arıyor; güncelliğini yitirmiş bir *aralık beyanı* ikisi de değil. `check_numbers.py`'ye `check_register_ranges()` eklendi — README'nin beyan ettiği D/Q aralığını sicilin gerçek en büyüğüyle karşılaştırıyor. Bayat hâle karşı negatif test edildi: iki bulguyla düşüyor. | KABUL | Kapının yakalayamadığı bir tutarlılık kuralı, kural değildir (D111 ile aynı ders, üçüncü kez). |
| **D119** ★ | **Q12 kapandı: `Φ` varyasyonel denklemlerden burada üretilir.** Arşiv `Φ`'yi saklamıyor (zorlanmış `δx` entegre ediyor, matris değil), dolayısıyla `tda.dynamics.propagate` 42 durumlu artırılmış sistemi tek geçişte koşuyor: `Φ̇=[[0,I],[G,0]]Φ`. **Maliyet ilan edildi:** merkezi-fark gradyanı çağrı başına 6 sentez ekliyor → varyasyonel yay ≈ **7× düz propagasyon** (eşit gradyan derecesinde); `estimate_syntheses_per_rhs` bunu manifest için hesaplıyor. Geri taşıma `Φ(t_0,t_i)` **matris tersi almadan**: akış Hamiltoniyen olduğu için `Φ⁻¹ = −JΦᵀJ` (blok devrik). Aynı özdeşliğin artığı **bedava doğruluk tanısı**; gradyan simetrikleştirildiği için artık entegratörü ölçüyor, fark gürültüsünü değil. | KABUL | Kampanyanın tek kritik girdisi hiç bütçelenmemişti. Tersini almak yerine simplektik yapıyı kullanmak hem çarpanlamayı hem koşullanma sorusunu kaldırıyor. |
| **D120** | **Paylaşılan bant yığını mevcut çekirdekte YOK.** `sh_accel_fixed_numba` (sabit derece) ve `_compute_sh_acceleration_dual_numba` (iki derece, tek Legendre geçişi) var; **kümülatif-dereceye-göre çıkış yok.** Makalenin §7.2'deki "bir yığın, bütün adaylar kısmi toplamla" maliyet modeli bu giriş noktasını gerektiriyor. Şimdilik `DifferencingBandStack` doğru sayıları derece başına bir sentezle veriyor ve `syntheses_per_point` bunu **dürüstçe** raporluyor, böylece ölçülen ek yük koşan kodun ek yükü oluyor. | AÇIK — yeni çekirdek işi | Çekirdeğin akümülatörü zaten derece derece topluyor, eklemek küçük; ama **yeni iş** ve doğrulanmalı (bantların toplamı sabit-derece sonucunu makine hassasiyetinde vermeli). O doğrulanmadan §7.2'nin maliyet modeli alıntılanamaz. |
| **D121** | `Δa(t,N)` çift-derece çekirdekle **tek geçişte** hesaplanıyor. Naif iki çağrı Legendre yinelemesini iki kez ödüyordu; yüksek truncation zaten düşüğün her terimini üretiyor. WP1'in tablo kurulumu M1'in en büyük tek hesabı olduğu için bu kozmetik değil. | KABUL | Arşivde hazır duran giriş noktası; kullanılmıyordu. |
| **D116** ★ | **`A-sign` ile `C-plan` aynı statüde değil; M2b dallanma aşaması eklendi.** `A-sign` yörünge-hatası fonksiyonelinden doğrudan türediği için savunulur; **`C-plan`'ın en iyi uçurulabilir çözüm olduğu iddiası geri çekildi** — ilk yorumlanabilir, denetlenebilir, düşük-karmaşıklıklı denemedir. WP21 (T1–T7, **0 yay**) hangi denetleyici ailesinin gerektiğini M3'ün 670 yayı harcanmadan ölçer; dallanma tablosu ROADMAP §M2b'de, eşikler PREREG OA-02'de. Yeni net iş yalnız T3/T4/T6 (T1/T2 WP7'nin, T5 WP6'nın, T7 `abl-probe`'un yan ürünü). | KABUL | İki bilinmeyen ayrı: (1) yörünge-farkında tahsiste ödül var mı, (2) uçurulabilir olarak nasıl yakalanır. Fantezi bir denetleyiciyle ikisini birden denemek, kaybedince "fikir mi yanlıştı yoksa uygulamam mı" sorusunu cevapsız bırakır. Ayrıca `C-plan`'ı kesin çözüm ilan etmek şu an **bilimsel olarak yanlış** olurdu. |
| **D117** | **T3 (eşleşmenin etkin rankı) planın kendi verisinden çıkan yeni bir olasılığı açıyor: `C-rank1`.** `M_j` rank-3; `A_i`'nin özdeğer dağılımı hiç ölçülmedi ve STM'in baskın along-track yönü nedeniyle izin çoğunun 1–2 özyönde toplanması muhtemel. `p ≤ 2` çıkarsa sadeleşme durumu `p` skalerdir, denetleyici dondurulmuş 6×6 `c` tablosu yerine `p` skaleri **çevrimiçi** günceller — IFBDA'dan ucuz, pilot kaymasına dayanıklı. Aynı ölçüm DP/Viterbi'nin durum patlaması sorusunu da sayıyla kapatır. | KABUL | Ölçüm bedava (`A_i` zaten var) ve iki ayrı mimari soruyu birden cevaplıyor. Uygulama kararı **ölçümden sonra**. |
| **D118** | **T4 (ufuk yeterliliği) ve T6 (politika ifade edilebilirliği) kampanyadan bağımsız olarak da gerekli.** T4 eşleşmenin ne kadar ileri uzandığını söylüyor — makalenin şu an eksik olan bir fizik sonucu ve MPC sorusunun cevabı. T6a makalenin **kendi tezinin yanlışlama testi**: düzgün-durum politikası kazancın %60'ından fazlasını koruyorsa kazanç doku değil geometridir ve bant probu gereksizdir. Bu sonuç çıkarsa yeniden yorumlanmadan, reddedildiği hâliyle raporlanır. | KABUL | Bir bilgi iddiası ancak yanlışlanabilir bir testle birlikte iddiadır. Tasarım A'da uydur / tasarım B'de ölç kuralı sızıntıyı kapatıyor. |
| **D113** | **İkinci senkronizasyon süpürmesi (D110/D106 kalıntıları).** Ölü skaler `w` manuskriptte üç yerde yaşıyordu (§4 FW karşılaştırması, Alg. 3 girdisi, §7.5 faz-indeksleme gerekçesi) → `K_i`. `\mathbf Q_{ii}` çevrimdışı plan tablosunda kalmıştı → `K_i`. `κ` §7'de hâlâ "the declared gate", ve §7.7'de `\dnote{Gate G4 lives here}` → tanı diline çevrildi. §3'te `B_1` "per-call mean" yazıyordu (§5 ile çelişik) → zaman ortalaması. §4'te yetim bir cümle parçası iki paragrafa bölünmüştü → birleştirildi. | KABUL | Bir turda tanım değişince kalıntı **manuskriptin de içinde** kalıyor; önceki tur yalnızca `.md` tarafını süpürmüştü. |
| **D114** | **`\mathbf A` sembol çakışması giderildi.** `A_i` sonek bloğu (§3) ile `A_{40}(t)` dinamik matrisi (§7) aynı harfti. Varyasyonel denklem artık nomenclature'da tanımlı ama hiç kullanılmayan `G` (gravite gradyanı) ile yazılıyor: `Φ̇ = [[0,I],[G_40,0]]Φ`. Böylece hem çakışma hem ölü sembol kapandı. | KABUL | `S_i`/`S_j` ve `w`/`ω_j` ile aynı sınıf hata; üçüncüsü. |
| **D115** | `check_stale.py`'ye üç desen daha: çıplak `$w$` / `w\,d(` / `(\mathbf c,w`; İngilizce kapı dili (`declared gate`, `Gate G4`); çıplak `\mathbf Q_{ii}`. Altı gerçek kalıntının altısını da yakaladığı, yamadan önceki metinlere karşı sınanarak doğrulandı. | KABUL | D111'in desenleri yalnızca `.md` biçimlerini yakalıyordu; manuskriptteki altı kalıntı kapıdan geçti. |
| **D112** | Nomenclature'da `B = M·B₁` → **`B = B₁·T = Σ_q W_q N_q²`**; `B₁` tanımı "çağrı başına" yerine "iş oranı, zaman ortalaması `⟨N²⟩_t`" olarak §5 ile eşitlendi. Ayrıca `K(τ,σ)` girdisindeki `\tau` bir betik turunda sekme karakterine dönüşmüştü (`K(<TAB>au,σ)` — LaTeX hatası vermeden yanlış basıyordu), düzeltildi. | KABUL | `M·B₁` düzgün ızgara kalıntısı; D55/D60 ile kısıt zaman integrali oldu. PREREG/PLAN zaten yeni sözleşmedeydi. |

---

### Kapatılmayan, bilerek bırakılan

| Konu | Durum | Neden bırakıldı |
|---|---|---|
| Manuskriptte 21 float'ın 13'ü metinden referanslanmıyor (`fig:ladder`, `tab:beta1`, `alg:online`, …) | **açık** | Referans cümleleri `\ph{one paragraph}` bölgelerinin içine gelecek. Şimdi eklemek içeriksiz dolgu yazmak olurdu. **WP19'un teslim listesine yazıldı**; gönderim öncesi sıfır olmalı. |
| `r_RHS(t)`'nin perilunda yükselmesi prob ek yükü integraline katılmadı | **açık** | Bu bir ölçüm, tahmin değil: adım geçmişi arşivde var ve WP5/WP16 hesaplayacak. Mevcut ~%12–19 üst taraftan sapıyor ve öyle etiketlendi. |
| Anahtarlama sayısının arşivden sayılması | **açık** | WP16'ya yazıldı; şu anki ~0.5 çağrı/anahtarlama kestirim olarak etiketli. |

---

## Açık sorular

Bunlar karar değil, **cevabı ölçümden gelecek** sorular. Her biri bir WP'ye
bağlı; WP koşulduğunda buraya karar olarak dönerler.

| # | Soru | Cevaplayan | Cevap gelmezse |
|---|---|---|---|
| **Q1** | Prob derinliği `k` kaç olmalı? | WP5 (`κ` vs `k` vs maliyet) | k=3 varsayılan |
| **Q2** | `T_coh` yay uzunluğundan kısa mı? | WP6 | Kısaysa WP11 zorunlu |
| **Q3** | Izgara kaçta yakınsıyor? Ölçüt **`n_s = τ_corr/Δt_i`** (D68), mutlak `Δt_acc` değil | WP4 | `n_s = 2` varsayılan, maliyeti kabul edilir |
| **Q4** | IFBDA kaç iterasyonda yakınsıyor? | WP9 (`abl-J`) | J=2 varsayılan |
| **Q5** | Medyan `g_E` 0.10'un altına iner mi? Kaç yörünge boş sertifika verir? | WP7 | İnmezse D1/D29 uyarınca uzun ad |
| ~~**Q6**~~ | ~~Vektör tabloları nasıl saklanacak?~~ **KAPANDI (D44/D55/D61):** uyarlanabilir ızgara ~87 → ~9 MB, `float64`, bellek-eşlemeli sıralı erişim. | WP1 | — |
| **Q7** | `A-sign` kazancı propagasyonda hayatta kalıyor mu? | WP17 + H1b | Kalmıyorsa 🟠 bandı, güçlü negatif sonuç |
| **Q8** | Geniş-eliptik popülasyonda da kazanıyor mu? | WP15 | Kazanmıyorsa iddia rejim-sınırlı yazılır |
| **Q9** | Hangi dergi — JGCD mi CMDA mı? | WP19 | JGCD varsayılan (mevcut takım kalibre) |
| **Q10** | `𝒩` derece kümesi ne kadar ince olmalı? | WP3 (yayılım) | `rev14_oracle.degree_grid` + politika dereceleri |
| ~~**Q12**~~ | ~~Kıyasın `Φ`'si nereden geliyor?~~ **KAPANDI (D119):** varyasyonel denklemlerden, `tda.dynamics`'te; maliyet ≈7× düz propagasyon, geri taşıma simplektik özdeşlikle. Eski soru metni: **`Φ`'si nereden geliyor, hangi gradyan derecesinde, hangi maliyetle?** `Q` tamamen `Φ`'ye bağlı ama hiçbir belge onu üretmeyi tarif etmiyor. Arşiv `Φ`'yi **saklamıyor**: `rev13_variational_check.py` zorlanmış `δx`'i entegre ediyor, 6×6 matrisi değil. Yani WP1 her yörünge için 42 durumlu bir varyasyonel entegrasyon koşmak zorunda — ROADMAP ise M1'i "**0 propagasyon**" diye etiketliyor. Ayrıca `u_i` için `Φ(t_0,t_i) = Φ(t_i,t_0)⁻¹` gerekiyor; koşullanma hiç tartışılmamış (`S_j` için tartışılmışken). | **koşu başlamadan** — WP1'in kapsamı ve M1'in maliyeti buna bağlı | plan yanlış maliyetle başlar |
| **Q13** ★ | **Gradyan derecesi kıyası bozar mı?** Önceki makalenin R21/R33'ü tam olarak bunu ölçtü: derece-120 gradyanı **31 km perilunlu yörüngelerde** yetersizdi, R33 gradyanı referans dereceye çıkarmak için koşuldu. Yeni kampanyada **düşük-perilun popülasyonu var** ve `abl-stm` yalnızca *denetleyicinin* Kepler-vs-düşük-derece STM'ini kapsıyor; kıyasın gradyan derecesi için ne ablasyon ne kontrol var. Ölçülmüş, devralınmış bir risk taşınmamış. | WP1 + yeni bir kontrol | düşük-perilun sonucu savunulamaz |
| **Q14** | **İleri prob, en ucuz rakibine karşı gerekçelendirilmemiş.** §7.4 ileri probu yalnızca *aralık içi geriye dönük* kullanıma karşı savunuyor. Bir hakemin soracağı asıl bedava alternatif: **bir önceki devirde aynı fazda ölçülen sapmayı yeniden kullan.** Aritmetik: LLO'da devir başına yer izi kayması ≈33 km; `πr/N` = 48 km (N=120), 19 km (N=300), 10 km (N=600). Yani yeniden kullanım yalnız düşük derecede kısmen tutarlı, sadeleşmenin önemli olduğu yüksek derecede 2–3 kat dekorele. **Argüman bir satır ve ileri probu tercih olmaktan çıkarıp zorunluluk yapıyor** — ama şu an metinde yok. | tek paragraf + arşivden bir sayı | tercih, zorunluluk gibi sunulmuş kalır |
| **Q15** | **Hızlı alan değerlendirme literatürüne karşı konum yok.** Alan enterpolasyonu / sıkıştırılmış gravite modelleri yüksek dereceyi çok daha ucuza değerlendirebiliyorsa, dereceyi bütçelemenin bütün öncülü çöker. §7.1'in çözünürlük argümanı bunu *zaten* kapatıyor ama yalnızca **prob** için yazılmış; aynı argüman **alanın kendisi** için de kurulmalı, yoksa §2 eksik görünür. | §2 + §7.1'e bir paragraf | öncül savunmasız |
| ~~**Q11**~~ | ~~M5'te kalibre aday 4 mü 3 mü?~~ **KAPANDI (D108): 3.** Kampanyanın %56'sı orada; 4→3 (`A-sign`'ı popülasyon düzeyinde bırakmak) ~920 yay tasarruf eder | kullanıcı kararı (D99) | **koşu başlamadan karara bağlanmalı**; plan tek başına almaz |

---

## Devralınan, tartışılmayacak konvansiyonlar

Önceki makaleden gelen ve yeniden tartışılmayacak olanlar — gerekçeleri
[README.md](README.md)'de:

- Çözünürlük kuralı `M_res > 1`; kararsızlar hiçbir tarafa sayılmaz.
- p-değeri yok; yön bağımsız tasarımlarda tekrar ile kurulur.
- Oran istatistiği = yörünge başına oranların medyanı.
- Ön-tescil + hash + post-hoc etiketleme.
- Kabul (admissibility) kontrolü.
- Manifest digest zinciri.
- Sansür: politika referans derecesine ulaşırsa hücre sansürlenir.
- Commit mesajlarına asistan imzası eklenmez.
