# oracle-adaptive

Yörünge-farkında (trajectory-aware) küresel harmonik derece tahsisi: `codebase`
makalesinin açıkça açık bıraktığı yönü bir yöntem hâline getirme çalışması.

## Neyin devamı

`codebase` makalesi (`../codebase`) şunu ölçtü: sabit bir gravite bütçesinde,
dereceyi yarıçapa göre dağıtmak kuvvet hatasını 4.6–5.0 kat düşürüyor ama
yedi günlük konum hatasını **büyütüyor**. Mekanizmayı da buldu: yörünge hatası
Φ'ye karşı **işaretli bir integral**, kuvvet defekt normu ise büyüklüklerin
toplamı; sıralamayı belirleyen gravite-gradyan eşleşmesini tam olarak
defekt normu atıyor.

O makale iki şey söyleyip bıraktı:

1. **Sonuç (§9):** yarıçaptan derece seçen bir kural, Φ(t_f,t) duyarlılığını,
   kuvvetin yönünü, yörünge fazını ve bir epokta işlenen hataya bağlı **kalan
   uçuş süresini** hesaba katamaz. Bunları birlikte kullanan bir kural
   tasarlamak ve doğrulamak "bu makalenin dışında".
2. **Tartışma (§8.x):** amaç fonksiyonu `E² = (1/T)∫‖H_r ∫Φ B Δa dτ‖² dt` (`H_r = [I₃ 0]`, D24)
   ayrıştırılabilir değil, dolayısıyla O26 tahsis kıyası (Lagrange gevşetmesi)
   basitçe Φ ile yeniden ağırlıklandırılıp yörünge seviyesine taşınamaz.

Bu klasör tam olarak o ikisini yapar: (1) yörünge seviyesinde bir tahsis
**kıyası** kurar, (2) ondan uçurulabilir bir **denetleyici** damıtmayı dener.

## Projenin kilitleri — altı kapı ve bir tanı, artan maliyet sırasıyla

Her kapı bir sonraki aşamanın masrafını riske atmadan önce sorulan bir soru.
Kapıda durursa maliyet kesildiği yerde kalır.

| Kapı | Soru | Aşama | Kümülatif maliyet | Düşerse |
|---|---|---|---|---|
| **G1** | Ödül var mı? `ρ̂(R-int, A-sign)` ≥ 1.5 | M1 | ~3 gece | 🟡 orta |
| **G2** | Enstrüman ölçebiliyor mu? kazanç/zarf ≥ 3 | M1 | ~3 gece | ⚫ enstrüman arızası |
| **G3** | Çizelge `Δt_acc` ile kararlı mı? | M1 | ~3 gece | Yakınsak ızgarada tekrar |
| *(tanı)* | Yön taşınabilir mi? `κ` — **kapı değil** (D106) | M1 | ~3 gece | uyarı; karar G5b'de |
| **G5b** | Bir aday ödülün %15'ini yakalıyor mu? `f̂ ≥ 0.15` | M2 | ~7 gece | 🟡, M3 harcanmaz |
| **G6** | Varyasyonel eleyici propagasyonla uyuşuyor mu? | M3 | ~10 gece | Teşhis; iki turdan sonra kaçış |
| **G7** | Kazanç β=1'e mi özgü? | **M4** | ~16 gece | M5+M6 iptal (~5700 yay), tek popülasyon |

Ayrıca `T_coh` (WP6) bir kapı değil ama **mimariyi belirler**: kısa çıkarsa
denetleyici alıcı-ufuklu hâle gelmek zorunda (WP11).

**M2b bir kapı değil, bir dallanmadır (D116).** `A-sign` yörünge-hatası
fonksiyonelinden doğrudan çıktığı için ona güveniyoruz; `C-plan` ise **ilk
yorumlanabilir denemedir**, en iyi uçurulabilir çözüm olduğu iddia edilmiyor.
WP21'in yedi testi (T1–T7, hepsi **0 yay**) hangi denetleyici ailesinin
gerektiğini M3'ün ~670 yayı harcanmadan **ölçer**: alıcı-ufuk MPC mi, düşük
ranklı `C-rank1` mi, vektör-vekil mi, yoksa `C-plan` mı. Tablo
[ROADMAP.md](ROADMAP.md#m2b) M3'ten önce okunur ve seçilen aile yazılı
olarak ilan edilir.

**G7 neden var:** önceki makalede iç üyenin bütün yapıcı sonucu β=0.5'te ters
dönmüştü. β=1'de kazanır kazanmaz popülasyon ve tam ızgara aşamalarına
**~5700 yay** harcamak aynı hatanın tekrarı olurdu; 576 yaylık prob bunu
önlüyor. Kampanyanın toplamı ~6230 yay (D108 ile M5 3 adaya indi); G7'de durursa ~1400.

Kapı olmayan ama ölçüm geçerliliğini belirleyen bir nokta daha: `Δa` yaklaşık
`πr/N` ölçeğinde dekorele oluyor (N=300 → ~18 km → ~11 s uçuş), dolayısıyla
120 s'lik ızgarada **işaretli** integral aliasing'e uğrayabilir. Bu yüzden
biriktirme ve karar ızgaraları ayrıldı (D9) ve yakınsama testi geç aşamadan
M1'e alındı (D8, G3).

## Üç araştırma sorusu

- **RQ1.** Yörünge-farkında bilgi altında ne kadar **erişilebilir kazanç** var?
- **RQ2.** Kazancın ne kadarı **duyarlılık ağırlığından**, ne kadarı **işaretli
  sadeleşmeden** geliyor?
- **RQ3.** Referans yay olmadan ve alanın yay boyunca **noktasal
  değerlendirmeleri** olmadan (yalnız `P_n` gibi çevrimdışı üst-veriyle) bunun
  anlamlı bir kısmı yakalanabilir mi?

İlk makale bu üçüyle sınırlıdır (karar D11).

## Dokümanlar — okuma sırası

| # | Dosya | İçerik |
|---|---|---|
| 1 | [PLAN.md](PLAN.md) | Bilimsel içerik: problem, algoritma, yöntem merdiveni L0–L4, açık problemler, makale iskeleti, şekil/tablo listesi |
| 2 | [WP.md](WP.md) | **İş paketi sicili — yapılacak her şey.** WP0–WP21 (WP17/WP20 a/b ayrık; WP21 = mimari seçim testleri): soru, girdi, çıktı, maliyet, kapı |
| 3 | [ROADMAP.md](ROADMAP.md) | Aşamalar M0–M9 (+ M2b dallanma), altı karar kapısı (G1, G2, G3, G5, G5b, G6, G7 — **G4 kaldırıldı, D106**), hipotezler H1a–H7, maliyet |
| 4 | [OUTCOMES.md](OUTCOMES.md) | Sonuç bantları, sahte-iyi artefakt tablosu, risk sicili, erken uyarı göstergeleri |
| 5 | [NOTATION.md](NOTATION.md) | **Dondurulmuş adlar ve tanımlar.** Politika kodları, maliyet düzeyleri, ızgaralar, terminoloji |
| 6 | [PREREG.md](PREREG.md) | Üç ön-tescilin tam metin taslağı (OA-01/02/03) + WP21 dallanma eşikleri |
| 7 | [DECISIONS.md](DECISIONS.md) | Karar günlüğü **D1–D140** (**D66 = kapsam donduruldu**) ve açık sorular **Q1–Q15** (Q6/Q11/Q12 kapalı) |
| 8 | [paper/](paper) | LaTeX manuskript taslağı — `main.tex`, `chapters/`, `references.bib` |
| 9 | [src/tda/](src/tda) | Kampanya kodu (aşağıdaki katman tablosu); [tests/](tests) altında birim testleri |

### Kod katmanları

Bağımlılıklar tek yönlü; hiçbir modül yukarı bakmıyor.

| Katman | Modül | Bağımlılık | Konu |
|---|---|---|---|
| 0 | `config` | — | donmuş sabitler + provenance hash |
| 0 | `stm` | — | `Φ` cebiri: simplektik ters, artık tanısı, satır/sütun dilimleri |
| 0 | `analytic` | — | kapalı-form alanlar (nokta kütle, `J₂`) |
| 1 | `field` | config | **korunan** alanı değerlendirmek: ivme, kesme kusuru, gradyan |
| 1 | `grids` | config | `τ_corr` eşdağıtımı, orta nokta/kenar şeması, karar ızgarası |
| 2 | `spectrum` | field | **ihmal edilen** kuyruk: derece varyansları, `γ`, bant yığını |
| 2 | `dynamics` | config, field, stm | referans yay + varyasyonel denklemler |

### Depo tek başına yeterli DEĞİL

`pip install -e ".[dev]" && pytest` **çalışır** — birim testleri bilerek
`tda.analytic`'in kapalı-form alanlarını kullanıyor, küresel harmonik
çekirdeğe ihtiyaç duymuyorlar. **Gerçek kampanya çalışmaz.** Üretim yolu
şunları istiyor ve hiçbiri burada değil:

| Bağımlılık | Ne için | Nerede |
|---|---|---|
| `lunaris` (`physics.spherical_harmonics`) | küresel harmonik sentez; **özel** `_compute_sh_acceleration_dual_numba` giriş noktası dahil | ayrı depo, **commit SHA pinlenmeli** |
| `rev3_common` (`load_model`, `kernel_args`, `warmup`, `OMEGA_MOON`) | model yükleme, çekirdek argümanları, gövde dönüşü | `../codebase/python_codes`, salt-okunur arşiv |
| GRAIL katsayı dosyası | alanın kendisi | arşivin çözdüğü yol |
| Arşiv `metrics/` kayıtları | WP0 kabul kontrolünün karşılaştırdığı sayılar | salt-okunur arşiv |

**Açık borç:** `lunaris` commit SHA'sı henüz pinlenmedi ve özel bir API'ye
bağımlıyız (`_`-önekli); ikisi de gönderim öncesi manifeste yazılmalı.
Lisans dosyası da henüz yok — public depo, karar bekliyor.

**Neden alt klasör yok.** Alt paket, ≥2 kardeş modülü olduğunda açılır; tek
dosyalık bir alt paket düz modülden kötüdür. Sırada olanlar zaten öyle
gelecek: `allocate/` (separable, descent, frankwolfe, rounding),
`controller/` (ifbda, plan, feedback), `analysis/` (WP21'in T1–T7'si).
Şu anki yedi modül tek bir taban katmanı ve düz durması doğru.

### Taslak makinesi (`paper/preamble.tex`)

| Makro | Anlamı | Bulma |
|---|---|---|
| `\ph{...}` | sayısal placeholder — kampanya `metrics/` makrosuyla değiştirilecek, elle yazılmayacak | `grep -rn '\\ph{' chapters/` |
| `\dnote{...}` | metnin **yönü** hangi sonuca bağlı — editoryal to-do değil | `grep -rn '\\dnote{' chapters/` |
| `\wpref{...}` | o pasajın sayılarını üreten iş paketi | — |

`\draftfalse` ile temiz okuma sürümü derlenir. Şu an: **28 sayfa**,
`\ph` sayısı 296, `\dnote` sayısı 11; tanımsız referans/atıf yok.
Bu üç sayı `check_numbers.py` tarafından **kaynaktan ve derleme çıktısından
yeniden hesaplanıp** buradaki beyanla karşılaştırılıyor; elle güncellenmeleri
gerekiyor ama bayat kalmaları mümkün değil.

> **Kapsam donduruldu (D66).** Bu noktadan sonra yeni kontrol, ablasyon, kapı
> veya WP eklenmez; yeni fikirler doğrudan PLAN §7 future work'e yazılır.
> Bundan sonrası ölçüm.

**Kural:** bir iş [WP.md](WP.md)'de yoksa yapılmaz; yapılması gerekiyorsa önce
oraya eklenir. Bir ad [NOTATION.md](NOTATION.md)'de yoksa kullanılmaz.

### Tutarlılık kapıları

Her değişiklikten sonra ikisi de çalıştırılır; çıkış kodu 0 olmalı.

```bash
python check_stale.py && python check_numbers.py
```

`check_stale.py` eskimiş **terimi** yakalar — on dört desen, her biri bir
kararla değiştirilmiş bir ifadeye karşılık geliyor (yarıçap sembolü, kuadratür
ağırlığı, prob maliyet modeli, genlik kaynağı, plan indeksi, depolama tipi,
eşleştirme sözleşmesi, aşama numaralandırması, sert dil). Desenlerin kendisi
betiğin içinde; burada tekrarlanmıyor, aksi hâlde tarayıcı kendi
dokümantasyonunu yakalar. `check_numbers.py` eskimiş **sayıyı** yakalar
(eşikler, panel büyüklükleri, aşama maliyetleri). İkisi ayrı kategori: bir eşik dosyalar arasında ayrışırsa
terim taraması göremez. `DECISIONS.md` ikisinde de hariç — karar günlüğü eski
ifadeyi **korumak zorunda**.

## Devralınan sözleşmeler

`codebase` kampanyalarının usul kurallarını **aynen** devralıyoruz; makalenin
tekrarlanabilirlik kredisinin çoğu buradan geliyor ve sıfırdan kurmanın anlamı
yok:

- **Çözünürlük kuralı:** bir karşılaştırma ancak iki hatanın farkı, toplanmış
  referans-dahil sayısal zarftan büyükse (`M_res > 1`) karara bağlanır.
  Karara bağlanmayanlar beraberlik değil, **kararsız**dır ve hiçbir tarafa
  sayılmaz.
- **İstatistik:** karıştırılmış Sobol noktaları bağımlı, çözünürlük kuralı
  örneklemi rastgele olmayan biçimde sansürlüyor. **p-değeri yok.** Yön,
  bağımsız tasarımlarda tekrar (replication) ile kurulur.
- **Maliyet:** dört basamak ayrı ayrı raporlanır (nominal çağrı başına ⟨N²⟩,
  gerçekten yapılan çağrı başına iş, RHS çağrı sayısı, gerçekleşen toplam
  kuadratik iş) + ölçülen çekirdek zamanı. Hiçbiri diğerinin yerine geçmez.
- **Ön-tescil:** hipotezler, ızgara, komparatör kuralı, sansür ve karar mantığı
  hiçbir toplu sonuca bakılmadan JSON'a yazılıp hash'lenir. Sonradan eklenen
  her şey **post-hoc** olarak etiketlenir.
- **Kabul (admissibility) kontrolü:** yeni bir koşu kabul edilmeden önce
  arşivlenmiş bir değeri aynı kod yolundan **birebir** yeniden üretmesi
  gerekir.
- **Manifest zinciri:** her kampanyanın çıktısı digest'lenir; sonraki kampanya
  öncekini yeniden hesaplamak yerine digest'e çivileyerek taşır.

## Kod ve veri

Bu klasör `../codebase` arşivini **salt-okunur** kullanır. Oradaki dosyalar
değiştirilmez; gerekli girdiler kopyalanır veya digest ile referanslanır.

Yeniden kullanılacak başlıca giriş noktaları:

- `../codebase/python_codes/rev14_oracle.py` — kuvvet seviyesi Lagrange tahsis kıyası (O26), `d(t,N)` tablolaması
- `../codebase/python_codes/rev42_variational_complete.py` — 128 yörüngelik zorlanmış varyasyonel çözüm (Φ, G)
- `../codebase/python_codes/rev14_budget_pareto.py` — bütçe kalibrasyonu, binleme, komparatör derecesi
- `../codebase/python_codes/rev18_span_sweep.py` — iç üye (k=0.5) ailesi
- `../codebase/python_codes/rev15_deployable_calibration.py` — pilot-yay ile öngörüsüz kalibrasyon (%8 maliyet)
- `../codebase/python_codes/population_registry.py` — tasarım A/B/C ve stratum kayıtları

## Adlandırma

Kampanya kodları `OA1, OA2, ...`; betikler `oa01_*.py` biçiminde. Çıktılar
`metrics/`, ön-tesciller `metrics/oaNN_preregistration.json`.

> Klasör adı boşluk yerine tire ile yazıldı (`oracle-adaptive`), çünkü içine
> PowerShell'den çağrılacak betikler girecek ve `../lunar-gravity-force-trajectory-gap`
> ile aynı düzeni izliyor. Boşluklu isim istenirse yeniden adlandırmak tek komut.
