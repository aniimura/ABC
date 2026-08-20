import ABC3.Skeleton.GaloisRep.WeilFunctionField
import ABC3.Found.GaloisRep.DvdCount
import ABC3.Found.GaloisRep.D2Bridge

/-!
# スケルトン —— **`div(f_P ∘ [n])` の因子計算**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★層 3 の `sorry` は **全部消えた**(2026-08-20)

| 節点 | 状態 |
|---|---|
| `dvd_count_pullback` (D1') | ✅ **第 149 ブロックで証明された** |
| `exists_fractionalIdeal_pow` (D1) | ✅ D1' + 第 142 から証明される |
| `exists_nthRoot_comp_mulByN`(`WeilRoot`) | ✅ D1 + D2 + 第 139 から証明される |
| `fractionalIdeal_isPrincipal` (D2) | ✅ **第 162-175 ブロック** |

★`#print axioms` で確かめた: `exists_nthRoot_comp_mulByN` は
`[propext, Classical.choice, Quot.sound]` の 3 つにしか依存していない。

### ★★D2 が閉じた道筋(第 162-175)

| ブロック | 内容 |
|---|---|
| 162-163 | 還元の核 `E₁` が部分群、**点の還元が群準同型** |
| 164 | **`P' = XYIdeal(n·Q_v)`**、台は `[n]⁻¹(P)` |
| 165-166 | `n` 乗根イデアルの明示、**`[J] = toClass(Σ e_v·Q_v)`** |
| 167 | **`hnn` は仮定でなく定理** |
| 168-171 | **`τ ∘ μ = μ`**、**素点の輸送**、**`e_v` の一定性** |
| 172-175 | 素点と点の対応、**台の総和が 0**、`n` 乗根の一意性 |

★★★**`[n]` の不分岐性(`e = 1`)を一度も使っていない**。
ファイバーと `E[n]∖{O}` の総和がそれぞれ独立に 0 だからである。

★D1' の証明は 2 つの場合の合流である:

| 場合 | 判定 | 内容 |
|---|---|---|
| A | `∀ r, w(μ r) ≤ 1` | 第 143(水準イデアルへの帰納法) |
| B | `∃ r, 1 < w(μ r)` | 第 144–149(無限遠、偶奇、超楕円対合) |

★★**分岐指数も、`deg[n] = n²` も、`#E[n] = n²` も、場所の分類定理も使っていない。**

## ★★★★★★不分岐性は要らない

一般の付値論で `ord_Q(μ f) = e_Q · ord_{[n]Q}(f)` であり、
`div(f_P) = n(P) − n(O)` から `ord_{[n]Q}(f_P) ∈ {n, −n, 0}`。
★**`e_Q` が何であっても `n` で割れる**。`deg[n] = n²` も `#E[n] = n²` も不要。

★★D2 でも `e_Q` は消える(平行移動不変性 + `Σ_{[n]Q=P} Q = 0` + `Σ_{T∈E[n]} T = 0`)。

## ★★★★D1' の中身は 2 つの場合に分かれる

`ν := ord_v ∘ μ̃` を `F(W)` 上の付値とし、`ν` が `F[W]` 上非負かどうかで分ける。

* **場合 A**(`[n]Q` がアフィン点): `ν` は `F[W]` 上非負なので
  イデアルの最小値が定義でき、`(f_P) = I_P^n` から `ν(f_P) = n·ν(I_P)`。
* **場合 B**(`[n]Q = O`): `ν(x) < 0`。★第 129 の `z² = Ψ₂Sq(x)` と
  `deg Ψ₂Sq = 3`(第 142、奇数)から `2ν(z) = 3ν(x)`、よって **`ν(x)` は偶数**。
  ★★`p(x)` の付値は偶数倍、`q(x)z` の付値は奇数倍で**決して一致しない**ので
  `ν(p + qz) = min` となり、`ν` は `F(x)` 上への制限で一意に決まる。
  ★★★したがって超楕円対合で不変になり、ノルムを使って `ν(f_P) = −m·n`。

## ★★★★因子を幾何で読む準備は済んでいる

| 材料 | ブロック |
|---|---|
| `IsDedekindDomain F[W]` | 第 137 |
| `HeightOneSpectrum F[W] ↔ 曲線上の点` | 第 138 |
| `J^n = (u)` かつ `J = (g)` ⟹ `n` 乗根 | 第 139 |
| `Σ_{T ∈ E[n]} T = 0`、類が自明 ⟹ 単項 | 第 140 |
| `f_P·f_{−P} = c(x − x_P)^n`、`μ` が `F`-代数射 | 第 141 |
| `∀ v, n ∣ count_v ⟹ n` 乗、`deg Ψ₂Sq = 3` | 第 142 |
| `(f_P) = XYIdeal_P^n` | 第 126 |
| `μ`(群法則で固定) | 第 118・125 |
| 平行移動 `τ_T` | 第 120・124 |
| 単元は定数 | 第 128 |
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve WeierstrassCurve.Affine
open IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★★★**D1' —— 各素点での指数は `n` で割れる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが層 3 の残り 2 つのうちの 1 つである。★★2 つの場合(アフィン/無限遠)に分かれる。 -/
theorem dvd_count_pullback [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    haveI := isDedekindDomain_coordinateRing h2 W
    ∀ v : HeightOneSpectrum W.CoordinateRing,
      (n : ℤ) ∣ FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) := by
  haveI := isDedekindDomain_coordinateRing h2 W
  exact ABC3.Found.GaloisRep.dvd_count_pullback W h2 μ hμinj hμF h n hP fP hfP

/-- ★★★★★**D1 —— `(μ f_P)` は `n` 乗である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★D1'(指数の `n` 可除性)と第 142 ブロックから**証明される**。 -/
theorem exists_fractionalIdeal_pow [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n) (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField,
      J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) := by
  haveI := isDedekindDomain_coordinateRing h2 W
  exact exists_pow_of_dvd_count (spanSingleton_mu_ne_zero W n fP hfP μ hμinj)
    (dvd_count_pullback h2 W h n hn hP μ hμinj hμF hns hμP hμx hμy fP hfP)

/-- ★★★★★**D2 —— その `J` は単項である**(Abel–Jacobi)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Point.toClass` の単射性(mathlib)で、類が 0 であることから単項性が出る。
★★類の計算では `Σ_{Q ∈ [n]⁻¹(P)} Q = 0` と `Σ_{T ∈ E[n]} T = 0`(第 140)を使う。 -/
theorem fractionalIdeal_isPrincipal [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
    (hJ : J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) :
    ∃ g : W.FunctionField, J = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g := by
  haveI := isDedekindDomain_coordinateRing h2 W
  exact ABC3.Found.GaloisRep.fractionalIdeal_isPrincipal_proof h2 W h n hn hchar hP μ hμinj
    hμF hns hμP hμx hμy fP hfP J hJ

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def dvd_count_pullback.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——div(f_P ∘ [n]) の各素点での指数が n で割れること)",
    sectionId := "genell-thm-3-8" }

def dvd_count_pullback.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**2026-08-20: 本節点は第 149 ブロックで証明された**(`Found/GaloisRep/DvdCount.lean` の `dvd_count_pullback`)。場合 A(第 143)と場合 B(第 144-149)の合流である。★**分岐指数も `deg[n] = n²` も `#E[n] = n²` も場所の分類定理も使っていない**(0 ブロック)" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8.1 の証明(div(f_P∘[n]) の計算)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .implicitStep
      "★★★★★★**不分岐性は要らない**(2026-08-20 の見通し変更)。`ord_Q(μ f) = e_Q · ord_{[n]Q}(f)` で右の因子は `div(f_P) = n(P) − n(O)` から `{n, −n, 0}` のいずれか。★`e_Q` が何であっても `n` で割れる。`deg[n] = n²` も `#E[n] = n²` も**不要**である" 19,
    .implicitStep
      "★★★★★★**加法的な指数表示は mathlib にあった**(2026-08-20 実測)——`FractionalIdeal.count : FractionalIdeal R⁰ K → ℤ` と `count_pow`・`count_finprod`・`finite_factors`・`finprod_heightOneSpectrum_factorization'`。付値(`WithZero (Multiplicative ℤ)`、乗法的)を経由せずに済む(0 ブロック)" 19,
    .implicitStep
      "★★★★★**最終組み立ては第 142 ブロックで済**——`∀ v, n ∣ count_v(I) ⟹ ∃ J, J^n = I`(`exists_pow_of_dvd_count`)(0 ブロック)" 19,
    .implicitStep
      "★★★★**場合 A(`F[W] ⊆ O_ν`、`[n]Q` がアフィン点)**: `ν := ord_v ∘ μ̃` は `F[W]` 上非負なので、イデアルの最小値 `ν(I) := min_{a ∈ I} ν(a)` が定義でき、付値だから `ν(I·J) = ν(I) + ν(J)`。★`(f_P) = I_P^n`(第 126)から **`ν(f_P) = n·ν(I_P)`** が直ちに出る。分岐指数を経由しない(8-20 ブロック)" 19,
    .implicitStep
      "★★★★★**場合 B(`F[W] ⊄ O_ν`、`[n]Q = O`)**: `ν(x) < 0` となる。★第 129 の `z² = Ψ₂Sq(x)` と `deg Ψ₂Sq = 3`(第 142、**奇数**)から `2ν(z) = 3ν(x)`、`gcd(2,3) = 1` より **`2 ∣ ν(x)`**。★★`p(x) ∈ F[x]` の付値は `deg(p)·ν(x)`(偶数倍)、`q(x)z` の付値は奇数倍なので**決して一致せず**、`ν(p + qz) = min` が成り立つ。★★★したがって `ν` は `F(x)` 上への制限で一意に決まり、超楕円対合 `ι` で不変(8-20 ブロック)" 19,
    .implicitStep
      "★★★★★場合 B の帰結: ノルム `N(a) = a·ι(a) ∈ F[x]` に `ν` を当てると `2ν(a) = deg N(a)·ν(x)`。★`ν(x) = −2m` として **`ν(a) = −m·deg N(a)`**。★★`(f_P) = I_P^n` から `deg N(f_P) = dim_F(F[W]/I_P^n) = n` なので **`ν(f_P) = −m·n`**、`n` で割れる(8-20 ブロック)" 19,
    .implicitStep
      "★★★★**`f_P · f_{−P} = c·(x − x_P)^n` は第 141 ブロックで済**(mathlib の `XYIdeal_neg_mul` + 第 128)。場合 B でノルムの代わりに使うこともできる(0 ブロック)" 19,
    .implicitStep
      "★★★★**付値側の API も実測できた**(2026-08-20)——`HeightOneSpectrum.valuation K v : Valuation K (WithZero (Multiplicative ℤ))`、`valuation_of_algebraMap`(`R` の元では `intValuation` に一致)、`intValuation_le_pow_iff_dvd : v.intValuation r ≤ WithZero.exp (-n) ↔ v.asIdeal^n ∣ span {r}`、`valuation_lt_one_iff_mem`。★**「位数 ≥ n ⟺ イデアル所属」の橋が既にある**ので、場合 A の超距離不等式は `Ideal.add_mem` に帰着できる(0 ブロック)" 19,
    .implicitStep
      "★★★★★場合 A の証明の形(2026-08-20 に確定): `(f_P) = I_P^n` の生成元を `a^i b^{n−i} = f_P · s_i` と書くと、`(s_0,…,s_n) = F[W]`(整域での約去)。★`ν ≥ 0` なら `min_i ν(s_i) = 0` となり、`ν(f_P) = min_i [i·ν(a) + (n−i)·ν(b)] = n·min(ν a, ν b)`。★★**分岐指数も素点の対応も使わない**(5-10 ブロック)" 19,
    .implicitStep
      "★★`μ` の単射性を仮定に置いた。第 118 の `pointHom_injective_of_transcendental` から出るが、`x([n]·)` の超越性を別途示す必要がある(3-8 ブロック)" 19,
    .implicitStep
      "★逸脱の記録: `[IsAlgClosed F]` と `IsUnit (2 : F)` を仮定に足した。第 137(Dedekind 性)と第 139(定数の `n` 乗根)が使う。消費側((G5) `det_cyclotomic`)は `[IsAlgClosed L]`・`[CharZero K]` の下で述べられているので後続に影響しない" 19 ]

def exists_fractionalIdeal_pow.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——(μ f_P) が n 乗であること)",
    sectionId := "genell-thm-3-8" }

def exists_fractionalIdeal_pow.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★**本節点の `sorry` は消えた**(2026-08-20)——`dvd_count_pullback`(D1'、各素点での指数が `n` で割れること)と第 142 ブロック(`exists_pow_of_dvd_count`)から証明される(0 ブロック)" 19,
    .implicitStep
      "★★`(μ f_P) ≠ 0` は第 142 ブロックで済(`spanSingleton_mu_ne_zero`)——`μ` の単射性から出る(0 ブロック)" 19,
    .implicitStep
      "★★★`IsDedekindDomain F[W]` は第 137 ブロックで済。`[IsAlgClosed F]` と `IsUnit (2 : F)` を仮定に置いているのはそのためである(0 ブロック)" 19 ]

def fractionalIdeal_isPrincipal.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の類が自明であること)",
    sectionId := "genell-thm-3-8" }

def fractionalIdeal_isPrincipal.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★★**2026-08-20: 本節点は証明された**（第 162-175 ブロック）。★道筋は次の通り: 第 162-163 で点の還元が群準同型になり、第 164 で台が `[n]⁻¹(P)` になり、第 165-166 で `[J] = toClass(Σ e_v·Q_v)` になり、第 167-171 で `e_v` がファイバー上一定になり、第 172-175 で総和が 0 になった。★★**`[n]` の不分岐性（`e = 1`）は一度も使っていない**（0 ブロック）" 19,
    .implicitStep
      "★逸脱の記録: `hchar`（`∀ k, 1 ≤ k → k ≤ n → (k : F) ≠ 0`）を仮定に足した。第 150（`Σ_{T ∈ E[n]} T = 0`）が (G1) の `E[n] ≅ (ℤ/n)²` を使うためである。★Weil 対 `e_n` はそもそも `char ∤ n` を要求するので、消費側（`det_cyclotomic`、`[CharZero K]`）に影響はない" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.3.5(Abel–Jacobi)"
      (.absent "mathlib に Weil 対およびその構成要素は 0 件(2026-08-20)") 19,
    .implicitStep
      "★★★★★★**分岐指数 `e_Q` はここでも消える**——平行移動 `τ_T`(`T ∈ E[n]`)は `[n]∘τ_T = [n]` により `μ` を保つので `e_Q` は各ファイバー上で一定。共通因子 `e` を括り出すと `Σ_{Q ∈ [n]⁻¹(P)} Q = n²Q₀ = n·P = 0` と `Σ_{T ∈ E[n]} T = 0` から類は 0 になる。★`e` の値を知る必要はない" 19,
    .implicitStep
      "★★★★`Point.toClass` の**単射性は mathlib にある**(`Point.toClass_injective`、2026-08-20 実測)。類が自明なら単項であることは **第 140 ブロックで済**(`isPrincipal_of_classGroup_eq_one`)(0 ブロック)" 19,
    .implicitStep
      "★★★平行移動 `τ_T` が体の自己同型であることは **第 120・124 ブロックで済**(`translateAlgEquiv`・`exists_translateAut_all`)(0 ブロック)" 19,
    .implicitStep
      "★★★★`Σ_{T ∈ E[n]} T = 0` は **第 140 ブロックで済**(`sum_univ_eq_zero_of_addEquiv`)——`E[n] ≅ (ℤ/n)²`((G1) で取得済)を当てるだけ(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**ファイバーの総和は第 150 ブロックで済**——`Σ_{T ∈ E[n]} T = 0`(`sum_torsion_eq_zero`)、`#E[n] = n²`(`card_torsion`)、`Σ_{T} (Q₀ + T) = n·(n·Q₀)`(`sum_fiber_eq`)。★(G1) の `E[n] ≅ (ℤ/n)²` に第 140 を当てるだけで出た。**G1 が G5 の中で初めて実際に消費された箇所である**(0 ブロック)" 19,
    .implicitStep
      "★★★★★★**台の同定は第 151 ブロックで済んだ**(`count_ne_zero_iff`)——引き戻した素イデアル `P' := {a | w(μ a) < 1}` を使うと `count_v(μ f_P) ≠ 0 ⟺ I_P = P'` が**イデアルの等式**として書ける。剰余体も点の還元も経由しない(0 ブロック)" 19,
    .implicitStep
      "★★★★★**D2 に残るのは 1 段だけ**——`P'` が定める点が `n·Q_v` であること(`Skeleton/GaloisRep/PointReduction.lean` の `pullbackPrime_eq_xyIdeal_nsmul`)。これは「点の還元が群準同型である」ことの言い換えであり、(G7) 半安定モデルでも要る(20-45 ブロック)" 19,
    .implicitStep
      "★★★★(旧記述)因子の台の同定——`count_v(μ f_P) ≠ 0` となる `v`(= 曲線上の点 `Q`)を決めること。場合 A では `count_v(μ f_P) = n·min(count_v(μ a), count_v(μ b))`(第 143)なので、`count ≠ 0 ⟺ μ a と μ b が共に `v` で消える ⟺ `[n]Q = P` である。★これを言うには**点の還元(特殊化)写像**が要る: `μ` を剰余体 `κ(v) ≅ F`(第 138)へ落として `([μ x], [μ y]) = [n]Q` を示す(8-20 ブロック)" 19,
    .implicitStep
      "★★★2026-08-20 実測: **mathlib の `AlgebraicGeometry/EllipticCurve/Reduction.lean` は曲線の還元**(極小モデル・good/multiplicative/additive reduction)であって、**点の還元写像ではない**。点の特殊化は自前で積む必要がある。★ただし `Point.map`(環準同型に沿った関手性)は mathlib にあるので、剰余体への環準同型を作れば流せる" 19,
    .implicitStep
      "★★★2026-08-20 実測: mathlib の `Point.toClass` には `toClass_injective` と **`toClass_eq_zero : toClass P = 0 ↔ P = 0`** があるが、**全射性は無い**。したがって `[J]` を `toClass(R)` の形で得るには、`J` を `XYIdeal'` の積として**明示的に構成する**必要がある(= 台の同定)。★`toClass_eq_zero` があるので、構成できれば結論は 1 行である" 19,
    .implicitStep
      "★★アフィン座標環の類群は `Pic⁰(E)` である(無限遠点 `O` を落とすと次数が消える)。`toClass` が `(R) − (O)` の類を与えることを使って、因子の類を点の和として読む(3-8 ブロック)" 19,
    .implicitStep
      "★★因子を幾何の言葉で読む対応は **第 138 ブロックで済**(`exists_point_of_heightOneSpectrum`)(0 ブロック)" 19 ]

end ABC3.Skeleton.GaloisRep
