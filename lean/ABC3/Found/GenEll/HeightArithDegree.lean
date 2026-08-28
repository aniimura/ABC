/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightProductFormula
import ABC3.Found.GenEll.HeightIdealNorm
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★超平面因子の算術次数は素朴高さである —— 段 C2c 後半（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★これは何か —— `degFin + degArch = log H`

段 C2c の後半は「**超平面因子の高さが素朴高さ `log max|x_i|` であること**」である。
★`§9-865`（段 C2c-1）で、点 `xF` に沿った超平面因子の引き戻しが
**`(x₀/x_{i₀})`** であることが取れた。★★本ファイルはその**算術次数**を計算する:

    `deg_fin + deg_arch = log H(x)`

ここで

| 項 | 中身 |
|---|---|
| `deg_fin` | `log N((x₀/x_{i₀}))` ——`§9-865` の引き戻しイデアルのノルム |
| `deg_arch` | `Σ_v mult(v) · log( sup_i v(x_i) / v(x₀) )` ——Fubini–Study 型の計量 |
| `log H(x)` | 素朴高さ（`Height.mulHeight`） |

## ★★★機構 —— 在庫 2 本の引き算

★`§9-853`: `H(x) · N(span{x_i}) = ∏_v (sup_i v(x_i))^{mult}`
★★`§9-855`: `∏_v v(x₀)^{mult} = N((x₀))`（主因子の算術次数は消える）

★★★点がチャート `D₊(x_{i₀})` を通ることは、代数の言葉では
「どの `x_i` も `x_{i₀}` の倍元である」＝ `span{x_i} = (x_{i₀})` である
（`span_range_eq_span_of_dvd`）。
★★★★あとは `N((x₀)) = N((x₀/x_{i₀})) · N((x_{i₀}))` で打ち消すだけである。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★**点が超平面の上にないこと**（`x₀ ≠ 0`）が要る。これは因子の高さの定義そのものが
台の外でしか意味を持たないからで、原典も同じ約束である。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★イデアルのノルムの補助 -/

/-- ★**単項イデアルのノルムは 0 でない**（`§9-853` の `ι = Unit` の場合）。 -/
theorem absNorm_span_singleton_ne_zero (K : Type) [Field K] [NumberField K]
    (z : 𝓞 K) (hz : z ≠ 0) : ((Ideal.span {z}).absNorm : ℝ) ≠ 0 := by
  have h := absNorm_span_range_ne_zero K (ι := Unit) (fun _ => z) (by
    intro h0
    exact hz (congrFun h0 ()))
  simpa [Set.range_const] using h

/-- ★**ノルムは積を積に送る**（単項イデアルの形）。 -/
theorem absNorm_span_mul (K : Type) [Field K] [NumberField K] (z w : 𝓞 K) :
    (Ideal.span {z * w}).absNorm = (Ideal.span {z}).absNorm * (Ideal.span {w}).absNorm := by
  rw [← Ideal.span_singleton_mul_span_singleton, map_mul]

/-! ## ★★チャートを通ることのイデアル的な言い換え -/

/-- ★★**点がチャート `D₊(x_{i₀})` を通ることの代数的な形**。

★「どの `x_i` も `x_{i₀}` の倍元である」——これはちょうど
`x_i/x_{i₀} ∈ 𝓞_K`、すなわち点が `Spec 𝓞_K` 上でチャートに入ることである。 -/
theorem span_range_eq_span_of_dvd (K : Type) [Field K] [NumberField K] {ι : Type}
    (x : ι → 𝓞 K) (i₀ : ι) (h : ∀ i, ∃ w : 𝓞 K, x i = w * x i₀) :
    Ideal.span (Set.range x) = Ideal.span {x i₀} := by
  apply le_antisymm
  · rw [Ideal.span_le]
    rintro _ ⟨i, rfl⟩
    obtain ⟨w, hw⟩ := h i
    rw [SetLike.mem_coe, hw]
    exact Ideal.mul_mem_left _ _ (Ideal.subset_span rfl)
  · rw [Ideal.span_le]
    rintro _ rfl
    exact Ideal.subset_span ⟨i₀, rfl⟩

/-! ## ★★★アルキメデス側は Fubini–Study の和である -/

/-- ★★★**アルキメデス側は Fubini–Study 型の和である**。

    `log ∏_v (sup_i v(x_i))^{mult} − log ∏_v v(x₀)^{mult}
       = Σ_v mult(v) · log( sup_i v(x_i) / v(x₀) )`

★これが `archDeg`（`§9-802` の計量つき因子の次数）が読む形である。 -/
theorem log_arch_ratio (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (j : ι) (hj : x j ≠ 0) :
    Real.log (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult)
      - Real.log (∏ v : InfinitePlace K, v ((x j : K)) ^ v.mult)
    = ∑ v : InfinitePlace K, (v.mult : ℝ) * Real.log ((⨆ i, v ((x i : K))) / v ((x j : K))) := by
  have hxj : ∀ v : InfinitePlace K, (0:ℝ) < v ((x j : K)) := fun v =>
    (AbsoluteValue.pos_iff _).2 (by simpa using hj)
  have hsup : ∀ v : InfinitePlace K, (0:ℝ) < ⨆ i, v ((x i : K)) := fun v =>
    lt_of_lt_of_le (hxj v) (le_ciSup (Finite.bddAbove_range (fun i => v ((x i : K)))) j)
  rw [Real.log_prod (fun v _ => (pow_pos (hsup v) v.mult).ne'),
    Real.log_prod (fun v _ => (pow_pos (hxj v) v.mult).ne'),
    ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun v _ => ?_)
  rw [Real.log_pow, Real.log_pow, Real.log_div (hsup v).ne' (hxj v).ne']
  ring

/-! ## ★★★★★★★★★★★段 C2c 後半の本体 -/

/-- ★★★★★★★★★★**有限側とアルキメデス側の和は素朴高さである**（生の形）。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★`§9-853` と `§9-855` の引き算だけである。 -/
theorem log_absNorm_add_archHeight (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) (i₀ j : ι) (z : 𝓞 K)
    (hz : x j = z * x i₀) (hj : x j ≠ 0)
    (hI : Ideal.span (Set.range x) = Ideal.span {x i₀}) :
    Real.log (((Ideal.span {z}).absNorm : ℝ))
        + (Real.log (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult)
          - Real.log (∏ v : InfinitePlace K, v ((x j : K)) ^ v.mult))
      = Real.log (Height.mulHeight (fun i => (x i : K))) := by
  have hzne : z ≠ 0 := by rintro rfl; exact hj (by rw [hz, zero_mul])
  have hi0ne : x i₀ ≠ 0 := by rintro h0; exact hj (by rw [hz, h0, mul_zero])
  have h1 : Real.log (∏ v : InfinitePlace K, v ((x j : K)) ^ v.mult)
      = Real.log (((Ideal.span {x j}).absNorm : ℝ)) := by
    rw [prod_infinitePlace_eq_absNorm]
  have h2 : ((Ideal.span {x j}).absNorm : ℝ)
      = ((Ideal.span {z}).absNorm : ℝ) * ((Ideal.span {x i₀}).absNorm : ℝ) := by
    rw [hz, absNorm_span_mul]; push_cast; ring
  have h3 := log_mulHeight_add_log_absNorm K x hx
  rw [hI] at h3
  rw [h1, h2, Real.log_mul (absNorm_span_singleton_ne_zero K z hzne)
    (absNorm_span_singleton_ne_zero K (x i₀) hi0ne)]
  linarith [h3]

/-- ★★★★★★★★★★★**段 C2c 後半** —— `deg_fin + deg_arch = log H(x)`。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `log N((x₀/x_{i₀})) + Σ_v mult(v)·log( sup_i v(x_i)/v(x₀) ) = log H(x)`

★左の第 1 項が**有限部分**（`§9-865` の引き戻しイデアルのノルム）、
第 2 項が**アルキメデス部分**（Fubini–Study 型の計量）である。
★★これが「超平面因子の高さは素朴高さ `log max|x_i|` である」ことの中身である。
★★★仮定は 3 つだけ:
点がチャート `D₊(x_{i₀})` を通ること（`hI`）、
`x_j = z·x_{i₀}`（`z = x₀/x_{i₀}` が整であること）、
点が超平面の上にないこと（`hj`）。 -/
theorem log_degFin_add_degArch_eq_log_mulHeight (K : Type) [Field K] [NumberField K]
    {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) (i₀ j : ι) (z : 𝓞 K)
    (hz : x j = z * x i₀) (hj : x j ≠ 0)
    (hI : Ideal.span (Set.range x) = Ideal.span {x i₀}) :
    Real.log (((Ideal.span {z}).absNorm : ℝ))
        + ∑ v : InfinitePlace K,
            (v.mult : ℝ) * Real.log ((⨆ i, v ((x i : K))) / v ((x j : K)))
      = Real.log (Height.mulHeight (fun i => (x i : K))) := by
  rw [← log_arch_ratio K x j hj]
  exact log_absNorm_add_archHeight K x hx i₀ j z hz hj hI

/-! ## ★出典の紐付け(`.src`) -/

def span_range_eq_span_of_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(点がチャートを通ることの代数的な形)",
    sectionId := "genell-prop-1-4" }

def log_arch_ratio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(アルキメデス側は Fubini–Study 型の和である)",
    sectionId := "genell-prop-1-4" }

def log_absNorm_add_archHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(有限側とアルキメデス側の和は素朴高さである——生の形)",
    sectionId := "genell-prop-1-4" }

def log_degFin_add_degArch_eq_log_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(段 C2c 後半——deg_fin + deg_arch = log H(x))",
    sectionId := "genell-prop-1-4" }

def log_degFin_add_degArch_eq_log_mulHeight.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "mulHeight_mul_absNorm(H(x)·N(I) = ∏_v (sup_i v(x_i))^mult、§9-853)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mulHeight_mul_absNorm") 2,
    .citation "[ABC3]" "prod_infinitePlace_eq_absNorm(∏_v v(x)^mult = N((x))、§9-855)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prod_infinitePlace_eq_absNorm") 2,
    .citation "[ABC3]" "pullbackIdeal_hyperplane_point(引き戻しは (x₀/x_i)、段 C2c-1、§9-865)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pullbackIdeal_hyperplane_point") 2,
    .implicitStep
      ("★逸脱: 点が超平面の上にないこと(x₀ ≠ 0)が要る。" ++
       "因子の高さの定義そのものが台の外でしか意味を持たないからで、原典も同じ約束である") 5,
    .implicitStep
      ("★★仮定 hI(点がチャート D₊(x_{i₀})を通ること)は " ++
       "span_range_eq_span_of_dvd で「どの x_i も x_{i₀} の倍元」から作る。" ++
       "★これは §9-C2b(exists_X_notMem)と §9-854(exists_integral_repr)から出る") 4,
    .implicitStep
      ("★★★残るのは段 C2c-3(Fubini–Study 型の計量が Definition 1.1 の意味の計量であること)" ++
       "であり、そこは §9-806 の htMetricAlg_sub_abs_le(Prop 1.4 (iii))が受ける") 4 ]

end ABC3.Found.GenEll
