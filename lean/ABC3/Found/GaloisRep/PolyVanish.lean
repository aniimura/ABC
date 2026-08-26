import ABC3.Found.GaloisRep.TateDefectBound

/-!
# Galois (G6) 第 236 ブロック —— **★★★★★★評価から係数の消滅へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★解析の評価を代数に戻す道具

第 235 で `‖tateDefectTrunc (n+1) u w (uw)‖ ≤ C_n‖q‖^{n+1}` が出た。
第 225 の橋でこれは**分子 `P ∈ ℤ[A,W]` の値の評価**になる。
そこから `Wⁿ⁺¹ ∣ P` を出すには、次の二つの一般論が要る:

| 主張 | 内容 |
|---|---|
| `X_pow_dvd_of_norm_le` | ★★★★★★`‖f(w)‖ ≤ C‖w‖ᵐ`(小さい `w ≠ 0` で)なら `Xᵐ ∣ f` |
| `poly_eq_zero_of_infinite_zeros` | ★★★★無限個の点で消える多項式は 0 |

## ★★★★`Xᵐ ∣ f` の証明は `m` の帰納法

| 段 | 内容 |
|---|---|
| `m = 0` | `X⁰ = 1 ∣ f` |
| `m+1` | `f(0) = 0` から `X ∣ f`、`f = X·g` と書き、`g` について `m` の場合を使う |

★`f(0) = 0` は、`w ≠ 0` での評価を `w → 0` に飛ばして出す(`eval_zero_le_of_bound`)。
**多項式の評価は連続なので、`𝓝[≠] 0` での極限が `f(0)` に一致する**。
★★ここだけ解析(連続性と極限)が要る。あとは純粋に代数である。

★★★仮定を `w ≠ 0` に限ってあるのが要点である——`w = 0` を含めると
`‖f(0)‖ ≤ C·0ᵐ` が仮定から直接出てしまい、帰納法の各段で同じ形を作れなくなる
(`g` の側で `w = 0` の評価が要るのに、それは仮定からは出ない)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eval_zero_le_of_bound` | ★★★★`w → 0` の極限で `w = 0` の評価を得る |
| `X_pow_dvd_of_norm_le` | ★★★★★★**`‖f(w)‖ ≤ C‖w‖ᵐ` なら `Xᵐ ∣ f`** |
| `coeff_eq_zero_of_norm_le` | ★★★★★係数の消滅の形 |
| `poly_eq_zero_of_infinite_zeros` | ★★★★無限個の零点をもつ多項式は 0 |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★★★`w → 0` の極限 -/

/-- ★★★★**`w ≠ 0` での評価を `w = 0` に延ばす**——多項式の評価は連続である。 -/
theorem eval_zero_le_of_bound (g : Polynomial ℂ) (C r : ℝ) (hr : 0 < r) (m : ℕ)
    (h : ∀ w : ℂ, w ≠ 0 → ‖w‖ < r → ‖g.eval w‖ ≤ C * ‖w‖ ^ m) :
    ‖g.eval 0‖ ≤ C * ‖(0 : ℂ)‖ ^ m := by
  have hc : ContinuousAt (fun x : ℂ => g.eval x) 0 := g.continuous.continuousAt
  have t1 : Filter.Tendsto (fun w : ℂ => ‖g.eval w‖) (nhdsWithin 0 {(0 : ℂ)}ᶜ)
      (nhds ‖g.eval 0‖) := (hc.continuousWithinAt).tendsto.norm
  have hcont2 : ContinuousAt (fun w : ℂ => C * ‖w‖ ^ m) 0 :=
    ((continuous_norm.pow m).continuousAt).const_mul C
  have t2 : Filter.Tendsto (fun w : ℂ => C * ‖w‖ ^ m) (nhdsWithin 0 {(0 : ℂ)}ᶜ)
      (nhds (C * ‖(0 : ℂ)‖ ^ m)) := (hcont2.continuousWithinAt).tendsto
  have hball : ∀ᶠ w : ℂ in nhds (0 : ℂ), ‖w‖ < r := by
    filter_upwards [Metric.ball_mem_nhds (0 : ℂ) hr] with w hw
    simpa [Complex.dist_eq] using hw
  refine le_of_tendsto_of_tendsto t1 t2 ?_
  filter_upwards [self_mem_nhdsWithin, hball.filter_mono nhdsWithin_le_nhds] with w hw0 hwr
  exact h w hw0 hwr

/-! ## ★★★★★★評価から整除性へ -/

/-- ★★★★★★**`‖f(w)‖ ≤ C‖w‖ᵐ`(小さい `w ≠ 0` で)なら `Xᵐ ∣ f`**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem X_pow_dvd_of_norm_le (C r : ℝ) (hr : 0 < r) :
    ∀ (m : ℕ) (f : Polynomial ℂ),
      (∀ w : ℂ, w ≠ 0 → ‖w‖ < r → ‖f.eval w‖ ≤ C * ‖w‖ ^ m) →
      Polynomial.X ^ m ∣ f := by
  intro m
  induction m with
  | zero =>
    intro f _
    simp
  | succ m ih =>
    intro f h
    have h0 : ‖f.eval 0‖ ≤ C * ‖(0 : ℂ)‖ ^ (m + 1) := eval_zero_le_of_bound f C r hr (m + 1) h
    have h0' : f.eval 0 = 0 := by
      have hz : C * ‖(0 : ℂ)‖ ^ (m + 1) = 0 := by simp
      rw [hz] at h0
      exact norm_eq_zero.1 (le_antisymm h0 (norm_nonneg _))
    have hXd : Polynomial.X ∣ f := by
      rw [Polynomial.X_dvd_iff, Polynomial.coeff_zero_eq_eval_zero]
      exact h0'
    obtain ⟨g, hg⟩ := hXd
    have hgb : ∀ w : ℂ, w ≠ 0 → ‖w‖ < r → ‖g.eval w‖ ≤ C * ‖w‖ ^ m := by
      intro w hw0 hwr
      have hf := h w hw0 hwr
      rw [hg] at hf
      simp only [Polynomial.eval_mul, Polynomial.eval_X, norm_mul] at hf
      have hwpos : (0 : ℝ) < ‖w‖ := norm_pos_iff.2 hw0
      have hcalc : ‖w‖ * ‖g.eval w‖ ≤ ‖w‖ * (C * ‖w‖ ^ m) := by
        calc ‖w‖ * ‖g.eval w‖ ≤ C * ‖w‖ ^ (m + 1) := hf
          _ = ‖w‖ * (C * ‖w‖ ^ m) := by ring
      exact le_of_mul_le_mul_left hcalc hwpos
    obtain ⟨c, hc⟩ := ih g hgb
    exact ⟨c, by rw [hg, hc]; ring⟩

/-- ★★★★★**係数の消滅の形**。 -/
theorem coeff_eq_zero_of_norm_le (f : Polynomial ℂ) (C r : ℝ) (hr : 0 < r) (m : ℕ)
    (h : ∀ w : ℂ, w ≠ 0 → ‖w‖ < r → ‖f.eval w‖ ≤ C * ‖w‖ ^ m) (j : ℕ) (hj : j < m) :
    f.coeff j = 0 :=
  Polynomial.X_pow_dvd_iff.1 (X_pow_dvd_of_norm_le C r hr m f h) j hj

/-! ## ★★★★無限個の零点 -/

/-- ★★★★**無限個の点で消える多項式は 0**(ℂ の上)。 -/
theorem poly_eq_zero_of_infinite_zeros (f : Polynomial ℂ) (s : Set ℂ) (hs : s.Infinite)
    (h : ∀ x ∈ s, f.eval x = 0) : f = 0 := by
  refine Polynomial.eq_zero_of_infinite_isRoot f ?_
  refine Set.Infinite.mono ?_ hs
  intro x hx
  exact h x hx

/-! ## ★出典の紐付け(`.src`) -/

def X_pow_dvd_of_norm_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——評価から整除性へ)",
    sectionId := "genell-def-3-3" }

def poly_eq_zero_of_infinite_zeros.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——無限個の零点をもつ多項式は 0)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
