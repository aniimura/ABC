import ABC3.Found.GenEll.ZetaSix

/-!
# GenEll 第 343 ブロック —— **★★★★★`Σ_{d>0} d⁻ᴷ` を `ζ(K)` と同定した**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★第 341 の右辺の第 1 因子を値にした

第 341 で

    Σ_{v ≠ 0} (v₁ + v₂τ)⁻ᴷ = (Σ_{d>0} d⁻ᴷ) · (Σ_{w 原始} (w₁ + w₂τ)⁻ᴷ)

を得た。★本ブロックは**第 1 因子を `ζ(K)` の値と同定**する:

    Σ_{d>0} d⁻⁴ = π⁴/90,    Σ_{d>0} d⁻⁶ = π⁶/945

## ★★★★段取り(2 つの橋)

1. **実 → 複素**: `hasSum_zeta_four`(実の級数)を `Complex.ofRealHom` で写す
   (`HasSum.map` + `Complex.continuous_ofReal`)。
2. **全体 → 正の部分**: `n = 0` の項は `(0^K)⁻¹ = 0` なので
   `hasSum_subtype_iff_of_support_subset` で `{d > 0}` に落とせる。

★★`K ≥ 1` があれば `0^K = 0` なので、台の包含はそのまま出る。

## ★★残り((i) の最後の 2 つ)

★`Σ_{w 原始} (w₁+w₂τ)⁻ᴷ = 2·E_K(τ)`(mathlib の `eisensteinSeries` との同定、索引の入れ替えつき)
★★`E₄³ − E₆² = 1728·Δ`(mathlib に無い。レベル 1 の次元公式と Sturm 境界から出る見込み)

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tsum_natPow_inv_eq` | ★★★★実の級数の値を複素の部分型和に移す |
| `tsum_natPow_inv_four` | ★★★★★**`Σ_{d>0} d⁻⁴ = π⁴/90`** |
| `tsum_natPow_inv_six` | ★★★★★**`Σ_{d>0} d⁻⁶ = π⁶/945`** |
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★実の級数から複素の部分型和へ -/

/-- ★★★★**実の級数の値を、複素数値の `{d > 0}` 上の和に移す**。

★`n = 0` の項は `(0^K)⁻¹ = 0` なので台の包含で部分型に落ちる。 -/
theorem tsum_natPow_inv_eq (K : ℕ) (hK : 1 ≤ K) (c : ℝ)
    (h : HasSum (fun n : ℕ => 1 / (n : ℝ) ^ K) c) :
    ∑' d : {d : ℕ // 0 < d}, (((d : ℕ) : ℂ) ^ K)⁻¹ = (c : ℂ) := by
  have hC : HasSum (fun n : ℕ => ((1 / (n : ℝ) ^ K : ℝ) : ℂ)) (c : ℂ) :=
    h.map (Complex.ofRealHom : ℝ →+* ℂ) Complex.continuous_ofReal
  have hC2 : HasSum (fun n : ℕ => (((n : ℂ)) ^ K)⁻¹) (c : ℂ) := by
    refine hC.congr_fun (fun n => ?_)
    push_cast
    rw [one_div]
  have hsupp : Function.support (fun n : ℕ => (((n : ℂ)) ^ K)⁻¹) ⊆ {d : ℕ | 0 < d} := by
    intro n hn
    simp only [Function.mem_support] at hn
    by_contra h0
    simp only [Set.mem_setOf_eq, not_lt, Nat.le_zero] at h0
    rw [h0] at hn
    simp [zero_pow (by omega : K ≠ 0)] at hn
  exact ((hasSum_subtype_iff_of_support_subset hsupp).2 hC2).tsum_eq

/-! ## ★★★★★`ζ(4)`・`ζ(6)` -/

/-- ★★★★★**`Σ_{d>0} d⁻⁴ = π⁴/90`**。 -/
theorem tsum_natPow_inv_four :
    ∑' d : {d : ℕ // 0 < d}, (((d : ℕ) : ℂ) ^ 4)⁻¹ = ((Real.pi ^ 4 / 90 : ℝ) : ℂ) :=
  tsum_natPow_inv_eq 4 (by norm_num) _ hasSum_zeta_four

/-- ★★★★★**`Σ_{d>0} d⁻⁶ = π⁶/945`**。 -/
theorem tsum_natPow_inv_six :
    ∑' d : {d : ℕ // 0 < d}, (((d : ℕ) : ℂ) ^ 6)⁻¹ = ((Real.pi ^ 6 / 945 : ℝ) : ℂ) :=
  tsum_natPow_inv_eq 6 (by norm_num) _ hasSum_zeta_six

/-! ## ★出典の紐付け(`.src`) -/

def tsum_natPow_inv_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def tsum_natPow_inv_four.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
